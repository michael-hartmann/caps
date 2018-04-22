import argparse
import numpy as np
from pyx import *
from sys import argv

parser = argparse.ArgumentParser(description="Create plot for plasma model at T=0")
parser.add_argument("filename", help="path to CSV file generated by data.py")
parser.add_argument("--fit_min", help="points used for fit: FIT_MIN ≤ L/R ≤ FIT_MAX", type=float, default=0)
parser.add_argument("--fit_max", help="points used for fit: FIT_MIN ≤ L/R ≤ FIT_MAX", type=float, default=0.001)
parser.add_argument("--xmin", help="minimum of x-axis", type=float, default=4e-4)
parser.add_argument("--xmax", help="maximum of x-axis", type=float, default=1e-1)
parser.add_argument("--ymin", help="minimum of y-axis", type=float, default=2e-5)
parser.add_argument("--ymax", help="maximum of y-axis", type=float, default=1e-1)
parser.add_argument("--xtitle", help="label of x-axis", type=str, default=r"$x = L/R$")
parser.add_argument("--ytitle", help="label of y-axis", type=str, default=r"$1-E/E_\mathrm{PFA}+\theta_1 x$")
parser.add_argument("-o", "--output", help="output filename", type=str, default="plot.pdf")
args = parser.parse_args()

data = np.loadtxt(args.filename, delimiter=",")

theta1 = float("nan")
inva = float("nan")
omegap = float("nan")
with open(args.filename) as f:
    for line in f:
        if line[0] == "#":
            if "theta1" in line and "=" in line:
                theta1 = float(line.split("=")[1])
            if "omega_p*L/c" in line and "=" in line:
                inva = float(line.split("=")[1])
            if "omega_p" in line and "=" in line and "[eV]" in line:
                omegap = float(line.split("=")[1][:-5])

LbyR = data[:,0]
corr = data[:,7]

pnts = zip(LbyR,corr)

LbyR_fit, corr_fit = [], []
LbyR_nofit, corr_nofit = [], []

for i,x in enumerate(LbyR):
    if args.fit_min <= x <= args.fit_max:
        LbyR_fit.append(x)
        corr_fit.append(corr[i])
    else:
        LbyR_nofit.append(x)
        corr_nofit.append(corr[i])

if len(LbyR_fit) < 2:
    print("not enough points for fit: adjust FIT_MIN and FIT_MAX")
    print()
    parser.print_help()
    exit(1)


LbyR_fit = np.array(LbyR_fit)
corr_fit = np.array(corr_fit)

y = corr_fit/LbyR_fit**1.5
x = np.sqrt(LbyR_fit)

theta3,theta2 = np.polyfit(x,y,1)

print("fit for θ_2, θ_3 with %d points %g ≤ L/R ≤ %g:" % (len(LbyR_fit), args.fit_min, args.fit_max))
print()
print("E ≈ E_PFA ( 1 + θ_1·x + θ_2·x^1.5 + θ_3·x² )")
print("θ_1 = %+.8g" % theta1)
print("θ_2 = %+.8g" % theta2)
print("θ_3 = %+.8g" % theta3)

text.set(text.LatexRunner)

g = graph.graphxy(
    key   = graph.key.key(pos="tl"),
    width = 8,
    x     = graph.axis.log(title=args.xtitle, min=args.xmin, max=args.xmax),
    y     = graph.axis.log(title=args.ytitle, min=args.ymin, max=args.ymax),
)

cmd = "y(x) = %.12g*x**1.5 + %.12g*x**2" % (theta2,theta3)
title=r"$\theta_2 x^{3/2} + \theta_3 x^2$, $\theta_2=%.3g$, $\theta_3=%.3g$" % (theta2,theta3)
g.plot(graph.data.function(cmd, title))

attrs = [graph.style.symbol(graph.style.symbol.triangle, symbolattrs=[color.cmyk.Maroon], size=0.07)]
g.plot(graph.data.values(x=LbyR_fit, y=corr_fit, title=r"multipole (used for fit)"), attrs)

attrs = [graph.style.symbol(graph.style.symbol.circle, symbolattrs=[color.cmyk.MidnightBlue], size=0.07)]
g.plot(graph.data.values(x=LbyR_nofit, y=corr_nofit, title=r"multipole"), attrs)

x, y = g.pos(8e-2, 6e-5)
deltay = 0.5
g.text(x,y, r"$\theta_1=%.5g$" % theta1, [text.halign.boxright])
g.text(x,y+deltay, r"$\omega_P L/c = %g$" % inva, [text.halign.boxright])
g.text(x,y+2*deltay, r"$\omega_P = %g \, \mathrm{eV}$" % omegap, [text.halign.boxright])

g.writePDFfile(args.output)
