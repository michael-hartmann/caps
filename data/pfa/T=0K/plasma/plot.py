import numpy as np
from pyx import *
from sys import argv

try:
    filename = argv[1]
except IndexError:
    print("Usage %s: filename" % argv[0])
    print("Example: python %s data_a10.csv" % argv[0])
    exit(1)

data = np.loadtxt(filename, delimiter=",")

theta1 = float("nan")
inva = float("nan")
omegap = float("nan")
with open(filename) as f:
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
    if x < 0.003:
        LbyR_fit.append(x)
        corr_fit.append(corr[i])
    else:
        LbyR_nofit.append(x)
        corr_nofit.append(corr[i])

LbyR_fit = np.array(LbyR_fit)
corr_fit = np.array(corr_fit)

y = corr_fit/LbyR_fit**1.5
x = np.sqrt(LbyR_fit)

theta3,theta2 = np.polyfit(x,y,1)

print("E ≈ E\_PFA ( 1 + θ_1·x + θ_2·x^1.5 + θ_3·x² )")
print("θ_1 = %+.8g" % theta1)
print("θ_2 = %+.8g" % theta2)
print("θ_3 = %+.8g" % theta3)

text.set(text.LatexRunner)

xmin,xmax = 5e-4, 1e-1
ymin,ymax = 3e-5, 1e-1
xtitle = r"$x = L/R$"
ytitle = r"$1-E/E_\mathrm{PFA}+\theta_1 x$"

g = graph.graphxy(
    key   = graph.key.key(pos="tl"),
    width = 8,
    x     = graph.axis.log(title=xtitle, min=xmin, max=xmax),
    y     = graph.axis.log(title=ytitle, min=ymin, max=ymax),
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

g.writePDFfile()
