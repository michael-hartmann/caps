import argparse
import numpy as np
from sys import stderr, argv
from pyx import *

def get_radius(filename):
    with open(filename) as f:
        for line in f:
            if "# R = " in line and "nm" in line:
                return int(line.split("=")[1].strip()[:-2])

    assert False

parser = argparse.ArgumentParser(description="plot data")
parser.add_argument("-e", "--experimental",  help="experimental data", type=str)
parser.add_argument("-n", "--numerical",     help="numerical data", type=str)
parser.add_argument("-o", "--output", type=str, default="plot.pdf")
args = parser.parse_args()

text.set(text.LatexRunner)

experiment = np.loadtxt(args.experimental, delimiter=",")
numerics   = np.loadtxt(args.numerical,    delimiter=",")

R_exp = get_radius(args.experimental)
R_num = get_radius(args.numerical)
if R_exp != R_num:
    print("radii of experimental and numerical data differ: %d (exp) vs %d (num)" % (R_exp, R_num))
    exit(1)

g = graph.graphxy(
    width = 11,
    key   = graph.key.key(pos="tl"),
    x     = graph.axis.lin(title=r"$L$ (nm)", min=150, max=800),
    y     = graph.axis.lin(title=r"$F/F_\mathrm{PR}$")
)

g.plot(graph.data.values(x=experiment[:,0], y=experiment[:,1], title=None), [graph.style.symbol(graph.style.symbol.circle, size=0.07)])

g.plot(graph.data.values(x=numerics[:,0], y=numerics[:,1], title="Drude (exact)"), [graph.style.line([color.cmyk.MidnightBlue])])
g.plot(graph.data.values(x=numerics[:,0], y=numerics[:,2], title="Drude (PFA)"), [graph.style.line([color.cmyk.MidnightBlue, style.linestyle.dashed])])

g.plot(graph.data.values(x=numerics[:,0], y=numerics[:,3], title="plasma (exact)"), [graph.style.line([color.cmyk.BrickRed])])
g.plot(graph.data.values(x=numerics[:,0], y=numerics[:,4], title="plasma (PFA)"), [graph.style.line([color.cmyk.BrickRed, style.linestyle.dashed])])

g.text(10.9,0.3, r"$R=%.1f\,\mu\mathrm{m}$" % (R_num/1000), [text.halign.right])

g.writePDFfile(args.output)
