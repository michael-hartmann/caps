from pyx import *
from sys import stderr
from glob import glob
import numpy as np
import scipy.interpolate

def cut_E(value):
    points = []
    for directory in glob("R*/"):
        try:
            filename = directory + "/E.csv"
            dataE = np.loadtxt(filename, delimiter=",")
            L     = dataE[:,1]
            ratio = dataE[:,4]

            if ratio[0] <= value <= ratio[-1]:
                interp = scipy.interpolate.interp1d(ratio, L)
                R = dataE[0,0]*dataE[0,1]
                L = float(interp(value))
                points.append((R*1e6, 1e9*L))
            else:
                print("E: %g not in %s" % (value, filename))
        except FileNotFoundError:
            print("Could not find %s" % filename, file=stderr)

    return np.array(sorted(points))


def cut_F(value):
    points = []
    for directory in glob("R*/"):
        try:
            filename = directory + "/F.csv"
            dataF = np.loadtxt(filename, delimiter=",")
            L     = dataF[:,1]
            ratio = dataF[:,4]

            if ratio[0] <= value <= ratio[-1]:
                interp = scipy.interpolate.interp1d(ratio, L)
                R = dataF[0,0]*dataF[0,1]
                L = float(interp(value))
                points.append((R*1e6, 1e9*L))
            else:
                print("F: %g not in %s" % (value, filename))
        except FileNotFoundError:
            print("Could not find %s" % filename, file=stderr)

    return np.array(sorted(points))


if __name__ == "__main__":
    text.set(text.LatexRunner)

    text.preamble(r"\usepackage{textcomp}")

    xmin,xmax = 1, 300
    ymin,ymax = 4e0, 5e3

    g = graph.graphxy(
        width = 8,
        key   = graph.key.key(pos="tl"),
        x     = graph.axis.log(title=r"$R$ in $\mu\mathrm{m}$", min=xmin, max=xmax),
        y     = graph.axis.log(title=r"$L$ in $\mathrm{nm}$", min=ymin, max=ymax)
    )

    plots = []
    for value in (0.0025, 0.005,0.0075):
        title = r"$1-F/F_\mathrm{PFA} = %g$\textperthousand" % (value*1000)
        data = cut_F(value)
        plots.append(graph.data.points(data, x=1, y=2, title=title))

    g.plot(plots, [graph.style.symbol(graph.style.symbol.changecircle, size=0.1), graph.style.line([style.linestyle.solid])])

    g.writePDFfile()
