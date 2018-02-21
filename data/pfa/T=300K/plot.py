from pyx import *
from sys import stderr
from glob import glob
import numpy as np
import scipy.interpolate
from bimonte import F_HT_cut, E_HT_cut

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


def plotF(out="pfaF.pdf"):
    values = (0.01, 0.005, 0.003) # 1%, 0.5%, 0.3%

    xmin,xmax = 1, 300
    ymin,ymax = 5e0, 1e4
    xtitle = r"$R$ in $\mu\mathrm{m}$"
    ytitle = r"$L$ in $\mathrm{nm}$"

    g = graph.graphxy(
        width = 8,
        key   = graph.key.key(pos="tl"),
        x     = graph.axis.log(title=xtitle, min=xmin, max=xmax),
        y     = graph.axis.log(title=ytitle, min=ymin, max=ymax)
    )

    plots = []
    for value in values:
        title = r"$1-F/F_\mathrm{PFA} = %g$" % (value*100) + r"\%"
        #title = r"$%g$" % (value*100) + r"\%"
        data = cut_F(value)
        plots.append(graph.data.points(data, x=1, y=2, title=title))

    g.plot(plots, [graph.style.symbol(graph.style.symbol.changecircle, size=0.1), graph.style.line([style.linestyle.solid])])

    plots_ht = []
    for value in values:
        LbyR = F_HT_cut(value)
        cmd = "y(x)=x*%.8g*1e3" % LbyR
        plots_ht.append(graph.data.function(cmd, min=20, title=None))

    g.plot(plots_ht, [graph.style.line([style.linestyle.dashed, color.cmyk.Maroon])])

    g.writePDFfile(out)


def plotE(out="pfaE.pdf"):
    values = (0.005,0.01,0.02) # 1%

    xmin,xmax = 1, 300
    ymin,ymax = 5e0, 1e4
    xtitle = r"$R$ in $\mu\mathrm{m}$"
    ytitle = r"$L$ in $\mathrm{nm}$"

    g = graph.graphxy(
        width = 8,
        key   = graph.key.key(pos="tl"),
        x     = graph.axis.log(title=xtitle, min=xmin, max=xmax),
        y     = graph.axis.log(title=ytitle, min=ymin, max=ymax)
    )

    plots = []
    for value in values:
        title = r"$1-E/E_\mathrm{PFA} = %g$" % (value*100) + r"\%"
        #title = r"$%g$" % (value*100) + r"\%"
        data = cut_E(value)
        plots.append(graph.data.points(data, x=1, y=2, title=title))

    g.plot(plots, [graph.style.symbol(graph.style.symbol.changecircle, size=0.1), graph.style.line([style.linestyle.solid])])

    plots_ht = []
    for value in values:
        LbyR = E_HT_cut(value)
        cmd = "y(x)=x*%.8g*1e3" % LbyR
        plots_ht.append(graph.data.function(cmd, min=20, title=None))

    g.plot(plots_ht, [graph.style.line([style.linestyle.dashed, color.cmyk.Maroon])])

    g.writePDFfile(out)


if __name__ == "__main__":
    text.set(text.LatexRunner)

    plotF()
    plotE()
