import numpy as np
from pyx import *

hbar = 6.582119514e-16 # eV*s/rad

olmonEV  = np.loadtxt("olmon_EV.csv", delimiter=",", usecols=(0,3))
olmonSC  = np.loadtxt("olmon_SC.csv", delimiter=",", usecols=(0,3))
olmonTS  = np.loadtxt("olmon_TS.csv", delimiter=",", usecols=(0,3))
palik    = np.loadtxt("palik.csv",    delimiter=",", usecols=(0,3))
brandli  = np.loadtxt("brandli.csv",  delimiter=",", usecols=(0,3))

xmin, xmax = 2e-3, 1e4
g = graph.graphxy(
    width = 8,
    key   = graph.key.key(pos="bl"),
    x     = graph.axis.log(title=r"$\omega$ (eV)", min=xmin, max=xmax),
    x2    = graph.axis.log(title=r"$\omega$ (rad/s)", min=xmin/hbar, max=xmax/hbar),
    y     = graph.axis.log(title=r"$\epsilon''(\omega)$")
)

g.plot(graph.data.points(palik, x=1, y=2, title="Palik"),
    [graph.style.symbol(size=0.07, symbolattrs=[color.cmyk.Peach])]
)

g.plot(graph.data.points(brandli, x=1, y=2, title=r"Br\"andli"),
    [graph.style.symbol(size=0.07, symbolattrs=[color.cmyk.OliveGreen])]
)

g.plot([
    graph.data.points(olmonEV, x=1, y=2, title="Olmon (EV)"),
    graph.data.points(olmonSC, x=1, y=2, title="Olmon (SC)"),
    graph.data.points(olmonTS, x=1, y=2, title="Olmon (TS)"),
    ], [graph.style.line([style.linestyle.solid, color.gradient.RedBlue])]
)

x2min, x2max = 0.04, 5
g2 = graph.graphxy(
    width = 2.5,
    xpos  = 5.3,
    ypos  = 3,
    x     = graph.axis.log(min=x2min, max=x2max),
    y     = graph.axis.log(min=1, max=2e4)
)
g2.plot([
    graph.data.points(olmonEV, x=1, y=2, title="Olmon (EV)"),
    graph.data.points(olmonSC, x=1, y=2, title="Olmon (SC)"),
    graph.data.points(olmonTS, x=1, y=2, title="Olmon (TS)"),
    ], [graph.style.line([style.linestyle.solid, color.gradient.RedBlue])]
)

g.insert(g2)
g.writePDFfile()
