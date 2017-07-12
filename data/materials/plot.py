import numpy as np
from pyx import *

hbar = 6.582119514e-16 # eV/s

# read data of gold
dalvit = np.loadtxt("GoldDalvit.dat")
palik  = np.loadtxt("GoldPalik.dat")
mixed  = np.loadtxt("GoldMixed.dat")

g = graph.graphxy(
    width = 8,
    key   = graph.key.key(),
    x     = graph.axis.log(title=r"$\omega$ (eV)", min=1e-4, max=1e3),
    y     = graph.axis.log(title=r"$\epsilon(i\omega)$")
)

g.plot([
    graph.data.values(x=dalvit[:,0]*hbar, y=dalvit[:,1], title="Dalvit"),
    graph.data.values(x=palik[:,0]*hbar, y=palik[:,1], title="Palik"),
    graph.data.values(x=mixed[:,0]*hbar, y=mixed[:,1], title="Olmon + Palik")
], [graph.style.line([color.gradient.RedBlue])])

g.writePDFfile()
