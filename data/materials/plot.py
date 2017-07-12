import numpy as np
from pyx import *

hbar = 6.582119514e-16 # eV/s

gradient = color.lineargradient_hsb(color.hsb(1, 1, 0.6), color.hsb(0.4, 1, 0.6))

# read data of gold
dalvit = np.loadtxt("GoldDalvit.dat")
mixed  = np.loadtxt("GoldMixed.dat")
#palik  = np.loadtxt("GoldPalik.dat")

g = graph.graphxy(
    width = 6,
    key   = graph.key.key(),
    x     = graph.axis.log(title=r"$\omega$ (eV)", min=1e-4, max=1e3),
    y     = graph.axis.log(title=r"$\epsilon(i\omega)$", min=1, max=1e8)
)

drude = "y(x) = 1+%g/(x*(x+%g))" % (9**2, 0.035)
attrs = graph.style.line([gradient, attr.changelist([style.linestyle.solid, style.linestyle.dashed, style.linestyle.dashdotted])])
g.plot([
    graph.data.values(x=dalvit[:,0]*hbar, y=dalvit[:,1], title="Dalvit"),
    graph.data.values(x=mixed[:,0]*hbar,  y=mixed[:,1], title="Olmon + Palik"),
    graph.data.function(drude, title="Drude")
], [attrs])


g.writePDFfile()
