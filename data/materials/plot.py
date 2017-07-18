import numpy as np
from pyx import *

hbar = 6.582119514e-16 # eV*s/rad

gradient = color.lineargradient_hsb(color.hsb(1, 1, 0.6), color.hsb(0.4, 1, 0.6))

dalvit = np.loadtxt("GoldDalvit.dat")
mixed  = np.loadtxt("GoldMixed.dat")
#palik  = np.loadtxt("GoldPalik.dat")

xmin, xmax = 1e11, 1e18
x2min, x2max = xmin*hbar, xmax*hbar

g = graph.graphxy(
    width = 8,
    key   = graph.key.key(),
    x2    = graph.axis.log(title=r"$\omega$ (rad/s)", min=xmin, max=xmax),
    x     = graph.axis.log(title=r"$\omega$ (eV)",    min=x2min, max=x2max),
    y     = graph.axis.log(title=r"$\epsilon(i\omega)$", min=1, max=1e8)
)

omegap, gamma = 9, 0.03
drude = "y(x) = 1+%g/(x*(x+%g))" % (omegap**2, gamma)

attrs = graph.style.line([gradient, style.linewidth.thick, attr.changelist([style.linestyle.solid, style.linestyle.dashed, style.linestyle.dashdotted])])

g.plot([
    graph.data.values(x=dalvit[:,0]*hbar, y=dalvit[:,1], title="Palik"),
    graph.data.values(x=mixed[:,0]*hbar,  y=mixed[:,1],  title="Olmon (SC) + Palik"),
    graph.data.function(drude, title=r"Drude")
], [attrs])


g.text(4.25,2.75, r"$\omega_P=9\,$eV, $\gamma=30\,$meV")


g.writePDFfile()
