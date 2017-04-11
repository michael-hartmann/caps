import numpy as np
from pyx import *

data = np.loadtxt("out", delimiter=",")

text.set(text.LatexRunner)

L = data[:,0]
P_drude      = data[:,1]
P_plasma     = data[:,2]
P_drude_pfa  = data[:,3]
P_plasma_pfa = data[:,4]
ratio_drude  = data[:,5]
ratio_plasma = data[:,6]

attrs = [color.gradient.RedBlue]


x = L/151300
y_drude  = 1-ratio_drude
y_plasma = 1-ratio_plasma

g = graph.graphxy(
    width = 8,
    key = graph.key.key(pos="tl"),
    x = graph.axis.lin(title=r"$d/R$"),
    y = graph.axis.lin(title=r"$1-F^\prime/F^\prime_\mathrm{PFA}$")
)

g.plot([
    graph.data.values(x=x, y=y_drude, title="Drude"),
    graph.data.values(x=x, y=y_plasma, title="Plasma")
    ], [graph.style.symbol(graph.style.symbol.changecircle, symbolattrs=attrs, size=0.02)])

g.plot(graph.data.function("y(x)=0.4*x"))

g.writePDFfile()
