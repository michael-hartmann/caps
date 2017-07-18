import numpy as np
from pyx import *

text.set(text.LatexRunner)

data = np.loadtxt("diff.csv")

g = graph.graphxy(
    width = 6,
    x     = graph.axis.lin(title=r"$L$ (nm)", min=150, max=800),
    y     = graph.axis.lin(title=r"$F_\mathrm{Palik}/F_\mathrm{Olmon}$")
)

g.plot(graph.data.points(data, x=1, y=2), [graph.style.line()])

g.writePDFfile()
