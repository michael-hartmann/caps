from pyx import *

text.set(mode="latex")

g = graph.graphxy(
    width=8,
    x = graph.axis.linear(title=r"$T$"),
    y = graph.axis.linear(title=r"$\mathcal{F}$")
)
g.plot(graph.data.file("S", x=3, y=4))
g.writePDFfile("F.pdf")

g = graph.graphxy(
    width=8,
    x = graph.axis.linear(title=r"$T$"),
    y = graph.axis.linear(title=r"$S$")
)
g.plot(graph.data.file("S", x=3, y=5))
g.writePDFfile("S.pdf")
