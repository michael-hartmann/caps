import numpy as np
from pyx import *
from sys import argv

text.set(text.LatexRunner)

filename = argv[1]
output = filename
if "." in output:
    output = filename[:output.rfind(".")]

data = np.loadtxt(filename, delimiter=",")

width  = 6
height = 6
gap = 0.2

xlabel = r"column number $j$"
ylabel = r"row number $i$"
zlabel = r"$\log_{10}\left|\mathcal{M}_{ij}^{(m)}(\xi)\right|$"

coloraxis = graph.axis.lin(title=zlabel)
kg1       = graph.graphx(ypos=width+gap, xpos=0, length=width, size=0.28*unit.v_cm, direction="horizontal", x2=coloraxis)
kg2       = graph.graphx(ypos=width+gap, xpos=8, length=width, size=0.28*unit.v_cm, direction="horizontal", x2=coloraxis)

g1 = graph.graphxy(
    height = height,
    width  = width,
    x      = graph.axis.lin(title=xlabel),
    y      = graph.axis.lin(title=ylabel)
)

g1.plot(graph.data.points(data, x=2, y=1, color=4),
       [graph.style.density(gradient=color.gradient.Jet, coloraxis=coloraxis, keygraph=kg1)])

g1.insert(kg1)

g2 = graph.graphxy(
    height = height,
    width  = width,
    xpos   = 8,
    x      = graph.axis.lin(title=xlabel),
    y      = graph.axis.lin(title=ylabel)
)

g2.plot(graph.data.points(data, x=2, y=1, color=3),
       [graph.style.density(gradient=color.gradient.Jet, coloraxis=coloraxis, keygraph=kg2)])

g2.insert(kg2)

g1.insert(g2)
g1.writePDFfile(output)
