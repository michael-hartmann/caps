import numpy as np
from sys import argv
from pyx import *

text.set(text.LatexRunner)

filename_exp = argv[1]
filename_num = argv[2]
R = argv[3]

experiment = np.loadtxt(filename_exp, delimiter=",")
numerics   = np.loadtxt(filename_num, delimiter=",")

g = graph.graphxy(
    width = 11,
    key   = graph.key.key(pos="tl"),
    x     = graph.axis.lin(title=r"$L$ (nm)", min=150, max=800),
    y     = graph.axis.lin(title=r"$F/F_\mathrm{PR}$")
)

g.plot(graph.data.values(x=experiment[:,0], y=experiment[:,1], title=None), [graph.style.symbol(graph.style.symbol.circle, size=0.07)])

g.plot(graph.data.values(x=numerics[:,0], y=numerics[:,1], title="Drude (exact)"), [graph.style.line([color.cmyk.MidnightBlue])])
g.plot(graph.data.values(x=numerics[:,0], y=numerics[:,2], title="Drude (PFA)"), [graph.style.line([color.cmyk.MidnightBlue, style.linestyle.dashed])])

g.plot(graph.data.values(x=numerics[:,0], y=numerics[:,3], title="plasma (exact)"), [graph.style.line([color.cmyk.BrickRed])])
g.plot(graph.data.values(x=numerics[:,0], y=numerics[:,4], title="plasma (PFA)"), [graph.style.line([color.cmyk.BrickRed, style.linestyle.dashed])])

g.text(10.9,0.3, r"$R=%s\,\mu\mathrm{m}$" % R, [text.halign.right])

g.writePDFfile()
