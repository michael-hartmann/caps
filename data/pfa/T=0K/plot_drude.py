import numpy as np
from pyx import *


if __name__ == "__main__":
    # read data
    drude = np.loadtxt("drude/data.csv", delimiter=",")
    pc    = np.loadtxt("pc/data_eta10.csv", delimiter=",")

    # use LaTeX for Pyx
    text.set(text.LatexRunner)

    g = graph.graphxy(
        width = 8,
        key   = graph.key.key(pos="tl"),
        x     = graph.axis.log(title=r"$L/R$", min=3e-4, max=1e-1),
        y     = graph.axis.log(title=r"$1-\mathcal{F}/\mathcal{F}_\mathrm{PFA}$")
    )

    LbyR_drude  = drude[:,0]
    ratio_drude = 1-drude[:,7]

    LbyR_pc  = pc[:,0]
    ratio_pc = 1-pc[:,5]
    g.plot([
       graph.data.values(x=LbyR_pc, y=ratio_pc, title=r"PR"),
       graph.data.values(x=LbyR_drude, y=ratio_drude, title="Drude"),
       ], [graph.style.symbol(graph.style.symbol.changecircle, symbolattrs=[color.gradient.RedBlue], size=0.07)]
    )

    g.writePDFfile()
