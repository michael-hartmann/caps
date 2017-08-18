import numpy as np
from pyx import *


if __name__ == "__main__":
    # read data
    drudeR10  = np.loadtxt("drude/data_R10mu.csv", delimiter=",")
    drudeR100 = np.loadtxt("drude/data_R100mu.csv", delimiter=",")
    pc    = np.loadtxt("pc/data_eta10.csv", delimiter=",")

    # use LaTeX for Pyx
    text.set(text.LatexRunner)

    g = graph.graphxy(
        width = 8,
        key   = graph.key.key(pos="tl"),
        x     = graph.axis.log(title=r"$L/R$", min=3e-4, max=1e-1),
        y     = graph.axis.log(title=r"$1-\mathcal{F}/\mathcal{F}_\mathrm{PFA}$", min=3e-4, max=2e-1)
    )

    LbyR_drude10  = drudeR10[:,0]
    ratio_drude10 = 1-drudeR10[:,7]

    LbyR_drude100  = drudeR100[:,0]
    ratio_drude100 = 1-drudeR100[:,7]

    LbyR_pc  = pc[:,0]
    ratio_pc = 1-pc[:,5]
    g.plot([
       graph.data.values(x=LbyR_pc, y=ratio_pc, title=r"PR"),
       graph.data.values(x=LbyR_drude10, y=ratio_drude10, title=r"Drude, $R=10\mu\mathrm{m}$"),
       graph.data.values(x=LbyR_drude100, y=ratio_drude100, title=r"Drude, $R=100\mu\mathrm{m}$"),
       ], [graph.style.symbol(graph.style.symbol.changecircle, symbolattrs=[color.gradient.RedBlue], size=0.07)]
    )

    g.writePDFfile()
