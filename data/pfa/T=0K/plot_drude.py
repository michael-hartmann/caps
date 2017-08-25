import numpy as np
from pyx import *


if __name__ == "__main__":
    # read data
    drudeR10_35   = np.loadtxt("drude/R10mu_gamma35meV.csv", delimiter=",", usecols=(0,8))
    drudeR100_35  = np.loadtxt("drude/R100mu_gamma35meV.csv", delimiter=",", usecols=(0,8))
    drudeR100_350 = np.loadtxt("drude/R100mu_gamma350meV.csv", delimiter=",", usecols=(0,8))
    pr = np.loadtxt("pr/data_eta10.csv", delimiter=",", usecols=(0,5))

    # use LaTeX for Pyx
    text.set(text.LatexRunner)

    g = graph.graphxy(
        width = 8,
        key   = graph.key.key(pos="tl"),
        x     = graph.axis.log(title=r"$L/R$", min=3e-4, max=1e-1),
        y     = graph.axis.log(title=r"$1-\mathcal{F}/\mathcal{F}_\mathrm{PFA}$", min=3e-4, max=2e-1)
    )

    g.plot([
       graph.data.values(x=drudeR10_35[:,0],   y=1-drudeR10_35[:,1],  title=r"$R=10\mu\mathrm{m}$, $\gamma=35$meV"),
       graph.data.values(x=drudeR100_35[:,0],  y=1-drudeR100_35[:,1], title=r"$R=100\mu\mathrm{m}$, $\gamma=35$meV"),
       graph.data.values(x=drudeR100_350[:,0], y=1-drudeR100_350[:,1], title=r"$R=100\mu\mathrm{m}$, $\gamma=350$meV"),
       graph.data.values(x=pr[:,0], y=1-pr[:,1], title=r"PR")
       ], [graph.style.symbol(graph.style.symbol.changecircle, symbolattrs=[color.gradient.RedBlue], size=0.07)]
    )

    g.writePDFfile()
