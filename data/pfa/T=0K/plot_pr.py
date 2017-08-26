from pyx import *
import numpy as np

# formula Bimonte: Bimonte, Emig, Jaffe, Kardar, Casimir forces beyond the proximity approximation
theta1 = 1/3 - 20/np.pi**2

# PFA formula
def pfa(LbyR):
    return -np.pi**3/720.*(LbyR+1)/LbyR**2


if __name__ == "__main__":
    # use LaTeX for Pyx
    text.set(text.LatexRunner)

    eta8  = np.loadtxt("pr/data_eta8.csv", delimiter=",")
    eta10 = np.loadtxt("pr/data_eta10.csv", delimiter=",")

    g = graph.graphxy(
        width = 10,
        key   = graph.key.key(pos="tl", dist=0.1),
        x     = graph.axis.log(title=r"$x=L/R$", min=3e-4),
        y     = graph.axis.log(title=r"$E/E_\mathrm{PFA}-1-\theta_1 x$")
    )

    LbyR, E = eta8[:,0], eta8[:,3]
    g.plot(
        graph.data.values(x=LbyR, y=E/pfa(LbyR)-1-theta1*LbyR, title=r"$\eta=8$"),
        [graph.style.symbol(graph.style.symbol.circle, size=0.05, symbolattrs=[color.cmyk.CadetBlue])]
    )

    LbyR, E = eta10[:,0], eta10[:,3]
    g.plot(
        graph.data.values(x=LbyR, y=E/pfa(LbyR)-1-theta1*LbyR, title=r"$\eta=10$"),
        [graph.style.symbol(graph.style.symbol.diamond, size=0.05, symbolattrs=[color.cmyk.BrickRed])]
    )

    g.writePDFfile()
