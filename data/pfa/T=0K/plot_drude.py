import numpy as np
from pyx import *
from glob import glob

def slurp(filenames):
    data = []
    for filename in filenames:
        with open(filename, "r") as f:
            for line in f:
                line = line.strip()
                if len(line) and line[0] != "#":
                    data.append(list(map(float, line.split(","))))

    return np.array(sorted(data))


if __name__ == "__main__":
    # read data
    data = slurp(glob("drude/eta10/*.out"))

    LbyR = data[:,0]
    F    = data[:,5]

    # use LaTeX for Pyx
    text.set(text.LatexRunner)

    g = graph.graphxy(
        width = 8,
        x     = graph.axis.log(title=r"$L/R$"),
        y     = graph.axis.log(title=r"$-\mathcal{F}$")
    )

    # L/R, ldim, F_PFA(T=0)*(L+R)/ħc), F(T=0)*(L+R)/(ħc), F/F_pfa
    g.plot(
       graph.data.values(x=LbyR, y=-F),
       [graph.style.line(), graph.style.symbol(graph.style.symbol.circle, size=0.07)]
    )

    g.writePDFfile()
