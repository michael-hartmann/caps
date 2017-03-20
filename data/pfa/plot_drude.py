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

    return list(sorted(data))


if __name__ == "__main__":
    desc = r"""$R=151.3\mu\mathrm{m}$,
    $\omega_p  = 8.9 \mathrm{eV}$,
    $\gamma    = 0.0357 \mathrm{eV}$
    """
    epsrel = r"$\varepsilon_\mathrm{rel} = 3\cdot 10^{-5}$"

    # read data
    data = slurp(glob("drude/*.out"))

    # use LaTeX for Pyx
    text.set(text.LatexRunner)

    g = graph.graphxy(
        width = 8,
        x     = graph.axis.log(title=r"$L/R$", min=0.00095),
        y     = graph.axis.lin(title=r"$\mathcal{F}/\mathcal{F}_\mathrm{PFA}$", max=1)
    )

    # L/R, ldim, F_PFA(T=0)*(L+R)/ħc), F(T=0)*(L+R)/(ħc), F/F_pfa
    g.plot(
       graph.data.points(data, x=1, y=5),
       [graph.style.line(), graph.style.symbol(graph.style.symbol.circle, size=0.07)]
    )

    g.finish()

    g.text(0.3, 0.8, desc,   [text.halign.flushright])
    g.text(0.3, 0.3, epsrel, [text.halign.flushright])

    g.writePDFfile()
