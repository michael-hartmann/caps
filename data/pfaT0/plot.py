#!/usr/bin/python3

from pyx import *
from glob import glob
from math import pi,log,exp
from numpy import polyfit


def slurp_alt(filename, data=[]):
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            empty = line == ""
            comment = line.startswith("#")
            if not(empty or comment):
                # L/R, lmax, order, alpha, F(T=0)
                LbyR,T,F,error = map(float, line.split(","))
                ratio = F/pfa(LbyR)
                data.append((LbyR, F, ratio, ratio/bimonte(LbyR)))

    return data


def slurp(filename, data=[]):
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            empty = line == ""
            comment = line.startswith("#")
            if not(empty or comment):
                # L/R, lmax, order, alpha, F(T=0)
                LbyR,lmax,order,alpha,F = map(float, line.split(","))
                ratio = F/pfa(LbyR)
                data.append((LbyR, F, ratio, ratio/bimonte(LbyR)))

    return data


if __name__ == "__main__":
    plotPFA     = True
    plotBimonte = True

    # use LaTeX for Pyx
    text.set(text.LatexRunner)

    # PFA formula
    pfa = lambda x: -pi**3/720.*(x+1)/x**2

    # formula Bimonte: Bimonte, Emig, Jaffe, Kardar, Casimir forces beyond the proximity approximation
    theta1,theta2 = 1./3 - 20/pi**2, -4.52
    bimonte = lambda x: 1+theta1*x+theta2*x**2*log(x)

    # read in data
    data_eta7 = []
    for filename in glob("eta7/slurm-*.out"):
        slurp(filename, data_eta7)
    data_eta7 = sorted(data_eta7, key=lambda x: x[0])

    data_eta8 = []
    for filename in glob("eta8/slurm-*.out"):
        slurp(filename, data_eta8)
    data_eta8 = sorted(data_eta8, key=lambda x: x[0])


    attrs = [color.gradient.RedBlue]

    if plotPFA:
        # plot F/F_PFA and Bimonte
        g = graph.graphxy(
            width = 10,
            x = graph.axis.log(title=r"$x=L/R$", max=0.2, min=0.001),
            y = graph.axis.lin(title=r"$\mathcal{F}/\mathcal{F}_\mathrm{PFA}(T=0)$", min=0.8),
            key=graph.key.key(pos="tr", dist=0.1)
        )

        g2 = g.insert(graph.graphxy(
            width = 4,
            xpos = 1.7,
            ypos = 1,
            x = graph.axis.lin(max=0.008, min=0.001),
            y = graph.axis.lin(max=0.997, min=0.9875)
        ))


        g.plot(
            graph.data.points(data_eta8, x=1, y=3, title="numerisch"),
            [graph.style.symbol(graph.style.symbol.circle, size=0.04, symbolattrs=attrs)]
        )

        f = "y(x)=1+%.15g*x+%.15g*x**2*log(x)" % (theta1,theta2)
        g.plot(graph.data.function(f, title="Bimonte et al."))

        g2.plot(
            graph.data.points(data_eta8, x=1, y=3, title="numerisch"),
            [graph.style.symbol(graph.style.symbol.circle, size=0.04, symbolattrs=attrs)]
        )

        g2.plot(graph.data.function(f))

        g.writePDFfile("pfaT0.pdf")


    if plotBimonte:
        # Plot (F/F_PFA) / Bimonte
        g = graph.graphxy(
            width = 10,
            key   = graph.key.key(pos="bl"),
            x     = graph.axis.log(title=r"$x=L/R$", max=0.02, min=0.0007),
            y     = graph.axis.lin(min=0.9998, max=1.0003, title=r"$\frac{\mathcal{F}(T=0)}{\mathcal{F}_\mathrm{PFA}(T=0)}/\left(1+\theta_1 x + \theta_2 x^2 \log x\right)$"),
        )

        g.plot(
            [
                #graph.data.points(data_eta7, x=1, y=4, title=r"$\eta=7$"),
                graph.data.points(data_eta8, x=1, y=4, title=r"$\eta=8$")
            ],
            #graph.data.points(data2, x=1, y=4, title=r"numerisch (unterschiedliche $\ell_\mathrm{max}$)"),
            #graph.data.points(data_alt, x=1, y=4, title=r"numerisch (alt, $\ell_\mathrm{max}=6.2 L/R$)")],
            [graph.style.symbol(graph.style.symbol.circle, size=0.06, symbolattrs=attrs), graph.style.line()]
        )

        g.finish()
        g.stroke(g.ygridpath(1), [style.linestyle.dashed])

        g.writePDFfile("ratio.pdf")
