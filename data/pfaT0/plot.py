#!/usr/bin/python

from __future__ import division

from pyx import *
from glob import glob
from math import pi,log,log


def slurp(filename, data=[]):
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if line == "" or line[0] == "#":
                continue
            # L/R, lmax, order, alpha, F(T=0)
            LbyR,lmax,order,alpha,F = map(float, line.split(","))
            ratio = F/pfa(LbyR)
            data.append((LbyR, F, ratio, ratio/bimonte(LbyR)))

    return data


if __name__ == "__main__":
    # use LaTeX for Pyx
    text.set(mode = "latex")

    # PFA formula
    pfa = lambda x: -pi**3/720.*(x+1)/x**2

    # formula Bimonte
    theta1,theta2 = 1./3 - 20/pi**2, -4.52
    bimonte = lambda x: 1+theta1*x+theta2*x**2*log(x)

    # read in data
    data = []
    for filename in glob("slurm-*.out"):
        slurp(filename, data)

    # sort data
    data = sorted(data, key=lambda x: x[0])

    # plot F/F_PFA and Bimonte
    g = graph.graphxy(
        width = 10,
        x = graph.axis.log(title=r"$x=L/R$", max=0.2, min=0.005),
        y = graph.axis.lin(title=r"$\mathcal{F}/\mathcal{F}_\mathrm{PFA}(T=0)$", min=0.8),
        key=graph.key.key(pos="tr", dist=0.1)
    )

    g2 = g.insert(graph.graphxy(
        width = 4,
        xpos = 1.5,
        ypos = 1,
        x = graph.axis.lin(max=0.008, min=0.005),
        y = graph.axis.lin(max=0.992, min=0.9875)
    ))


    attrs = [color.gradient.RedBlue]
    g.plot(
        graph.data.points(data, x=1, y=3, title="numerisch vs Bordag"),
        [graph.style.symbol(graph.style.symbol.circle, size=0.04, symbolattrs=attrs)]
    )

    f = "y(x)=1+%.15g*x+%.15g*x**2*log(x)" % (theta1,theta2)
    g.plot(graph.data.function(f, title="Bimonte et al."))

    g2.plot(
        graph.data.points(data, x=1, y=3, title="numerisch"),
        [graph.style.symbol(graph.style.symbol.circle, size=0.04, symbolattrs=attrs)]
    )

    g2.plot(graph.data.function(f))

    g.writePDFfile("pfaT0.pdf")


    # Plot (F/F_PFA) / Bimonte
    g = graph.graphxy(
        width = 10,
        x = graph.axis.log(title=r"$x=L/R$", max=0.02, min=0.004),
        y = graph.axis.lin(min=0.9995, max=1.0005, title=r"$\frac{\mathcal{F}(T=0)}{\mathcal{F}_\mathrm{PFA}(T=0)}/\left(1+\theta_1 x + \theta_2 x^2 \log x\right)$"),
    )

    g.plot(
        graph.data.points(data, x=1, y=4),
        [graph.style.symbol(graph.style.symbol.circle, size=0.06, symbolattrs=attrs), graph.style.line()]
    )

    g.finish()
    g.stroke(g.ygridpath(1), [style.linestyle.dashed])

    g.writePDFfile("ratio.pdf")
