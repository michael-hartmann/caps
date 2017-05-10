#!/usr/bin/python3

from pyx import *
import numpy as np
from glob import glob
from math import pi,log,exp

def slurp(filenames):
    data = []
    for filename in filenames:
        with open(filename, "r") as f:
            for line in f:
                line = line.strip()
                empty = line == ""
                comment = line.startswith("#")
                if not(empty or comment):
                    # support old and new format
                    try:
                        # L/R, lmax, order, alpha, F(T=0)
                        LbyR,lmax,order,alpha,F = map(float, line.split(","))
                    except ValueError:
                        # L/R, L, R, T, ldim, F*(L+R)/(ħc)
                        LbyR,L,R,T,ldim,F = map(float, line.split(","))

                    ratio = F/pfa(LbyR)
                    data.append((LbyR, F, ratio, ratio/bimonte(LbyR)))

    return np.array(sorted(data, key=lambda x: x[0]))


if __name__ == "__main__":
    plotPFA     = True
    plotBimonte = False

    # use LaTeX for Pyx
    text.set(text.LatexRunner)

    # PFA formula
    pfa = lambda x: -pi**3/720.*(x+1)/x**2

    # formula Bimonte: Bimonte, Emig, Jaffe, Kardar, Casimir forces beyond the proximity approximation
    theta1,theta2 = 1./3 - 20/pi**2, -4.52
    bimonte = lambda x: 1+theta1*x+theta2*x**2*log(x)

    data = slurp(glob("pc_gk/slurm-*.out"))

    attrs = [color.gradient.RedBlue]

    if plotPFA:
        # plot F/F_PFA and Bimonte
        x = data[:,0]
        y = (1-data[:,2])/x
        g = graph.graphxy(
            width = 10,
            x = graph.axis.log(title=r"$x=L/R$", min=4e-4),
            y = graph.axis.lin(title=r"$\mathcal{F}/\mathcal{F}_\mathrm{PFA}(T=0)$"),
            key=graph.key.key(pos="tr", dist=0.1)
        )

        #g2 = g.insert(graph.graphxy(
        #    width = 4,
        #    xpos = 1.7,
        #    ypos = 1,
        #    x = graph.axis.lin(max=0.0045, min=0.00045),
        #    y = graph.axis.lin(max=0.9998, min=0.9925)
        #))


        g.plot(
            graph.data.values(x=x, y=y, title="numerisch"),
            [graph.style.symbol(graph.style.symbol.circle, size=0.04, symbolattrs=attrs)]
        )

        #f = "y(x)=1+%.15g*x+%.15g*x**2*log(x)" % (theta1,theta2)
        #g.plot(graph.data.function(f, title="Bimonte et al."))

        #g2.plot(
        #    graph.data.points(data, x=1, y=3, title="numerisch"),
        #    [graph.style.symbol(graph.style.symbol.circle, size=0.04, symbolattrs=attrs)]
        #)

        #g2.plot(graph.data.function(f))

        g.writePDFfile("pfaT0.pdf")


    if plotBimonte:
        # Plot (F/F_PFA) / Bimonte
        g = graph.graphxy(
            width = 10,
            key   = graph.key.key(pos="bl"),
            x     = graph.axis.log(title=r"$x=L/R$", max=0.02, min=0.00045),
            y     = graph.axis.lin(min=0.9998, max=1.0003, title=r"$\frac{\mathcal{F}(T=0)}{\mathcal{F}_\mathrm{PFA}(T=0)}/\left(1+\theta_1 x + \theta_2 x^2 \log x\right)$"),
        )

        g.plot(graph.data.points(data, x=1, y=4, title=r"$\eta=8$"), [graph.style.symbol(graph.style.symbol.circle, size=0.06, symbolattrs=attrs), graph.style.line()])

        g.finish()
        g.stroke(g.ygridpath(1), [style.linestyle.dashed])

        g.writePDFfile("ratio.pdf")
