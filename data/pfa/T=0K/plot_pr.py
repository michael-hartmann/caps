#!/usr/bin/python3

from pyx import *
import numpy as np
from glob import glob

from scipy.optimize import curve_fit

# formula Bimonte: Bimonte, Emig, Jaffe, Kardar, Casimir forces beyond the proximity approximation
theta1 = 1/3 - 20/np.pi**2

# PFA formula
def pfa(x):
    return -np.pi**3/720.*(x+1)/x**2


def slurp(filenames):
    data = []
    for filename in filenames:
        with open(filename, "r") as f:
            absrel = float("nan")

            for line in f:
                line = line.strip()
                empty = line == ""
                comment = line.startswith("#")
                if comment:
                    index = line.find("absrel=")
                    if index > 0:
                        line = line[index+7:]
                        absrel = float(line)
                elif not empty:
                    # support old and new format
                    try:
                        # L/R, lmax, order, alpha, F(T=0)
                        LbyR,lmax,order,alpha,F = map(float, line.split(","))
                    except ValueError:
                        # L/R, L, R, T, ldim, F*(L+R)/(ħc)
                        LbyR,L,R,T,ldim,F = map(float, line.split(","))

                    ratio      = F/pfa(LbyR)
                    ratio_up   = ratio*(1+absrel)
                    ratio_down = ratio*(1-absrel)
                    data.append((LbyR, F, ratio, ratio-1, ratio-1-theta1*LbyR, ratio_up-1-theta1*LbyR, ratio_down-1-theta1*LbyR))

    return np.array(sorted(data, key=lambda x: x[0]))



if __name__ == "__main__":
    # use LaTeX for Pyx
    text.set(text.LatexRunner)

    data_eta8  = slurp(glob("pc/eta8/slurm-*.out"))
    data_eta10 = slurp(glob("pc/eta10/slurm-*.out"))

    for entry in data_eta10:
        print(entry[0], -100*entry[3])

    g = graph.graphxy(
        width = 10,
        key   = graph.key.key(pos="tl", dist=0.1),
        x     = graph.axis.log(title=r"$x=L/R$", min=3e-4),
        y     = graph.axis.log(title=r"$1-\mathcal{F}/\mathcal{F}_\mathrm{PFA}-\theta_1 L/R$")
    )

    g.plot(
        graph.data.points(data_eta8,  x=1, y=5, ymax=6, ymin=7, title=r"$\eta=8$"),
        [graph.style.symbol(graph.style.symbol.circle, size=0.05, symbolattrs=[color.cmyk.CadetBlue]), graph.style.errorbar()]
    )

    g.plot(
        graph.data.points(data_eta10, x=1, y=5, ymax=6, ymin=7, title=r"$\eta=10$"),
        [graph.style.symbol(graph.style.symbol.diamond, size=0.05, symbolattrs=[color.cmyk.BrickRed]), graph.style.errorbar()]
    )

    x,y = data_eta8[:,0], data_eta8[:,4]

    fit_left, fit_right = 2e-3, 8e-3
    fitx, fity = [], []
    for i,LbyR in enumerate(x):
        if fit_left <= LbyR <= fit_right:
            fitx.append(x[i])
            fity.append(y[i])


    fitfunc = lambda x, theta2: theta2*x**1.5
    params,pcov = curve_fit(fitfunc,fitx,fity,p0=(2.57,))
    theta2 = params
    cmd = "y(x) = %.10g*x**1.5" % (theta2)
    g.plot(graph.data.function(cmd, title=r"$\theta_2 x^{1.5}$, $\theta_2=%.3f$" % theta2))

    g.finish()
    g.stroke(g.xgridpath(fit_left), [style.linestyle.dotted])
    g.stroke(g.xgridpath(fit_right), [style.linestyle.dotted])

    g.writePDFfile()
