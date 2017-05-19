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
                    data.append((LbyR, F, ratio, ratio-1, ratio-1-theta1*LbyR))

    return np.array(sorted(data, key=lambda x: x[0]))



if __name__ == "__main__":
    # use LaTeX for Pyx
    text.set(text.LatexRunner)

    data = slurp(glob("pc_gk/slurm-*.out"))

    attrs = [color.gradient.RedBlue]

    x,y = data[:,0], data[:,4]

    g = graph.graphxy(
        width = 10,
        key   = graph.key.key(pos="tl", dist=0.1),
        x     = graph.axis.log(title=r"$x=L/R$", min=3e-4),
        y     = graph.axis.log(title=r"$1-\mathcal{F}/\mathcal{F}_\mathrm{PFA}-\theta_1 L/R$")
    )

    g.plot(
        graph.data.values(x=x, y=y, title=None),
        [graph.style.symbol(graph.style.symbol.circle, size=0.04, symbolattrs=attrs)]
    )

    fit_left, fit_right = 2e-3, 8e-3
    fitx, fity = [], []
    for i,LbyR in enumerate(x):
        if fit_left <= LbyR <= fit_right:
            fitx.append(x[i])
            fity.append(y[i])


    fitfunc = lambda x, theta2, theta3: theta2*x**1.5 + theta3*x**2
    params,pcov = curve_fit(fitfunc,fitx,fity,p0=(2.57,-3))
    theta2,theta3 = params
    cmd = "y(x) = %.10g*x**1.5 %+.10g*x**2" % (theta2,theta3)
    g.plot(graph.data.function(cmd, title=r"$\theta_2 x^{1.5} + \theta_3 x^2$, $\theta_2=%.3f$, $\theta_3=%.3f$" % (theta2, theta3)))

    g.finish()
    g.stroke(g.xgridpath(fit_left), [style.linestyle.dotted])
    g.stroke(g.xgridpath(fit_right), [style.linestyle.dotted])

    g.writePDFfile()
