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

    # read in data
    data = np.loadtxt("data_eta10.csv", delimiter=",")
    LbyR, E = data[:,0], data[:,3]
    correction = E/pfa(LbyR)-1-theta1*LbyR

    xmin, xmax = 2e-4, 0.1
    ymin, ymax = 8e-6, 1e-1

    xtitle = r"$x=L/R$"
    ytitle = r"$1-E/E_\mathrm{PFA}-\theta_1 x$"

    g = graph.graphxy(
        width = 8,
        key   = graph.key.key(pos="tl", dist=0.1),
        x     = graph.axis.log(title=xtitle, min=xmin, max=xmax),
        y     = graph.axis.log(title=ytitle, min=ymin, max=ymax)
    )

    g.plot([
        graph.data.values(x=LbyR, y=correction, title="multipole data"),
        ], [graph.style.symbol(graph.style.symbol.changecircle, size=0.05, symbolattrs=[color.cmyk.BrickRed])]
    )


    left, right = 0.0015, 0.007
    fit_x, fit_y = [], []
    for i,LbyR in enumerate(LbyR):
        if left <= LbyR <= right:
            fit_x.append(LbyR)
            fit_y.append(correction[i])

    # parameters obtained from fit
    theta2, theta3 = 2.59928774, -3.54471423
    cmd = "y(x)=%.8g*x**1.5+%.8g*x**2" % (theta2,theta3)
    g.plot(graph.data.function(cmd, title=r"$\theta_2 x^{3/2}+\theta_3 x^2$, $\theta_2=2.6$, $\theta_3=-3.54$"))

    g.plot([
        graph.data.values(x=fit_x, y=fit_y, title=None),
        ], [graph.style.symbol(graph.style.symbol.changecircle, size=0.05, symbolattrs=[color.cmyk.NavyBlue])]
    )

    # plot correction according to Bimonte et al.
    # Bimonte, Emig, Jaffe, Kardar, Casimir forces beyond the proximity approximation
    g.plot(graph.data.function("y(x)=-4.52*x**2*log(x)", title=r"Bimonte \textit{et al.}"))

    g.writePDFfile()
