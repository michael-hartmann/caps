from pyx import *
import numpy as np

# formula Bimonte: Bimonte, Emig, Jaffe, Kardar, Casimir forces beyond the proximity approximation
theta1 = 1/3 - 20/np.pi**2

# PFA formula
def pfa(LbyR):
    return -np.pi**3/720.*(LbyR+1)/LbyR**2


def fit_theta23(LbyR, corr):
    LbyR_fit = np.array(LbyR)
    corr_fit = np.array(corr)

    x = np.sqrt(LbyR_fit)
    y = corr_fit/LbyR_fit**1.5
    theta3,theta2 = np.polyfit(x,y,1)

    return theta2, theta3


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

    left, right = 0, 0.02
    fit_x, fit_y = [], []
    nofit_x, nofit_y = [], []
    for i,LbyR_ in enumerate(LbyR):
        if left <= LbyR_ <= right:
            fit_x.append(LbyR_)
            fit_y.append(correction[i])
        else:
            nofit_x.append(LbyR_)
            nofit_y.append(correction[i])

    theta2, theta3 = fit_theta23(fit_x, fit_y)
    cmd = "y(x)=%.8g*x**1.5 %+.8g*x**2" % (theta2,theta3)
    title = r"$\theta_2 x^{3/2}+\theta_3 x^2$, $\theta_2=%.2f$, $\theta_3=%.2f$" % (theta2,theta3)
    g.plot(graph.data.function(cmd, title=title))

    g.plot(
        graph.data.values(x=fit_x, y=fit_y, title="multipole data (used for fit)"),
        [graph.style.symbol(graph.style.symbol.triangle, size=0.05, symbolattrs=[color.cmyk.Maroon])]
    )

    g.plot(
        graph.data.values(x=nofit_x, y=nofit_y, title="multipole data"),
        [graph.style.symbol(graph.style.symbol.circle, size=0.05, symbolattrs=[color.cmyk.MidnightBlue])]
    )

    # plot correction according to Bimonte et al.
    # Bimonte, Emig, Jaffe, Kardar, Casimir forces beyond the proximity approximation
    g.plot(graph.data.function("y(x)=-4.52*x**2*log(x)", title=r"Bimonte \textit{et al.}"))

    g.writePDFfile()
