#!/usr/bin/python

from __future__ import division
from scipy.special import zeta
from sys import exit
from math import *
from pyx import *
from glob import glob
import numpy as np
from fit import *

def emig_perf(LbyR):
    zeta3 = 1.202056903159594
    gamma0 = 0.174897
    mu  = log1p(LbyR+sqrt(LbyR*(2+LbyR)))
    Ep  = -zeta3/(2*mu**2) + log(mu)/12 + 1/8 - gamma0 - 7/2880*mu**2 - 31/725760*mu**4
    return Ep/2


def emig_drude(LbyR):
    gamma1 = 1.270362
    gamma2 = 1.35369
    gamma3 = 1.59409
    gamma4 = 2.51153

    mu  = log1p(LbyR+sqrt(LbyR*(2+LbyR)))
    Ep = emig_perf(LbyR)
    return Ep + (log(gamma1-log(mu)) + 1/6*(-gamma2+log(mu))/(-gamma1+log(mu))*mu**2 - 1/180*(gamma3-gamma4*log(mu)+log(mu)**2)/((-gamma1+log(mu))**2)*mu**4)/2


output = ""

fitbereich = [-5, -3] # punkte mit -5 <= ln(x) <= -3 fuer fit verwenden
fit_x   = []
fit_lnx = []
fit_y_p = []
fit_y_d = []

CD = zeta(3,1)/8
CP = zeta(3,1)/4

text.set(mode="latex")

plot = []
special = []
filenames = glob("slurm-*.out")
# die dateinamen sortieren, damit die datenpunkte geordnet sind
filenames.sort()

for filename in filenames:
    fh = open(filename, "r")
    data = ""
    for line in fh:
        line = line.strip()
        if line == "" or line[0] == "#":
            continue
        else:
            x,FP,FD,time = map(float, line.split(","))
            logdetD_perf  = FP*2
            logdetD_drude = FD*2

            rhoP = -4*x*logdetD_perf/zeta(3,1)
            rhoD = -8*x*logdetD_drude/zeta(3,1)

            FP *= 2 # casimi_hiT spuckt 0.5*log det D(n=0) aus
            FD *= 2
            logx  = log(x)

            FP_emig = emig_perf(x)
            FD_emig = emig_drude(x)

            #print x, rhoP, rhoD

            PhiP      = -FP
            PhiD      = -FD
            PhiP_emig = -FP_emig
            PhiD_emig = -FD_emig

            rhoP_emig  = x*PhiP_emig/CP
            rhoD_emig  = x*PhiD_emig/CD

            betaP      = (rhoP-1)/x
            betaD      = (rhoD-1)/x
            betaP_emig = (rhoP_emig-1)/x
            betaD_emig = (rhoD_emig-1)/x

            #            1, 2,    3,  4,  5,       6,       7,    8,    9,         10,        11,   12,   13,        14,        15,    16,    17,         18
            plot.append((x, logx, -FP, -FD, -FP_emig, -FD_emig, PhiP, PhiD, PhiP_emig, PhiD_emig, rhoP, rhoD, rhoP_emig, rhoD_emig, betaP, betaD, betaP_emig, betaD_emig, log10(x), logdetD_perf/logdetD_drude))
            #print x,FD

            if logx >= fitbereich[0] and logx <= fitbereich[1]:
                fit_x.append(x)
                fit_lnx.append(logx)
                fit_y_p.append(betaP)
                fit_y_d.append(betaD)
    fh.close()


attrs = [color.gradient.RedBlue]

print len(fit_x)
print len(fit_lnx)

# 1) rho vs x
g_inset = graph.graphxy(
    width = 4,
    #key   = graph.key.key(pos="tr"),
    xpos = 1.5,
    ypos = 1.3,
    x     = graph.axis.lin(min=0,max=1),
    y     = graph.axis.lin(min=0.2)
)
g_inset.plot(
    [graph.data.points(plot, x=1, y=11, title="perfect"), graph.data.points(plot, x=1, y=12, title="Drude")],
    [
        graph.style.line([color.gradient.RedBlue])
    ]
)
#g.writePDFfile("hiT_1_rho_vs_x")

# 2) rho vs x
g = graph.graphxy(
    width = 10,
    key   = graph.key.key(pos="tr"),
    #x     = graph.axis.log(title=r"$x=L/R$", min=4e-4,max=0.1),
    #y     = graph.axis.lin(title=r"$\rho(x)$",min=0.7,max=1.02)
    x     = graph.axis.log(title=r"$x=L/R$", min=4.5e-4, max=0.25),
    y     = graph.axis.lin(title=r"$\varrho(x)$",max=1.03,min=0.6)
)
g.plot(
    [graph.data.points(plot, x=1, y=11, title="perfect"), graph.data.points(plot, x=1, y=12, title="Drude")],
    [
        graph.style.line([color.gradient.RedBlue]),
        graph.style.symbol(graph.style.symbol.changecircle, symbolattrs=attrs, size=0.08)
    ]
)
g.insert(g_inset)
g.writePDFfile("rho_vs_lnx")



# 3) beta_p and beta_d vs lnx
g = graph.graphxy(
    width = 8,
    key   = graph.key.key(pos="br"),
    x     = graph.axis.lin(title=r"$\mathrm{ln}(x=L/R)$", min=-7.6, max=2),
    y2    = graph.axis.lin(title=r"$\beta(x)$", min=-16, max=0)
)
g.plot(
    [graph.data.points(plot, x=2, y=15, title="perfect"), graph.data.points(plot, x=2, y=16, title="Drude")],
    [
        graph.style.line(attrs)
    ]
)
g.writePDFfile("hiT_3_beta_vs_lnx")


# do the fits
# beta_l
c = np.polyfit(fit_lnx, fit_y_p, 3)
beta_l_p = "y(x) = %f*x**3 %+f*x**2 %+f*x %+f" % (c[0], c[1], c[2], c[3])
output += "perfect: a0=%f, a1=%f, a2=%f, a3=%f\n" % (c[3], c[2], c[1], c[0])

c = np.polyfit(fit_lnx, fit_y_d, 3)
beta_l_d = "y(x) = %f*x**3 %+f*x**2 %+f*x %+f" % (c[0], c[1], c[2], c[3])
output += "drude:   a0=%f, a1=%f, a2=%f, a3=%f\n" % (c[3], c[2], c[1], c[0])

# beta_x
c = np.polyfit(fit_x, fit_y_p, 3)
beta_x_p = "y(x) = %f*exp(x)**3 %+f*exp(x)**2 %+f*exp(x) %+f" % (c[0], c[1], c[2], c[3])
output += "perfect: b0=%f, b1=%f, b2=%f, b3=%f\n" % (c[3], c[2], c[1], c[0])

c = np.polyfit(fit_x, fit_y_d, 3)
beta_x_d = "y(x) = %f*exp(x)**3 %+f*exp(x)**2 %+f*exp(x) %+f" % (c[0], c[1], c[2], c[3])
output += "drude:   b0=%f, b1=%f, b2=%f, b3=%f\n" % (c[3], c[2], c[1], c[0])

# beta m
# startpunkte fuer c0,c1,c2 und c3 finden idee: wir nehmen die ersten vier
# punkte in fit_y_X und bestimmen die koefficienzen c0,c1,c2 und c3, so dass
# die fit-funktion durch diese punkte laeuft. Das benutzen wir dann als
# initialwerte fuer den fit.
M = np.zeros((4,4))
yp = np.zeros(4)
yd = np.zeros(4)
for i in range(4):
    x       = fit_x[i]
    lnx     = fit_lnx[i]
    yp[i]   = fit_y_p[i]
    yd[i]   = fit_y_d[i]
    M[i][0] = 1
    M[i][1] = lnx
    M[i][2] = x
    M[i][3] = x**2

cp = np.linalg.solve(M,yp)
cd = np.linalg.solve(M,yd)

c0 = Parameter(cp[0])
c1 = Parameter(cp[1])
c2 = Parameter(cp[2])
c3 = Parameter(cp[3])

def betam(x):
    return c0() + c1()*x + c2()*exp(x) + c3()*exp(x)**2

fit(betam, [c0, c1, c2, c3], np.array(fit_y_p), x=np.array(fit_lnx))
beta_m_p = "y(x) = %f %+f*x %+f*exp(x) %+f*exp(x)**2" % (c0(), c1(), c2(), c3())
output += "perfect: c0=%f, c1=%f, c2=%f, c3=%f\n" % (c0(), c1(), c2(), c3())

c0 = Parameter(cd[0])
c1 = Parameter(cd[1])
c2 = Parameter(cd[2])
c3 = Parameter(cd[3])

fit(betam, [c0, c1, c2, c3], np.array(fit_y_d), x=np.array(fit_lnx))
beta_m_d = "y(x) = %f %+f*x %+f*exp(x) %+f*exp(x)**2" % (c0(), c1(), c2(), c3())
output += "drude:   c0=%f, c1=%f, c2=%f, c3=%f\n" % (c0(), c1(), c2(), c3())

output += "\n"
output += "(x entspricht ln(x) im plot und exp(x) entspricht x\n"
output += "beta_l,perfect: %s\n" % beta_l_p
output += "beta_l,drude:   %s\n" % beta_l_d
output += "beta_x,perfect: %s\n" % beta_x_p
output += "beta_x,drude:   %s\n" % beta_x_d
output += "beta_m,perfect: %s\n" % beta_m_p
output += "beta_m,drude:   %s\n" % beta_m_d

# 4) beta vs lnx for perfect with fits
g = graph.graphxy(
    width = 8,
    key   = graph.key.key(pos="br"),
    x     = graph.axis.lin(title=r"$\mathrm{ln}(x=L/R)$", min=-7.6, max=-2),
    y     = graph.axis.lin(title=r"$\beta(x)$", min=-16, max=0)
)
g.plot(
    graph.data.points(plot, x=2, y=15, title="perfect"),
    [graph.style.symbol(graph.style.symbol.circle, symbolattrs=[color.gradient.RedBlue], size=0.07)]
)
g.plot([
    graph.data.function(beta_l_p, title=r"$\beta_l$"),
    graph.data.function(beta_x_p, title=r"$\beta_x$"),
    graph.data.function(beta_m_p, title=r"$\beta_m$")
    ])
g.writePDFfile("hiT_4_beta_vs_lnx_perfect_fits")


## 5) beta vs lnx for perfect with emig bimonte
#g = graph.graphxy(
#    width = 10,
#    key   = graph.key.key(pos="br"),
#    x     = graph.axis.lin(title=r"$\mathrm{ln}(x=L/R)$", min=-7.6, max=0),
#    y     = graph.axis.lin(title=r"$\beta(x)$", min=-16, max=0)
#)
#g.plot(
#    graph.data.points(plot, x=2, y=15, title="perfect"),
#    [graph.style.symbol(graph.style.symbol.circle, symbolattrs=[color.gradient.RedBlue], size=0.07)]
#)
#g.plot(
#    graph.data.points(plot, x=2, y=17, title="Emig, Bimonte"),
#    styles=[graph.style.line([color.cmyk.Blue])]
#)
#g.writePDFfile("hiT_5_beta_vs_lnx_perfect_bimonte")

# 6) beta vs lnx for Drude with fits
g = graph.graphxy(
    width = 8,
    key   = graph.key.key(pos="br"),
    x     = graph.axis.lin(title=r"$\mathrm{ln}(x=L/R)$", min=-7.6, max=-2),
    y2    = graph.axis.lin(title=r"$\beta(x)$", min=-4.1, max=-1.8)
)
g.plot(
    graph.data.points(plot, x=2, y=16, title="Drude"),
    [graph.style.symbol(graph.style.symbol.circle, symbolattrs=[color.gradient.RedBlue], size=0.07)]
)
g.plot([
    graph.data.function(beta_l_d, title=r"$\beta_l$"),
    graph.data.function(beta_x_d, title=r"$\beta_x$"),
    graph.data.function(beta_m_d, title=r"$\beta_m$")
    ])
g.writePDFfile("hiT_6_beta_vs_lnx_drude_fit")

# 7) beta vs lnx for Drude with emig,bimonte
g = graph.graphxy(
    width = 8,
    key   = graph.key.key(pos="br"),
    x     = graph.axis.lin(title=r"$\mathrm{ln}(x=L/R)$", min=-7.6, max=2),
    y     = graph.axis.lin(title=r"$\beta(x)$", min=-4.1, max=0)
)
g.plot(
    graph.data.points(plot, x=2, y=16, title="Drude"),
    [graph.style.symbol(graph.style.symbol.circle, symbolattrs=[color.gradient.RedBlue], size=0.07)]
)
g.plot(
    graph.data.points(plot, x=2, y=18, title="Bimonte, Emig"),
    styles=[graph.style.line()]
)
g.writePDFfile("hiT_7_beta_vs_lnx_drude_bimonte")


# inset)
g_inset = graph.graphxy(
    width = 3.5,
    xpos = 0.4,
    ypos = 1.5,
    x     = graph.axis.log(title=r"$x=L/R$", max=10, min=4.5e-4),
    y2    = graph.axis.lin(title=r"$\Phi^\mathrm{P}/\Phi^\mathrm{D}$", max=2.04, min=1.46, density=0.5)
)
g_inset.plot(
    graph.data.points(plot, x=1, y=20),
    [
        graph.style.line([color.gradient.RedBlue, attr.changelist([style.linestyle.solid, style.linestyle.dashed, style.linestyle.dashdotted])])
    ]
)

# 6)
g = graph.graphxy(
    width = 8,
    key   = graph.key.key(pos="bl"),
    x     = graph.axis.log(title=r"$x=L/R$", max=10, min=4.5e-4),
    y     = graph.axis.log(title=r"$-\Phi(x)$", max=1e3)
)
g.plot(
    [graph.data.points(plot, x=1, y=4, title="Drude"),
     graph.data.points(plot, x=1, y=3, title="perfect"),
     graph.data.points(plot, x=1, y=6, title=r"Bimonte, Emig")],
    [
        graph.style.line([color.gradient.RedBlue, attr.changelist([style.linestyle.solid, style.linestyle.dashed, style.linestyle.dashdotted])])
    ]
)
#g.insert(g_inset)
g.writePDFfile("hiT_perf_drude_emig_bimonte")

print output,
