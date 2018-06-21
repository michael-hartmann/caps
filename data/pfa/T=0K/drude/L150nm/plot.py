import numpy as np
from sys import argv
from pyx import *

text.set(text.LatexRunner)

filename = "data.csv"
if len(argv) > 1:
    filename = argv[1]
data = np.loadtxt(filename)

LbyR = data[:,0]
corr = data[:,2]

y = corr/LbyR**1.5
rho = np.sqrt(LbyR)

fit_x = rho[-2:]
fit_y = y[-2:]
m,t = np.polyfit(fit_x, fit_y, 1)
print(t,m)

xtitle = r"$\rho = \sqrt{L/R}$"
ytitle = r"$\left(E/E_\mathrm{PFA}-1-\theta_1\rho^2\right)/\rho^3$"

g = graph.graphxy(
    width = 8,
    x     = graph.axis.lin(title=xtitle, max=0.3),
    y     = graph.axis.lin(title=ytitle, max=2.5)
)

cmd = "y(x)=%g*x+%g" % (m,t)
g.plot(graph.data.function(cmd, min=0, max=0.04))

g.plot(graph.data.values(x=rho, y=y), [graph.style.symbol(graph.style.symbol.circle, symbolattrs=[color.cmyk.MidnightBlue], size=0.07)])

g.writePDFfile()
