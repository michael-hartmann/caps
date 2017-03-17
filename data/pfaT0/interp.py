from math import log,e
from pyx import *
import numpy as np
from scipy import interpolate


def slurp(filename):
    data = []
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if "# k=" in line:
                line = line[2:]
                line = line.replace("k=", "")
                line = line.replace("logdetD=", "")
                line = line.replace("t=", "")
                line = line.replace("xi=", "")
                data.append(list(map(float, line.split(","))))
    return data


def interp(x,y, eps=1e-7):
    f_linear = interpolate.interp1d(x, y, kind="linear")
    f_cubic  = interpolate.interp1d(x, y, kind="cubic")

    def g(x):
        lin = f_linear(x)
        if abs(lin) < eps:
            return lin

        cub = f_cubic(x)
        if cub > 0:
            return lin
        else:
            return cub

    return np.vectorize(g)


# k, xi, logdetD, t
data = slurp("drude/0014.out")
x = []
y = []
for k,xi,logdetD,t in data:
    x.append(xi)
    y.append(logdetD)

f = interp(x,y)

# use LaTeX for Pyx
text.set(text.LatexRunner)

d_inter = [(0,f(0))]
for x in np.logspace(log(1e-3),log(10),500,base=e):
    d_inter.append((x,f(x)))
for x in np.linspace(log(10),500,200):
    d_inter.append((x,f(x)))

d_inter = list(sorted(d_inter))

g = graph.graphxy(
    width = 8,
    x     = graph.axis.lin(title=r"$\xi (L+R)/\mathrm{c}$", max=500),
    y     = graph.axis.lin(title=r"$\log\det\mathcal{D}(\xi)$")
)

g.plot(
   graph.data.points(data, x=2, y=3),
   [graph.style.symbol(graph.style.symbol.circle, size=0.07)]
)
g.plot(graph.data.points(d_inter, x=1,y=2), [graph.style.line([color.rgb.red])])

g2 = graph.graphxy(
    width = 4, xpos = 3.5, ypos = 1,
    x     = graph.axis.log(min=1e-2,max=25),
    y     = graph.axis.lin()
)
g2.plot(graph.data.points(d_inter, x=1,y=2), [graph.style.line([color.rgb.red])])
g2.plot(
   graph.data.points(data, x=2, y=3),
   [graph.style.symbol(graph.style.symbol.circle, size=0.07)]
)
g.insert(g2)


g.text(0.3,4.5, "$L/R=0.00336$")

g.writePDFfile()
