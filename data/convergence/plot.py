from pyx import *
import numpy as np

text.set(text.LatexRunner)

RbyL_10_0    = np.loadtxt("conv_10_0.csv", delimiter=",")
RbyL_200_0   = np.loadtxt("conv_200_0.csv", delimiter=",")
RbyL_200_300 = np.loadtxt("conv_200_300.csv", delimiter=",")

# fit
x = RbyL_200_0[:,0]
y = np.log(RbyL_200_0[:,3])
m,t = np.polyfit(x,y,1)
print("From fit:")
print("Δ=%g*exp(%g*η)" % (np.exp(t), m))
print("log(Δ)=%g%+gη" % (t, m))
print("η=(%g-log(Δ))/%g" % (-t, -m))

g = graph.graphxy(
    width = 8,
    key   = graph.key.key(pos="bl"),
    x     = graph.axis.lin(title=r"$\eta$",min=0.5, max=10.5),
    y     = graph.axis.log(title=r"$\Delta$")
)

attrs = [graph.style.line([style.linestyle.solid, color.gradient.RedBlue]), graph.style.symbol(graph.style.symbol.changecircle, size=0.08, symbolattrs=[deco.filled([color.grey(1)])])]
g.plot([
    graph.data.values(x=RbyL_10_0[:,0], y=RbyL_10_0[:,3], title=r"$R/L=10$, $T=0$"),
    graph.data.values(x=RbyL_200_0[:,0], y=RbyL_200_0[:,3], title=r"$R/L=200$, $T=0$"),
    graph.data.values(x=RbyL_200_300[:,0], y=RbyL_200_300[:,3], title=r"$R/L=200$, $T=300\,\mathrm{K}$")
    ], attrs)

g.writePDFfile()
