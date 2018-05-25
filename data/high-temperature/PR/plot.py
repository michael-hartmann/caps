import numpy as np
from pyx import *

text.set(text.LatexRunner)

zeta3 = 1.202056903159594

def pfa_drude(RbyL):
    return -zeta3*RbyL/8

data  = np.loadtxt("pr.csv", delimiter=",")
RbyL  = data[:,0]
LbyR  = 1/RbyL
drude = data[:,1]
perf  = data[:,2]
TE    = perf-drude-pfa_drude(RbyL)-np.log(LbyR)**2/16

x = []
y = []
for a,b in zip(LbyR,TE):
    if 2e-4 < a < 5e-4:
        x.append(a)
        y.append(b)

m,t = np.polyfit(np.log(np.array(x)),np.array(y),1)
print(m,t)

g = graph.graphxy(
    width = 8,
    x     = graph.axis.log(title=r"$L/R$"),
    y     = graph.axis.lin(title=r"foo")
)

g.plot(graph.data.values(x=LbyR, y=TE), [graph.style.symbol(graph.style.symbol.changecircle, size=0.05)])
g.plot(graph.data.values(x=x, y=y), [graph.style.symbol(graph.style.symbol.changecircle, symbolattrs=[color.gradient.RedBlue], size=0.05)])
cmd = "y(x)=%g*log(x)+%g" % (m,t)
print(cmd)
g.plot(graph.data.function(cmd))
g.writePDFfile()
