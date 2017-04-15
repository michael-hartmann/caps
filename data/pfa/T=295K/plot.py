import numpy as np
from pyx import *

data_R10 = np.loadtxt("R=10e-6/data.csv", delimiter=",")
data_R30 = np.loadtxt("R=30e-6/data.csv", delimiter=",")
data_R50 = np.loadtxt("R=50e-6/data.csv", delimiter=",")
data_R151 = np.loadtxt("R=151.3e-6/data.csv", delimiter=",")

x_R10        = data_R10[:,0]
y_R10_drude  = data_R10[0,1]*(1-data_R10[:,3]/data_R10[:,5])
y_R10_plasma = data_R10[0,1]*(1-data_R10[:,4]/data_R10[:,6])

x_R30        = data_R30[:,0]
y_R30_drude  = data_R30[0,1]*(1-data_R30[:,3]/data_R30[:,5])
y_R30_plasma = data_R30[0,1]*(1-data_R30[:,4]/data_R30[:,6])

x_R50        = data_R50[:,0]
y_R50_drude  = data_R50[0,1]*(1-data_R50[:,3]/data_R50[:,5])
y_R50_plasma = data_R50[0,1]*(1-data_R50[:,4]/data_R50[:,6])

x_R151        = data_R151[:,0]
y_R151_drude  = data_R151[0,1]*(1-data_R151[:,3]/data_R151[:,5])
y_R151_plasma = data_R151[0,1]*(1-data_R151[:,4]/data_R151[:,6])

attrs = [color.gradient.RedBlue]

text.set(text.LatexRunner)

g_drude = graph.graphxy(
    width = 8,
    key = graph.key.key(pos="tr"),
    x = graph.axis.lin(title=r"$d$ (nm)", divisor=1e-9),
    y = graph.axis.lin(title=r"$\frac{R}{d} \left(1-F^\prime/F^\prime_\mathrm{PFA}\right)$")
)

g_drude.plot([
    graph.data.values(x=x_R10,  y=y_R10_drude/x_R10,  title=r"$R=10\mu\mathrm{m}$"),
    graph.data.values(x=x_R30,  y=y_R30_drude/x_R30,  title=r"$R=30\mu\mathrm{m}$"),
    graph.data.values(x=x_R50,  y=y_R50_drude/x_R50,  title=r"$R=50\mu\mathrm{m}$"),
    graph.data.values(x=x_R151, y=y_R151_drude/x_R151, title=r"$R=151.3\mu\mathrm{m}$"),
    ], [graph.style.symbol(graph.style.symbol.changecircle, symbolattrs=attrs, size=0.03)])

g_drude.writePDFfile("drude")


g_plasma = graph.graphxy(
    width = 8,
    key = graph.key.key(pos="tl"),
    x = graph.axis.lin(title=r"$d$ (nm)", divisor=1e-9),
    y = graph.axis.lin(title=r"$\frac{R}{d} \left(1-F^\prime/F^\prime_\mathrm{PFA}\right)$")
)

g_plasma.plot([
    graph.data.values(x=x_R10,  y=y_R10_plasma/x_R10,  title=r"$R=10\mu\mathrm{m}$"),
    graph.data.values(x=x_R30,  y=y_R30_plasma/x_R30,  title=r"$R=30\mu\mathrm{m}$"),
    graph.data.values(x=x_R50,  y=y_R50_plasma/x_R50,  title=r"$R=50\mu\mathrm{m}$"),
    graph.data.values(x=x_R151, y=y_R151_plasma/x_R151, title=r"$R=151.3\mu\mathrm{m}$"),
    ], [graph.style.symbol(graph.style.symbol.changecircle, symbolattrs=attrs, size=0.03)])

g_plasma.writePDFfile("plasma")
