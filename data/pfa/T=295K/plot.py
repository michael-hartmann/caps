import numpy as np
from pyx import *

pfa = np.loadtxt("data_pfa", delimiter=",")
pfa_L        = pfa[:,0]
pfa_P_drude  = pfa[:,1]
pfa_P_plasma = pfa[:,2]


numerics = np.loadtxt("data_numerics", delimiter=",")
numerics_L        = numerics[:,0]
numerics_P_drude  = numerics[:,1]
numerics_P_plasma = numerics[:,2]

attrs = [color.gradient.RedBlue]

g = graph.graphxy(
    width = 12,
    key = graph.key.key(pos="br"),
    x = graph.axis.lin(title="$L$ (nm)", min=270, max=750),
    y = graph.axis.lin(title="$P$ (mPa)", max=0, min=-155)
)

g.plot([
    graph.data.values(x=numerics_L, y=numerics_P_drude, title="Drude"),
    graph.data.values(x=numerics_L, y=numerics_P_plasma, title="Plasma")
    ], [graph.style.symbol(graph.style.symbol.changecircle, symbolattrs=attrs, size=0.1)])

g.plot([
    graph.data.values(x=pfa_L, y=pfa_P_drude, title="Drude (PFA)"),
    graph.data.values(x=pfa_L, y=pfa_P_plasma, title="Plasma (PFA)")
    ], [graph.style.line(attrs)])


g.writePDFfile()
