import numpy as np
from pyx import *

def fit(RbyL,t,min=200):
    x,y = [],[]
    for xi,ti in zip(RbyL,t):
        if xi >= min:
            x.append(np.log(xi))
            y.append(np.log(ti))

    return np.polyfit(x,y,1)
    

if __name__ == "__main__":
	text.set(text.LatexRunner)

	eta = 5 # ldim = eta*R/L
	xmin, xmax = 50, 2000
	ymin, ymax = 0.2, 1000

	data = np.loadtxt("detalg5.csv", delimiter=",", usecols=(0,3,4))

	RbyL = data[:,0] # R/L
	timing_cholesky = data[:,1] # average runtime for Cholesky (s)
	timing_hodlr    = data[:,2] # average runtime for HODLR (s)

	mytrianglesymbol = lambda c, x_pt, y_pt, size_pt, attrs: graph.style._trianglesymbol(c, x_pt, y_pt, 0.8*size_pt, attrs)
	triangle = [graph.style.symbol(mytrianglesymbol, size=0.09, symbolattrs=[deco.filled([color.grey(0)])])]
	square   = [graph.style.symbol(graph.style._squaresymbol, size=0.09, symbolattrs=[deco.filled([color.grey(1)])])]

	g = graph.graphxy(
		width = 8,
		key   = graph.key.key(pos="tl"),
		x     = graph.axis.log(title=r"$R/L$",min=xmin, max=xmax),
		x2    = graph.axis.log(title=r"$\ell_\mathrm{dim}$", min=xmin*eta, max=xmax*eta, texter=graph.axis.texter.decimal()),
		y     = graph.axis.log(title=r"average time (s)", min=ymin, max=ymax)
	)

	# cholesky
	a,b = fit(RbyL, timing_cholesky)
	g.plot(graph.data.function("y(x)=%g*x**%g" % (np.exp(b),a), title=None, min=200))
	x,y = g.pos(240,40)
	cmd = r"$\propto \left(R/L\right)^{%.2f}$" % a
	g.text(x,y, cmd)

	# hodlr
	a,b = fit(RbyL, timing_hodlr)
	g.plot(graph.data.function("y(x)=%g*x**%g" % (np.exp(b),a), title=None, min=200), [graph.style.line([style.linestyle.solid])])
	x,y = g.pos(700,3.3)
	cmd = r"$\propto \left(R/L\right)^{%.2f}$" % a
	g.text(x,y, cmd)

	# hodlr and cholesky
	g.plot(graph.data.values(x=RbyL, y=timing_cholesky, title=r"Cholesky"), square)
	g.plot(graph.data.values(x=RbyL, y=timing_hodlr,    title=r"HODLR"),    triangle)

	g.writePDFfile()
