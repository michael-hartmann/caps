from numpy import linspace
from functools import partial
from math import log10, cos, sin, radians, sqrt
from pyx import *

def crgb(r,g,b, v=255):
    return color.rgb(r/v, g/v, b/v)

def geometry():
    d = 0.2         # Dicke der Platte
    R = 1           # Radius der Kugel
    b = 2*R         # Breite der Platte
    LbyR = 0.5      # aspect ratio R/L
    phi = 30        # Winkel für Pfeil für R
    dx = 0.1        # Abstand für Pfeile für L
    lend = 0.2      # Länge Strich für L
    cfill1  = crgb(250,200,10) # Farbe
    cfill2  = crgb(178,178,255) # Farbe
    cfill2  = crgb(178,178,255) # Farbe
    cborder = crgb(0,0,178)

    c = canvas.canvas()

    # sphere
    c.fill(path.circle(0, 0, R), [cfill1])

    # arrow
    c.stroke(path.line(0, 0, R*cos(radians(phi)), R*sin(radians(phi))), [deco.earrow])

    # label R
    x,y = 0.2*R*cos(radians(phi)), 0.2*R*sin(radians(phi))+0.2
    c.text(x,y, r"$R$")

    # plate
    c.fill(path.rect(-b/2, -R*(1+LbyR), b, -d), [cfill1])

    # label L, arrow and lines
    c.stroke(path.line(b/2+dx, -R,   b/2+dx+lend, -R))
    c.stroke(path.line(b/2+dx, -R*(1+LbyR), b/2+dx+lend, -R*(1+LbyR)))
    c.stroke(path.line(b/2+dx+lend/2, -R, b/2+dx+lend/2, -R*(1+LbyR)), [deco.earrow, deco.barrow])
    c.text(b/2+dx+lend+0.1, -R*(1+LbyR/2), r"$L$", [text.valign.middle])

    c.stroke(path.line(-3*dx, -R, b/2+dx, -R), [style.linestyle.dashed])

    return c

def scaling(x, xmin, xmax, length):
    return (x-xmin)/(xmax-xmin)*length

experiments = [
    {"R": (20.7e-2,), "Lmin": 0.48e-6, "Lmax": 6.5e-6,
     "label": r"Masuda, Sasaki (2009)"},

    {"R": (15.6e-2,), "Lmin": 0.7e-6, "Lmax": 6e-6,
     "label": r"Sushkov \textit{et al.}~(2011)"},
    
    {"R": (11.3e-2,), "Lmin": 0.6e-6,  "Lmax": 6e-6, "label": r"Lamoreaux (1997)"},
    
    {"R": (4e-3,), "Lmin": 0.1e-6,  "Lmax": 2e-6,
     "label": r"Garcia-Sanchez \textit{et al.}~(2012)" },
    
    {"R": (296e-6,), "Lmin": 0.2e-6,  "Lmax": 2e-6,
     "label": r"Decca \textit{et al.}~(2003)"},
    
    {"R": (151.3e-6,), "Lmin": 160e-9, "Lmax": 750e-9,
     "label": r"Decca \textit{et al.}~(2007)"},
    
    {"R": (98e-6,), "Lmin": 0.1e-6, "Lmax": 0.9e-6, "label": r"Mohideen, Roy (1998)" },
    
    {"R": (100e-6,), "Lmin": 50e-9, "Lmax": 195e-9,
     "label": r"Man \textit{et al.}~(2009)"},
    
    {"R": (100e-6,), "Lmin": 0.09e-6, "Lmax": 2.1e-6,
     "label": r"Chan \textit{et al.}~(2001)"},
    
    {"R": (50e-6,), "Lmin": 12e-9, "Lmax": 200e-9,
     "label": r"van Zwol \textit{et al.}~(2008)"},
    
    {"R": (41.3e-6,), "Lmin": 235e-9, "Lmax": 420e-9,
     "label": r"Chang \textit{et al.}~(2012)"},
    
    {"R": (20e-6,), "Lmin": 100e-9, "Lmax": 600e-9,
     "label": r"Jourdan \textit{et al.}~(2009)"},
    
    {"R": (10.5e-6, 31.4e-6, 52.3e-6, 102.8e-6, 148.2e-6), "Lmin": 164e-9, "Lmax": 986e-9,
     "label": r"Krause \textit{et al.}~(2007)"},
    
    {"R": (10e-6,), "Lmin": 60e-9, "Lmax": 300e-9,
     "label": r"Torricelli \textit{et al.}~(2011)"},
    
    {"R": (61.71e-6,), "Lmin": 223e-9, "Lmax": 1000e-9,
     "label": r"Banishev \textit{et al.}~(2013)"},
    
    {"R": (150e-6,), "Lmin": 200e-9, "Lmax": 700e-9,
     "label": r"Bimonte \textit{et al.}~(2016)"},
    
    {"R": (50e-6,), "Lmin": 20e-9, "Lmax": 120e-9,
     "label": r"Zwol \textit{et al.}~(2008)"},
    
    {"R": (20e-6,), "Lmin": 20e-9, "Lmax": 100e-9,
     "label": r"Munday \textit{et al.}~(2008)"},
    
    {"R": (20e-6,), "Lmin": 20e-9, "Lmax": 100e-9,
     "label": r"Munday \textit{et al.}~(2009)"},
]
experiments.sort(key=lambda x: -x["R"][-1]/x["Lmin"])


ratio = 1.333
width = 8
height = width/ratio
axislens = (width, height)
epsinvmin = 1
epsinvmax = 6
ticklen = 0.1

text.set(text.LatexRunner)
c = canvas.canvas()
nrscaling = partial(scaling, xmin=0, xmax=len(experiments), length=axislens[0])
epsscaling = partial(scaling, xmin=epsinvmin, xmax=epsinvmax, length=axislens[1])

xmax = axislens[0]
ymin = epsscaling(log10(10))
ymax = epsscaling(log10(5000))
c.fill(path.rect(0, ymin, xmax, ymax-ymin), [color.rgb(0.66, 0.82, 1)])

c.stroke(path.rect(0, 0, *axislens))

for y_ in (5e3,):
    y = epsscaling(log10(y_))
    c.stroke(path.line(0, y, width, y), [style.linestyle.dashed])

for ny in range(epsinvmin, epsinvmax+1):
    y = epsscaling(ny)
    c.stroke(path.line(0, y, -ticklen, y))
    c.text(-1.5*ticklen, y, '$10^{%i}$' % ny, [text.halign.right, text.valign.middle])
    if ny < epsinvmax:
        for n in range(1, 10):
            y = epsscaling(ny+log10(n))
            c.stroke(path.line(0, y, -0.5*ticklen, y))
c.text(-0.8, axislens[1]/2, r'\Large $\frac{R}{L}$', [text.halign.right, text.valign.middle])
c.text(axislens[0]/2, -0.6, 'experiment', [text.halign.center, text.valign.top])
for nrexp, expdata in enumerate(experiments):
    radii = expdata["R"]
    lmin = expdata["Lmin"]
    lmax = expdata["Lmax"]
    ref = nrexp+1 #expdata["ref"]
    label = expdata["label"]
    print(ref, label)
    nrexp_sc = nrscaling(nrexp+0.5)

    aspect_min = []
    aspect_max = []
    for nr_r, r in enumerate(radii):
        aspect_min.append(epsscaling(log10(r/lmax)))
        aspect_max.append(epsscaling(log10(r/lmin)))

    aspectratiomin_sc = min(aspect_min)
    aspectratiomax_sc = max(aspect_max)

    dy = 0.12
    c.stroke(path.rect(nrexp_sc-dy, aspectratiomin_sc, 2*dy, aspectratiomax_sc-aspectratiomin_sc), [deco.filled([color.cmyk.Maroon])])

    c.stroke(path.line(nrexp_sc, 0, nrexp_sc, -ticklen))
    c.text(nrexp_sc, -1.7*ticklen, ref, [text.size.footnotesize, text.halign.center, text.valign.top])

for nrexp, expdata in enumerate(experiments, start=1):
    label = expdata["label"]

    dy = 0.37
    c.text(width+0.8, 0.3+height-nrexp*dy, "[%d]" % nrexp, [text.size.footnotesize, text.halign.right])
    c.text(width+0.9, 0.3+height-nrexp*dy, label, [text.size.footnotesize])

c.insert(geometry(), [trafo.translate(6,4.97)])

c.writePDFfile()
