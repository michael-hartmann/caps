from pyx import *

text.set(text.LatexRunner)

length = 2

def cs(xm=-2,xp=2,ym=-2,yp=2,gap=0.2):
    c = canvas.canvas()

    c.stroke(path.line(xm, 0, xp, 0), [deco.earrow()]) # x-axis
    c.stroke(path.line(0, ym, 0, yp), [deco.earrow()]) # y-axis

    c.text(gap, yp, r"$\mathrm{Im}\,x$")
    c.text(xp-gap, gap, r"$\mathrm{Re}\,x$")

    return c


c1 = cs(ym=-1.3)
c1.text(-2,2, "a)")
c1.fill(path.circle(1.0, 0, 0.05), [color.cmyk.Maroon])
c1.text(1, -0.3, r"$\omega$", [text.halign.boxcenter])

c1.stroke(path.line(-1.7, 0.07, 1.7, 0.07), [deco.earrow(pos=0.3), color.cmyk.OliveGreen])
c1.stroke(path.path(path.arc(0,0.07, 1.69, 0, 180)), [color.cmyk.OliveGreen, deco.earrow(pos=0.3)])

c1.text(-1.5, 0.5, r"$\mathcal{C}_1$", [color.cmyk.OliveGreen])


c2 = cs(ym=-1.3)
c2.text(-2,2, "b)")
c2.fill(path.circle(0, +1, 0.05), [color.cmyk.Maroon])
c2.text(0.1, 1, r"$i\omega$")
c2.fill(path.circle(0, -1, 0.05), [color.cmyk.Maroon])
c2.text(0.1, -1, r"$-i\omega$")

c2.stroke(path.line(-1.7, 0, 1.7, 0), [deco.earrow(pos=0.3), color.cmyk.OliveGreen])
c2.stroke(path.path(path.arc(0,0, 1.69, 0, 180)), [color.cmyk.OliveGreen, deco.earrow(pos=0.3)])

c2.text(-1.5, 0.5, r"$\mathcal{C}_2$", [color.cmyk.OliveGreen])

c1.insert(c2, [trafo.translate(5,0)])

c1.writePDFfile()
