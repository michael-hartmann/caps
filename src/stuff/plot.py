#!/usr/bin/python3

import mmap
from pyx import *
from sys import argv
from math import log

def slurp(filename, modulo=8):
    data = []
    log10 = log(10)
    with open(filename, "r+b") as f:
        m = mmap.mmap(f.fileno(), 0, prot=mmap.PROT_READ)
        while True:
            line = m.readline().decode()
            if line == "":
                break

            i,j,Mij = line.split(",")
            i = int(i)
            if i % modulo == 0:
                j = int(j)
                if j % modulo == 0:
                    data.append((i,j,float(Mij)/log10))

    return data


if __name__ == "__main__":
    filename = argv[1]

    print("Reading...")
    data = slurp(filename, modulo=1)

    print("Plotting...")
    g = graph.graphxy(
        width = 8,
        ratio = 1,
        x     = graph.axis.lin(title="$j$ (column)", min=0, max=800),
        y     = graph.axis.lin(title="$i$ (row)",    min=0, max=800, reverse=True)
    )

    g.plot(graph.data.points(data, x=2, y=1, color=3, title=r"$\log_{10}|M_{ij}|$"),
           [graph.style.density(gradient=color.gradient.Jet)])

    g.writePDFfile(filename)
