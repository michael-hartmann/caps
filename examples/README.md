examples
========

This directory contains examples that show how to use libcasimir. Every example
is kept simple and well documented. You can compile the programs by running
`make`.


AuAl.c
------

Compute the Casimir free energy at T=300K for a sphere of radius R=50µm and a
separation between sphere and plate of L=1µm. The sphere is consists of gold
(Drude model) and the plate consists of aluminium (Drude model). The program
computes the free energy in units of k_B T.

Output (runtime: ~5mins):
```
$ ./AuAl 
# We assume Drude metals; sphere is gold, plate is aluminium
#
# L/R    = 0.02
# L      = 1e-06
# R      = 5e-05
# ldim   = 320
# epsrel = 1.0e-08
# detalg = HODLR
free energy: -12.584229 kb*T
```
