Optical data
============

eps.py
------

The script `eps.py` allows you to compute the dielectric function ε(iξ) at the
imaginary axis from optical data. `python eps.py --help` will print a usage.
The following example computes the dielectric function at the imaginary axis
from Palik's optical data:
```
$ python eps.py gold/palik.csv --omegap_low 9 --gamma_low 0.035 --description "Gold (Palik)"
```
Here, we have used a Drude extrapolation of the optical data for small
frequencies with plasma frequency omegap=9eV/hbar and relaxation frequency
gamma=35meV/hbar. The string given by `--description` will be used for a
comment in the output.

gold/
-----

This directory contains the optical data of Gold (real part and imaginary part
of the dielectric function for various frequencies) from different sources as
csv files.

The script `plot.py` plots the imaginary part of the dielectric function for
the different optical data sets. Similarly, `plot_pyx.py` does the same but
uses PyX to generate a more beautiful output.
