The Casimir effect
==================

In 1948 Hendrik Casimir considered two parallel, perfectly conducting plates in
vacuum at temperature T=0 and predicted an attracting force. This force arises
due to vacuum fluctuations and was experimentally verified in 1956 by
Derjaguin, Abrikosova and Lifshitz, as well as in 1958 by Sparnaay.

Since the Casimir effect is a manifestation of vacuum fluctuations in the
mesoscopic world, it has relations to many open physical questions. As the
Casimir force is the dominant force between electically neutral bodies at
micron or sub-micron distances, the Casimir effect plays an important role in
the search for new hypothetical forces predicted by unified theories. The
Casimir effect is also linked to the theory of gravitation and the puzzle of the
cosmological constant. All energy gravitates and thus zero point fluctuations
are expected to contribute to the stress–energy tensor in Einstein's field
equations. In fact, several cosmological observations like the discovery that
the universe expands with an increasing rate suggest a non-zero energy density of
vacuum. However, estimations of the cosmological constant and measurements
disagree by about 120 orders of magnitude.


libcasimir
==========

What is libcasimir?
-------------------
libcasimir implements the numerics for the Casimir effect in the plane-sphere
geometry for arbitrary materials at zero and finite temperature using the
scattering approach. A sphere of radius R is separated by a distance of L from
a plane. The plane is infinite in the xy-direction.

Features
--------
 - Calculate the free energy F(T,L,R) for different separations and temperatures
 - Calculate the free energy F(T→∞,L,R) in the high temperature limit
 - Full support for perfect conductors, Drude metals, and generic metals
   described by a user-defined dielectric function
 - libcasimir is fast and reliable
 - ready to use programs: you don't have to modify the code
 - libcasimir is free software – you may use it or even modify it

Installation
------------
If you use Linux or Unix, you need the gcc and development libraries and header
files for the standard C library, and MPI. On a Debian-like Linux the command
```
$ sudo apt-get install gcc g++ libc6-dev libc++-dev make libopenmpi-dev openmpi-bin liblapack-dev
```
will install all dependencies. You can compile the sources with:
```
$ cd src/
$ make
```
This will build the shared objects `libhodlr.so` and `libcasimir.so`,
and the executables `casimir` and `casimir_logdetD`.

If you want to run the programs, make sure that `libcasimir.so` and
`libhodlr.so` are in the search path or you will get an error similar to:
```
./casimir_logdetD: error while loading shared libraries: libcasimir.so: cannot open shared object file: No such file or directory
```
If the shared libraries are not in the search path, you can still run the
programs by specifying the directories that contain the shared libraries in
`LD_LIBRARY_PATH`:
```
$ LD_LIBRARY_PATH=/path/to/the/libraries ./casimir_logdetD -R 1 -L 0.01 -m 1 --nT 1
# ./casimir_logdetD, -R, 1, -L, 0.01, -m, 1, --nT, 1
# L/R    = 0.01
# ldim   = 500
# epsrel = 1.0e-08
# detalg = HODLR
#
# L, R, ξ*(L+R)/c, m, logdet(Id-M), ldim, time
0.01, 1, 1, 1, -6.46397198784309, 500, 0.524163
```

Besides gcc, the sources may also be compiled with icc or clang (llvm). You can
compile the sources for example with clang with:
```
$ CC=clang make
```

Usage
-----
To compute the Casimir free energy between a sphere of radius R=150µm and a
plane separated by a distance L=2µm at room temperature T=300K, use the
command:
```
mpirun -c 7 ./casimir -L 2e-6 -R 150e-6 -T 300
# L      = 2e-06
# R      = 0.00015
# LbyR   = 0.0133333333333333
# T      = 300
# cutoff = 1e-09
# epsrel = 1e-06
# ldim   = 525
# cores  = 7
# alpha  = 0.0263157894736842
# k=1, xi=0, logdetD=-20.5414371366319, t=0.893615
# k=2, xi=125.121224504, logdetD=-0.697593023726068, t=47.042
# k=3, xi=250.242449008, logdetD=-0.0258391836330834, t=47.5505
# k=4, xi=375.363673512, logdetD=-0.000959951381331692, t=44.3509
# k=5, xi=500.484898016, logdetD=-3.56248208858549e-05, t=39.9074
# k=6, xi=625.60612252, logdetD=-1.31970422831906e-06, t=29.2294
#
# L/R, L, R, T, ldim, F*(L+R)/(ħc)
0.01333333333333333, 2e-06, 0.00015, 300, 525, -1049.959261174529
```

Bugs, developing and contributing
---------------------------------

The latest version of libcasimir is available at
[github](https://github.com/michael-hartmann/libcasimir-dev). If you find a bug, please
create an [issue](https://github.com/michael-hartmann/libcasimir-dev/issues). If you have
improvements, create a pull request.

Authors, license and credits
----------------------------

 * [Michael Hartmann](https://myweb.rz.uni-augsburg.de/~hartmmic/), michael.hartmann@physik.uni-augsburg.de  
   main developer

 * [Gert-Ludwig Ingold](http://www.physik.uni-augsburg.de/theo1/ingold/), gert.ingold@physik.uni-augsburg.de  
   contributions to documentation, speedup of PFA calculation, bugfixes

For a full list, see CREDITS.

The code is licensed under GPLv2, see LICENSE.

Also, libcasimir uses some third-party software:
 * [HODLR](https://github.com/sivaramambikasaran/HODLR) A fast direct solver
   and determinant computation for dense linear systems. (MPL2)
 * [libeigen](http://eigen.tuxfamily.org) HODLR depends on libeigen. However,
   only header files are needed. (MPL2)
 * [cquadpack](https://github.com/ESSS/cquadpack) for integration. cquadpack is
   a a C port of the QUADPACK software originally written in Fortran for
   solving integrals. (public domain)
 * [cephes](http://www.netlib.org/cephes/) is a software collection with
   special functions. libcasimir uses the implementation for the modified
   Bessel functions [I0(x)](http://www.netlib.org/cephes/doubldoc.html#i0) and
   [I1(x)](http://www.netlib.org/cephes/doubldoc.html#i1). Boths files have
   been slightly modified, see besselI.c. (No license, probably BSD licensed.)
 * [LAPACK](http://www.netlib.org/lapack/) Linear algebra library. LAPACK may
   be used to calculate the determinant of the scattering matrices. However,
   for small separations using HODLR is much faster. (Modified BSD)
 * [buf](https://github.com/skeeto/growable-buf) growable memory buffers for
   C99. (public domain)


Publications
------------
 * [Advancing numerics for the Casimir effect to experimentally relevant aspect ratios](https://arxiv.org/abs/1803.05791) (on arxiv)
   Michael Hartmann, Gert-Ludwig Ingold, Paulo A. Maia Neto,
   accepted by Physica Scripta

 * [Plasma vs Drude modelling of the Casimir force: beyond the proximity force approximation](https://doi.org/10.1103/PhysRevLett.119.043901) (on [arxiv](https://arxiv.org/abs/1705.04196))  
   Michael Hartmann, Gert-Ludwig Ingold, and Paulo A. Maia Neto,  
   Phys. Rev. Lett. **119**, 043901 (2017)

 * [Disentangling geometric and dissipative origins of negative Casimir entropies](https://dx.doi.org/10.1103/PhysRevE.92.042125) (on [arxiv](http://arxiv.org/abs/1507.05891))  
   Stefan Umrath, Michael Hartmann, Gert-Ludwig Ingold, and Paulo A. Maia Neto,  
   Phys. Rev. E **92**, 042125 (2015)

 * [Geometric origin of negative Casimir entropies: A scattering-channel analysis](https://dx.doi.org/10.1103/PhysRevE.91.033203) (on [arxiv](http://arxiv.org/abs/1411.1866))  
   Gert-Ludwig Ingold, Stefan Umrath, Michael Hartmann, Romain Guérout, Astrid Lambrecht, Serge Reynaud, and Kimball A. Milton,  
   Phys. Rev. E **91**, 033203 (2015)
