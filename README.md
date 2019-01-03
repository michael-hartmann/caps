libcasimir
==========

What is libcasimir?
-------------------

libcasimir is an implementation of the Casimir effect in the plane-sphere
geometry. The geometry consists of a sphere of radius R separated by a
distance L from an infinite plate.

With libcasimir you can compute the free Casimir energy in the plane-sphere
geometry for arbitrary materials at zero and finite temperature. The library is
highly optimized and allows you - depending on parameters and your hardware - to
compute the free energy for aspect ratios of R/L~10'000.

Features
--------
 - Calculate the free energy for different separations and temperatures
 - Calculate the free energy in the high temperature limit
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
$ sudo apt-get install gcc g++ libc6-dev libc++-dev make libopenmpi-dev openmpi-bin liblapack-dev libgfortran-7-dev gfortran-7
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
If the shared libraries are not in the search path, you can run the programs
by specifying the directories that contain the shared libraries in
`LD_LIBRARY_PATH`:
```
$ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/libcasimir-dev/src/:/path/to/libcasimir-dev/src/libhodlr
```

Usage
-----
To compute the Casimir free energy between a sphere of radius R=150µm and a
plane separated by a distance L=2µm at room temperature T=300K assuming
that both objects are perfect reflectors, use the command:
```
$ mpirun -c 7 ./casimir -L 2e-6 -R 150e-6 -T 300
# compiler: gcc
# compile time: Dec 17 2018 13:39:58
# compiled on: Linux jonas.physik.uni-augsburg.de 4.9.0-8-amd64 #1 SMP Debian 4.9.130-2 (2018-10-27) x86_64 GNU/Linux
# git HEAD: f62659b0a419a153c517c54ad1657f7d5df4e2fe
# git branch: master
# pid: 7649
# start time: Mon Dec 17 14:30:58 2018
#
# LbyR = 0.0133333333333333
# RbyL = 75
# L = 2e-06
# R = 0.00015
# T = 300
# using Matsubara spectrum decomposition (MSD)
# cutoff = 1e-09
# epsrel = 1e-06
# iepsrel = 1e-08
# ldim = 525
# cores = 7
# model = perfect reflectors
#
# xi*(L+R)/c=0, logdetD=-20.5414434056605, t=0.124541
# xi*(L+R)/c=125.121224504047, logdetD=-0.697593024441683, t=10.5891
# xi*(L+R)/c=250.242449008094, logdetD=-0.0258391838200906, t=10.7015
# xi*(L+R)/c=375.363673512142, logdetD=-0.000959951407663542, t=10.0898
# xi*(L+R)/c=500.484898016189, logdetD=-3.56248238853923e-05, t=8.87825
# xi*(L+R)/c=625.606122520236, logdetD=-1.31970448193816e-06, t=7.02438
#
# 486 determinants computed
# stop time: Mon Dec 17 14:31:45 2018
#
# L/R, L, R, T, ldim, E*(L+R)/(hbar*c)
0.01333333333333333, 2e-06, 0.00015, 300, 525, -437.9074196681792
```

Documentation
-------------

A user manual is available in `manual/`, the API is documented using doxygen.


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

 * [Gert-Ludwig Ingold](https://www.physik.uni-augsburg.de/theo1/ingold/), gert.ingold@physik.uni-augsburg.de
   contributions to documentation, bugfixes

For a full list, see CREDITS.

The code is licensed under GPLv2, see LICENSE.

Also, libcasimir uses some third-party software:
 * [HODLR](https://github.com/sivaramambikasaran/HODLR) A fast direct solver
   and determinant computation for dense linear systems. (MPL2)
 * [libeigen](https://eigen.tuxfamily.org) Eigen is a C++ template library for
   linear algebra (MPL2)
 * [cquadpack](https://github.com/ESSS/cquadpack) cquadpack is
   a a C port of the QUADPACK software originally written in Fortran for
   solving integrals. (public domain)
 * [cephes](https://www.netlib.org/cephes/) is a software collection with
   special functions. libcasimir uses the implementation for the modified
   Bessel functions [I0(x)](https://www.netlib.org/cephes/doubldoc.html#i0) and
   [I1(x)](https://www.netlib.org/cephes/doubldoc.html#i1). Boths files have
   been slightly modified, see besselI.c. (No license, probably BSD licensed.)
 * [LAPACK](https://www.netlib.org/lapack/) Linear algebra library. LAPACK may
   be used to calculate the determinant of the scattering matrices. However,
   for small separations using HODLR is much faster. (Modified BSD)
 * [buf](https://github.com/skeeto/growable-buf) growable memory buffers for
   C99. (public domain)
 * [argparse](https://github.com/cofyc/argparse) command line arguments parsing
   library in C (MIT)


Publications
------------
 * [Casimir effect in the plane-sphere geometry: Beyond the proximity force approximation](https://opus.bibliothek.uni-augsburg.de/opus4/44798)  
   Michael Hartmann, PhD thesis (Universität Augsburg, 2018)

 * [Advancing numerics for the Casimir effect to experimentally relevant aspect ratios](https://doi.org/10.1088/1402-4896/aae34e) (on [arxiv](https://arxiv.org/abs/1803.05791))  
   Michael Hartmann, Gert-Ludwig Ingold, Paulo A. Maia Neto,  
   Phys. Scr. 93, 114003 (2018), DOI: 10.1088/1402-4896/aae34e

 * [Plasma vs Drude modelling of the Casimir force: beyond the proximity force approximation](https://doi.org/10.1103/PhysRevLett.119.043901) (on [arxiv](https://arxiv.org/abs/1705.04196))  
   Michael Hartmann, Gert-Ludwig Ingold, and Paulo A. Maia Neto,  
   Phys. Rev. Lett. **119**, 043901 (2017), DOI: 10.1103/PhysRevLett.119.04390

 * [Disentangling geometric and dissipative origins of negative Casimir entropies](https://doi.org/10.1103/PhysRevE.92.042125) (on [arxiv](https://arxiv.org/abs/1507.05891))  
   Stefan Umrath, Michael Hartmann, Gert-Ludwig Ingold, and Paulo A. Maia Neto,  
   Phys. Rev. E **92**, 042125 (2015), DOI: 10.1103/PhysRevE.92.042125

 * [Geometric origin of negative Casimir entropies: A scattering-channel analysis](https://doi.org/10.1103/PhysRevE.91.033203) (on [arxiv](https://arxiv.org/abs/1411.1866))  
   Gert-Ludwig Ingold, Stefan Umrath, Michael Hartmann, Romain Guérout, Astrid Lambrecht, Serge Reynaud, and Kimball A. Milton,  
   Phys. Rev. E **91**, 033203 (2015). DOI: 10.1103/PhysRevE.91.033203
