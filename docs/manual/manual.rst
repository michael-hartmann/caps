.. CaPS documentation master file, created by
   sphinx-quickstart on Mon Dec 17 15:08:14 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

CaPS user manual
================

.. toctree::
   :maxdepth: 2


Overview and Features
=====================

CaPS provides software to describe the Casimir effect in the plane-sphere
geometry for arbitrary temperatures and arbitrary non-magnetic materials
constituting sphere and plate. Both objects are assumed to be in vacuum. The
plane-sphere geometry is sketched in the inset of :numref:`aspect_ratios` and is
characterized by the sphere radius :math:`R` and the distance :math:`L` between
sphere and infinite plane.

The main goal of the library and the associated programs is to compute the free
energy :math:`\mathcal{F}` as a function of the radius :math:`R` of the sphere,
the separation :math:`L` between sphere and plate, the temperature :math:`T`,
and the dielectric properties of the sphere and the plane. The code is highly
optimized and -- depending on parameters and the available resources -- allows
to compute the free energy for aspect ratios up to :math:`R/L\sim 5\,000`.
An idea of typical aspect ratios used in Casimir experiments in the plane-sphere
geometry can be obtained from :numref:`aspect_ratios`.

.. _aspect_ratios:
.. figure:: images/overview.*
   :scale: 100 %

   The aspect ratio :math:`R/L` used in experiments in the plane-sphere geometry
   is shown by the red stripes. The blue area indicates the aspect ratios that
   are accessible using CaPS. The inset depicts the plane-sphere geometry where
   a sphere of radius :math:`R` is placed at a distance :math:`L` from a plane.

Features
--------

CaPS provides the following main features:

 - Computation of the free energy for aspect ratios used in typical experiments.
 - Full support for perfect reflectors, metals described by the Drude and plasma model, and generic materials described by a user-defined dielectric function.
 - Support for parallelization using MPI.
 - Computation of the free energy in the high-temperature limit for perfect reflectors and metals described by the Drude or plasma model.

The computation of the high-temperature limit for the Drude model is based on
`G. Bimonte, T. Emig, "Exact Results for Classical Casimir Interactions: Dirichlet
and Drude Model in the Sphere-Sphere and Sphere-Plane Geometry", Phys. Rev. Lett. 109,
160403 (2012) <https://doi.org/10.1103/PhysRevLett.109.160403>`_.

Basic support for further geometries is provided for the special case of zero
temperature and perfect reflectors:

 - Computation of the free energy in the plane-cylinder geometry.
 - Computation of the free energy for two spheres with equal radii.

The implementation for the plane-cylinder geometry is based on a symmetrized
version of the matrix elements given in `T. Emig, R. L. Jaffe, M. Kardar, and
A. Scardicchio, "Casimir Interaction between a Plate and a Cylinder", Phys.
Rev. Lett. 96, 080403 (2006) <https://doi.org/10.1103/PhysRevLett.96.080403>`_.


Further reading
---------------

Some of the numerical ideas used in this library are described in `M. Hartmann,
G.-L. Ingold, P. A. Maia Neto, "Advancing numerics for the Casimir effect to experimentally
relevant aspect ratios", Phys. Scr. 93, 114003 (2018)
<https://doi.org/10.1088/1402-4896/aae34e>`_. A more detailed description can
be found in `M. Hartmann, "Casimir effect in the plane-sphere geometry: Beyond the
proximity force approximation", Ph.D. thesis (Universität Augsburg, 2018)
<https://opus.bibliothek.uni-augsburg.de/opus4/44798>`_.


Installation
============

In the following, we assume the operating system to be Ubuntu 18.04. The
commands should also work on other Debian-like systems.

Compilation
-----------

The easiest way to obtain the source code of the CaPS package is by cloning it
from Github. If not already present, install git by running

.. code-block:: none

    $ sudo apt install git

in a terminal. Here, the dollar sign indicates the shell prompt. Once git is
installed, the command

.. code-block:: none

    $ git clone https://github.com/michael-hartmann/caps.git

will download the complete CaPS repository and store it in the directory
``caps/``. As an alternative, you can also download and extract the zip- or
tar.gz-archive of the latest release by navigating to
``https://github.com/michael-hartmann/caps/releases``.

The CaPS library and the programs are written in C and C++ using LAPACK and
MPI. The requirements to compile the source code are

 * a C and C++ compiler,
 * the development files for LAPACK and MPI,
 * the build tools make and cmake.

These dependencies can be installed with:

.. code-block:: none

    $ sudo apt install gcc g++ libc6-dev libc++-dev cmake make libopenmpi-dev openmpi-bin liblapack-dev


In order to compile the code, create a directory ``bin`` in the ``caps/``
directory and run ``cmake`` followed by ``make``:

.. code-block:: none

    $ cd caps/
    $ mkdir bin
    $ cd bin/
    $ cmake ..
    $ make

The last command compiles the `HODLR library
<https://github.com/sivaramambikasaran/HODLR/>`_, the libcaps library, and
builds the shared objects ``libhodlr.so`` and ``libcaps.so``. Then, the
programs ``caps`` and ``caps_logdetD`` are built.

Note that alternative compilers can be specified by setting the variables
``CC`` and ``CXX``. For the Intel C and C++ compilers, the last two commands
above should be replaced by:

.. code-block:: none

    $ CC=icc CXX=icpc cmake ..
    $ make

In order to run the programs, the system must be able to find the libraries
``libhodlr.so`` and ``libcaps.so``. As the corresponding directories are not in
the default search path, they need to be added to ``LD_LIBRARY_PATH``

.. code-block:: none

    $ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/hendrik/caps/bin

where we have assumed that the CaPS repository is in the directory
``/home/hendrik`` [#hendrik]_ .

After the previous steps, the programs ``capc`` (plane-cylinder geometry) and
``cass`` (sphere-sphere geometry with equal radii) to compute the Casimir free
energy for perfect reflectors at :math:`T=0` can be built with:

.. code-block:: none

    $ make capc
    $ make cass

Under Ubuntu 18.10 we encountered problems linking to `OpenBLAS
<https://www.openblas.net/>`_ resulting in error messages similar to:

.. code-block:: none

    undefined reference to 'dgemm_'
    undefined reference to 'dstemr_'
    undefined reference to 'dpotrf_'
    undefined reference to 'dgemm_'
    undefined reference to 'dgetrf_'
    undefined reference to 'dgeqrf_'
    undefined reference to 'ddot_'

In this case, we recommend using `Atlas <http://math-atlas.sourceforge.net/>`_
as a `BLAS <https://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms>`_
implementation. Make sure that Atlas is installed

.. code-block:: none

    $ sudo apt install libatlas-dev libatlas3-base

and compile the code using:

.. code-block:: none

    $ cd caps/
    $ mkdir bin
    $ cd bin/
    $ cmake .. -DBLA_VENDOR=ATLAS
    $ make


Testing
-------

In order to verify that the compilation was successful, build and run
the unit tests in ``bin/``:

.. code-block:: none

    $ make tests
    $ ./caps_tests

All tests should pass. Running the tests takes (depending on your hardware)
about 7 minutes.


Improving performance
---------------------

In order to improve performance, it might be necessary to tweak some compiler
options. By default, the optimization level is ``-O3`` which usually yields
reasonable performance.

If you either run and compile the code on the same machine, or the target
machine supports the same instruction set, the option ``-march=native`` may
increase performance. On an Intel Core i7-2600 machine, a performance boost
of about 5% was found. To compile the code with this option, run:

.. code-block:: none

    $ cmake .. -DOPT="-O3 -march=native"


Similarly, link time optimization ``-flto`` might also increase performance.
However, a test on an Intel Core i7-2600 machine showed basically no
performance gain.


Programs
========

.. _caps:

caps
----

The program ``caps`` computes the Casimir free energy :math:`\mathcal{F}` for
the plane-sphere geometry as a sum

.. math::
  \mathcal{F} = \frac{k_\mathrm{B}T}{2} \sum_{n=-\infty}^\infty \sum_{m=-\infty}^\infty \log\mathrm{det}\left(1-\mathcal{M}^{(m)}(\xi_n)\right)
  :label: matsubara_sum

over the Matsubara frequencies :math:`\xi_n=2\pi n k_\mathrm{B} T /\hbar`. For
zero temperature :math:`T=0`, the sum over the Matsubara frequencies becomes an
integration. :math:`\mathcal{M}^{(m)}` denotes the round-trip operator associated with
the scattering of an electromagnetic wave propagating from the sphere to the
plane and back. Due to the axial symmetry of the plane-sphere geometry, in the
multipole basis the round-trip operator becomes block-diagonal in the
eigenvalues :math:`m` of the :math:`z`-component of the angular momentum.

The program supports a wide variety of options. A summary of all options can be
obtained with ``caps --help``. By default, the temperature is set to
:math:`T=0`, and the sphere and plane are assumed to be perfect reflectors.

Please note that due to parallelization, the number of terms computed in the
summation over :math:`m` may vary from run to run. As a consequence, the
numerical value obtained for the free energy also may vary from run to run
while respecting the prescribed relative error. The abort criterion for the
summation over :math:`m` can be changed by the option ``--cutoff`` described in
more detail below.

Mandatory options
^^^^^^^^^^^^^^^^^

There are two mandatory parameters:

* separation :math:`L` between sphere and plane
* radius :math:`R` of the sphere.

The program expects the lengths to be given in units of meters. As an example,
the following command computes the Casimir interaction at :math:`T=0` for
perfect reflectors for a sphere of radius :math:`R=50\mu\mathrm{m}` and a
separation :math:`L=500\,\mathrm{nm}`:

.. code-block:: none

    $ mpirun -n 8 ./caps -R 50e-6 -L 500e-9

The command ``mpirun`` will set up the environment for MPI and the flag ``-n``
specifies how many processes the program should use. The first process is the
master process. It delegates works to the other minion processes and collects
the results once they are available. Due to this design, ``caps`` needs at
least two processes to be able to start. As the master process does not consume
much CPU time, set ``-n`` to N+1 to fully utilize the computational power of a
computer with N processor cores.

The output of the above command looks similar to:

.. code-block:: none

    $ mpirun -n 8 ./caps -R 50e-6 -L 500e-9
    # version: 0.5
    # compiler: gcc
    # compile time: Nov 18 2019 06:58:55
    # compiled on: Linux host.name 5.0.0-36-generic x86_64
    # git HEAD: a7e8abd
    # git branch: master
    # pid: 8816
    # start time: Mon Nov 18 07:11:26 2019
    #
    # LbyR = 0.009999999999999998
    # RbyL = 100
    # L = 5e-07
    # R = 5e-05
    # T = 0
    # cutoff = 1e-09
    # epsrel = 1e-06
    # iepsrel = 1e-08
    # ldim = 701
    # cores = 8
    # quad = adaptive Gauss-Kronrod
    #
    # xi*(L+R)/c=50.5, logdetD=-9.606867042703382, t=16.748
    # xi*(L+R)/c=11769.79101884239, logdetD=-5.195558844576977e-101, t=2.54115
    # xi*(L+R)/c=0.216677594013123, logdetD=-27.84616165794131, t=14.2098
    # xi*(L+R)/c=1934.091409969965, logdetD=-5.351518755017293e-16, t=5.30664
    # xi*(L+R)/c=1.318577801883526, logdetD=-27.490829731263, t=14.5494
       ...
    # xi*(L+R)/c=1.947625694054103, logdetD=-27.20161097974516, t=14.6707
    # xi*(L+R)/c=4.123327036713686, logdetD=-26.06136838981697, t=14.9208
    # xi*(L+R)/c=2.630682896323236, logdetD=-26.85839019124925, t=14.8417
    #
    # ier=0, integral=-26.66082541075548, neval=165, epsrel=4.22789e-07
    #
    # 13687 determinants computed
    # stop time: Mon Nov 18 07:53:31 2019
    #
    # L/R, L, R, T, ldim, E*(L+R)/(hbar*c)
    0.009999999999999998, 5e-07, 5e-05, 0, 701, -428.5634172542063

The output adopts the CSV format and additional comments start with a number
sign (#). The first comment section contains information on the compilation
like the version of CaPS, time of compilation, name of compiler, machine where
it was compiled and so on. A second comment section gives information about the
geometry (radius :math:`R`, minimal separation :math:`L`, aspect ratio
:math:`R/L`, and inverse aspect ratio :math:`L/R`) as well as numerical
parameters and the employed integration technique (cutoff, epsrel, iepsrel,
ldim, cores, quad). We will discuss the latter in more detail below. The value
of the cores parameter is the number of MPI processes that were used for the
computation.

The following section lists the numerical results for the determinant of the
scattering matrix printed for the different Matsubara frequencies at which it
was evaluated. The comment starting with ``ier`` gives the result of the
integration. If the integration was successful, the value is ``ier=0``, see
also the description of `dqags <http://www.netlib.org/quadpack/dqags.f>`_ of
`QUADPACK <http://www.netlib.org/quadpack/>`_. The program ends by printing
the result of the computation. The free energy is given in units of
:math:`(L+R)/\hbar c`. For the present example the free energy is

.. math::
  \mathcal{F}\approx \frac{-428.6 \hbar c}{50\mu\mathrm{m}+500\mathrm{nm}} \approx -2.68\times10^{-19} \mathrm{J}.

The PFA result in this case is :math:`\mathcal{F}_\mathrm{PFA} = \hbar c \pi^3 R/720 L^2 \approx -2.72\times10^{-19} \mathrm{J}`.

The desired relative accuracy of the integration over the Matsubara frequencies
can be set using ``--epsrel``. By default, ``EPSREL`` is :math:`10^{-6}`. Note
that the integrand needs to be sufficiently smooth. In particular, for very low
values of ``EPSREL`` you might need to decrease the value of ``CUTOFF`` using
``--cutoff``. The value of ``CUTOFF`` determines when the summation over
:math:`m` is stopped. The abort criterion is:

.. math::
    \frac{\log\mathrm{det}\left(1-\mathcal{M}^{(m)}(\xi)\right)}{\log\mathrm{det}\left(1-\mathcal{M}^{(0)}(\xi)\right)} < \mathrm{CUTOFF}

The default value of ``CUTOFF`` is :math:`10^{-9}`. As a rule of thumb, in
order for the integrand to be sufficiently smooth for the integration routine,
``CUTOFF`` should be at least two orders of magnitude smaller than ``EPSREL``.

By default, the integration routine uses an adaptive Gauss-Kronrod method
provided by `CQUADPACK <https://github.com/ESSS/cquadpack>`_. For perfect
reflectors it is sometimes faster to use an adaptive exponentially convergent
Fourier-Chebyshev quadrature scheme (FCQS), see `J. P. Boyd, "Exponentially
convergent Fourier-Chebshev quadrature schemes on bounded and infinite
intervals", J. Sci. Comput. 2, 99 (1987) <https://doi.org/10.1007/BF01061480>`_.
FCQS can be enabled by the flag ``--fcqs``. Since the adaptive algorithm for
FCQS is not well tested, this option is considered experimental. Moreover, it
is not recommended to use FCQS for any other materials than perfect reflectors.

Temperature
^^^^^^^^^^^

The temperature in Kelvin can be set by using the flag ``-T``. The following
call computes the free energy just like in the last example but at room
temperature :math:`T=300\mathrm{K}`:

.. code-block:: none

    $ mpirun -n 8 ./caps -R 50e-6 -L 500e-9 -T 300
    # version: 0.5
    # compiler: gcc
    # compile time: Nov 18 2019 06:58:55
    # compiled on: Linux host.name 5.0.0-36-generic x86_64
    # git HEAD: a7e8abd
    # git branch: master
    # pid: 10689
    # start time: Mon Nov 18 07:54:55 2019
    #
    # LbyR = 0.009999999999999998
    # RbyL = 100
    # L = 5e-07
    # R = 5e-05
    # T = 300
    # using Matsubara spectrum decomposition (MSD)
    # cutoff = 1e-09
    # epsrel = 1e-06
    # iepsrel = 1e-08
    # ldim = 701
    # cores = 8
    # model = perfect reflectors
    #
    # xi*(L+R)/c=0, logdetD=-27.86199079699379, t=0.217287
    # xi*(L+R)/c=41.56988050956833, logdetD=-11.57945164872656, t=16.3904
    # xi*(L+R)/c=83.13976101913666, logdetD=-4.919484480176409, t=17.4001
    # xi*(L+R)/c=124.709641528705, logdetD=-2.131755366198671, t=18.1506
    # xi*(L+R)/c=166.2795220382733, logdetD=-0.9308974790223026, t=18.6908
    # xi*(L+R)/c=207.8494025478416, logdetD=-0.4078030842242257, t=18.8725
    # xi*(L+R)/c=249.41928305741, logdetD=-0.1788875087017671, t=19.0952
    # xi*(L+R)/c=290.9891635669783, logdetD=-0.07851443395562838, t=18.956
    # xi*(L+R)/c=332.5590440765466, logdetD=-0.03446767502060717, t=19.1627
    # xi*(L+R)/c=374.128924586115, logdetD=-0.01513223879283446, t=18.9359
    # xi*(L+R)/c=415.6988050956833, logdetD=-0.006643455732743601, t=18.6485
    # xi*(L+R)/c=457.2686856052516, logdetD=-0.00291655594485459, t=18.444
    # xi*(L+R)/c=498.83856611482, logdetD=-0.001280335351840918, t=18.162
    # xi*(L+R)/c=540.4084466243883, logdetD=-0.0005620156141146819, t=17.8208
    # xi*(L+R)/c=581.9783271339566, logdetD=-0.0002466829395738677, t=17.121
    # xi*(L+R)/c=623.548207643525, logdetD=-0.0001082657488926631, t=17.108
    # xi*(L+R)/c=665.1180881530933, logdetD=-4.751159181908004e-05, t=15.9656
    # xi*(L+R)/c=706.6879686626615, logdetD=-2.084779251432432e-05, t=15.0527
    #
    # 1683 determinants computed
    # stop time: Mon Nov 18 07:59:59 2019
    #
    # L/R, L, R, T, ldim, E*(L+R)/(hbar*c)
    0.009999999999999998, 5e-07, 5e-05, 300, 701, -452.7922091974657

For finite temperatures, the free energy is no longer given as an integral, but
as the sum :eq:`matsubara_sum` over Matsubara frequencies :math:`\xi_n`. The
summation over :math:`n` is stopped once

.. math::

    \frac{\log\mathrm{det}\left( 1-\mathcal{M}(\xi_n) \right)}{\log\mathrm{det}\left( 1-\mathcal{M}(0) \right)} < \mathrm{EPSREL} .

By default, ``EPSREL`` is :math:`10^{-6}`. Its value can be modified by means
of the option ``--epsrel``.

By default, the free energy is computed by means of :eq:`matsubara_sum`
referred to as Matsubara spectrum decomposition (MSD). An alternative approach
is the Padé spectrum decomposition (PSD). PSD is an optimal sum-over-poles
expansion scheme explained in `J. Hu, M. Luo, F. Jiang, R.-X. Xu, Y. J. Yan,
"Padé spectrum decompositions of quantum distribution functions and optimal
hierarchical equations of motion construction for quantum open systems", J.
Chem. Phys. 134, 244106 (2010) <https://doi.org/10.1063/1.3602466>`_. The PSD
requires the evaluation of fewer terms as compared to the MSD. PSD can be
enabled with the flag ``--psd``. The order is determined
automatically to achieve a relative error of the order specified by
``--epsrel``, but can also be set manually using ``--psd-order``. Since
the automatic determination of the order is not well tested, PSD is
considered experimental.

The high-temperature limit of the Casimir free energy requires only the
evaluation of :math:`\log\mathrm{det}(1-\mathcal{M}(0))` and can be obtained by
means of the flag ``--ht``. It will be determined for Drude metals and perfect
reflectors. In contrast to the cases of zero and finite temperatures, the
result is given in units of :math:`k_\mathrm{B} T`. While for perfect
reflectors, no analytical result is known, for Drude metals the analytical
formula given in `G. Bimonte, T. Emig, PRL 109, 160403 (2012)
<https://doi.org/10.1103/PhysRevLett.109.160403>`_ is used.

For example, for a sphere of radius :math:`R=100\mu\mathrm{m}` and a smallest
separation :math:`L=100\mathrm{nm}` one finds:

.. code-block:: none

    $ mpirun -n 8 ./caps -R 100e-6 -L 100e-9 --ht
    # version: 0.5
    # compiler: gcc
    # compile time: Nov 18 2019 06:58:55
    # compiled on: Linux host.name 5.0.0-36-generic x86_64
    # git HEAD: a7e8abd
    # git branch: master
    # pid: 10901
    # start time: Mon Nov 18 08:01:34 2019
    #
    # LbyR = 0.0009999999999999998
    # RbyL = 1000
    # L = 1e-07
    # R = 0.0001
    # high-temperature limit
    # cutoff = 1e-09
    # epsrel = 1e-06
    # iepsrel = 1e-08
    # ldim = 7001
    # cores = 8
    #
    # L/R, L, R, ldim, E_Drude/(kB*T), E_PR/(kB*T)
    0.0009999999999999998, 1e-07, 0.0001, 7001, -149.6981411829862, -296.4343145093178

Compared to zero or finite temperatures, it is considerably less demanding to
compute the high-temperature limit of the Casimir free energy. For the above
example, the aspect ratio is :math:`R/L=1000`, but the computation time on a
standard desktop computer using 8 cores is only about 13 seconds.

Material parameters
^^^^^^^^^^^^^^^^^^^

The examples presented so far were mostly for sphere and plate made of
perfect reflectors. In addition, it is possible to do calculations for the
plasma model with the imaginary-frequency dielectric function

.. math::
  \epsilon(\mathrm{i}\xi) = 1+\frac{\omega_\mathrm{P}^2}{\xi^2},
  :label: plasma_model

for the Drude model with

.. math::
  \epsilon(\mathrm{i}\xi) = 1+\frac{\omega_\mathrm{P}^2}{\xi(\xi+\gamma)},
  :label: drude_model

and for materials with a user-defined dielectric function. Both objects, sphere
and plate, are assumed to consist of the same material.

The plasma model allows to account for high-frequency transparency of the
material and is characterized by the plasma frequency :math:`\omega_\mathrm{P}`
appearing in :eq:`plasma_model`. This parameter can be set using ``--omegap``.
and its value is expected to be given in units of :math:`\mathrm{eV}/\hbar`. For example,
for :math:`R=50\mu\mathrm{m}`, :math:`L=500\mathrm{nm}`,
:math:`T=300\mathrm{K}`, and plasma frequency
:math:`\omega_\mathrm{P}=9\mathrm{eV}/\hbar` appropriate for gold, the Casimir
free energy assuming the plasma model is:

.. code-block:: none

    $ mpirun -n 8 ./caps -R 50e-6 -L 500e-9 -T 300 --omegap 9
    # version: 0.5
    # compiler: gcc
    # compile time: Nov 18 2019 06:58:55
    # compiled on: Linux host.name 5.0.0-36-generic x86_64
    # git HEAD: a7e8abd
    # git branch: master
    # pid: 10988
    # start time: Mon Nov 18 08:03:02 2019
    #
    # LbyR = 0.009999999999999998
    # RbyL = 100
    # L = 5e-07
    # R = 5e-05
    # T = 300
    # using Matsubara spectrum decomposition (MSD)
    # cutoff = 1e-09
    # epsrel = 1e-06
    # iepsrel = 1e-08
    # ldim = 701
    # cores = 8
    # omegap = 9
    # gamma = 0
    # model = plasma
    #
    # xi*(L+R)/c=0, logdetD=-26.69763442281914, t=0.997097
    # xi*(L+R)/c=41.56988050956833, logdetD=-10.46875933145602, t=28.0664
    # xi*(L+R)/c=83.13976101913666, logdetD=-4.171274450117993, t=29.1975
    # xi*(L+R)/c=124.709641528705, logdetD=-1.691108794704161, t=30.2682
    # xi*(L+R)/c=166.2795220382733, logdetD=-0.6897717724345741, t=30.5765
    # xi*(L+R)/c=207.8494025478416, logdetD=-0.2819480734885185, t=31.4948
    # xi*(L+R)/c=249.41928305741, logdetD=-0.1153300107335112, t=31.8762
    # xi*(L+R)/c=290.9891635669783, logdetD=-0.04718515917377312, t=31.8926
    # xi*(L+R)/c=332.5590440765466, logdetD=-0.01930591989009091, t=31.7388
    # xi*(L+R)/c=374.128924586115, logdetD=-0.007899245650576751, t=31.072
    # xi*(L+R)/c=415.6988050956833, logdetD=-0.00323219686668101, t=30.9696
    # xi*(L+R)/c=457.2686856052516, logdetD=-0.001322632063466512, t=30.2518
    # xi*(L+R)/c=498.83856611482, logdetD=-0.0005412792023858141, t=29.5332
    # xi*(L+R)/c=540.4084466243883, logdetD=-0.0002215415508830627, t=28.231
    # xi*(L+R)/c=581.9783271339566, logdetD=-9.068791275557157e-05, t=27.8686
    # xi*(L+R)/c=623.548207643525, logdetD=-3.712883414740327e-05, t=26.6126
    # xi*(L+R)/c=665.1180881530933, logdetD=-1.5203622898573e-05, t=25.0134
    #
    # 1580 determinants computed
    # stop time: Mon Nov 18 08:10:57 2019
    #
    # L/R, L, R, T, ldim, E*(L+R)/(hbar*c)
    0.009999999999999998, 5e-07, 5e-05, 300, 701, -408.1688660043022

The Drude model :eq:`drude_model` not only accounts for the high-frequency
cutoff but also a finite zero-frequency conductivity :math:`\sigma_0 =
\omega_\mathrm{P}^2/\gamma`. The additional parameter :math:`\gamma` can be
specified by the flag ``--gamma`` with a value given in units of
:math:`\mathrm{eV}/\hbar`. For gold, :math:`\gamma=35\mathrm{meV}/\hbar` and
extending the previous example to the Drude model yields:

.. code-block:: none

    $ mpirun -n 8 ./caps -R 50e-6 -L 500e-9 -T 300 --omegap 9 --gamma 0.035
    # version: 0.5
    # compiler: gcc
    # compile time: Nov 18 2019 06:58:55
    # compiled on: Linux host.name 5.0.0-36-generic x86_64
    # git HEAD: a7e8abd
    # git branch: master
    # pid: 12151
    # start time: Mon Nov 18 10:14:33 2019
    #
    # LbyR = 0.009999999999999998
    # RbyL = 100
    # L = 5e-07
    # R = 5e-05
    # T = 300
    # using Matsubara spectrum decomposition (MSD)
    # cutoff = 1e-09
    # epsrel = 1e-06
    # iepsrel = 1e-08
    # ldim = 701
    # cores = 8
    # omegap = 9
    # gamma = 0.035
    # model = drude
    #
    # xi*(L+R)/c=0, logdetD=-14.56972271677286, t=5.91278e-05
    # xi*(L+R)/c=41.56988050956833, logdetD=-10.36538156526362, t=27.6398
    # xi*(L+R)/c=83.13976101913666, logdetD=-4.136280953711821, t=32.2039
    # xi*(L+R)/c=124.709641528705, logdetD=-1.677636571797743, t=32.745
    # xi*(L+R)/c=166.2795220382733, logdetD=-0.6843948068144456, t=32.3712
    # xi*(L+R)/c=207.8494025478416, logdetD=-0.2797738555246083, t=33.3351
    # xi*(L+R)/c=249.41928305741, logdetD=-0.1144462419581209, t=33.8076
    # xi*(L+R)/c=290.9891635669783, logdetD=-0.04682514492128204, t=32.6435
    # xi*(L+R)/c=332.5590440765466, logdetD=-0.01915913156991313, t=32.6817
    # xi*(L+R)/c=374.128924586115, logdetD=-0.007839375398166228, t=32.2618
    # xi*(L+R)/c=415.6988050956833, logdetD=-0.003207775547051069, t=32.0442
    # xi*(L+R)/c=457.2686856052516, logdetD=-0.001312670760241527, t=31.6889
    # xi*(L+R)/c=498.83856611482, logdetD=-0.0005372163538879504, t=30.4185
    # xi*(L+R)/c=540.4084466243883, logdetD=-0.0002198846198504241, t=29.254
    # xi*(L+R)/c=581.9783271339566, logdetD=-9.001224030793225e-05, t=28.6272
    # xi*(L+R)/c=623.548207643525, logdetD=-3.685333004404979e-05, t=26.779
    # xi*(L+R)/c=665.1180881530933, logdetD=-1.50912958984194e-05, t=25.6237
    # xi*(L+R)/c=706.6879686626615, logdetD=-6.180965431945956e-06, t=24.0002
    #
    # 1605 determinants computed
    # stop time: Mon Nov 18 10:23:11 2019
    #
    # L/R, L, R, T, ldim, E*(L+R)/(hbar*c)
    0.009999999999999998, 5e-07, 5e-05, 300, 701, -325.8011897616507

General materials can be defined by the user and specified by the flag
``--material``. Its parameter should be a path to a file describing
the material in the following format:

.. code-block:: none

    # Drude parameters for low frequencies
    # omegap_low = 9.0eV
    # gamma_low  = 0.03eV
    #
    # Drude parameters for high frequencies
    # omegap_high = 54.475279eV
    # gamma_high  = 211.48855eV
    #
    # xi in rad/s            epsilon(i*xi)
    151900000000.0000        27202177.31278406
    167090000000.0000        24722324.84701070
    182280000000.0000        22655566.74077545
    ...
    1.4886200000000000E+018   1.002550948190463
    1.5038100000000000E+018   1.002504731170786
    1.5190000000000000E+018   1.002459773692494

Each line either starts with a number sign (#) or contains a frequency
:math:`\xi` in units of :math:`\mathrm{rad/s}` and the corresponding value of
the dielectric function :math:`\epsilon(\mathrm{i}\xi)` separated by either tabs or spaces.
The frequencies have to be in ascending order. The dielectric function for an
arbitrary frequency is then computed using linear interpolation. For
frequencies smaller than the smallest frequency provided in the file, the
dielectric function is computed using the Drude model :eq:`drude_model`
with the plasma frequency given by ``omegap_low`` and the relaxation frequency
given by ``gamma_low``. If ``omegap_low`` and ``gamma_low`` are not given in
the file, the dielectric function is assumed to be 1. The behavior for
frequencies larger than the largest provided frequency is analogous using the
parameters given by ``omegap_high`` and ``gamma_high``. More details can be
found in the directory ``materials/``.

The following example computes the Casimir energy for a sphere of
:math:`R=50\mu\mathrm{m}` at separation :math:`L=500\mathrm{nm}` at room
temperature :math:`T=300\mathrm{K}` for real gold:

.. code-block:: none

    $ mpirun -n 8 ./caps -R 50e-6 -L 500e-9 -T 300 --material ../materials/gold.csv
    # version: 0.5
    # compiler: gcc
    # compile time: Nov 18 2019 06:58:55
    # compiled on: Linux host.name 5.0.0-36-generic x86_64
    # git HEAD: a7e8abd
    # git branch: master
    # pid: 12528
    # start time: Mon Nov 18 10:29:46 2019
    #
    # LbyR = 0.009999999999999998
    # RbyL = 100
    # L = 5e-07
    # R = 5e-05
    # T = 300
    # using Matsubara spectrum decomposition (MSD)
    # cutoff = 1e-09
    # epsrel = 1e-06
    # iepsrel = 1e-08
    # ldim = 701
    # cores = 8
    # filename = ../materials/gold.csv
    # model = optical data (xi=0: Drude)
    # plasma = -26.69763442234945 (logdetD(xi=0) for plasma model with omegap=9eV)
    #
    # xi*(L+R)/c=0, logdetD=-14.56972271677286, t=0.986975
    # xi*(L+R)/c=41.56988050956833, logdetD=-10.39087965476027, t=28.3849
    # xi*(L+R)/c=83.13976101913666, logdetD=-4.156258938103592, t=29.4192
    # xi*(L+R)/c=124.709641528705, logdetD=-1.691910799997694, t=30.6499
    # xi*(L+R)/c=166.2795220382733, logdetD=-0.6936539818518235, t=31.1129
    # xi*(L+R)/c=207.8494025478416, logdetD=-0.2853780306974359, t=32.1065
    # xi*(L+R)/c=249.41928305741, logdetD=-0.1176725558423893, t=32.0842
    # xi*(L+R)/c=290.9891635669783, logdetD=-0.04860885045324241, t=32.2575
    # xi*(L+R)/c=332.5590440765466, logdetD=-0.02011518916812526, t=32.1692
    # xi*(L+R)/c=374.128924586115, logdetD=-0.008338175617856349, t=31.7953
    # xi*(L+R)/c=415.6988050956833, logdetD=-0.003462411566035189, t=31.8069
    # xi*(L+R)/c=457.2686856052516, logdetD=-0.001440113625772619, t=30.6612
    # xi*(L+R)/c=498.83856611482, logdetD=-0.0006000317807945225, t=30.6105
    # xi*(L+R)/c=540.4084466243883, logdetD=-0.0002504008864329687, t=29.151
    # xi*(L+R)/c=581.9783271339566, logdetD=-0.0001046571681066221, t=28.2663
    # xi*(L+R)/c=623.548207643525, logdetD=-4.380432428803017e-05, t=27.5338
    # xi*(L+R)/c=665.1180881530933, logdetD=-1.836011835885276e-05, t=25.7766
    # xi*(L+R)/c=706.6879686626615, logdetD=-7.705107816276758e-06, t=24.594
    #
    # 1680 determinants computed
    # stop time: Mon Nov 18 10:38:16 2019
    #
    # L/R, L, R, T, ldim, E*(L+R)/(hbar*c)
    0.009999999999999998, 5e-07, 5e-05, 300, 701, -326.8806691538903

As indicated in the comment line referring to the model used, the zeroth
Matsubara frequency is evaluated by means of a Drude model.  In order to obtain
the corresponding result where the plasma model is used instead for the zeroth
Matsubara frequency, the value given in the comment line starting with ``#
plasma =`` can be used. The value given there, i.e. -26.69... in our example,
is given in units of :math:`k_\mathrm{B}T/2` and corresponds to the additional
contribution in the high-temperature limit to the energy in the plasma model.
In our example, the free energy using the Drude model at zero frequency is

.. math::
  \mathcal{F}_\mathrm{Drude} \approx -326.88 \frac{\hbar c}{L+R} \approx -2.0464\times10^{-19}\mathrm{J},

while assuming the plasma model for zero frequency we find

.. math::
  \mathcal{F}_\mathrm{plasma} \approx \mathcal{F}_\mathrm{Drude} - 26.69763 \frac{k_\mathrm{B}T}{2} \approx -2.5993 \times 10^{-19} \mathrm{J} .

Truncation of the vector space
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The required dimension of the vector space

.. math::
   :label: eta

   \ell_\mathrm{dim} = \eta\frac{R}{L}

scales linearly with the aspect ratio :math:`R/L` with a prefactor
:math:`\eta` related to the estimated relative error according to the following table:

+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| relative error  | :math:`10^{-2}` | :math:`10^{-3}` | :math:`10^{-4}` | :math:`10^{-5}` | :math:`10^{-6}` | :math:`10^{-7}` | :math:`10^{-8}` |
+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| :math:`\eta`    | 2.8             | 4               | 5.2             | 6.4             | 7.6             | 8.8             | 10              |
+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+

The truncation of the vector space is discussed in more detail in `M. Hartmann,
"Casimir effect in the plane-sphere geometry: Beyond the proximity force
approximation", PhD thesis (Universität Augsburg, 2018)
<https://opus.bibliothek.uni-augsburg.de/opus4/44798>`_.

The dimension of the vector space can be specified by means of the flag
``--ldim`` fixing :math:`\ell_\mathrm{dim}`. Alternatively, the flag ``--eta`` can
be used from which the dimension is obtained according to

.. math::
    \ell_\mathrm{dim} = \mathrm{max}\left(20, \left\lceil \eta R/L \right\rceil\right)

where :math:`\lceil x\rceil` denotes the smallest integer larger than :math:`x`.
By default, the dimension of the vector space is determined from :eq:`eta` with
:math:`\eta=7`.


Other options
^^^^^^^^^^^^^

The computation of the matrix elements of the round-trip operator involves an
integration. The desired relative error for this integration can be set using
``--iepsrel``. The default value of :math:`10^{-8}` should be sufficient for
almost all purposes. If the Casimir energy needs to be determined to very high
accuracy with a relative error of :math:`10^{-7}` or smaller, it is recommended
to decrease ``IEPSREL`` accordingly.

If ``caps`` was interrupted, e.g. when the time limit on a compute cluster was
exceeded, the ``--resume`` option can be used to resume the computation on the
basis of the partial output created so far. If it is given the option
``--resume FILENAME``, ``caps`` reads the content of ``FILENAME`` and re-uses
the computed values. It is the responsibility of the user to make sure that all
other parameters given to ``caps`` exactly match the parameters used to
generate ``FILENAME`` in a previous run.

caps_logdetD
------------

The program ``caps_logdetD`` computes the building block of :eq:`matsubara_sum`

.. math::
    \log\mathrm{det}\left(1-\mathcal{M}^{(m)}(\xi)\right)

which, in addition to the parameters introduced above, depends on :math:`m` and
:math:`\xi`. Accordingly, there exist two additional mandatory options besides
``-L`` and ``-R``, namely ``-m`` and ``--xi``. The frequency given by ``--xi``
is expected in units of :math:`c/(L+R)`.


The options ``-L``, ``-R``, ``--ldim``, ``--material``, and ``--iepsrel`` are
the same as described in :numref:`caps` for the program ``caps``. In addition,
the algorithm used to compute the determinant can be specified with ``--detalg``.
Valid values are HODLR, QR, LU, and Cholesky.

A typical output looks like

.. code-block:: none

    $ ./caps_logdetD -R 100e-6 -L 1e-6 -m 1 --xi 1
    # ./caps_logdetD, -R, 100e-6, -L, 1e-6, -m, 1, --xi, 1
    # L/R    = 0.009999999999999998
    # L      = 1e-06
    # R      = 0.0001
    # ldim   = 501
    # epsrel = 1.0e-08
    # detalg = HODLR
    #
    # L, R, ξ*(L+R)/c, m, logdet(Id-M), ldim, time
    1e-06, 0.0001, 1, 1, -6.463973151333012, 501, 0.500394

Sometimes, it is useful to dump the round-trip matrix in Numpy format. If the
environment variable ``CAPS_DUMP`` is set and ``detalg`` is not HODLR, the
round-trip matrix will be saved to the filename contained in ``CAPS_DUMP``.
Also note that if ``detalg`` is Cholesky, only the upper half of the matrix is
computed.

The following example demonstrates how to generate and save a round-trip matrix:

.. code-block:: none

    $ CAPS_DUMP=M.npy ./caps_logdetD -R 100e-6 -L 1e-6 -m 1 --xi 2 --detalg LU
    # ./caps_logdetD, -R, 100e-6, -L, 1e-6, -m, 1, --xi, 2, --detalg, LU
    # L/R    = 0.009999999999999998
    # L      = 1e-06
    # R      = 0.0001
    # ldim   = 501
    # epsrel = 1.0e-08
    # detalg = LU
    #
    # L, R, ξ*(L+R)/c, m, logdet(Id-M), ldim, time
    1e-06, 0.0001, 2, 1, -6.349155228127988, 501, 0.758627

The newly created file ``M.npy`` containing the NumPy dump of the round-trip matrix
can be read back into a Python session as follows:

.. code-block:: none

    $ python
    Python 3.7.3 (default, Mar 27 2019, 22:11:17)
    [GCC 7.3.0] :: Anaconda, Inc. on linux
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import numpy as np
    >>> M = np.load("M.npy")    # load matrix
    >>> dim, _ = M.shape        # get dimension
    >>> Id = np.eye(dim)        # identity matrix
    >>> np.linalg.slogdet(Id-M) # compute log of determinant
    (1.0, -6.349155228128016)


capc
----

The program ``capc`` computes the Casimir interaction in the cylinder-plane
geometry for perfect reflectors at zero temperature. The radius of the cylinder
is specified by means of the flag ``-R`` and the smallest separation between
cylinder and plate is specified by ``-d``. Both lengths are expected to be
given in meters. The length :math:`L` of the cylinder is assumed to be large
compared to the cylinder radius :math:`R` and the cylinder-plate distance
:math:`d` so that the cylinder can be assumed to be infinitely long.

The following example computes the Casimir energy per unit length for an
infinitely long cylinder of radius :math:`R=100\mu\mathrm{m}` and a separation
of :math:`d=100\mathrm{nm}`:

.. code-block:: none

    ./capc -R 100e-6 -d 100e-9
    # R/d = 1000
    # d = 1e-07
    # R = 0.0001
    # T = 0
    # lmax = 6000
    # epsrel = 1e-08
    #
    # d/R, d, R, T, lmax, E_PFA/(L*hbar*c), E_D/E_PFA, E_N/E_PFA, E_EM/E_PFA
    0.001, 1e-07, 0.0001, 0, 6000, -72220981652413.5, 0.500089151077938, 0.499432943796738, 0.999522094874677

``E_D`` and ``E_N`` correspond to Dirichlet and Neumann boundary conditions, respectively,
and ``E_EM`` is the energy for the electromagnetic field with ``E_EM = E_D + E_N``.
All results are given as ratios with respect to the free energy per unit length calculated
within the proximity-force approximation.

A full list of options accepted by ``capc`` can be obtained by ``capc --help``.


cass
----

The program ``cass`` computes the Casimir energy in the sphere-sphere geometry
for perfect reflectors at zero temperature. It is assumed that both spheres have the
same radius :math:`R=R_1=R_2`. The radius is specified by means of the flag ``-R`` and
the smallest distance between the two spheres is specified by ``-d``. Both lengths are
expected to be given in units of meters.

The round-trip operator in the sphere-sphere geometry is given as

.. math::

     \mathcal{M}_\mathrm{SS} = \mathcal{R}_1 \mathcal{T}_{12} \mathcal{R}_2 \mathcal{T}_{21} \,.

Here, :math:`\mathcal{R}_j` denotes the reflection operator at sphere
:math:`j`, and :math:`\mathcal{T}_{ij}` is the translation operator from sphere
:math:`j` to sphere :math:`i`. After symmetrization and for equal radii, the
round-trip operator can be written as:

.. math::

     \widehat{\mathcal{M}}_\mathrm{SS} = \underbrace{\sqrt{\mathcal{R}} \mathcal{T} \sqrt{\mathcal{R}}}_{=A} \underbrace{\sqrt{\mathcal{R}} \mathcal{T} \sqrt{\mathcal{R}}}_{=A} = A^2 \,.

Since for perfect reflectors the Fresnel coefficients are
:math:`r_\mathrm{TM}=-r_\mathrm{TE}=1`, the operator :math:`A` can be expressed
as :math:`A=\mathcal{M}_\mathrm{PS}` where :math:`\mathcal{M}_\mathrm{PS}` is
the symmetrized round-trip operator in the plane-sphere geometry. This idea
basically amounts to using the method of image charges.

The following example computes the Casimir energy in the sphere-sphere geometry
for :math:`R_1=R_2=100\mu\mathrm{m}` and :math:`d=10\mu\mathrm{m}`:

.. code-block:: none

    $ ./cass -R 100e-6 -d 10e-6
    # version: 0.5
    # compiler: gcc
    # compile time: Nov 18 2019 06:58:55
    # compiled on: Linux host.name 5.0.0-36-generic x86_64
    # git HEAD: a7e8abd
    # git branch: master
    #
    # sphere-sphere geometry
    # model: perfect reflectors
    # R1 = R2 = 0.0001
    # d = 1e-05
    # T = 0
    # ldim = 50
    # epsrel = 1e-06
    # iepsrel = 1e-09
    # cutoff = 1e-10
    #
    # xi*d/c=1, logdet=-0.303363189004
    # xi*d/c=233.06516869, logdet=-3.13138917509e-204
    # xi*d/c=0.004290645426, logdet=-1.68510617725
    ...
    # xi*d/c=0.0385668454268, logdet=-1.65999004318
    # xi*d/c=0.081650040331, logdet=-1.60322260646
    # xi*d/c=0.0520927306203, logdet=-1.64443634206
    #
    # ier = 0, neval = 135
    #
    # R1, R2, L, E*d/(hbar*c), relative error (due to integration)
    0.0001, 0.0001, 1e-05, -0.166155548334, 5.80089e-07

The runtime of this program is about 4 minutes. For a full list of options see
``cass --help``.

The program ``cass`` carries out a full matrix-matrix multiplication using LAPACK to
compute :math:`\widehat{\mathcal{M}}_\mathrm{SS}` and a LU decomposition to compute
the determinant of the scattering matrix. Note that this program does not support
parallelization. Therefore, the program is useful only for not too large
aspect ratios :math:`R/d`.



API Documentation
=================

The API documentation of the API is available online as
`html <https://www.speicherleck.de/michael/caps/api/index.html>`_ or
`PDF <https://www.speicherleck.de/michael/caps/api.pdf>`_. You can also
generate the API documentation running

.. code-block:: none

  $ doxygen doxygen.conf

in the directory ``src/``. The documentation will be generated in ``docs/api/``.
You need doxygen installed on your computer.


Examples
========

In the directory ``examples/`` examples of simple programs can be found that
demonstrate how to use the CaPS API. The examples are kept simple and well
documented.

.. [#hendrik] In honor of Hendrik Brugt Gerhard Casimir (1909-2000).
