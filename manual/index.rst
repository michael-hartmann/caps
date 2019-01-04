.. libcasimir documentation master file, created by
   sphinx-quickstart on Mon Dec 17 15:08:14 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to libcasimir's documentation!
======================================

.. toctree::
   :maxdepth: 2


Overview and Features
=====================

libcasimir is an implementation of the Casimir effect in the plane-sphere
geometry. The geometry consists of a sphere of radius :math:`R` separated by a
distance :math:`L` from an infinite plate. The main goal of the library and the
associated programs is to compute the free energy :math:`\mathcal{F}` depending
on the radius of the sphere, the separation between sphere and plate, the
temperature, and the dielectric properties of the objects. The library is
highly optimized and allows you - depending on parameters and your hardware -
to compute the free energy for aspect ratios up to :math:`R/L\sim 10\,000`.

.. figure:: geometry.png
   :scale: 33 %

   plane-sphere geometry

Features
--------

 - Calculate the free energy for different separations and temperatures
 - Calculate the free energy in the high temperature limit
 - Full support for perfect conductors, Drude metals, and generic metals described by a user-defined dielectric function
 - libcasimir is fast and reliable
 - ready to use programs: you don't have to modify the code
 - libcasimir is free software – you may use it or even modify it


Further reading
---------------

Some of the numerical ideas used in this library are described in `Hartmann,
Ingold, Maia Neto, "Advancing numerics for the Casimir effect to experimentally
relevant aspect ratios", Phys. Scr. 93, 114003 (2018)
<https://doi.org/10.1088/1402-4896/aae34e>`_. A more detailed description can
be found in `Hartmann, "Casimir effect in the plane-sphere geometry: Beyond the
proximity force approximation", phd thesis (2018)
<https://opus.bibliothek.uni-augsburg.de/opus4/44798>`_.


Installation
============

In the following, we assume the operating system to be Ubuntu 18.04. The
commands should also work on other Debian-like systems.

Compilation
-----------

The easiest way to get the source code is to use git. To install git, run

.. code-block:: console

    $ sudo apt install git

in a terminal. Once git is installed, the command

.. code-block:: console

    $ git clone https://github.com/michael-hartmann/libcasimir-dev.git

will get you the complete libcasimir repository and stores it in the directory
``libcasimir-dev/``.

The libcasimir library and the programs are written in C and C++ using LAPACK
and MPI. In order to compile the source files, you need a C and C++ compiler,
the development files for LAPACK and MPI, and the build tool make. You can install
all dependencies with:

.. code-block:: console

	$ sudo apt install gcc g++ libc6-dev libc++-dev make libopenmpi-dev openmpi-bin liblapack-dev libgfortran-7-dev gfortran-7


In order to compile the code, run ``make`` in the directory ``libcasimir-dev/src/``:

.. code-block:: console

    $ cd libcasimir-src/dev
    $ make

This command compiles the HODLR library and builds the shared object
``libhodlr.so`` in the directory ``libhodlr/``. Then, the Makefile compiles the
libcasimir library ``libcasimir.so`` and builds the programs ``casimir`` and
``casimir_logdetD``.

In order to run the programs, the system must be able to find the libraries
``libhodlr.so`` and ``libcasimir.so``. As both libraries are not in the default
search path, you have to add the directories to ``LD_LIBRARY_PATH``

.. code-block:: console

    $ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/hendrik/libcasimir-dev/src:/home/hendrik/libcasimir-dev/src/libhodlr

where we have assumed that the libcasimir repository is in the directory
``/home/hendrik``.

It is also possible to use different compilers than gcc. To compile the code with clang and clang++ run

.. code-block:: console

    $ make clean
    $ CC=clang CXX=clang++ make

and for Intel's compiler:

.. code-block:: console

    $ make clean
    $ CC=icc CXX=icpc make

Please make sure to run ``make clean`` first.

Testing
-------

In order to check whether the compilation was successful, you can build and run
the unit tests in ``src/``:

.. code-block:: console

    $ make casimir_tests

All tests should pass. Running the tests takes up (depending on your hardware)
about 7 minutes.


Adapting the Makefile
---------------------

In order to improve performance, it might be neccessary to tweak some compiler
options. The options are described in the file ``src/Makefile``. The most
interesting option is ``-mmarch=native`` which tells the compiler to optimize
the code for the architecture the compiler is running on. This might improve
performance by ~5%.


Programs
========

casimir
-------

The program ``casimir`` computes the free Casimir energy :math:`\mathcal{F}`
for the plane-sphere geometry as a sum

.. math::
  \mathcal{F} = \frac{k_\mathrm{B}T}{2} \sum_{n=-\infty}^\infty \sum_{m=-\infty}^\infty \log\mathrm{det}\left(1-\mathcal{M}^m(\xi_n)\right)

over the Matsubara frequencies :math:`\xi_n=2\pi n k_\mathrm{B} T /\hbar`. For
zero temperature :math:`T=0`, the sum over the Matsubara frequencies becomes an
integration.

The program supports a wide variety of options. You can get a summary of all
options using ``casimir --help``. By default, the temperature is set to
:math:`T=0`, and the sphere and plane are assumed to be perfect reflectors.

Mandatory options
^^^^^^^^^^^^^^^^^

There are two mandatory options: the separation :math:`L` between sphere and
plane, and the radius of the sphere :math:`R`. The program expects the lengths
given in units of meters. As an example, the following command computes the
Casimir interaction at :math:`T=0` for perfect reflectors for a sphere of
radius :math:`R=50\mu\mathrm{m}` and a separation :math:`L=2\mu\mathrm{m}`:

.. code-block:: console

	$ mpirun -n 8 ./casimir -R 50e-6 -L 2e-6

The command ``mpirun`` will set up the environment for MPI and the flag ``-n``
specifies how many processes the program should use. If you want to utilize the
full capacity of your machine and your processor has N cores, set ``-n`` to
N+1. The output of this command looks similar to:

.. code-block:: console

    # compiler: gcc
    # compile time: Jan  2 2019 12:35:13
    # compiled on: Linux jonas.physik.uni-augsburg.de 4.9.0-8-amd64 #1 SMP Debian 4.9.130-2 (2018-10-27) x86_64 GNU/Linux
    # git HEAD: 51b8df47618d15747e59321aee77e1f86fab32a8
    # git branch: master
    # pid: 23336
    # start time: Thu Jan  3 12:27:38 2019
    #
    # LbyR = 0.04
    # RbyL = 25
    # L = 2e-06
    # R = 5e-05
    # T = 0
    # cutoff = 1e-09
    # epsrel = 1e-06
    # iepsrel = 1e-08
    # ldim = 176
    # cores = 8
    # quad = adaptive Gauss-Kronrod
    #
    # xi*(L+R)/c=13, logdetD=-2.35165297953047, t=1.14694
    # xi*(L+R)/c=3029.84719296933, logdetD=0, t=9.53674e-07
    	...
    # xi*(L+R)/c=61.4928297600332, logdetD=-0.0548151501755969, t=1.37253
    # xi*(L+R)/c=51.8422518708047, logdetD=-0.115164831511918, t=1.38072
    #
    # ier=0, integral=-6.41428305795366, neval=135, epsrel=2.91062e-07
    #
    # 6154 determinants computed
    # stop time: Thu Jan  3 12:30:01 2019
    #
    # L/R, L, R, T, ldim, E*(L+R)/(hbar*c)
    0.03999999999999999, 2e-06, 5e-05, 0, 176, -26.54248623166202

The output is in the format of a CSV file and additional comments start with a
pound (#). The program first outputs some information on the compilation, i.e.,
time of compilation, name of compiler, machine where it was compiled and so on.
Then, information about the geometry (radius, seperation, aspect ratio and
inverse of aspect ratio), numerical parameters (cutoff, epsrel, iepsrel, ldim,
cores) are printed. We will discuss the numerical parameters in more detail
later. The value of cores is the number of MPI processes that are used for the
computation. Then, the determinant of the scattering matrix for different
Matsubara frequencies are printed. The comment starting with ``ier`` gives the
result of the integration and is 0 if the integration was successful. The
program ends by printing the result of the computation. The free energy is
outputed in units of :math:`(L+R)/\hbar c`, i.e., in this case, the free energy
is

.. math::
  \mathcal{F}\approx \frac{-26.54 \hbar c}{50\mu\mathrm{m}+2\mu\mathrm{m}} \approx -1.61\times10^{-20} \mathrm{J}.

The PFA result in this case is :math:`\mathcal{F}_\mathrm{PFA}\approx-1.64\times10^{-20} \mathrm{J}`.

The desired relative accuracy of the integration over the Matsubara frequencies
can be set using ``--epsrel``. By default, ``EPSREL`` is :math:`10^{-6}`.  Note
that the integrand needs to be sufficiently smooth. In particular, for very low
values of ``EPSREL`` you might need to decrease the value of ``CUTOFF`` using
``--cutoff``. The value of ``CUTOFF`` determines when the summation over
:math:`m` is stopped:

.. math::
    \frac{\log\mathrm{det}\left(1-\mathcal{M}^m(\xi)\right)}{\log\mathrm{det}\left(1-\mathcal{M}^0(\xi)\right)} < \mathrm{CUTOFF}

The default value of ``CUTOFF`` is :math:`10^{-9}`. As a rule of thumb, in
order that the integrand is sufficiently smooth for the integration routine,
``CUTOFF`` should be at least two orders of magnitude smaller than ``EPSREL``.

By default, the integration routine uses an adaptive Gauss-Kronrod method
provided by CQUADPACK. For perfect reflectors it is sometimes faster to use an
adaptive exponentially convergent Fourier-Chebshev quadrature scheme (FCQS),
see `Boyd, "Exponentially convergent Fourier-Chebshev quadrature schemes on
bounded and infinite intervals", JOSC 2, 2 (1987)
<https://doi.org/10.1007/BF01061480>`_. You can use FCQS using the flag
``--fcqs``. Since the adaptive algorithm is not well tested, this option is
considered experimental.

Temperature
^^^^^^^^^^^

You can set the temperature using ``-T``. The following program computes the
free energy just like in the last example but at room temperature
:math:`T=300\mathrm{K}`:

.. code-block:: console

    $ mpirun -n 8 ./casimir -R 50e-6 -L 2e-6 -T 300
    # compiler: gcc
    # compile time: Jan  2 2019 12:35:13
    # compiled on: Linux jonas.physik.uni-augsburg.de 4.9.0-8-amd64 #1 SMP Debian 4.9.130-2 (2018-10-27) x86_64 GNU/Linux
    # git HEAD: 51b8df47618d15747e59321aee77e1f86fab32a8
    # git branch: master
    # pid: 24111
    # start time: Thu Jan  3 12:51:15 2019
    #
    # LbyR = 0.04
    # RbyL = 25
    # L = 2e-06
    # R = 5e-05
    # T = 300
    # using Matsubara spectrum decomposition (MSD)
    # cutoff = 1e-09
    # epsrel = 1e-06
    # iepsrel = 1e-08
    # ldim = 176
    # cores = 8
    # model = perfect reflectors
    #
    # xi*(L+R)/c=0, logdetD=-6.16165739556559, t=0.025461
    # xi*(L+R)/c=42.8046294355951, logdetD=-0.230987113168236, t=1.22125
    # xi*(L+R)/c=85.6092588711902, logdetD=-0.00858232692932167, t=1.186
    # xi*(L+R)/c=128.413888306785, logdetD=-0.000319205031747869, t=1.10549
    # xi*(L+R)/c=171.21851774238, logdetD=-1.18538336130793e-05, t=0.949226
    # xi*(L+R)/c=214.023147177976, logdetD=-4.39335203111237e-07, t=0.77597
    #
    # 302 determinants computed
    # stop time: Thu Jan  3 12:51:21 2019
    #
    # L/R, L, R, T, ldim, E*(L+R)/(hbar*c)
    0.03999999999999999, 2e-06, 5e-05, 300, 176, -45.24539531432269

For finite temperature the free energy is no longer given as an integral, but
as a sum over the Matsubara frequencies :math:`\xi_n`.  The summation is
stopped once

.. math::

	\frac{\log\mathrm{det}\left( 1-\mathcal{M}(\xi_n) \right)}{\log\mathrm{det}\left( 1-\mathcal{M}(0) \right)} < \mathrm{EPSREL} .

By default, ``EPSREL`` is :math:`10^{-6}`. You can change the value of
``EPSREL`` using the option ``--epsrel``.

By default, the free energy is computed summing over the Matsubara frequencies
:math:`\xi_n`, also called Matsubara spectrum decomposition (MSD). Another
option is to compute the free energy using Padé spectrum decomposition (PSD).
PSD is an optimal sum-over-poles expansion scheme, more information can be
found in `Hu, Xu, Yan, "Padé spectrum decompositions of quantum distribution
functions and optimal hierarchical equations of motion construction for quantum
open systems", J. Chem. Phys. 133, 101106 (2010)
<https://doi.org/10.1063/1.3602466>`_. The PSD requires less terms to be
computed compared to the MSD. You can tell the program to use PSD with the flag
``--psd``. The order is determined automatically to achieve a relative error of
the order specified by ``--epsrel``. You can also manually set the order using
``--psd-order``. Since the automatic determination of the order is not well
tested, this option is considered experimental.

If you are only interested in the high-temperature limit, the flag ``--ht``
will only compute :math:`\log\mathrm{det}(1-\mathcal{M}(0))` and output the
Casimir energy in the limit :math:`T\to\infty` in units of :math:`k_\mathrm{B}
T`.


Material parameters
^^^^^^^^^^^^^^^^^^^

Up to this point, we have assumed that the sphere and the plate are perfect
reflectors.  If you want to model the sphere and the plate using the plasma
model, you can set the plasma frequency using ``--omegap``. The plasma
frequency is expected in units of :math:`\mathrm{eV}/\hbar`. For example, for
:math:`R=50\mu\mathrm{m}`, :math:`L=800\mathrm{nm}`, :math:`T=300\mathrm{K}`,
and gold (plasma frequency :math:`\omega_\mathrm{P}=9\mathrm{eV}/\hbar`), the
Casimir free energy using the plasma model is:

.. code-block:: console

    mpirun -n 8 ./casimir -R 50e-6 -L 1e-6 -T 300 --omegap 9
    # compiler: gcc
    # compile time: Jan  2 2019 12:35:13
    # compiled on: Linux jonas.physik.uni-augsburg.de 4.9.0-8-amd64 #1 SMP Debian 4.9.130-2 (2018-10-27) x86_64 GNU/Linux
    # git HEAD: 51b8df47618d15747e59321aee77e1f86fab32a8
    # git branch: master
    # pid: 25590
    # start time: Thu Jan  3 13:45:32 2019
    #
    # LbyR = 0.02
    # RbyL = 50
    # L = 1e-06
    # R = 5e-05
    # T = 300
    # using Matsubara spectrum decomposition (MSD)
    # cutoff = 1e-09
    # epsrel = 1e-06
    # iepsrel = 1e-08
    # ldim = 351
    # cores = 8
    # omegap = 9
    # gamma = 0
    # model = plasma
    #
    # xi*(L+R)/c=0, logdetD=-12.9879117945873, t=0.386934
    # xi*(L+R)/c=41.9814634849106, logdetD=-2.25212799211675, t=8.56138
    # xi*(L+R)/c=83.9629269698212, logdetD=-0.399130249157291, t=8.79863
    # xi*(L+R)/c=125.944390454732, logdetD=-0.0715772034739604, t=8.86131
    # xi*(L+R)/c=167.925853939642, logdetD=-0.0128539503490175, t=8.90524
    # xi*(L+R)/c=209.907317424553, logdetD=-0.00230798551800234, t=8.56306
    # xi*(L+R)/c=251.888780909464, logdetD=-0.000414281644513791, t=8.18385
    # xi*(L+R)/c=293.870244394374, logdetD=-7.43406394982715e-05, t=7.62512
    # xi*(L+R)/c=335.851707879285, logdetD=-1.33358593258147e-05, t=6.89742
    # xi*(L+R)/c=377.833171364195, logdetD=-2.39148466183249e-06, t=6.48929
    #
    # 685 determinants computed
    # stop time: Thu Jan  3 13:46:45 2019
    #
    # L/R, L, R, T, ldim, E*(L+R)/(hbar*c)
    0.02, 1e-06, 5e-05, 300, 351, -123.3743917511159

To describe the objects using the Drude model, you can additional specify the
relaxation frequency :math:`\gamma` (also in units of
:math:`\mathrm{eV}/\hbar`). For gold, :math:`\gamma=35\mathrm{meV}/\hbar`, so
the same example as above for Drude gives:

.. code-block:: console

    mpirun -n 8 ./casimir -R 50e-6 -L 1e-6 -T 300 --omegap 9 --gamma 0.035
    # compiler: gcc
    # compile time: Jan  2 2019 12:35:13
    # compiled on: Linux jonas.physik.uni-augsburg.de 4.9.0-8-amd64 #1 SMP Debian 4.9.130-2 (2018-10-27) x86_64 GNU/Linux
    # git HEAD: 51b8df47618d15747e59321aee77e1f86fab32a8
    # git branch: master
    # pid: 25643
    # start time: Thu Jan  3 13:48:50 2019
    #
    # LbyR = 0.02
    # RbyL = 50
    # L = 1e-06
    # R = 5e-05
    # T = 300
    # using Matsubara spectrum decomposition (MSD)
    # cutoff = 1e-09
    # epsrel = 1e-06
    # iepsrel = 1e-08
    # ldim = 351
    # cores = 8
    # omegap = 9
    # gamma = 0.035
    # model = drude
    #
    # xi*(L+R)/c=0, logdetD=-7.09741176750412, t=0.000463963
    # xi*(L+R)/c=41.9814634849106, logdetD=-2.23320258322643, t=8.40539
    # xi*(L+R)/c=83.9629269698212, logdetD=-0.396003476948632, t=8.85682
    # xi*(L+R)/c=125.944390454732, logdetD=-0.071024911146045, t=9.09324
    # xi*(L+R)/c=167.925853939642, logdetD=-0.0127552647899639, t=8.74376
    # xi*(L+R)/c=209.907317424553, logdetD=-0.00229031081195299, t=8.5974
    # xi*(L+R)/c=251.888780909464, logdetD=-0.000411115126231591, t=8.13146
    # xi*(L+R)/c=293.870244394374, logdetD=-7.3773492717613e-05, t=7.64662
    # xi*(L+R)/c=335.851707879285, logdetD=-1.32343284458642e-05, t=7.20287
    # xi*(L+R)/c=377.833171364195, logdetD=-2.37331945823411e-06, t=6.58019
    #
    # 627 determinants computed
    # stop time: Thu Jan  3 13:50:03 2019
    #
    # L/R, L, R, T, ldim, E*(L+R)/(hbar*c)
    0.02, 1e-06, 5e-05, 300, 351, -83.71300491448063

The Casimir energy in the high-temperature limit for the Drude and the plasma
model differ by a factor of 2. This is the reason why in this example the
Casimir energy using the plasma model is considerably larger than using the
Drude model.

General materials can be described using ``--material`` which expects the path
to a material file. A material file has the following format:

.. code-block:: console

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

Each line either starts with a pound (#) or contains a frequency :math:`\xi` in
units of :math:`\mathrm{rad/s}` and the corresponding value of the dielectric
function :math:`\epsilon(\mathrm{i}\xi)` separated by spaces. The frequencies
have to be in ascending order. The dielectric function for an arbitrary
frequency is then computed using linear interpolation. For frequencies smaller
than the smallest frequency provided in the file, the dielectric function is
computed using the Drude model

.. math::
	\epsilon(\mathrm{i}\xi) = 1+\frac{\omega_\mathrm{P}^2}{\xi(\xi+\gamma)}

with the plasma frequency given by ``omegap_low`` and the relaxation frequency
given by ``gamma_low``. If ``omegap_low`` and ``gamma_low`` are not given in
the file, the dielectric function is assumed to be 1. The behaviour for
frequencies larger than the largest provided frequency is analougous.

Here is an example that computes the Casimir energy for a sphere of
:math:`R=50\mu\mathrm{m}` at separation :math:`L=1\mu\mathrm{m}` at room
temperature :math:`T=300\mathrm{K}` for real gold:

.. code-block:: console

  $ mpirun -n 8 ./casimir -R 50e-6 -L 1e-6 -T 300 --material ../data/materials/GoldDalvit.dat 
  # compiler: gcc
  # compile time: Jan  2 2019 12:35:13
  # compiled on: Linux jonas.physik.uni-augsburg.de 4.9.0-8-amd64 #1 SMP Debian 4.9.130-2 (2018-10-27) x86_64 GNU/Linux
  # git HEAD: 51b8df47618d15747e59321aee77e1f86fab32a8
  # git branch: master
  # pid: 26164
  # start time: Thu Jan  3 14:08:10 2019
  #
  # LbyR = 0.02
  # RbyL = 50
  # L = 1e-06
  # R = 5e-05
  # T = 300
  # using Matsubara spectrum decomposition (MSD)
  # cutoff = 1e-09
  # epsrel = 1e-06
  # iepsrel = 1e-08
  # ldim = 351
  # cores = 8
  # filename = ../data/materials/GoldDalvit.dat
  # model = optical data (xi=0: Drude)
  # plasma = -12.9879117939843 (logdetD(xi=0) for plasma model with omegap=9eV)
  #
  # xi*(L+R)/c=0, logdetD=-7.09741176750412, t=0.347423
  # xi*(L+R)/c=41.9814634849106, logdetD=-2.23671565687706, t=8.46615
  # xi*(L+R)/c=83.9629269698212, logdetD=-0.397198209143505, t=8.72867
  # xi*(L+R)/c=125.944390454732, logdetD=-0.0714429961178252, t=8.97971
  # xi*(L+R)/c=167.925853939642, logdetD=-0.0128860192197612, t=8.94925
  # xi*(L+R)/c=209.907317424553, logdetD=-0.00232743121335635, t=8.76096
  # xi*(L+R)/c=251.888780909464, logdetD=-0.0004208797235257, t=8.81326
  # xi*(L+R)/c=293.870244394374, logdetD=-7.62850679760849e-05, t=8.17316
  # xi*(L+R)/c=335.851707879285, logdetD=-1.38130965079446e-05, t=7.38953
  # xi*(L+R)/c=377.833171364195, logdetD=-2.50883086124893e-06, t=6.88574
  #
  # 684 determinants computed
  # stop time: Thu Jan  3 14:09:25 2019
  #
  # L/R, L, R, T, ldim, E*(L+R)/(hbar*c)
  0.02, 1e-06, 5e-05, 300, 351, -83.78391970580533

The energy printed in the last line assumes a Drude model for the zero-th
Matsubara frequency. If you want to use the plasma model for the zero-th
Matsubara frequency, you can use the value given by ``# plasma =``. This
number, i.e., -12.98... is given in units of :math:`k_\mathrm{B}T/2` and
corresponds to the additional contribution in the high-temperature limit to the
energy in the plasma model. In this example, the free energy using the Drude
model for zero-frequency is

.. math::
  \mathcal{F}_\mathrm{Drude} \approx -83.8 \frac{\hbar c}{L+R} \approx -5.19\times10^{-20}\mathrm{J},

and assuming the plasma model for zero frequency

.. math::
  \mathcal{F}_\mathrm{plasma} \approx \mathcal{F}_\mathrm{Drude} + \frac{-12.98 k_\mathrm{B}T}{2} \approx -7.88 \times 10^{-20} \mathrm{J} .

Truncation of the vector space
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The truncation of the vector space is described in more detail in `Hartmann,
"Casimir effect in the plane-sphere geometry: Beyond the proximity force
approximation", phd thesis (2018)
<https://opus.bibliothek.uni-augsburg.de/opus4/44798>`_.  You can either
specify the dimension of the vector space using ``--ldim``, or you choose the
vector space using the parameter ``--eta``:

.. math::
	\ell_\mathrm{dim} = \mathrm{max}\left(20, \left\lceil \eta R/L \right\rceil\right) .

The estimated error due to the truncation of the vector space depending on eta
is given by:

+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| numerical error | :math:`10^{-2}` | :math:`10^{-3}` | :math:`10^{-4}` | :math:`10^{-5}` | :math:`10^{-6}` | :math:`10^{-7}` | :math:`10^{-8}` |
+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+
| eta             | 2.8             | 4               | 5.2             | 6.4             | 7.6             | 8.8             | 10              | 
+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+


Other options
^^^^^^^^^^^^^

The computation of the matrix elements of the round-trip operator contain an
integration. The desired relative error for this integration can be set using
``--iepsrel``. The default value of :math:`10^{-8}` should be sufficient for
almost all purposes. If you want to compute the Casimir energy to very high
accuracy, to :math:`10^{-7}` or better, you might want to set a smaller value.

casimir_logdetD
---------------

cylinder
--------

API Documentation
=================

The documentation of the API is available at ``manual/api.pdf`` or can be
generated running

.. code-block:: console

  $ make doc

in the directory ``src/``. You need doxygen installed on your computer.
