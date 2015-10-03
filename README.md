libcasimir
==========

The Casimir effect
------------------

In 1984 Hendrik Casimir considered two parallel, perfectly conducting plates in vacuum at temperature T=0 and predicted an attracting force [1]. This force arises due to vacuum fluctuations and was experimentally verified in 1956 by Derjaguin, Abrikosova and Lifshitz [2], as well as in 1958 by Sparnaay [3].

Since the Casimir effect is a manifestation of vacuum fluctuations in the mesoscopic world, it has relations to many open physical questions. As the Casimir force is the dominant force between electically neutral bodies at micron or sub-micron distances, the Casimir effect plays an important role in the search for new hypothetical forces predicted by unified theories [4]. The Casimir effect is also linked to the theory geometry and the puzzle of the cosmological constant. All energy gravitates and thus zero point fluctuations are expected to contribute to the stress–energy tensor in Einstein's field equations [5]. In fact, several cosmological observations like the discovery that the universe expands in an increasing rate [6] suggest a non-zero energy density of vacuum. However, estimations of the cosmological constant and meassurements disagree by about 120 orders of magnitude.

Moreover, negative entropies are found for some geometries and parameters in the Casimir effect. Negative entropies, for example, occur in the plane–plane geometry for metals described by the Drude model. In addition, this effect also occurs in the plane–sphere geometry even for perfect reflectors, thus suggesting a geometrical origin of negative entropies [7,8]. This is in general not a problem, since the Casimir free energy is an interaction energy and does not describe the entire physical system. However, the origin of negative entropies is not understood very well.


What is libcasimir?
-------------------
geometry libcasimir implements the numerics for the Casimir effect in the plane-sphere geometry with perfect spheres using a scattering approach [8, 9]. We use the same approach, but derive slightly different formulas for the matrix element of the round-trip operator that don't need Wigned-d-symbols.

A sphere of radius R is separated by a distance of L from a plane. The plane is infinite in the xy -direction and both plane and sphere are assumed to be perfect reflectors. The programs calculate the free energy F(T,R/L) in scaled quantities:
F scaled = L+R ℏc F T scaled = 2π k (L+R) ℏc T
Using this scaling the free energy only depends on the temperature T and the geometric ratio L/R .

This code is part of my master thesis. This thesis is also a good reference for libcasimir. It costed much time and much work to write a stable and a fast implementation. So, if you find this piece of code useful and you use it for plots, please consider to cite libcasimir.

libcasimir
----------

Features:
 - Calculate the free energy F(T,L/R) for different separations and temperatures
 - Calculate the free energy F(T→∞,L/R) in the high temperature limit
 - Support for 80-bit extended floating point precision and quadrupole precision
 - libcasimir is fast and reliable
 - ready to use programs: you don't have to modify the code
 - libcasimir is free software – you may use it or even modify it
