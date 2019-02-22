---
title: 'libcasimir: The Casimir effect in the plane-sphere geometry'
tags:
  - Casimir effect
  - electromanetic scattering
  - multipole basis
  - plane-sphere geometry
  - hierarchical matrices
authors:
  - name: Michael Hartmann
    affiliation: 1
  - name: Gert-Ludwig Ingold
    affiliation: 1
affiliations:
  - name: Institut für Physik, Universität Augsburg, 86135 Augsburg, Germany
    index: 1
date: 27th February 2019
bibliography: paper.bib
---

----------------------

# Summary

In 1984 Hendrik Casimir considered two parallel and perfectly reflecting planes
at zero temperature and found an attractive force pushing the planes together
[@casimir_pr_1948]. While the Casimir effect has been measured in a number of
historic experiments [@bordag_advances], accurate measurements are only
available since the late 1990ies. In contrast to Casimir's original setup, to
avoid misalignments between the planes, most modern experiments measure the
Casimir effect in the plane-sphere geometry, as sketeched in Figure 1.
Moreover, an accurate description of current experiments requires to take
finite temperature and the material properties of the plane and the sphere into
account.

![Geometry of the plane-sphere setup. A sphere with radius $R$ is separated by
the smallest distance $L$ from an infintely extended plane.](geometry.pdf)

``libcasimir`` allows to compute the Casimir free energy $\mathcal{F}$ in the
plane-sphere geometry. The free energy $\mathcal{F}$ depends on the radius $R$
of the sphere, the smallest separation $L$ between plane and sphere, the
temperature $T$, and the optical properties of the plane and sphere. We assume
the plane and the sphere to be in vacuum. The optical properties of the plane
and the sphere can be modelled as perfect reflectors, using the Drude or plasma
model, or by a user-defined dielectric function. The library also contains
functions to compute the Casimir effect in the high-temperature limit for
perfect reflectors, and the Drude and plasma model. Finally, there is support
to compute the Casimir free energy in the plane-cylinder geometry for perfect
reflectors at $T=0$.

``libcasimir`` computes the free energy $\mathcal{F}$ using the scattering
approach [@lambrecht_njp_2006; @emig_prl_2007; @rahi_prd_2009]. In contrast to
other methods, the scattering approach is in exact. In particular, this means
that the numerical error can be controlled.

performance improvements Legendre [@bogaert_siam_2012] and HODLR
[@ambikasaran_josc_2013; @ambikasaran_arxiv_2014] using
[@ambikasaran_joss_2019]

``libcasimir`` has been used to analyze negative Casimir entropies
[@ingold_pre_2015; @umrath_pre_2015], and to study corrections to the widely
used proximity force approximation for experimentally relevant parameters
[@hartmann_prl_2017].

In summary, ``libcasimir`` has the following features:

 - Computation of the free energy for aspect ratios used in typical experiments
 - Computation of the free energy in the high-temperature limit for perfect reflectors, and metals described by the Drude or plasma model
 - Full support for perfect reflectors, metals described by the Drude and plasma model, and generic materials described by a user-defined dielectric function
 - Computation of the free energy at zero temperature for perfect reflectors in the plane-cylinder geometry
 - Support for parallelization using MPI

``libcasimir`` is released under the GPLv2 license.

# Acknowledgements

The authors would like to thank Sivaram Ambikasaran and Shyam Sundar Sankaran
for their help and support.

# References
