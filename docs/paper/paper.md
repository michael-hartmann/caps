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

This package contains functionality to compute the Casimir free energy
$\mathcal{F}$ for the geometry of a sphere with radius $R$ separated by a
minimum distance $L$ from a plane as shown in Figure 1. We assume the plane and
the sphere to be in vacuum. This geometry is important since the most sensitive
experiments of the Casimir effect typically employ this configuration. In these
experiments the aspect ratio $R/L$ is typically of the order $R/L\sim1000$.

![Geometry of the plane-sphere setup. A sphere with radius $R$ is separated by
the smallest distance $L$ from an infintely extended plane.](geometry.pdf)

Within the scattering approach [@lambrecht_njp_2006; @emig_prl_2007], the
Casimir free energy is given as a sum
$$\mathcal{F} = \frac{k_\mathrm{B}T}{2} \sum_{n=-\infty}^\infty \sum_{m=-\infty}^\infty \log\det\left(1-\mathcal{M}^{(m)}(|\xi_n|)\right)$$
over the Matsubara frequencies $\xi_n=2\pi n k_\mathrm{B}T/\hbar$. The
round-trip operator $\mathcal{M}$ represents a complete round-trip on an
electromagnetic wave between the sphere and the plane. Explicit expressions for
the matrix elements of $\mathcal{M}^{(m)}$ are given in [@hartmann_phscr_2018;
@hartmann_phd_2018].

The free energy $\mathcal{F}$ depends on the radius $R$ of the sphere, the
separation $L$ between plane and sphere, the temperature $T$, and the optical
properties of the plane and sphere. The optical properties of the plane and the
sphere can be modelled as perfect reflectors, using the Drude or plasma model,
or by a user-defined dielectric function. The library also contains functions
to compute the Casimir effect in the high-temperature limit for perfect
reflectors, and the Drude and plasma model. In the Drude case, an analytical
expression is used instead of a numerical evaluation [@bimonte_prl_2012].

``libcasimir`` allows a fast evaluation of the determinants appearing in
$\mathcal{F}$ In contrast to other methods, the scattering approach is in
exact. In particular, this means that the numerical error can be controlled.

performance improvements due to symmetrization [@hartmann_phscr_2018], fast Legendre polynomials [@bogaert_siam_2012] and HODLR
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
