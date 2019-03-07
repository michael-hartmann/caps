---
title: 'CaPS: The Casimir effect in the plane-sphere geometry'
tags:
  - Casimir effect
  - electromanetic scattering
  - multipole basis
  - plane-sphere geometry
  - hierarchical matrices
authors:
  - name: Michael Hartmann
    orcid: 0000-0002-9245-8582
    affiliation: 1
  - name: Gert-Ludwig Ingold
    orcid: 0000-0003-4851-2198
    affiliation: 1
affiliations:
  - name: Institut für Physik, Universität Augsburg, 86135 Augsburg, Germany
    index: 1
date: 27th February 2019
bibliography: paper.bib
---

----------------------

# Summary

CaPS allows you to compute the Casimir free energy $\mathcal{F}$ for the
geometry of a sphere with radius $R$ separated by a minimum distance $L$ from a
plane as sketched in Figure 1. We assume that the plane and the sphere are
surrounded by vacuum. The plane-sphere geometry is important since the most
sensitive measurements of the Casimir effect employ this configuration. In
such experiments the aspect ratio $R/L$ is typically of the order
$R/L\sim1000$.

![Geometry of the plane-sphere setup. A sphere with radius $R$ is separated by
the distance $L$ from an infintely extended plane. The aspect ratio $R/L=2$ in
this Figure is about three orders of magnitudes smaller than in typical
experiments.](geometry.pdf)

Within the scattering approach [@lambrecht_njp_2006; @emig_prl_2007], the
Casimir free energy is given as a sum $$\mathcal{F} = \frac{k_\mathrm{B}T}{2}
\sum_{n=-\infty}^\infty \log\det\left(1-\mathcal{M}(|\xi_n|)\right)$$ over the
Matsubara frequencies $\xi_n=2\pi n k_\mathrm{B}T/\hbar$. The round-trip
operator $\mathcal{M}$ represents a complete round-trip of an electromagnetic
wave between the sphere and the plane. Commonly, the round-trip operator is
expressed in the multipole basis. The dimension of $\mathcal{M}$ in this basis
is infinite and a numerical evaluation of the determinants requires a
truncation of the vector space. In order to achieve a sufficient accuracy, one
must scale the dimension of $\mathcal{M}$ with the aspect ratio $R/L$. Hence,
the matrices for aspect ratios used in typical experiments are huge. Moreover,
the matrices are ill-conditioned rendering a numerical evaluation difficult.
For these reasons, numerical evaluations were limited to aspect ratios
$R/L\lesssim100$ [@Durand_phd].

To overcome this problem, we use a symmetrized version of the round-trip
operator $\mathcal{M}$ as described in [@hartmann_phscr_2018;
@hartmann_phd_2018]. Moreover, it turns out that the matrices of the
symmetrized round-trip operator are hierarchical off-diagonal low-rank (HODLR)
matrices. HODLR matrices allow a fast computation of the determinant
[@ambikasaran_josc_2013; @ambikasaran_arxiv_2014]. In particular, we use
HODLRlib [@ambikasaran_joss_2019] to compute the determinants. Explicit
expressions for the matrix elements of the symmetrized round-trip operator and
further information can be found in [@hartmann_phscr_2018; @hartmann_phd_2018].

The free energy $\mathcal{F}$ depends on the radius $R$ of the sphere, the
separation $L$ between plane and sphere, the temperature $T$, and the optical
properties of the plane and sphere. The optical properties of the plane and the
sphere can be modelled as perfect reflectors, using the Drude or plasma model,
or by a user-defined dielectric function. The library also contains functions
to compute the Casimir effect in the high-temperature limit for perfect
reflectors, and the Drude and plasma model. In the Drude case, an analytical
expression is used instead of a numerical evaluation [@bimonte_prl_2012].

In summary, ``CaPS`` has the following features:

 - Computation of the free energy for aspect ratios used in typical experiments
 - Computation of the free energy in the high-temperature limit for perfect reflectors, and metals described by the Drude or plasma model
 - Full support for perfect reflectors, metals described by the Drude and plasma model, and generic materials described by a user-defined dielectric function
 - Computation of the free energy at zero temperature for perfect reflectors in the plane-cylinder geometry
 - Support for parallelization using MPI


``CaPS`` has been used to analyze negative Casimir entropies [@ingold_pre_2015;
@umrath_pre_2015], and to study corrections to the widely used proximity force
approximation for experimentally relevant parameters [@hartmann_prl_2017]. The
package is released under the GPLv2 license.

# Acknowledgements

The authors would like to thank Sivaram Ambikasaran and Shyam Sundar Sankaran
for their help and support.

# References
