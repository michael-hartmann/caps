---
title: 'CaPS: Casimir Effect in the Plane-Sphere Geometry'
tags:
  - Casimir effect
  - electromagnetic scattering
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
date: 23th June 2019
bibliography: paper.bib
---

----------------------

# Summary

CaPS is a package for the analysis of the Casimir effect in the plane-sphere
geometry. This geometry depicted in Fig. 1 is most commonly used in precision
measurements of the Casimir force arising from quantum and thermal fluctuations
of the electromagnetic field [@bordag_advances]. Specifically, CaPS allows to
compute the Casimir free energy as a function of the sphere radius $R$, the
minimal separation $L$ between sphere and plane, the temperature $T$, and the
material properties of plane and sphere. It is assumed that both objects are
non-magnetic and placed in vacuum. In typical experiments the aspect ratio
$R/L$ is of the order of 1000. The main purpose of this package is to
make aspect ratios as large as $R/L\sim5000$ accessible. Higher aspect ratios
are usually well covered by the proximity force approximation.

![Geometry of the plane-sphere setup: A sphere with radius $R$ is separated by
the distance $L$ from an infinitely extended plate. In typical experiments,
the aspect ratio $R/L$ is about three orders of magnitude larger than
shown here.](geometry.pdf)

Within the scattering approach [@lambrecht_njp_2006; @emig_prl_2007], the
Casimir free energy is given as a Matsubara sum
$$\mathcal{F} = \frac{k_\mathrm{B}T}{2} \sum_{n=-\infty}^\infty \log\det\left(1-\mathcal{M}(|\xi_n|)\right)$$
with $\xi_n=2\pi n k_\mathrm{B}T/\hbar$. $k_\mathrm{B}$ and $\hbar$ are the Boltzmann
constant and Planck constant, respectively. The round-trip operator $\mathcal{M}$
represents a complete round trip of an electromagnetic wave between the sphere
and the plane. Commonly, the round-trip operator is expressed in the multipole
basis. The dimension of $\mathcal{M}$ in this basis is infinite and a numerical
evaluation of the determinants requires a truncation of the vector space. In
order to achieve a sufficient accuracy, one must scale the dimension of
$\mathcal{M}$ with the aspect ratio $R/L$. Hence, the matrices for aspect
ratios used in typical experiments are huge. Moreover, the matrices are
ill-conditioned rendering a numerical evaluation difficult. These problems
limit the aspect ratios accessible in standard implementations to
$R/L\lesssim100$ [@Durand_phd].

CaPS addresses these issues by using a symmetrized version of the round-trip
operator $\mathcal{M}$ as described in [@hartmann_phscr_2018;
@hartmann_phd_2018]. The matrix representation of the symmetrized round-trip
operator yields hierarchical off-diagonal low-rank (HODLR) matrices which allow
for a fast computation of determinants [@ambikasaran_josc_2013;
@ambikasaran_arxiv_2014]. Specifically, we use HODLRlib
[@ambikasaran_joss_2019] for this purpose. Further information including
explicit expressions for the matrix elements of the symmetrized round-trip
operator can be found in [@hartmann_phscr_2018; @hartmann_phd_2018].

``CaPS`` provides the following main features:

 - Computation of the free energy for aspect ratios used in typical experiments.
 - Full support for perfect reflectors, metals described by the Drude and plasma model, and generic materials described by a user-defined dielectric function.
 - Support for parallelization using MPI.
 - Computation of the free energy in the high-temperature limit for perfect reflectors and metals described by the Drude or plasma model.

The computation of the high-temperature limit for the Drude model is based on
[@bimonte_prl_2012].

Basic support for further geometries is provided for the special case of zero
temperature and perfect reflectors:

 - Computation of the free energy in the plane-cylinder geometry.
 - Computation of the free energy for two spheres with equal radii.

The implementation for the plane-cylinder geometry is based on [@emig_prl_2006].

``CaPS`` has been used to analyze negative Casimir entropies [@ingold_pre_2015;
@umrath_pre_2015], and to study corrections to the widely used proximity force
approximation for experimentally relevant parameters [@hartmann_prl_2017]. The
package is released under the GPLv2 license.

# Acknowledgements

The authors thank the DAAD for financial support through the PROBRAL program.
Fruitful discussions with Paulo Maia Neto and the help provided by 
Sivaram Ambikasaran and Shyam Sundar Sankaran are gratefully acknowledged.

# References
