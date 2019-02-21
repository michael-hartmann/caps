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

general description

derivative expansion / asymptotic aproaches vs scattering theory

performance improvements

parallelization

``libcasimir`` has been used to analyze negative Casimir entropies
[@ingold_pre_2015; @umrath_pre_2015], and to study corrections to the commonly
used proximity force approximation [@hartmann_prl_2017].


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
