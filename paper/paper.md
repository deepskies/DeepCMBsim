---
title: 'Title Title Title Title'
tags:
  - Python
  - astronomy
  - cosmology
  - cosmic microwave background
  - gravitational lensing
authors:
  - name: Samuel McDermott
    orcid: 0000-0001-5513-1938
    equal-contrib: true
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Camille Avestruz
    orcid: 0000-0000-0000-0000
    equal-contrib: true
    affiliation: "1, 2" # (Multiple affiliations must be quoted
  - name: Brian Nord
    orcid: 0000-0000-0000-0000
    equal-contrib: true
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Humna Awan
    orcid: 0000-0000-0000-0000
    equal-contrib: true
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Maggie Voetberg
    orcid: 0000-0000-0000-0000
    equal-contrib: true
    affiliation: "1, 2" # (Multiple affiliations must be quoted)

affiliations:
 - name: University of Chicago Department of Astronomy and Astrophysics
   index: 1
 - name: Affil 2
   index: 2
 - name: Affil 3
   index: 3
date: 13 August 2023
bibliography: paper.bib


---

[//]: # (based on: https://joss.theoj.org/papers/10.21105/joss.00388)

# Summary

[//]: # (1. Cosmic microwave background and lensing)
[//]: # (2. CAMB)
[//]: # (3. Fast, so that it can make lots and lots at multiple levels of fidelity.)
[//]: # (4. Useful for computational experiments with machine learning and SBI)
[//]: # (5. We create a simple interface for building simulations around camb and namaster.)

The cosmic microwave background (CMB) radiation is a direct link to the earliest moments after the birth of the Universe.
Patterns in the CMB can tell us about the contents of the Universe at otherwise-inaccessible times and at otherwise-unacheivable energies.

This information is encoded in the two-point correlation function of the hotspots of CMB, or its "power spectrum".
However, this information can be difficult to access.
Sources of confusion which can interfere with our ability to observe the CMB include: _lensing_ from structures that grow during the evolution of the universe; _noise_ from foregrounds; and _beam_ artifacts from the instruments we use to perform the measurements.

The `simcmb` package combines these sources of noise in a straightforward and accessible framework that enables fast and realistic simulation of the CMB with lensing and noise.


# Statement of need

[//]: # (1. Most CMB lensing simulators are not use-friendly and don't easily )
[//]: # (2. needed to study B-modes and the effects of additional physics, like different kinds of noise)
[//]: # (3. a bridge to more complex simulations and provides a benchmark for those more high-fidelity simulations)
[//]: # (4. How does this compare in middle level of detail with pixell and others)

Most CMB lensing simulators are not user-friendly.
The `simcmb` package emphasizes user-friendliness by enabling simple variable specification in a `yaml` file, including noise, beam, and lensing.
The `yaml` file contents are used to specify a `CAMBparams` instance which is used by the `CAMB` package `[@Lewis:1999bs; @Howlett:2012mh]` to generate a power spectrum.

The noise level and beam shape are specified in the `yaml` file.
The power spectrum of the noise follows the form in `[@Hu:2001kj]` and relies on the assumption of statistical independence in the Stokes parameters `[@Knox:1995dq; @Zaldarriaga:1996xe]`.
This noisy, lensed power spectrum output is made available to the user.

Optionally, the user can make maps from the power spectra, which is done internally by `simcmb` by calling `namaster` `[@Alonso:2018jzx]`.
We provide intuitive functionality for single maps or for a set of (T,E,B) or (T,Q,U) maps, assuming parity conservation on the sky.


# Workflow

![Example workflow for the `simcmb` package.\label{fig:workflow}](ex_workflow.png)
The package workflow is demonstrated in \autoref{fig:workflow}.


# Citations
1. deeplenstronomy
2. pixell
3. camb
4. webcmb?
5. namaster
6. 

# Acknowledgements

We acknowledge the Deep Skies Lab as a community of multi-domain experts and collaborators who’ve facilitated an environment of open discussion, idea-generation, and collaboration. This community was important for the development of this project.

# References