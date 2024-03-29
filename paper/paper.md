---
title: 'DeepCMBSim'

tags:
  - Python
  - astronomy
  - cosmology
  - cosmic microwave background
  - gravitational lensing

authors:
  - name: Samuel D. McDermott
    orcid: 0000-0001-5513-1938
    affiliation: 1
  - name: M. Voetberg
    orcid: 0009-0005-2715-4709
    affiliation: 2
  - name: Humna Awan
    orcid: 0000-0003-2296-7717
    affiliation: 3
  - name: Camille Avestruz
    orcid: 0000-0001-8868-0810
    affiliation: 3
  - name: Samantha Usman
    orcid: 0000-0003-0918-7185
    affiliation: 1
  - name: Ashia Livaudais-Lewis
    orcid: 0000-0003-3734-335X
    affiliation: 2
  - name: Brian Nord
    orcid: 0000-0001-6706-8972
    affiliation: "1, 2, 4, 5" # (Multiple affiliations must be quoted)

affiliations:
 - name: University of Chicago Department of Astronomy and Astrophysics
   index: 1
 - name: Fermilab Scientific Computing Division
   index: 2
 - name: University of Michigan
   index: 3
 - name: University of Chicago Kavli Institute for Cosmological Physics
   index: 4
 - name: MIT Laboratory for Nuclear Science
   index: 5
date: 14 June 2023
bibliography: paper.bib


---

# Summary

The cosmic microwave background (CMB) radiation is a direct link to the earliest moments after the birth of the Universe.
The CMB is nearly, but not exactly, homogeneous, with small temperature anistropies that deviate from the average temperature at the level of a part in a million.
Patterns in the anisotropies in the CMB can tell us about the contents of the Universe at otherwise-inaccessible times and at otherwise-unacheivable energies.

Because the patterns in the anisotropy field are derived from underlying physical mechanisms that are described by an almost perfectly Gaussian random process, most of the information content of the CMB is encoded in the two-point correlation function of the hotspots of CMB, or its "power spectrum".
However, because the magnitude of temperature anisotropies are at the level of a part in a million, as mentioned above, the information encoded in the CMB can be difficult to access.
Sources of confusion which can interfere with our ability to observe the CMB include: _lensing_ from structures that grow during the evolution of the universe; _noise_ from foregrounds; and _beam_ artifacts from the instruments we use to perform the measurements.

The `DeepCMBsim` package combines these physical process and these sources of noise in a straightforward and accessible framework that enables fast and realistic simulation of the CMB with lensing and noise.


# Statement of need

We present a user-friendly interface to create fast CMB lensing simulations, potentially for applications with Simulation Based Inference analysis approaches.
The `DeepCMBsim` package emphasizes user-friendliness by enabling simple variable specification in a `yaml` file, including baseline examples that incorporate noise, beam-size, and lensing effects.
This will be the first part of a simulation-based inference module for estimation of cosmological parameters.

# Workflow

![Example workflow for the `DeepCMBsim` package.\label{fig:workflow}](ex_workflow.pdf)
The package workflow is demonstrated in \autoref{fig:workflow}. Two `yaml` files are provided in `simcmb/settings`.
One of these files (`base_config.yaml`) contains a number of parameters that allow reproduction of the Planck 2018 cosmology `[@Planck:2018vyg]`.
The user does not need to interact with or edit it for basic functionality, but it specifies a number of parameters necessary for correct `CAMB` functionality.
The `user_config.yaml` file is the main interface for the user, into which different cosmological parameters and experimental descriptions can be entered.
These `yaml` files contents are used to specify a `CAMBparams` instance, which is accessed as an attribute of a `config_obj` class in the `params_io` module.

The noise level and beam shape are also specified in the `yaml` file.
The power spectrum of the noise follows the form in `[@Hu:2001kj]` and relies on the assumption of statistical independence in the Stokes parameters `[@Knox:1995dq; @Zaldarriaga:1996xe]`.

The primary physics module is `camb_power_spectrum`, which defines the `CAMBPowerSpectrum` class.
This calls `CAMB` `[@Lewis:1999bs; @Howlett:2012mh]`
This noisy, lensed power spectrum output is made available to the user.

Optionally, the user can make maps from the power spectra, which is done internally by `DeepCMBsim` by calling `namaster` `[@Alonso:2018jzx]`.
We provide intuitive functionality for single maps or for a set of (T,E,B) or (T,Q,U) maps, assuming parity conservation on the sky.

We provide an example notebook in `notebooks/simcmb_example.ipynb` which demonstrates core package functionality.


# Acknowledgements

We acknowledge the Deep Skies Lab as a community of multi-domain experts and collaborators who’ve facilitated an environment of open discussion, idea-generation, and collaboration. This community was important for the development of this project.

# References
