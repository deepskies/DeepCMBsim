# simcmb

Code for delensing the CMB.

This code relies on `camb` to generate power spectra and `namaster` to simulate CMB temperature maps. Then we use the simulation-based inference code `sbi` to train a machine learning algorithm to extract fundamental physics parameters from (noisy, lensed) CMB maps.