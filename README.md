# `simcmb`

Code for producing realistic simulations of the CMB with noise and lensing capabilities.

## installation

This code relies on `camb` to generate power spectra and (optionally) `namaster` to simulate CMB temperature maps. We provide an environment specification file for `conda` or `mamba` users at `conda-env.yml`.

If you have a newer Mac with Apple Silicon (eg, M1 or M2 chip) you may have issues with `namaster`. If you use `conda` or `mamba` for managing packages, you will need to follow the instructions [here](https://conda-forge.org/docs/user/tipsandtricks.html#installing-apple-intel-packages-on-apple-silicon).

From the top-level directory, you will do `pip install .` like normal.

## usage

The core functions are `PS_maker`, which relies on an instance of the `Yobj` class

The usage of the code is documented in `notebooks/simcmb_example.ipynb`, and a simple bash script that you can modify for your own purposes is given in `simcmb/simcmb.py`


## citation

If you use this code in your research, please cite our JOSS paper. Please also make use of the citation instructions for `camb` provided [here](https://camb.info).


## contributing

If you would like to contribute, please be in touch with the [authors](#contact)

## contact

The code was developed by [Samuel D. McDermott](https://samueldmcdermott.github.io) and is maintained by the [DeepSkies lab](https://deepskieslab.com)
