# simcmb

Code for producing noisy, lensed simulations of the CMB.

This code relies on `camb` to generate power spectra and (optionally) `namaster` to simulate CMB temperature maps. If you have a newer Mac with Apple Silicon (eg, M1 or M2 chip) you may have issues with `namaster`. If you use `conda` or `mamba` for managing packages, you will need to follow the instructions [here](https://conda-forge.org/docs/user/tipsandtricks.html#installing-apple-intel-packages-on-apple-silicon). If you use `poetry` to manage software, you can try to implement the solution [here](https://github.com/python-poetry/poetry/issues/3458#issuecomment-766195534).

From the top-level directory, you will do `pip install .` like normal.

The usage of the code is documented in `notebooks/simcmb_example.ipynb`, and a simple bash script that you can modify for your own purposes is given in `simcmb/simcmb.py`