"""
tests camb_ps_maker.py
"""

import numpy as np
import deepcmbsim.noise
from pytest import approx


def test_max_multipole():
    fwhm_arcmin = 5
    additional_factor = 3
    max_multipole_5_3 = deepcmbsim.noise.max_multipole(fwhm_arcmin, additional_factor=additional_factor)
    assert max_multipole_5_3 == 180*60*additional_factor/fwhm_arcmin


def test_detector_white_noise():
    arcmin_to_rad = np.pi / 180 / 60
    ell = int(1e4)
    noise_uK_arcmin = 10
    fwhm_arcmin = 1
    out = (noise_uK_arcmin * arcmin_to_rad)**2 * np.exp( ell * (ell + 1) * (fwhm_arcmin * arcmin_to_rad)**2 / (8 * np.log(2)) )
    final_wn_10_1_1e4 = deepcmbsim.noise.detector_white_noise(noise_uK_arcmin, fwhm_arcmin, ell, TT=True, units_uK = True)[-1]
    assert final_wn_10_1_1e4 == approx(out)

