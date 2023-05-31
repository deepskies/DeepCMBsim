"""
tests camb_ps_maker.py
"""

import numpy as np
from simcmb.params_io import config_obj
from simcmb.camb_power_spectrum import CAMBPowerSpectrum


def test_get_noise():
    baseline_config_obj = config_obj()
    new_max_l_use = 100
    baseline_config_obj.update_val('max_l_use', new_max_l_use)
    new_extra_l = 0
    baseline_config_obj.update_val('extra_l', new_extra_l)
    assert np.array(CAMBPowerSpectrum(baseline_config_obj).get_noise()).shape == (2, new_max_l_use+new_extra_l+1)


def test_get_cls():
    baseline_config_obj = config_obj()
    baseline_config_obj.update_val("verbose", False)
    baseline_config_obj.update_val("noise_type", None)
    clTT = CAMBPowerSpectrum(baseline_config_obj).get_cls()['clTT']
    assert len(clTT) == baseline_config_obj.UserParams['max_l_use']+1
