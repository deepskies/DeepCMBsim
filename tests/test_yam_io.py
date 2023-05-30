"""
tests yam_io.py
"""

from simcmb.params_io import config_obj


def test_update_val():
    co = config_obj()
    co.update_val("Alens", 10)
    co.update_val("InitPower.r", 10)
    assert [co.CAMBparams.Alens, co.CAMBparams.InitPower.r] == [10, 10]
