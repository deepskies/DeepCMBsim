"""
code for delensing the CMB and extracting r and Alens even in the presence of systematics
"""

from simcmb.params_io import config_obj
from simcmb.camb_power_spectrum import CAMBPowerSpectrum

try:
    from simcmb.cl_plotting import flatmap
except ModuleNotFoundError:
    pass