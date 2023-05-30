"""
code for delensing the CMB and extracting r and Alens even in the presence of systematics
"""

from CAMBPowerSpectrum.params_io import config_obj
from CAMBPowerSpectrum.camb_power_spectrum import CAMBPowerSpectrum
from CAMBPowerSpectrum.cl_plotting import flatmap