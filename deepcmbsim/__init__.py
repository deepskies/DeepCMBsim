"""
code for delensing the CMB and extracting r and Alens even in the presence of systematics
"""

from deepcmbsim.params_io import config_obj
from deepcmbsim.camb_power_spectrum import CAMBPowerSpectrum

try:
    from deepcmbsim.cl_plotting import flatmap
except ModuleNotFoundError:
    pass