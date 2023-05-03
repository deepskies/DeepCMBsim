"""
code for delensing the CMB and extracting r and Alens even in the presence of systematics
"""

from simcmb.yam_io import Ydict
from simcmb.camb_ps_maker import PS_Maker
from simcmb.clplotting import flatmap
from simcmb.clplotting import power_spectrum_plot