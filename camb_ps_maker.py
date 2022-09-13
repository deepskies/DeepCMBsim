import sys, platform, os
import numpy as np

import camb
from camb import model, initialpower


camb_dir = os.path.dirname(camb.__file__)



print('Using CAMB %s installed at %s'%(camb.__version__, camb_dir))



pars=camb.read_ini(os.path.join(camb_dir,'inifiles', sys.argv[1]))



results = camb.get_results(pars)



results.save_cmb_power_spectra("./"+sys.argv[2])