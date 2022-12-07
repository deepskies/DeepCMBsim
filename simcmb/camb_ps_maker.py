import os
import camb
import sys
import json
import numpy as np
from datetime import datetime as dt

"""
Code to create a single power spectrum or map from CAMB and/or namaster
Uses the .ini file in inifiles/ as a baseline, changes r and A_lens
"""

# _parentdir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
_basedir = os.getcwd()

base_pars = camb.read_ini(os.path.join(_basedir, 'inifiles', 'planck_2018_1e4.ini'))

max_l_calc, max_l_use = 10100, 10000

base_pars.max_l, base_pars.max_l_tensor = max_l_calc, max_l_calc

flags = ''.join([x for x in sys.argv if '-' in x])

# the following line determines the normalization of the C_ell's
# if this evaluates to False, which is the default, then the output
# power spectra include ell*(ell+1)/2/π for TT, EE, BB, TE and
# [ell*(ell+1)]**2/2/π for the PP, PT, and PE
# Tried to make it user-friendly so that you can write anything
# as long as it has both "raw" and "cl"
cls_raw = True if (('raw' in flags) and ('cl' in flags)) else False

# the following line determines the units of the TT and TE columns
# if this evaluates to 'muK', which is the default, then the output
# power spectra for TT carry units of µK^2 and for TE carry µK.
# This does not affect PT or PE, which are kept unitless.
units = None if (('dimensionless' in flags) and ('TT' in flags)) else 'muK'

# the following line determines the normalization of the C_ell's
# if this evaluates to False, which is the default, then the output
# power spectra include ell*(ell+1)/2/π for TT, EE, BB, TE and
# [ell*(ell+1)]**2/2/π for the PP, PT, and PE
# Tried to make it user-friendly so that you can write anything
# as long as it has both "raw" and "cl"
saveflatmap = True if (('save' in flags) and (('flat' in flags) or ('fm' in flags))) else False

ta = dt.now()

with open(os.path.join(_basedir, 'inifiles/config.json'), 'r') as j:
    j_data = json.load(j)

rr, aa = j_data['log10_r'], j_data['Alens']

base_pars.InitPower.r, base_pars.Alens = 10**rr, aa

results = camb.get_results(base_pars)

tt, ee, bb, te = results.get_total_cls(raw_cl=cls_raw, CMB_unit=units).T
pp, pt, pe = results.get_lens_potential_cls(raw_cl=cls_raw)[:max_l_use + 1].T
lvals = range(max_l_use + 1)

outarr = np.array([lvals, tt, ee, bb, te, pp, pt, pe]).T

namestr = "lr" + f'{rr:0.2f}' + "_A" + f'{aa:0.2f}' + "_d" + dt.strftime(dt.now(), '%y%m%d')
if cls_raw:
    namestr += "_rawCl"

if saveflatmap:
    import flatmaps as fm
    from astropy.wcs import WCS
    import pymaster as nmt
    pixels = 192.  # 192 pixels on each side
    side = 5  # 5 degrees on each side
    reso = side / pixels
    reso_arcmin = reso * 60
    dx = reso * np.pi / 180.0  # Resolution in radians
    lmax = 180. / reso  # Maximum l-mode achievable given these parameters
    lstep = lmax * 2 / pixels
    tfac = dx / pixels  # converts from pixels to radians
    w = WCS(naxis=2)
    nx = int(pixels)
    ny = int(pixels)
    w.wcs.crpix = [nx / 2, ny / 2]  # Center pixel X, Y
    w.wcs.cdelt = np.array([-reso, reso])
    w.wcs.crval = [0, 0]  # Center coordinates RA, DEC at 0,0
    w.wcs.ctype = ["RA---AIR", "DEC--AIR"]  # Airy projection; can be adjusted. Previous used Azimuthal equal-area
    fmi = fm.FlatMapInfo(w, nx=nx, ny=ny, lx=side, ly=side)
    sim_map = nmt.synfast_flat(int(fmi.nx), int(fmi.ny), fmi.lx_rad, fmi.ly_rad, [tt], [0], seed=0)
#    final step here: figure out how to save this image
else:
    np.savetxt(os.path.join(_basedir, j_data['outfiles'], namestr + '.txt'), outarr,
               header="d = " + dt.strftime(dt.now(), '%a, %b %d %Y, %I:%M:%S.%f %p') +
                      "\nusing log_10(r) = " + f'{rr:0.2f}' + ", A = " + f'{aa:0.2f}' +
                      '\nconfigured with json file ' + str(j_data)
               )

tb = dt.now()

if 'v' in flags:
    print('from', dt.strftime(ta, '%H:%M:%S.%f %P'), 'to', dt.strftime(tb, '%H:%M:%S.%f %P'), 'or', end=" ")
    print(str((tb - ta).seconds) + '.' + str((tb - ta).microseconds), 'seconds total')
