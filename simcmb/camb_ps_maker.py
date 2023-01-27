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

max_l_use = int(base_pars.max_l - 100) #according to the CAMB documentation, errors affect the last "100 or so" multipoles

with open(os.path.join(_basedir, 'inifiles/config.json'), 'r') as j:
    j_data = json.load(j)

cls_raw, units = bool(j_data['cls_raw']), bool(j_data['TT_dimensionless'])

rr, aa = j_data['log10_r'], j_data['Alens']

base_pars.InitPower.r, base_pars.Alens = 10**rr, aa

if bool(j_data["verbose"]):
    ta = dt.now()
results = camb.get_results(base_pars)
if bool(j_data["verbose"]):
    tb = dt.now()
    print('from', dt.strftime(ta, '%H:%M:%S.%f %P'), 'to', dt.strftime(tb, '%H:%M:%S.%f %P'), 'or', end=" ")
    print(str((tb - ta).seconds) + '.' + str((tb - ta).microseconds), 'seconds total')

tt, ee, bb, te = results.get_total_cls(raw_cl=cls_raw, CMB_unit=units).T
pp, pt, pe = results.get_lens_potential_cls(raw_cl=cls_raw)[:max_l_use + 1].T
lvals = range(max_l_use + 1)

outarr = np.array([lvals, tt, ee, bb, te, pp, pt, pe]).T

namestr = "lr" + f'{rr:0.2f}' + "_A" + f'{aa:0.2f}' + "_d" + dt.strftime(dt.now(), '%y%m%d')
if cls_raw:
    namestr += "_rawCl"

outfilename = os.path.join(_basedir, j_data['outfiles'], namestr)

def saveflatmap():
    if bool(j_data["verbose"]):
        ta = dt.now()
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
    np.save(outfilename, sim_map)
    if bool(j_data["verbose"]):
        tb = dt.now()
    if bool(j_data["verbose"]):
        print('from', dt.strftime(ta, '%H:%M:%S.%f %P'), 'to', dt.strftime(tb, '%H:%M:%S.%f %P'), 'or', end=" ")
        print(str((tb - ta).seconds) + '.' + str((tb - ta).microseconds), 'seconds total')

def savecls():
    np.savetxt(outfilename + '.txt', outarr,
               header="originally written at " + dt.strftime(dt.now(), '%a, %b %d %Y, %I:%M:%S.%f %p') +
                      # "\nusing log10(r) = " + f'{rr:0.2f}' + ", A = " + f'{aa:0.2f}' +
                      '\nconfigured with json file ' + str(j_data)
               )

