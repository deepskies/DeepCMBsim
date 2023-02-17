import os
import re

import camb
import json
import numpy as np
from datetime import datetime as dt

"""
Code to create a single power spectrum or map from CAMB and/or namaster
Uses the .ini file in inifiles/ as a baseline, changes r and A_lens
"""

# _parentdir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
_basedir = os.getcwd()
if 'inifiles' in os.listdir():
    _inidir = os.path.join(_basedir, 'inifiles')
else:
    _basedir = os.path.join(_basedir, 'simcmb')
    _inidir = os.path.join(_basedir, 'inifiles')


class PS_Maker(object):
    def __init__(self, param_file = 'config.json', base_cosmo_param_file = 'planck_2018_1e4.ini'):
        self.base_pars = camb.read_ini(os.path.join(_inidir, base_cosmo_param_file))
        # according to the CAMB documentation, errors affect the last "100 or so" multipoles
        self.max_l_use = int(self.base_pars.max_l - 100)

        with open(os.path.join(_inidir, param_file), 'r') as j:
            self.j_data = json.load(j)
        self._outdir = os.path.join(_basedir, self.j_data['outfiles'])

        self.cls_raw, self.units = bool(self.j_data['cls_raw']), self.j_data['TT_dimension']

        # self.rr, self.aa = self.j_data['log10_r'], self.j_data['Alens']
        #
        # self.base_pars.InitPower.r, self.base_pars.Alens = 10**self.rr, self.aa
        for X in self.j_data.keys():
            sX = re.split("\.", X)
            try:
                if len(sX) == 1:
                    getattr(self.base_pars, X)
                    setattr(self.base_pars, X, self.j_data[X])
                    print(X, getattr(self.base_pars, X))
                elif len(sX) == 2:
                    getattr(getattr(self.base_pars, sX[0]), sX[1])
                    setattr(getattr(self.base_pars, sX[0]), sX[1], self.j_data[X])
                    print(X, getattr(getattr(self.base_pars, sX[0]), sX[1]))
                # for k in range(len(sX)):
                    #recursively getattr: seems like a while loop should work, but it's a bit tricky because you don't want to overwrite self.base_pars
            except Exception:
                continue

        if bool(self.j_data["verbose"]):
            self.ta = dt.now()
        self.results = camb.get_results(self.base_pars)
        self.tt, self.ee, self.bb, self.te = self.results.get_total_cls(raw_cl=self.cls_raw, CMB_unit=self.units).T
        self.pp, self.pt, self.pe = self.results.get_lens_potential_cls(raw_cl=self.cls_raw)[:self.max_l_use + 1].T
        self.lvals = range(self.max_l_use + 1)
        self.outarr = np.array([self.lvals, self.tt, self.ee, self.bb, self.te, self.pp, self.pt, self.pe]).T
        if bool(self.j_data["verbose"]):
            self.tb = dt.now()
            print('from', dt.strftime(self.ta, '%H:%M:%S.%f %P'), 'to', dt.strftime(self.tb, '%H:%M:%S.%f %P'), end=" ")
            print('or', str((self.tb - self.ta).seconds) + '.' + str((self.tb - self.ta).microseconds), 'seconds total')

        self.namestr = "lr" + f'{self.base_pars.InitPower.r:0.2f}' + "_A" + f'{self.base_pars.Alens:0.2f}' + "_d" + dt.strftime(dt.now(), '%y%m%d')
        if self.cls_raw:
            self.namestr += "_rawCl"

        self.outfilename = os.path.join(self._outdir, self.namestr)

    def saveflatmap(self):
        if bool(self.j_data["verbose"]):
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
        sim_map = nmt.synfast_flat(int(fmi.nx), int(fmi.ny), fmi.lx_rad, fmi.ly_rad, [self.tt], [0], seed=0)
        np.save(self.outfilename, sim_map)
        if bool(self.j_data["verbose"]):
            tb = dt.now()
            print('from', dt.strftime(ta, '%H:%M:%S.%f %P'), 'to', dt.strftime(tb, '%H:%M:%S.%f %P'), 'or', end=" ")
            print(str((tb - ta).seconds) + '.' + str((tb - ta).microseconds), 'seconds total')

    def savecls(self):
        np.savetxt(self.outfilename + '.txt', self.outarr,
                   header="originally written at " + dt.strftime(dt.now(), '%a, %b %d %Y, %I:%M:%S.%f %p') +
                          # "\nusing log10(r) = " + f'{rr:0.2f}' + ", A = " + f'{aa:0.2f}' +
                          '\nconfigured with json file ' + str(self.j_data)
                   )
