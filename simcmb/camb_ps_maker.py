import os
import re
import camb
import json
import numpy as np
from datetime import datetime as dt
import h5py
import yam_in as yi
import noise

"""
Code to create a single power spectrum or map from CAMB and/or namaster
Uses the .ini file in inifiles/ as a baseline, changes r and A_lens
"""
class PS_Maker(object):
    def __init__(self, infile = "planck_2018_1e4.yaml"):
        self.infile = infile
        self.Ydict = yi.Ydict(infile = self.infile)

        # according to the CAMB documentation, errors affect the last "100 or so" multipoles
        self.max_l_use = min(self.Ydict.user_params['lmax_use'], noise.max_multipole(self.Ydict.user_params['beam_fwhm']) )

        self._outdir = self.Ydict.user_params['outfile_dir']

        self.cls_raw, self.units = bool(self.Ydict.user_params['cls_raw']), self.Ydict.user_params['TT_dimension']

    def get_cls(self, cpars):
        if bool(self.Ydict.user_params["verbose"]):
            ta = dt.now()
        results = camb.get_results(cpars)
        tt, ee, bb, te = results.get_total_cls(raw_cl=self.cls_raw, CMB_unit=self.units).T
        pp, pt, pe = results.get_lens_potential_cls(raw_cl=self.cls_raw)[:self.max_l_use + 1].T
        lvals = range(self.max_l_use + 1)
        outarr = np.array([lvals, tt, ee, bb, te, pp, pt, pe])
        outlabs = ['l', 'clTT', 'clEE', 'clBB', 'clTE', 'clPP', 'clPT', 'clPE']
        outdict = {}
        for i in range(len(outlabs)):
            outdict[outlabs[i]] = outarr[i]
        outdict['lensed_CLs'] = self.Ydict.user_params['FORCAMB']['DoLensing']

        if bool(self.Ydict.user_params["verbose"]):
            tb = dt.now()
            print('from', dt.strftime(ta, '%H:%M:%S.%f %P'), 'to', dt.strftime(tb, '%H:%M:%S.%f %P'), end=" ")
            print('or', str((tb - ta).seconds) + '.' + str((tb - ta).microseconds), 'seconds total')

        outdict['r'] = cpars.InitPower.r
        outdict['Alens'] = cpars.Alens

        return outdict

    def get_namestr(self, cpars):
        namestr = "r" + f'{cpars.InitPower.r:0.2f}' + "_A" + f'{cpars.Alens:0.2f}' + "_lmax" + str(self.max_l_use) + "_noise" + str(self.Ydict.user_params['noise_level']) +"." + str(self.Ydict.user_params['beam_fwhm'])
        if self.cls_raw:
            namestr += "_rawCl"
        outfilename = os.path.join(self._outdir, namestr)

        return outfilename

    def loop_cls_rA(self, overwrite=False):
        for rr in self.Ydict.rs:
            self.Ydict.pars.InitPower.r = rr
            for aa in self.Ydict.As:
                self.Ydict.pars.Alens = aa
                cpars_cur = self.Ydict.pars
                if overwrite:
                    outdict, namestr = self.get_cls(cpars_cur), self.get_namestr(cpars_cur)
                    self.savecls(outdict, namestr)
                else:
                    namestr = self.get_namestr(cpars_cur)
                    if sum([namestr in x for x in os.listdir(self._outdir)])==0:
                        outdict, namestr = self.get_cls(cpars_cur), self.get_namestr(cpars_cur)
                        self.savecls(outdict, namestr)

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

    def savecls(self, odict, oname):
        # np.savetxt(self.outfilename + '.txt', self.outarr,
        #            header="originally written at " + dt.strftime(dt.now(), '%a, %b %d %Y, %I:%M:%S.%f %p') +
        #                   # "\nusing log10(r) = " + f'{rr:0.2f}' + ", A = " + f'{aa:0.2f}' +
        #                   '\nconfigured with json file ' + str(self.j_data)
        #            )
        with h5py.File(oname + '.h5', 'w') as f:
            for k, v in odict.items():
                f.create_dataset(k, data = v)