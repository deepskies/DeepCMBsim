import os
import camb
import numpy as np
from datetime import datetime as dt
import h5py
import noise

"""
Code to create an array of power spectra from CAMB based on a yaml file
"""

class PS_Maker(object):
    def __init__(self, in_dict): #read in the yaml such that in_dict is yi.Ydict(infile) (this will be done in the ipynb)
        self.Ydict = in_dict
        self.Ydictu = self.Ydict.user_params #just an alias

        # according to the CAMB documentation, errors affect the last "100 or so" multipoles
        self.max_l_use = int( min(self.Ydictu['max_l_use'], noise.max_multipole(self.Ydictu['beam_fwhm'])) )
        self.Ydict.pars.max_l, self.Ydict.pars.max_l_tensor = int(self.max_l_use + self.Ydictu['extra_l']), int(self.max_l_use  + self.Ydictu['extra_l'])

        self._outdir = self.Ydictu['outfile_dir']

        self.cls_raw, self.units = bool(self.Ydictu['cls_raw']), self.Ydictu['TT_dimension']

    def get_noise(self):
        if self.Ydictu['noise_type'] == 'white':
            return noise.white_noise(self.Ydictu['noise_level'], self.Ydictu['beam_fwhm'], self.max_l_use, TT=True, units_uK=True), noise.white_noise(self.Ydictu['noise_level'], self.Ydictu['beam_fwhm'], self.max_l_use, TT=False, units_uK=False)
        elif self.Ydictu['noise_type'] is None:
            return np.zeros(self.max_l_use)
        else:
            print("only white noise is currently implemented")
            return np.zeros(self.max_l_use)

    def get_namestr(self, cpars):
        namestr = f"cls_camb_r{cpars.InitPower.r:0.2f}_A{cpars.Alens:0.2f}_lmax{self.max_l_use}_noise" + str(self.Ydictu['noise_level']) +"." + str(self.Ydictu['beam_fwhm'])
        if self.cls_raw:
            namestr += "_rawCl"
        return namestr

    def get_cls(self, cpars):
        if bool(self.Ydictu["verbose"]):
            ta = dt.now()
        results = camb.get_results(cpars)
        tt, ee, bb, te = results.get_total_cls(raw_cl=self.cls_raw, CMB_unit=self.units)[:self.max_l_use + 1].T
        if self.Ydictu['noise_type'] is not None:
            noise = self.get_noise()
            tt += noise[0]
            ee += noise[1]
            bb += noise[1]
            te += noise[1]
        pp, pt, pe = results.get_lens_potential_cls(raw_cl=self.cls_raw)[:self.max_l_use + 1].T
        lvals = range(self.max_l_use + 1)
        outarr = np.array([lvals, tt, ee, bb, te, pp, pt, pe])
        outlabs = ['l', 'clTT', 'clEE', 'clBB', 'clTE', 'clPP', 'clPT', 'clPE']
        outdict = {}
        for i in range(len(outlabs)):
            outdict[outlabs[i]] = outarr[i]
        outdict['lensed_CLs'] = self.Ydictu['FORCAMB']['DoLensing']

        if bool(self.Ydictu["verbose"]):
            tb = dt.now()
            print('from', dt.strftime(ta, '%H:%M:%S.%f %P'), 'to', dt.strftime(tb, '%H:%M:%S.%f %P'), end=" ")
            print('or', str((tb - ta).seconds) + '.' + str((tb - ta).microseconds), 'seconds total')

        outdict['r'] = cpars.InitPower.r
        outdict['Alens'] = cpars.Alens

        return outdict

    def loop_cls_rA(self, overwrite=False):
        for rr in self.Ydict.rs: #for par in self.Ydict['pass parameter']
            self.Ydict.pars.InitPower.r = rr
            for aa in self.Ydict.As:
                self.Ydict.pars.Alens = aa
                cpars_cur = self.Ydict.pars
                if overwrite:
                    outdict, namestr = self.get_cls(cpars_cur), self.get_namestr(cpars_cur)
                    savecls(outdict, os.path.join(self._outdir, namestr))
                else:
                    namestr = self.get_namestr(cpars_cur)
                    if sum([namestr in x for x in os.listdir(self._outdir)])==0:
                        outdict, namestr = self.get_cls(cpars_cur), self.get_namestr(cpars_cur)
                        savecls(outdict, os.path.join(self._outdir, namestr))

    def update_vals(self, attr, new_val, incamb = False):
        if incamb:
            camb_attr = "Initpower."+attr if attr=="r" else attr
            setattr(self.Ydict.pars, camb_attr, new_val)
        else:
            self.Ydictu[attr] = new_val

def savecls(out_dict, out_name): #move this to yam_in.py
    with h5py.File(out_name + '.h5', 'w') as f:
        for k, v in out_dict.items():
            f.create_dataset(k, data = v)

