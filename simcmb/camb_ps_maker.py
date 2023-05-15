import itertools
import camb
import numpy as np
from datetime import datetime as dt
from simcmb import noise
from collections import namedtuple

"""
Code to create an array of power spectra from CAMB based on a yaml file
includes noise, loops, and option to update parameters
"""

class PS_Maker(object):
    def __init__(self, in_dict): #read in the yaml such that in_dict is yi.Ydict(infile) (this will be done in the ipynb)
        self.Ydict = in_dict
        self.Ydictu = self.Ydict.all_params_dict['USERPARAMS']  # just an alias

        # according to the CAMB documentation, errors affect the last "100 or so" multipoles
        self.max_l_use = int( min(self.Ydictu['max_l_use'], noise.max_multipole(self.Ydictu['beam_fwhm'])) )
        self.Ydict.CAMBparams.max_l, self.Ydict.CAMBparams.max_l_tensor = int(self.max_l_use + self.Ydictu['extra_l']), int(self.max_l_use + self.Ydictu['extra_l'])
        # could make this^ its own function so that it gets recalculated even when you update params

        self._outdir = self.Ydictu['outfile_dir']

        self.cls_raw, self.units = bool(self.Ydictu['cls_raw']), self.Ydictu['TT_dimension']

        self.results = {}  # initialize an empty dictionary that can be filled in later, if desired

    def get_noise(self):
        if self.Ydictu['noise_type'] == 'white':
            return noise.white_noise(self.Ydictu['noise_level'], self.Ydictu['beam_fwhm'], self.max_l_use, TT=True), noise.white_noise(self.Ydictu['noise_level'], self.Ydictu['beam_fwhm'], self.max_l_use, TT=False)
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

    def get_cls(self, cpars, save_to_dict = None):
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

        if bool(self.Ydictu["verbose"]):
            tb = dt.now()
            print('from', dt.strftime(ta, '%H:%M:%S.%f %P'), 'to', dt.strftime(tb, '%H:%M:%S.%f %P'), end=" ")
            print('or', str((tb - ta).seconds) + '.' + str((tb - ta).microseconds), 'seconds total')

        if save_to_dict is not None:
            self.results[save_to_dict] = outdict
        else:
            return outdict

    def loop_sims(self, overwrite=False):#make this for general arguments!
        # old version, only for r and Alens:
        # for rr in self.Ydict.rs: #for par in self.Ydict['pass parameter'] - use itertools.product!
        #     self.Ydict.CAMBparams.InitPower.r = rr
        #     for aa in self.Ydict.As:
        #         self.Ydict.CAMBparams.Alens = aa
        #         cpars_cur = self.Ydict.CAMBparams
        #         if overwrite:
        #             outdict, namestr = self.get_cls(cpars_cur), self.get_namestr(cpars_cur)
        #             savecls(outdict, os.path.join(self._outdir, namestr))#replace this with something other than save!
        #         else:
        #             namestr = self.get_namestr(cpars_cur)
        #             if sum([namestr in x for x in os.listdir(self._outdir)])==0:
        #                 outdict, namestr = self.get_cls(cpars_cur), self.get_namestr(cpars_cur)
        #                 savecls(outdict, os.path.join(self._outdir, namestr))
        keys, values = list(self.Ydict.dict_iterables.keys()), self.Ydict.dict_iterables.values()
        for vector in itertools.product(*values):
            for i in range(len(vector)):
                self.Ydict.update_val(keys[i], vector[i])
            self.get_cls(self.Ydict.CAMBparams, save_to_dict=zip(keys, vector))

        return self.results



def savecls(all_sims, out_name, sims_to_save_start=None, sims_to_save_end=None, permission='r+', overwrite=False):
    if (sims_to_save_end is None) and (sims_to_save_start is None):
        sims_to_save_start = 0
        sims_to_save_end = len(all_sims)
    with h5py.File(out_name + '.h5', permission) as f:
        for i in range(sims_to_save_start, sims_to_save_end):
            out_dict = all_sims[i]
            for k, v in out_dict.items():
                try:
                    f.create_dataset(f"r{out_dict['r']}/Alens{out_dict['Alens']}/{k}", data=v)
                except ValueError:
                    if overwrite:
                        f[f"r{out_dict['r']}/Alens{out_dict['Alens']}/{k}"] = v
                    else:
                        print(
                            f"skipping because r{out_dict['r']}/Alens{out_dict['Alens']}/{k} already exists and overwrite set to False")
