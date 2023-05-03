import itertools
import camb
import numpy as np
from datetime import datetime as dt
from simcmb import noise

"""
Code to create an array of power spectra from CAMB based on a yaml file
includes noise, loops, and option to update parameters
"""


class PS_Maker:
    """
    main object for getting power spectra for set parameters, or looped over values of arbitrary
    numbers of parameters
    """
    def __init__(self, in_Ydict):  # read in the yaml such that in_Ydict is Ydict(infile)
        self.Ydict = in_Ydict
        self.CAMBparams = self.Ydict.CAMBparams
        self.UserParams = self.Ydict.UserParams

        # according to the CAMB documentation, errors affect the last "100 or so" multipoles
        self.max_l_use = min(self.UserParams['max_l_use'], noise.max_multipole(self.UserParams['beam_fwhm']))
        self.max_l_calc = int(self.max_l_use + self.UserParams['extra_l'])
        self.CAMBparams.max_l = self.max_l_calc
        self.CAMBparams.max_l_tensor = self.max_l_calc
        # could make this^ its own function so that it gets recalculated even when you update params

        self._outdir = self.UserParams['outfile_dir']

        self.normalize_cls = bool(self.UserParams['normalize_cls'])
        self.TT_units = self.UserParams['TT_units']

        self.results = {}  # initialize an empty dictionary that can be filled in later, if desired

    def get_noise(self):
        if self.UserParams['noise_type'] == 'white':
            t_noise = noise.white_noise(self.UserParams['noise_level'], self.UserParams['beam_fwhm'], self.max_l_use,
                                        TT=True)
            eb_noise = noise.white_noise(self.UserParams['noise_level'], self.UserParams['beam_fwhm'], self.max_l_use,
                                         TT=False)
            return t_noise, eb_noise
        elif self.UserParams['noise_type'] is None:
            return np.zeros((2, self.max_l_use))
        else:
            print("only white noise is currently implemented")
            return np.zeros((2, self.max_l_use))

    def get_namestr(self, cpars):
        namestr = f"cls_camb_r{cpars.InitPower.r:0.2f}_A{cpars.Alens:0.2f}_lmax{self.max_l_use}_noise" + str(
            self.UserParams['noise_level']) + "." + str(self.UserParams['beam_fwhm'])
        if self.normalize_cls:
            namestr += "_rawCl"
        return namestr

    def get_cls(self, save_to_dict=None):
        if bool(self.UserParams["verbose"]):
            time_start = dt.now()

        # main calculation: https://camb.readthedocs.io/en/latest/camb.html#camb.get_results
        results = camb.get_results(self.CAMBparams)

        # https://camb.readthedocs.io/en/latest/results.html#camb.results.CAMBdata.get_total_cls
        tt, ee, bb, te = results.get_total_cls(raw_cl=self.normalize_cls, CMB_unit=self.TT_units)[:self.max_l_use + 1].T
        if self.UserParams['noise_type'] is not None:
            noise = self.get_noise()
            tt += noise[0]
            ee += noise[1]
            bb += noise[1]
            te += noise[1]

        #https://camb.readthedocs.io/en/latest/results.html#camb.results.CAMBdata.get_lens_potential_cls
        pp, pt, pe = results.get_lens_potential_cls(raw_cl=self.normalize_cls)[:self.max_l_use + 1].T
        lvals = range(self.max_l_use + 1)
        outdict = {
            'l': lvals,
            'clTT': tt,
            'clEE': ee,
            'clBB': bb,
            'clTE': te,
            'clPP': pp,
            'clPT': pt,
            'clPE': pe
        }

        if bool(self.UserParams["verbose"]):
            time_end = dt.now()
            print('from', dt.strftime(time_start, '%H:%M:%S.%f %P'), 'to', dt.strftime(time_end, '%H:%M:%S.%f %P'),
                  end=" ")
            print('or', str((time_end - time_start).seconds) + '.' + str((time_end - time_start).microseconds),
                  'seconds total')

        if save_to_dict is not None:
            self.results[save_to_dict] = outdict
        else:
            return outdict

    def loop_sims(self, overwrite=False):
        # old version, only for r and Alens:
        # for rr in self.Ydict.rs: #for par in self.Ydict['pass parameter'] - use itertools.product!
        #     self.CAMBparams.InitPower.r = rr
        #     for aa in self.Ydict.As:
        #         self.CAMBparams.Alens = aa
        #         cpars_cur = self.CAMBparams
        #         if overwrite:
        #             outdict, namestr = self.get_cls(cpars_cur), self.get_namestr(cpars_cur)
        #             savecls(outdict, os.path.join(self._outdir, namestr))#replace this with something other than save!
        #         else:
        #             namestr = self.get_namestr(cpars_cur)
        #             if sum([namestr in x for x in os.listdir(self._outdir)])==0:
        #                 outdict, namestr = self.get_cls(cpars_cur), self.get_namestr(cpars_cur)
        #                 savecls(outdict, os.path.join(self._outdir, namestr))
        iterables = self.UserParams['ITERABLES']
        keys, values = list(iterables.keys()), iterables.values()
        for vector in itertools.product(*values):
            for i in range(len(vector)):
                self.Ydict.update_val(keys[i], vector[i])
            self.get_cls(save_to_dict=zip(keys, vector))

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
