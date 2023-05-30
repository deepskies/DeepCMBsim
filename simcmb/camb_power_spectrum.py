import itertools
import camb
import numpy as np
from datetime import datetime as dt
from simcmb import noise
import h5py
import json
from collections.abc import Iterable
import os

"""
Code to create an array of power spectra from CAMB based on a yaml file
includes noise, loops, and option to update parameters
"""


class CAMBPowerSpectrum:
    """
    main object for getting power spectra for set parameters, or looped over values of arbitrary
    numbers of parameters
    """
    def __init__(self, in_config_obj):  # read in the yaml such that in_Ydict is Ydict(infile)
        self.camb_params_to_dict = lambda user_params: in_config_obj.camb_params_to_dict(user_params=user_params)
        self.update_val = lambda k, v: in_config_obj.update_val(k, v)
        self.CAMBparams = in_config_obj.CAMBparams
        self.UserParams = in_config_obj.UserParams

        # according to the CAMB documentation, errors affect the last "100 or so" multipoles
        self.max_l_use = min(self.UserParams['max_l_use'], noise.max_multipole(self.UserParams['beam_fwhm']))
        self.max_l_calc = int(self.max_l_use + self.UserParams['extra_l'])
        self.CAMBparams.max_l = self.max_l_calc
        self.CAMBparams.max_l_tensor = self.max_l_calc
        # could make this^ its own function so that it gets recalculated even when you update params

        self._outdir = self.UserParams['outfile_dir']

        self.normalize_cls = bool(self.UserParams['normalize_cls'])
        self.TT_units = self.UserParams['TT_units']

        # initialize empty dictionaries to be filled in later with results and their params
        self.loop_runids = []
        self.results = {}
        self.result_parameters = {}

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

    def get_cls(self, save_to_dict=None, user_params=True):
        if bool(self.UserParams["verbose"]):
            time_start = dt.now()

        # main calculation: https://camb.readthedocs.io/en/latest/camb.html#camb.get_results
        results = camb.get_results(self.CAMBparams)

        # https://camb.readthedocs.io/en/latest/results.html#camb.results.CAMBdata.get_total_cls
        tt, ee, bb, te = results.get_total_cls(raw_cl=self.normalize_cls, CMB_unit=self.TT_units)[:self.max_l_use + 1].T
        if self.UserParams['noise_type'] is not None:
            _noise = self.get_noise()
            tt += _noise[0]
            ee += _noise[1]
            bb += _noise[1]
            te += _noise[1]

        # https://camb.readthedocs.io/en/latest/results.html#camb.results.CAMBdata.get_lens_potential_cls
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
            self.result_parameters[save_to_dict] = self.camb_params_to_dict(user_params=user_params).copy()
        else:
            return outdict

    def loop_sims(self, user_params=True):
        iterables = self.UserParams['ITERABLES']
        keys, values = list(iterables.keys()), iterables.values()
        for vector in itertools.product(*values):
            for i in range(len(vector)):
                self.update_val(keys[i], vector[i])
            _single_param_id = _generate_run_id()
            self.loop_runids.append(_single_param_id)
            self.get_cls(save_to_dict=_single_param_id, user_params=user_params)

    def savecls(self, savedir=os.path.join(os.getcwd(), "outfiles"),
                saveids=None, randomids=False, permission='w', overwrite=False):
        if not os.path.exists(os.path.join(os.getcwd(), "outfiles")):
            print("making local directory `./outfiles`")
            os.mkdir(os.path.join(os.getcwd(), "outfiles"))
        if saveids is not None:
            if type(saveids) == int:
                saveids = np.random.choice(range(len(self.loop_runids)), saveids, replace=False) if randomids else self.loop_runids[:saveids]
            elif isinstance(saveids, Iterable):
                saveids = [self.loop_runids[x] for x in saveids]
            else:
                saveids = saveids
        else:
            saveids = self.loop_runids

        for run_id in saveids:
            if True or overwrite:  # todo check to see if results with these parameters have already been run
                with h5py.File(os.path.join(savedir, f"{run_id}_results.h5"), permission) as f:
                    for k, v in self.results[run_id].items():
                        f.create_dataset(k, data=v)
                with open(os.path.join(savedir, f"{run_id}_params.yaml"), permission) as f:
                    json.dump(self.result_parameters[run_id], f, default=lambda x: x.tolist())
            else:
                print(f"skipping because {run_id}/parameters already exists and overwrite set to False")


def _generate_run_id(random_digits=6):
    _rint = np.random.randint(10**random_digits)
    return 'runid_'+dt.now().strftime('%y%m%d%H%M%S%f_')+str(_rint).zfill(random_digits)
