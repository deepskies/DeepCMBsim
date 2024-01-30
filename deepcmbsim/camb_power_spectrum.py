import itertools
import camb
import numpy as np
from datetime import datetime as dt
from deepcmbsim import noise
import h5py
import json
from collections.abc import Iterable
import os

"""
Code to create an array of power spectra from CAMB based on a yaml file.
Requires a config_obj instance.
"""


class CAMBPowerSpectrum:
    """
    main object for getting power spectra for set parameters, or looped over values of arbitrary
    numbers of parameters
    """
    def __init__(self, in_config_obj):
        """
        Parameters
        ----------
        in_config_obj: config_obj instance
            Object that posesses a CAMBparams instance, a UserParams dict, a camb_params_to_dict method,
            and an update_val method
        """
        self.camb_params_to_dict = lambda user_params: in_config_obj.camb_params_to_dict(user_params=user_params)

        self.CAMBparams = in_config_obj.CAMBparams
        self.UserParams = in_config_obj.UserParams

        self.update_val = lambda k, v: in_config_obj.update_val(k, v, verbose = self.UserParams["verbose"])

        # according to the CAMB documentation, errors affect the last "100 or so" multipoles
        self.max_l_use = min(self.UserParams['max_l_use'], noise.max_multipole(self.UserParams['beamfwhm_arcmin']))
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
        """
        Parameters
        ----------
        None

        Returns
        -------
        np.ndarray
            provides noise for the TT power spectrum and the polarization power spectra;
            shape is (2, max_l_use)
        """
        if self.UserParams['noise_type'] == 'detector-white':
            t_noise = noise.detector_white_noise(self.UserParams['noise_uKarcmin'], self.UserParams['beamfwhm_arcmin'], self.max_l_use,
                                                 TT=True)
            eb_noise = noise.detector_white_noise(self.UserParams['noise_uKarcmin'], self.UserParams['beamfwhm_arcmin'], self.max_l_use,
                                                  TT=False)
            return t_noise, eb_noise
        elif self.UserParams['noise_type'] is None:
            return np.zeros((2, self.max_l_use))
        else:
            print("only detector white noise is currently implemented, via `noise_type = 'detector-white'` in `user_config.yaml`")
            return np.zeros((2, self.max_l_use))

    def get_cls(self, save_to_dict=None, user_params=True):
        """
        central method of CAMBPowerSpectrum, used to obtain power spectra for a given config file

        Parameters
        ----------
        save_to_dict : str, optional
            if specified, will populate the self.results dictionary with the output dictionary
        user_params : bool, default True
            if save_to_dict is not None, then user_params=True saves only the CAMBparams that differ
            from the base_config in the result_parameters dict

        Returns
        -------
        dict
            dictionary of values of l, clTT, clEE, clBB, clTE, clPP, clPT, clPE
        """
        if bool(self.UserParams["verbose"]):
            time_start = dt.now()

        # main calculation: https://camb.readthedocs.io/en/latest/camb.html#camb.get_results
        results = camb.get_results(self.CAMBparams)

        outdict = { 'l': np.arange(self.max_l_use+1) }
        if self.UserParams['cls_to_output'] == 'all':
            cls_needed = ['clTT', 'clEE', 'clBB', 'clTE', 'clPP', 'clPT', 'clPE']
        else:
            cls_needed = self.UserParams['cls_to_output']
            if 'EB' in cls_needed:
                raise ValueError('camb does NOT output EB => it cant be in cls_to_output.')

        # https://camb.readthedocs.io/en/latest/results.html#camb.results.CAMBdata.get_total_cls
        if ('clTT' in cls_needed) or ('clEE' in cls_needed) or ('clBB' in cls_needed) or ('clTE' in cls_needed):
            # need to run things to get one/some/all of tt, ee, bb, te
            tt, ee, bb, te = results.get_total_cls(raw_cl=self.normalize_cls, CMB_unit=self.TT_units)[:self.max_l_use + 1].T
            # add noise
            if self.UserParams['noise_type'] is not None:
                _noise = self.get_noise()
                tt += _noise[0]
                ee += _noise[1]
                bb += _noise[1]
                te += _noise[1]
            # now add to outdict
            for key in ['clTT', 'clEE', 'clBB', 'clTE']:
                if key in cls_needed:
                    if key == 'clTT':
                        outdict[key] = tt
                    elif key == 'clEE':
                        outdict[key] = ee
                    elif key == 'clBB':
                        outdict[key] = bb
                    elif key == 'clTE':
                        outdict[key] = te
                    else:
                        raise ValueError('somethings wrong.')

        # https://camb.readthedocs.io/en/latest/results.html#camb.results.CAMBdata.get_lens_potential_cls
        if ('clPP' in cls_needed) or ('clPT' in cls_needed) or ('clPE' in cls_needed):
            pp, pt, pe = results.get_lens_potential_cls(raw_cl=self.normalize_cls)[:self.max_l_use + 1].T
            # now add to outdict
            for key in ['clPP', 'clPT', 'clPE']:
                if key in cls_needed:
                    if key == 'clPP':
                        outdict[key] = pp
                    elif key == 'clPT':
                        outdict[key] = pt
                    elif key == 'clPE':
                        outdict[key] = pe
                    else:
                        raise ValueError('somethings wrong.')

        # see if lmin is specified - if so, discard the ells < lmin
        if self.UserParams['lmin'] != 0:
            ell_inds = np.where(outdict['l'] >= self.UserParams['lmin'])[0]
            if bool(self.UserParams["verbose"]):
                print(f'## discarding {len(outdict["l"]) - len(ell_inds)} ells given lmin specification.')
            for key in outdict.keys():
                outdict[key] = outdict[key][ell_inds]

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
        """
        method for looping get_cls() over a range of values specified in the user_config.yaml file
        automatically saves results to self.results and parameters to self.result_parameters using
        unique identifiers

        Parameters
        ----------
        user_params : bool, default True
            if save_to_dict is not None, then user_params=True saves only the CAMBparams that differ
            from the base_config in the result_parameters dict

        Returns
        -------
        None
        """
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
        """
        method for saving power spectra

        Parameters
        ----------
        savedir : str, default CWD/outfiles
            directory into which results will be saved
        saveids : array_like, optional
            IDs to save. Can specify total number to save, specific runs to save, or (default) will save all
        randomids : bool, default False
            if saveids is an int, this will choose int random runs to save
        permission : {'w', 'a', 'w-', 'r'}, default 'w'
            permission settings for the file. The default behavior allows to overwrite an existing file.
        overwrite : bool, default False
            currently not being used

        Returns
        -------
        None
        """
        if not os.path.exists(os.path.join(os.getcwd(), "outfiles")):
            print("making local directory `./outfiles`")
            os.mkdir(os.path.join(os.getcwd(), "outfiles"))
        if saveids is not None:
            if type(saveids) == int:
                saveids = np.random.choice(range(len(self.loop_runids)), saveids, replace=False) if randomids else self.loop_runids[:saveids]
            elif isinstance(saveids, Iterable):
                if len(self.loop_runids) > 0:
                    saveids = [self.loop_runids[x] for x in saveids]
                else:
                    saveids = list(self.results.keys())
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
    """
    generate unique run ID including a random number whose length can be specified

    Parameters
    ----------
    random_digits : int, default 6
        number of random digits (including leading zeros) in the random number

    Returns
    -------
    str
        'runid_' plues year, month, day, hour, minute, second, microsecond, plus random number
    """
    _rint = np.random.randint(10**random_digits)
    return 'runid_'+dt.now().strftime('%y%m%d%H%M%S%f_')+str(_rint).zfill(random_digits)
