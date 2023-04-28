import camb
import yaml
import numpy as np
import h5py
from collections.abc import Iterable

class Ydict(object):
    def __init__(self, infile = "example_config.yaml"):  # put .yaml into a settings folder
        with open(infile, "r") as f:
            self.myyam = yaml.safe_load(f)

        self.pars = camb.CAMBparams() # creates a base CAMBparams instance

        self.basepars = self.myyam['BASECAMBPARAMS']
        for x, y in self.basepars.items():
            try:
                setattr(self.pars, x, y)
            except AttributeError:  # this is possible because some CAMBparams attributes have depth 2
                for a, b in y.items():
                    try:
                        setattr( getattr(self.pars, x), a, b)
                    except AttributeError:
                        print("this attribute doesn't exist")
                        continue

        self.myCAMBpars = self.myyam['USERPARAMS']['FORCAMB']
        for x, y in self.myCAMBpars.items():
            try:
                setattr(self.pars, x, y)
            except Exception:
                for a, b in y.items():
                    try:
                        setattr( getattr(self.pars, x), a, b)
                    except Exception:
                        continue

        rv = self.myyam['USERPARAMS']['rvals'] #more descriptive name for r
        if (isinstance(rv, Iterable) and (len(rv)<=3)):
            rv[:2] = np.log10(rv[:2]) if (rv[0]>0) else rv[:2]
            self.rs = np.logspace(*rv)
        else:
            self.rs = np.array(rv)

        Av = self.myyam['USERPARAMS']['Avals'] #more descriptive name for A
        if (isinstance(Av, Iterable) and (len(Av)<=3)):
            self.As = np.linspace(*Av)
        else:
            self.As = np.array(Av)

        self.user_params = self.myyam['USERPARAMS']

    def update_vals(self, attr, new_val, incamb = False):
        if incamb:
            if attr=="r":
                setattr(self.Ydict.pars.InitPower, attr, new_val)
            else:
                setattr(self.Ydict.pars, attr, new_val)
        else:
            self.Ydictu[attr] = new_val
        if "fwhm" in attr:
            self.max_l_use = int(min(self.Ydictu['max_l_use'], noise.max_multipole(self.Ydictu['beam_fwhm'])))
            self.Ydict.pars.max_l, self.Ydict.pars.max_l_tensor = int(self.max_l_use + self.Ydictu['extra_l']), int(self.max_l_use + self.Ydictu['extra_l'])



def savecls(all_sims, out_name, sims_to_save_start=None, sims_to_save_end=None, permission = 'r+', overwrite=False):
    if (sims_to_save_end is None) and (sims_to_save_start is None):
        sims_to_save_start = 0
        sims_to_save_end = len(all_sims)
    with h5py.File(out_name + '.h5', permission) as f:
        for i in range(sims_to_save_start, sims_to_save_end):
            out_dict = all_sims[i]
            for k, v in out_dict.items():
                try:
                    f.create_dataset(f"r{out_dict['r']}/Alens{out_dict['Alens']}/{k}", data = v)
                except ValueError:
                    if overwrite:
                        f[f"r{out_dict['r']}/Alens{out_dict['Alens']}/{k}"] = v
                    else:
                        print(f"skipping because r{out_dict['r']}/Alens{out_dict['Alens']}/{k} already exists and overwrite set to False")