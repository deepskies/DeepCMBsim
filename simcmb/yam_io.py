import camb
import yaml
import numpy as np
import re
from collections.abc import Iterable


def set_camb_attr(cambparams_instance, x, y):
    try:
        setattr(cambparams_instance, x, y)
    except TypeError:  # this is possible because some CAMBparams attributes have depth 2
        for a, b in y.items():
            try:
                setattr(getattr(cambparams_instance, x), a, b)
            except Exception:  # this is only necessary because of some weird type error with CustomSources.custom_source_ell_scales
                continue


class Ydict(object):
    def __init__(self, infile="example_config.yaml"):  # put .yaml into a settings folder
        with open(infile, "r") as f:
            self.all_params_dict = yaml.safe_load(f)

        self.CAMBparams = camb.CAMBparams()  # creates a base CAMBparams instance

        for x, y in self.all_params_dict['BASECAMBPARAMS'].items():  # these get set first and are potentially overwritten later
            set_camb_attr(self.CAMBparams, x, y)

        for x, y in self.all_params_dict['USERPARAMS']['FORCAMB'].items():
            set_camb_attr(self.CAMBparams, x, y)

        for x, y in self.all_params_dict['USERPARAMS']['ITERABLES'].items():
            if (isinstance(y, Iterable) and (len(y) <= 3)):
                self.all_params_dict['USERPARAMS']['ITERABLES'][x] = np.linspace(*y)
            else:
                if not isinstance(y, Iterable):
                    print(x, "is not iterable; are you sure it should be in ITERABLES?")
                self.all_params_dict['USERPARAMS']['ITERABLES'][x] = np.array(y)

        self.dict_iterables = self.all_params_dict['USERPARAMS']['ITERABLES']  # want to make this more easily accessible, for loops

    def update_val(self, attr, new_val):
        attr_split = re.split("\.", attr)
        if (len(attr_split) == 1) and (hasattr(self.CAMBparams, attr)):
            setattr(self.CAMBparams, attr, new_val)
        elif (len(attr_split) == 2) and (hasattr( getattr(self.CAMBparams, attr_split[0]), attr_split[1])):
            setattr(getattr(self.CAMBparams, attr_split[0]), attr_split[1], new_val)
        elif attr in self.all_params_dict['USERPARAMS']:
            self.all_params_dict['USERPARAMS'][attr] = new_val
        else:
            print("not a valid attribute")
