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
            except Exception:  # necessary because of error with CustomSources.custom_source_ell_scales
                continue


class Ydict:
    def __init__(self, infile="example_config.yaml"):  # put .yaml into a settings folder
        with open(infile, "r") as f:
            self._all_params_dict = yaml.safe_load(f)

        self.CAMBparams = camb.CAMBparams()  # creates a base CAMBparams instance

        for x, y in self._all_params_dict['BASECAMBPARAMS'].items():  # get set first (potentially overwritten later)
            set_camb_attr(self.CAMBparams, x, y)

        for x, y in self._all_params_dict['USERPARAMS']['FORCAMB'].items():
            set_camb_attr(self.CAMBparams, x, y)

        for x, y in self._all_params_dict['USERPARAMS']['ITERABLES'].items():
            if isinstance(y, Iterable) and len(y) <= 3:
                self._all_params_dict['USERPARAMS']['ITERABLES'][x] = np.linspace(*y)
            else:
                if not isinstance(y, Iterable):
                    print(x, "is not iterable; are you sure it should be in ITERABLES?")
                self._all_params_dict['USERPARAMS']['ITERABLES'][x] = np.array(y)

        self.UserParams = self._all_params_dict['USERPARAMS']

        self.dict_iterables = self._all_params_dict['USERPARAMS']['ITERABLES']  # make this more easily accessible

        self.out_dict = {}

    def update_val(self, attr, new_val):
        attr_split = re.split("\.", attr)
        if (len(attr_split) == 1) and (hasattr(self.CAMBparams, attr)):
            setattr(self.CAMBparams, attr, new_val)
        elif (len(attr_split) == 2) and (hasattr(getattr(self.CAMBparams, attr_split[0]), attr_split[1])):
            setattr(getattr(self.CAMBparams, attr_split[0]), attr_split[1], new_val)
        elif attr in self._all_params_dict['USERPARAMS']:
            self._all_params_dict['USERPARAMS'][attr] = new_val
        else:
            print("not a valid attribute")

    def cpars_to_dict(self):
        _long_str = str(self.CAMBparams)
        _long_str_lines = re.split("\\n", _long_str)[1:]  # first entry is just 'class: <CAMBparams>'
        outer_dict = {}
        i = 0
        while i < len(_long_str_lines):
            _line = _long_str_lines[i]
            _line2 = re.split("=", _line)
            if len(_line2)>1:
                outer_dict[_line2[0].replace(" ", "")] = strconvert(_line2[1])
            elif len(_line)>1:
                _line = re.split(":", _line)[0]
                inner_dict = {}
                while i < len(_long_str_lines):
                    i += 1
                    if _long_str_lines[i][:3] == '   ':
                        _line3 = re.split("=", _long_str_lines[i])
                        inner_dict[_line3[0].replace(" ", "")] = strconvert(_line3[1]) if "None" not in _line3[1] else "~"
                    else:
                        break
                outer_dict[_line.replace(" ", "")] = inner_dict
            i += 1

        return outer_dict


    def save_cpars(self, filename):
        _long_str = str(self.CAMBparams)
        _long_str_lines = re.split("\\n", _long_str)[1:]  # first entry is just 'class: <CAMBparams>'
        with open(filename, "w") as f:
            for _line in _long_str_lines:
                _line = re.split("<", _line.replace("=", ":"))[0]
                if len(re.split(":", _line))==2 and 'None' in re.split(":", _line)[1]:
                    _line = _line.replace("None", " ")
                f.write(_line)
                f.write("\n")


def strconvert(x):
    try:
        return eval(x)
    except NameError:
        return x.replace(" ","")