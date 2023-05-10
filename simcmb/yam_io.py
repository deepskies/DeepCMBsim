import camb
import yaml
import numpy as np
import re
from collections.abc import Iterable
import os


def _set_camb_attr(cambparams_instance, x, y):
    try:
        setattr(cambparams_instance, x, y)
    except TypeError:  # this is possible because some CAMBparams attributes have depth 2
        for a, b in y.items():
            try:
                setattr(getattr(cambparams_instance, x), a, b)
            except Exception:  # necessary because of error with CustomSources.custom_source_ell_scales
                continue


def _quick_yaml_load(infile=None):
    if infile is None:
        return {}
    else:
        with open(infile, "r") as f:
            return yaml.safe_load(f)


class Yobj:
    def __init__(
            self,
            user_config=os.path.join(os.path.dirname(__file__), "settings", "user_config.yaml"),
            base_config=os.path.join(os.path.dirname(__file__), "settings", "base_config.yaml")
    ):

        self._all_params_dict = {
            'USERPARAMS': _quick_yaml_load(user_config),
            'BASECAMBPARAMS': _quick_yaml_load(base_config)
        }

        self.CAMBparams = camb.CAMBparams()  # creates a base CAMBparams instance

        for x, y in self._all_params_dict['BASECAMBPARAMS'].items():  # get set first (potentially overwritten later)
            _set_camb_attr(self.CAMBparams, x, y)

        if len(self._all_params_dict['USERPARAMS']) > 0:
            try:
                for x, y in self._all_params_dict['USERPARAMS']['FORCAMB'].items():
                    _set_camb_attr(self.CAMBparams, x, y)
            except KeyError:
                print("nothing overwriting `base_config.yaml`")
            try:
                for x, y in self._all_params_dict['USERPARAMS']['ITERABLES'].items():
                    if isinstance(y, Iterable) and len(y) <= 3 and type(y[-1]) == int:
                        self._all_params_dict['USERPARAMS']['ITERABLES'][x] = np.linspace(*y)
                    else:
                        if not isinstance(y, Iterable):
                            print(x, "is not iterable; are you sure it should be in ITERABLES?")
                        self._all_params_dict['USERPARAMS']['ITERABLES'][x] = np.array(y)
            except KeyError:
                print("no iterables specified")
                self._all_params_dict['USERPARAMS']['ITERABLES'] = {}

        self.UserParams = self._all_params_dict['USERPARAMS']

        self.dict_iterables = self._all_params_dict['USERPARAMS']['ITERABLES']  # make this more easily accessible

        self.out_dict = {}

    def update_val(self, attr, new_val):
        attr_split = re.split("\.", attr)
        if (len(attr_split) == 1) and (hasattr(self.CAMBparams, attr)):
            setattr(self.CAMBparams, attr, new_val)
        elif (len(attr_split) == 2) and (hasattr(getattr(self.CAMBparams, attr_split[0]), attr_split[1])):
            setattr(getattr(self.CAMBparams, attr_split[0]), attr_split[1], new_val)
        # elif attr in self._all_params_dict['USERPARAMS']:
        self._all_params_dict['USERPARAMS'][attr] = new_val
        # else:
        #     print("not a valid attribute")

    def cpars_to_dict(self, user_params=True):
        cpd = _cpars_to_dict(self.CAMBparams)
        if user_params:
            self._all_params_dict['USERPARAMS']['FORCAMB'] = _nested_dict_diff(cpd, self._all_params_dict["BASECAMBPARAMS"])
            return self._all_params_dict['USERPARAMS']
        else:
            return cpd


def _strconvert(x):
    try:
        return eval(x)
    except NameError:
        return x.replace(" ", "")


def _cpars_to_dict(cambparams_instance):
    _long_str = str(cambparams_instance)
    _long_str_lines = re.split("\\n", _long_str)[1:]  # first entry is just 'class: <CAMBparams>'
    outer_dict = {}
    i = 0
    while i < len(_long_str_lines):
        _line = _long_str_lines[i]
        _line2 = re.split("=", _line)
        if len(_line2) > 1:
            outer_dict[_line2[0].replace(" ", "")] = _strconvert(_line2[1])
        elif len(_line) > 1:
            _line = re.split(":", _line)[0]
            inner_dict = {}
            while i < len(_long_str_lines):
                i += 1
                if _long_str_lines[i][:3] == '   ':
                    _line3 = re.split("=", _long_str_lines[i])
                    inner_dict[_line3[0].replace(" ", "")] = _strconvert(_line3[1]) if "None" not in _line3[1] else "~"
                else:
                    i -= 1
                    break
            outer_dict[_line.replace(" ", "")] = inner_dict
        i += 1
    return outer_dict


def _nested_dict_diff(d1, d2):
    diff_dict = {}
    for k, v in d1.items():
        if type(v) == dict:
            for k_inner, v_inner in v.items():
                if d2[k][k_inner] != v_inner:
                    diff_dict[k] = {k_inner: v_inner}
        else:
            if d2[k] != v:
                diff_dict[k] = v
    return diff_dict
