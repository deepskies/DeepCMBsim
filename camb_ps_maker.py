import os, camb, json, sys
import numpy as np
from datetime import datetime as dt

flags = ''.join([x for x in sys.argv if '-' in x])

# _parentdir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
_basedir = os.path.join(os.getcwd())

base_pars = camb.read_ini(os.path.join(_basedir, 'inifiles', 'planck_2018_1e4.ini'))

max_l_calc, max_l_use = 10100, 10000

base_pars.max_l, base_pars.max_l_tensor = max_l_calc, max_l_calc

with open(os.path.join(_basedir, 'sim_ranges.json'), 'r') as j:
    j_data = json.load(j)

# the following line determines the normalization of the C_ell's
# if this evaluates to False, which is the default, then the output
# power spectra include ell*(ell+1)/2/π for TT, EE, BB, TE and
# [ell*(ell+1)]**2/2/π for the PP, PT, and PE
# Tried to make it user-friendly so that you can write anything
# as long as it has both "raw" and "cl"
cls_raw = True if (('raw' in flags) and ('cl' in flags)) else False

# the following line determines the units of the TT and TE columns
# if this evaluates to 'muK', which is the default, then the output
# power spectra for TT carry units of µK^2 and for TE carry µK.
# This does not affect PT or PE, which are kept unitless.
units = None if (('dimensionless' in flags) and ('TT' in flags)) else 'muK'

ta = dt.now()

for rr in np.logspace(j_data['log10_r']['low'], j_data['log10_r']['high'], j_data['log10_r']['steps']):
    for aa in np.linspace(j_data['Alens']['low'], j_data['Alens']['high'], j_data['Alens']['steps']):
        base_pars.InitPower.r, base_pars.Alens = rr, aa

        results = camb.get_results(base_pars)

        tt, ee, bb, te = results.get_total_cls(raw_cl=cls_raw, CMB_unit=units).T
        pp, pt, pe = results.get_lens_potential_cls(raw_cl=cls_raw)[:max_l_use + 1].T
        lvals = range(max_l_use + 1)

        outarr = np.array([lvals, tt, ee, bb, te, pp, pt, pe]).T

        namestr = "lr" + f'{np.log10(rr):0.2f}' + "_A" + f'{aa:0.2f}' + "_d" + dt.strftime(dt.now(), '%y%m%d')
        if cls_raw:
            namestr += "_rawCl"

        np.save(os.path.join(_basedir, "psfiles", namestr), outarr)

tb = dt.now()

if 'v' in flags:
    print('from', dt.strftime(ta, '%H:%M:%S.%f %P'), 'to', dt.strftime(tb, '%H:%M:%S.%f %P'), 'or', end=" ")
    print(str((tb - ta).seconds) + '.' + str((tb - ta).microseconds), 'seconds total')
