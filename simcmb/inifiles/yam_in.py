with open("inifiles/planck_2018_1e4.yaml", "r") as f:
    myyam = yaml.safe_load(f)

basepars = myyam['BASECAMBPARAMS']
for x, y in basepars.items():
    try:
        setattr(pars, x, y)
    except Exception:
        for a, b in y.items():
            try:
                setattr( getattr(pars, x), a, b)
            except Exception:
                continue

mypars = myyam['USERPARAMS']['FORCAMB']
for x, y in mypars.items():
    try:
        setattr(pars, x, y)
    except Exception:
        for a, b in y.items():
            try:
                setattr( getattr(pars, x), a, b)
            except Exception:
                continue