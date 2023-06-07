import camb_ps_maker
import yam_io
import sys

infile = sys.argv[1]

Ydict = yam_io.config_obj(infile)

PS_maker = camb_ps_maker.CAMBPowerSpectrum(Ydict)