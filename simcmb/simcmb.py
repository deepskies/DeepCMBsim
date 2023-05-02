import camb_ps_maker
import yam_io
import sys

infile = sys.argv[1]

Ydict = yam_io.Ydict(infile)

PS_maker = camb_ps_maker.PS_Maker(Ydict)