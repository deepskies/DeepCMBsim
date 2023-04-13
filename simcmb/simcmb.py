import camb_ps_maker
import yam_in
import sys

infile = sys.argv[1]

Ydict = yam_in.Ydict(infile)

PS_maker = camb_ps_maker.PS_Maker(Ydict).loop_cls_rA()