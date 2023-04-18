import pymaster as nmt
from astropy.wcs import WCS
import numpy as np

class flatmap(object):
    def __init__(self, pixels, degrees, seed = -1, cl_dict = None):
        self.pixels, self.degrees = pixels, degrees  # pixels on each side, degrees on each side
        self.seed, self.cl_dict = seed, cl_dict # random seed if applicable, dictionary of CLs if you want to include them
        self.reso = self.degrees / self.pixels
        self.w = WCS(naxis=2)
        self.nx, self.lx, self.lx_rad, self.ny, self.ly, self.ly_rad = [int(self.pixels), int(self.degrees), self.degrees*np.pi/180]*2
        self.w.wcs.crpix = [self.nx // 2, self.ny // 2]  # Center pixel X, Y
        self.w.wcs.cdelt = np.array([-self.reso, self.reso])
        self.w.wcs.crval = [0, 0]  # Center coordinates RA, DEC at 0,0
        self.w.wcs.ctype = ["RA---AIR", "DEC--AIR"]  # Airy projection; can be adjusted. Previous used Azimuthal equal-area
        self.xlabel, self.ylabel = "longitude", "latitude"

    def flatmap(self, cl_array = None, cl_name = None, spin = [0], seed = None):
        seed = seed if seed is not None else self.seed
        if cl_array is not None:
            return nmt.synfast_flat(self.nx, self.ny, self.lx_rad, self.ly_rad, [cl_array], spin, seed = seed)
        elif (cl_name is not None) and (self.cl_dict is not None):
            return nmt.synfast_flat(self.nx, self.ny, self.lx_rad, self.ly_rad, [self.cl_dict[x] for x in cl_name], spin, seed=seed)
        else:
            print("need to specify Cls either with an array here or by passing a dictionary to the class and a name here")