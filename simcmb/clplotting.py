import pymaster as nmt
import numpy as np

class flatmap(object):
    def __init__(self, pixels, degrees, seed = -1, cl_dict = None):
        self.pixels, self.degrees = pixels, degrees  # pixels on each side, degrees on each side
        self.seed, self.cl_dict = seed, cl_dict # random seed if applicable, dictionary of CLs if you want to include them
        self.reso = self.degrees / self.pixels
        self.nx, self.lx, self.lx_rad, self.ny, self.ly, self.ly_rad = [int(self.pixels), int(self.degrees), self.degrees*np.pi/180]*2

    def _flatmap(self, cl_array, spin_array = [0], seed = None):
        seed = seed if seed is not None else self.seed
        return nmt.synfast_flat(self.nx, self.ny, self.lx_rad, self.ly_rad, np.array([x for x in cl_array]), spin_array, seed = seed)

    def flatmap(self, maps_out, seed = None):
        if self.cl_dict is not None:
            if len(maps_out)==1:
                cl_arr, spin_arr = np.array([self.cl_dict["cl"+maps_out*2]]), [0]
            elif maps_out=='TEB': # the cl_arr has to be in a specific order; we enforce by hand that the two other cross spectra TB and EB are zero
                cl_arr, spin_arr = np.array([self.cl_dict["clTT"], self.cl_dict["clTE"], np.zeros_like(self.cl_dict["clTE"]), self.cl_dict["clEE"], np.zeros_like(self.cl_dict["clTE"]), self.cl_dict["clBB"]]), [0, 0, 0]
            elif maps_out=='TQU': # the cl_arr has to be in a specific order; we enforce by hand that the two other cross spectra TB and EB are zero
                cl_arr, spin_arr = np.array([self.cl_dict["clTT"], self.cl_dict["clTE"], np.zeros_like(self.cl_dict["clTE"]), self.cl_dict["clEE"], np.zeros_like(self.cl_dict["clTE"]), self.cl_dict["clBB"]]), [0, 2]
            return self._flatmap(cl_arr, spin_arr, seed=seed)
        else:
            print("if you don't want to restrict to a `cl_dict` dictionary, use `self._flatmap` instead")