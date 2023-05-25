import pymaster as nmt
import numpy as np
from matplotlib import pyplot as plt


class flatmap:
    def __init__(self, pixels, degrees, seed = -1, cl_dict = None):
        self.pixels, self.degrees = pixels, degrees  # pixels on each side, degrees on each side
        self.seed, self.cl_dict = seed, cl_dict # random seed if applicable, dictionary of CLs if you want to include them
        self.reso = self.degrees / self.pixels
        self.nx, self.lx, self.lx_rad, self.ny, self.ly, self.ly_rad = [int(self.pixels), int(self.degrees), self.degrees*np.pi/180]*2
        self.ticks, self.labels = np.linspace(0, self.pixels, self.degrees+1), np.arange(0, self.degrees+1, dtype=float)

    def _flatmap(self, cl_array, spin_array = [0], seed = None):
        seed = seed if seed is not None else self.seed
        return nmt.synfast_flat(self.nx, self.ny, self.lx_rad, self.ly_rad, np.array([x for x in cl_array]), spin_array, seed = seed)

    def flatmap(self, maps_out, seed = None):
        if self.cl_dict is not None:
            if len(maps_out)==1 and maps_out in ["T", "E", "B", "P"]:
                cl_arr, spin_arr = np.array([self.cl_dict["cl"+maps_out*2]]), [0]
            elif (len(maps_out)==2) and maps_out in ["TT", "EE", "BB", "TE", "PP", "PT", "PE"]:
                cl_arr, spin_arr = np.array([self.cl_dict["cl"+maps_out]]), [0]
            elif (len(maps_out)==4) and maps_out[:2]=='cl' and (maps_out[2:] in ["TT", "EE", "BB", "TE", "PP", "PT", "PE"]):
                cl_arr, spin_arr = np.array([self.cl_dict[maps_out]]), [0]
            elif maps_out=='TEB': # the cl_arr has to be in a specific order; we enforce by hand that the two other cross spectra TB and EB are zero
                cl_arr, spin_arr = np.array([self.cl_dict["clTT"], self.cl_dict["clTE"], np.zeros_like(self.cl_dict["clTE"]), self.cl_dict["clEE"], np.zeros_like(self.cl_dict["clTE"]), self.cl_dict["clBB"]]), [0, 0, 0]
            elif maps_out=='TQU': # the cl_arr has to be in a specific order; we enforce by hand that the two other cross spectra TB and EB are zero
                cl_arr, spin_arr = np.array([self.cl_dict["clTT"], self.cl_dict["clTE"], np.zeros_like(self.cl_dict["clTE"]), self.cl_dict["clEE"], np.zeros_like(self.cl_dict["clTE"]), self.cl_dict["clBB"]]), [0, 2]
            else:
                print("not a valid map specification")
                return None
            return self._flatmap(cl_arr, spin_arr, seed=seed)
        else:
            print("if you don't want to restrict to a `cl_dict` dictionary, use `self._flatmap` instead")


def power_spectrum_plot(cldict = None, cls = None):
    if cldict is not None:
        if len(cls)==1 and cls in ["T", "E", "B", "P"]:
            cl_arr, ells = np.array([cldict["cl"+cls*2]]), cldict['l']
        elif (len(cls)==2) and cls in ["TT", "EE", "BB", "TE", "PP", "PT", "PE"]:
            cl_arr, ells = np.array([cldict["cl"+cls]]), cldict['l']
        elif (len(cls)==4) and cls[:2]=='cl' and (cls[2:] in ["TT", "EE", "BB", "TE", "PP", "PT", "PE"]):
            cl_arr, ells = np.array([cldict[cls]]), cldict['l']
        for cvals in cl_arr:
            plt.plot(ells, cvals)
    else:
        if len(cls) == 2:
            print("assuming that ells and Cls are provided")
            plt.plot(cls[0], cls[1])
        else:
            print("assuming that only Cls are provided")
            plt.plot(range(len(cls)), cls)