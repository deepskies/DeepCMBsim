"""
optional module for plotting flat-sky map realizations of power spectra
"""

import pymaster as nmt
import numpy as np


class flatmap:
    """
    map-making class
    Attributes
    ----------
    pixels : int
        number of pixels on each side of a flat map
    degrees : float
        number of degrees on each side of a flat map
    namaster_seed : int, default -1
        random namaster_seed to provide `namaster` (potential bug in `namaster`, does not always have the
        behavior that is described at https://namaster.readthedocs.io/en/latest/pymaster.html
    cl_dict : dict, optional
        dictionary of power spectra. Necessary to use the primary method of this class.
        If not specified, can still use the hidden mapping function.
    reso : float
        resolution of the maps, given by degrees/pixels
    nx, ny : int
        number of pixels on the x and y axes (fixed to be the same for now, but generalizable)
    lx, ly:
        number of degrees on x and y axes (fixed to be the same for now, but generalizable)
    lx_rad, ly_rad : float
        number of radians on x and y axes (fixed to be the same for now, but generalizable)
    ticks : np.ndarray
        tickmark locations for annotating maps
    lables : np.ndarray
        tickmark labels for annotating maps

    Methods
    -------
    flatmap
        function for plotting a realization of a flat-sky map with the given number of pixels
        representing the given size of sky patch based on the provided power spectra
    """
    def __init__(self, pixels, degrees, namaster_seed = -1, cl_dict = None):
        """
        initialize a map-making class
        Parameters
        ----------
        pixels : int
            number of pixels on each side of a flat map
        degrees : float
            number of degrees on each side of a flat map
        namaster_seed : int, default -1
            random namaster_seed to provide `namaster` (potential bug in `namaster`, does not always have the
            behavior that is described at https://namaster.readthedocs.io/en/latest/pymaster.html
        cl_dict : dict, optional
            dictionary of power spectra. Necessary to use the primary method of this class.
            If not specified, can still use the hidden mapping function.
        """
        self.pixels, self.degrees = pixels, degrees  # pixels on each side, degrees on each side
        self.namaster_seed, self.cl_dict = namaster_seed, cl_dict # random namaster_seed if applicable, dictionary of CLs if you want to include them
        self.reso = self.degrees / self.pixels
        self.nx, self.lx, self.lx_rad, self.ny, self.ly, self.ly_rad = [int(self.pixels), int(self.degrees), self.degrees*np.pi/180]*2
        self.ticks, self.labels = np.linspace(0, self.pixels, self.degrees+1), np.arange(0, self.degrees+1, dtype=float)

    def _flatmap(self, cl_array, spin_array = [0], namaster_seed = None):
        namaster_seed = namaster_seed if namaster_seed is not None else self.namaster_seed
        return nmt.synfast_flat(self.nx, self.ny, self.lx_rad, self.ly_rad, np.array([x for x in cl_array]), spin_array, namaster_seed = namaster_seed)

    def flatmap(self, maps_out, namaster_seed = None):
        """

        Parameters
        ----------
        maps_out : ["T", "E", "B", "P", "TT", "EE", "BB", "TE", "PP", "PT", "PE", "clTT", "clEE", "clBB", "clTE", "clPP", "clPT", "clPE", "TEB", "TQU"]
            map(s) that you would like to produce
        namaster_seed : int, optional
            random namaster_seed to provide `namaster` (potential bug in `namaster`, does not always have the
            behavior that is described at https://namaster.readthedocs.io/en/latest/pymaster.html
        Returns
        -------
        np.ndarray
            returns entries of flat maps in a (Nm, Np, Np) array, where Nm is the number of maps
            (equal to 1 if maps_out is in ["T", "E", "B", "P", "TT", "EE", "BB", "TE", "PP", "PT",
            "PE", "clTT", "clEE", "clBB", "clTE", "clPP", "clPT", "clPE"] and equal to 3 if
            maps_out is "TEB" or "TQU"]) and Np is the number of pixels
        """
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
            return self._flatmap(cl_arr, spin_arr, namaster_seed=namaster_seed)
        else:
            print("if you don't want to restrict to a `cl_dict` dictionary, use `self._flatmap` instead")
