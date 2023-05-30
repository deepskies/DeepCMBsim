import numpy as np
"""
noise module that implements Eq 8 of astro-ph/0111606 (Hu and Okamoto Astrophys.J. 574 (2002) 566-574)
also provides a maximum multipole given a beam size 
"""

def max_multipole(fwhm_arcmin, additional_factor=3):
    """
    this gives the maximum multipole advised to use given a certain beam size
    the additional_factor is set because there is not a sharp cutoff at a certain multipole
    for additional_factor = 3 (2) the noise has increased by roughly 7 (3) orders of magnitude from its low-l value
    """
    return 180*60*additional_factor/fwhm_arcmin

def white_noise(noise_uK_arcmin, fwhm_arcmin, lmax, TT=True, units_uK = True):
    """
    implements Eq 8 of astro-ph/0111606 (Hu and Okamoto Astrophys.J. 574 (2002) 566-574)
    first in Seljak and Zaldarriaga 1996 astro-ph/9609170, which is based on Knox 1995 astro-ph/9504054
    TT is True if the noise level is the noise of the temperature map, which is standard
    """
    noise_uK_arcmin = noise_uK_arcmin if TT else noise_uK_arcmin*np.sqrt(2)
    noise_uK_arcmin = noise_uK_arcmin if units_uK else noise_uK_arcmin/2.72548e6
    ells = np.arange(lmax+1)
    arcmin_to_rad = np.pi / 180 / 60
    return (noise_uK_arcmin * arcmin_to_rad)**2 * np.exp( ells * (ells + 1) * (fwhm_arcmin * arcmin_to_rad)**2 / (8 * np.log(2)) )