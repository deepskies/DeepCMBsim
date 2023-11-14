import numpy as np
"""
noise module that implements Eq 8 of astro-ph/0111606 (Hu and Okamoto Astrophys.J. 574 (2002) 566-574)
also provides a maximum multipole given a beam size 
"""

arcmin_to_rad = np.pi / 180 / 60

def max_multipole(beamfwhm_arcmin, additional_factor=3):
    """
    this gives the maximum multipole advised to use given a certain beam size

    Parameters
    ----------
    beamfwhm_arcmin : float
        the full width at half maximum of the beam in the Gaussian approximation, in units of arcminutes
    additional_factor : float
        this is a fudge factor to prevent exponentially large noise levels,
        because there is not a sharp cutoff at a certain multipole
        for additional_factor = 3 (2) the noise has increased by roughly 7 (3) orders of magnitude from its low-l value

    Returns
    -------
    float
        the maximum multipole to which we should calculate given a beam whose size is beamfwhm_arcmin in arcminutes
    """
    return 180 * 60 * additional_factor / beamfwhm_arcmin

def detector_white_noise(noise_uKarcmin, beamfwhm_arcmin, lmax, TT=True, units_uK = True):
    """
    describes white (no angular scale) noise from a detector
    implements Eq 8 of astro-ph/0111606 (Hu and Okamoto Astrophys.J. 574 (2002) 566-574)
    related ideas first in Knox 1995 astro-ph/9504054 for T only
    and in Seljak and Zaldarriaga 1996 astro-ph/9609170 for polarization

    Parameters
    ----------
    noise_uKarcmin : float
        noise level in units of microKelvin*arcminutes
    beamfwhm_arcmin : float
        the full width at half maximum of the beam in the Gaussian approximation, in units of arcminutes
    lmax : int
        maximum multipole to which to calculate
    TT : bool, default True
        TT is True if the noise level is the noise of the temperature map, which is standard
    units_uK : bool, default True
        whether or not to return the units in microKelvin (default) or dimensionless

    Returns
    -------
    np.ndarray
        provides noise for the TT power spectrum or polarization power spectra in an array of length lmax+1
    """
    noise_uKarcmin = noise_uKarcmin if TT else noise_uKarcmin * np.sqrt(2)
    noise_uKarcmin = noise_uKarcmin if units_uK else noise_uKarcmin / 2.72548e6
    ells = np.arange(lmax+1)
    return (noise_uKarcmin * arcmin_to_rad) ** 2 * np.exp(ells * (ells + 1) * (beamfwhm_arcmin * arcmin_to_rad) ** 2 / (8 * np.log(2)))