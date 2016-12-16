"""
Simulation and analysis tools for use with products from imSim and
the LSST Stack.
"""
from __future__ import absolute_import, print_function
import os
import copy
import lsst.sims.photUtils as photUtils
import lsst.utils as lsstUtils

__all__ = ['ApparentMagnitude']

class ApparentMagnitude(object):
    """
    Class to compute apparent magnitudes for a given rest frame SED.
    The SED normalization, internal extinction, redshift, and Galactic
    extinction are applied given the parameters in an instance catalog
    object line to produce the apparent magnitude in the desired band.

    Attributes
    ----------
    bps : dict
        Dictionary of LSST bandpasses.
    control_bandpass : lsst.photUtils.Bandpass instance
        The "imsim bandpass" which is used to set magnorm of an object's
        spectrum.
    sed_unnormed : lsst.photUtils.Sed object
        The un-normalized SED.
    max_mag : float
        Sentinal value for underflows of Sed.calcMag
    """
    def __init__(self, sed_name, max_mag=1000.):
        """
        Set up the LSST bandpasses and un-normalized SED.
        """
        self.bps = dict()
        throughput_dir = lsstUtils.getPackageDir('throughputs')
        for band in 'ugrizy':
            self.bps[band] = photUtils.Bandpass()
            self.bps[band].readThroughput(os.path.join(throughput_dir,
                                                       'baseline',
                                                       'total_%s.dat' % band))

        self.control_bandpass = photUtils.Bandpass()
        self.control_bandpass.imsimBandpass()

        sed_dir = lsstUtils.getPackageDir('sims_sed_library')
        self.sed_unnormed = photUtils.Sed()
        self.sed_unnormed.readSED_flambda(os.path.join(sed_dir, sed_name))
        self.max_mag = max_mag

    def _sed_copy(self):
        """
        Return a copy of the unnormalized SED.
        """
        return copy.deepcopy(self.sed_unnormed)

    def __call__(self, pars, band):
        """
        Compute the object's SED in the observer frame.

        Parameters
        ----------
        pars : pandas.Series
            This contains the parameters (ra, dec, magnorm, etc.) from
            an "object" entry in an instance catalog.
        band : str
            The LSST band ('u', 'g', 'r', 'i', 'z', or 'y') to use for
            the apparent magnitude calculation.

        Returns
        -------
        float
            The apparent magnitude in the desired band.
        """
        # Normalize the spectrum to magnorm.
        spectrum = self._sed_copy()
        fnorm = spectrum.calcFluxNorm(pars.magNorm, self.control_bandpass)
        spectrum.multiplyFluxNorm(fnorm)

        iA_v, iR_v = pars.internalAv, pars.internalRv
        gA_v, gR_v = pars.galacticAv, pars.galacticRv

        if iA_v != 0 or iR_v != 0:
            # Apply internal dust extinction.
            a_int, b_int = spectrum.setupCCMab()
            spectrum.addCCMDust(a_int, b_int, A_v=iA_v, R_v=iR_v)

        if pars.redshift > 0:
            spectrum.redshiftSED(pars.redshift, dimming=True)

        # Apply Galactic extinction.
        if gA_v != 0 or gR_v != 0:
            a_int, b_int = spectrum.setupCCMab()
            spectrum.addCCMDust(a_int, b_int, A_v=gA_v, R_v=gR_v)


        try:
            mag = spectrum.calcMag(self.bps[band])
        except Exception as eObj:
            if str(eObj).startswith("This SED has no flux"):
                mag = self.max_mag
            else:
                raise eObj

        return mag
