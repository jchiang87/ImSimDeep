"""
Simulation and analysis tools for use with products from imSim and
the LSST Stack.
"""
from __future__ import absolute_import, print_function
import os
import lsst.sims.photUtils as photUtils
import lsst.utils as lsstUtils

__all__ = ['ObservedSEDs']

class SED(object):
    """
    Thin wrapper over lsst.photUtils.Sed class to provide more natural
    syntax for accessing the magnitudes in a requested band.

    Attributes
    ----------
    sed_instance : lsst.photUtils.Sed instance
        The underlying Sed instance.

    bandpasses : dict
        A dictionary of lsst.photUtils.Bandpass instances.
    """
    def __init__(self, sed_instance, magnorm, iA_v, iR_v, gA_v, gR_v,
                 bandpasses):
        """
        Parameters
        ----------
        sed_instance : lsst.photUtils.Sed instance
            The underlying Sed instance.

        magnorm : float

        iA_v : float

        iR_v : float

        gA_v : float

        gR_v : float

        bandpasses : dict
            A dictionary of lsst.photUtils.Bandpass instances.
        """
        self.sed_instance = sed_instance
        self.magnorm = magnorm
        self.iA_v = iA_v
        self.iR_v = iR_v
        self.gA_v = gA_v
        self.gR_v = gR_v
        self.bandpasses = bandpasses

    def calcMag(self, band):
        """
        Compute the magnitude in the specified band.  This effectively
        overloads lsst.photUtils.Sed.calcMag.

        Parameters
        ----------
        band : str
            The desired band.  For LSST, it would be 'u', 'g', 'r', 'i',
            'z', or 'y'.

        Returns
        -------
        float
            The magnitude in the desired band.
        """
        return self.sed_instance.calcMag(self.bandpasses[band])

    def __getattr__(self, attrname):
        "Delegate access to all other attributes to underlying Sed instance."
        return getattr(self.sed_instance, attrname)

class ObservedSEDs(object):
    """
    Class to compute observer frame SEDs.  SED normalization, internal
    extinction, redshift, and Galactic extinction are applied given the
    parameters in an instance catalog object line.

    Attributes
    ----------
    bps : dict
        Dictionary of LSST bandpasses.
    control_bandpass : lsst.photUtils.Bandpass instance
        The "imsim bandpass" which is used to set magnorm of an object's
        spectrum.
    """
    def __init__(self):
        """
        Set up the LSST bandpasses.
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

    def get_SED(self, object_line):
        """
        Compute the object's SED in the observer frame.

        Parameters
        ----------
        object_line : str
            The object line entry from a phosim instance catalog.

        Returns
        -------
        SED instance
            SED is a thin wrapper class around the lsst.photUtils.Sed
            class.  It can be used to compute the object's apparent
            magnitude in a given band.
        """
        obj_tokens = object_line.strip().split()

        # Read the SED for the object.
        sed_dir = lsstUtils.getPackageDir('sims_sed_library')
        spectrum = photUtils.Sed()
        spectrum.readSED_flambda(os.path.join(sed_dir, obj_tokens[5]))

        # Normalize the spectrum to magnorm.
        magnorm = float(obj_tokens[4])
        fnorm = spectrum.calcFluxNorm(magnorm, self.control_bandpass)
        spectrum.multiplyFluxNorm(fnorm)

        if len(obj_tokens) > 17:
            # Object has internal extinction.
            iA_v = float(obj_tokens[-5])
            iR_v = float(obj_tokens[-4])
            gA_v = float(obj_tokens[-2])
            gR_v = float(obj_tokens[-1])
        else:
            # Object only should only have Galactic extinction applied.
            iA_v = 0.
            iR_v = 0.
            gA_v = float(obj_tokens[14])
            gR_v = float(obj_tokens[15])

        if iA_v != 0 or iR_v != 0:
            # Apply internal dust extinction.
            a_int, b_int = spectrum.setupCCMab()
            spectrum.addCCMDust(a_int, b_int, A_v=iA_v, R_v=iR_v)

        redshift = float(obj_tokens[6])
        if redshift > 0:
            spectrum.redshiftSED(redshift, dimming=True)

        # Apply Galactic extinction.
        a_int, b_int = spectrum.setupCCMab()
        spectrum.addCCMDust(a_int, b_int, A_v=gA_v, R_v=gR_v)

        return SED(spectrum, magnorm, iA_v, iR_v, gA_v, gR_v, self.bps)
