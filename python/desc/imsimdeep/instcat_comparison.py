"""
Using apparent magnitudes, coordinates, etc., based on the instance
catalog used to produce the simulated images, create a pandas data
frame combining those info with the positionally associated sources
from the forced source catalogs produced by the Stack.

The columns are
objectId, ra_meas, dec_meas, mag_meas, magerr_meas, ra_true, dec_true,
mag_true, offset, galSimType
"""
from __future__ import absolute_import, print_function, division
import os
import numpy as np
import astropy.io.fits as fits
import pandas as pd
import sklearn.neighbors
import lsst.daf.persistence as dp
import desc.imsim

__all__ = ['instcat_comparison']

logger = desc.imsim.get_logger('INFO')

class MagFromAdu(object):
    def __init__(self, fluxmag0):
        self.fluxmag0 = fluxmag0
    def __call__(self, counts):
        return -2.5*np.log10(counts/self.fluxmag0)
    def error(self, counts, counts_err):
        return 2.5*np.log10(self.fluxmag0)*counts_err/counts

def instcat_comparison(app_mag_file, repo, visit, raft, sensor,
                       flux_col='base_PsfFlux_flux', tract='0'):
    """
    Do a positional association between the forced source objects and
    the instance catalog coordintaes.  Return a data frame with
    coordinates and magnitudes to check astrometric and photometric
    results from the Stack.

    Parameters
    ----------
    app_mag_file : str
        Filename of a pickled pandas data frame of in-band apparent
        magnitudes computed from an instance catalog via
        compute_apparent_mags.py.
    repo : str
        Output repo containing the Stack forced source catalogs.
    visit : int
        Visit number to use.
    raft : str
        Raft id to use, e.g., '2,2'.
    sensor : str
        Sensor id to use, e.g., '1,1'
    flux_col : str, optional
        Name of the column in the forced source catalog to use for the
        flux measurment.  Default: 'base_PsfFlux_flux'
    tract : str, optional
        Tract to use.  Default: '0'

    Returns
    -------
    tuple(pandas.DataFrame, str) :
        The first item is the data frame containing the measured
        positions and magnitudes and associated true positions and
        magnitudes, positional offsets (in arcsec), and input object
        galSimType. The second item is the visit name derived from the
        visit number and retrieved band, e.g., 'v230-fr'.

    """
    butler = dp.Butler(repo)
    dataId = dict(visit=visit, raft=raft, sensor=sensor)
    calexp = butler.get('calexp', dataId=dataId)
    mag_from_adu = MagFromAdu(calexp.getCalib().getFluxMag0()[0])
    band = calexp.getFilter().getName()
    visit_name = 'v%i-f%s' % (visit, band)

    # The data butler doesn't work properly for forced_src catalogs
    # (It complains with "OperationalError: no such column: tract"
    # when executing butler.get('forced_src').), so access the data
    # using astropy.io.fits.
    forced_catalog = os.path.join(repo, 'forced', tract, visit_name,
                                  'R%s' % raft[::2], 'S%s.fits' % sensor[::2])
    forced = fits.open(forced_catalog)

    instcat = pd.read_pickle(app_mag_file)

    # Use a KD tree to find nearest instance catalog object for each
    # detected and measured object.
    tree = sklearn.neighbors.KDTree(
        np.array(((instcat['raICRS'].values*np.pi/180.,
                   instcat['decICRS'].values*np.pi/180.))).transpose())
    candidates = np.array(((forced[1].data['coord_ra'],
                            forced[1].data['coord_dec']))).transpose()
    offset, index = tree.query(candidates, k=1)

    # Since the KDTree is Euclidean for 2 dimensions by default,
    # account for the effect of a non-zero latitude and convert from
    # radians to arcsec.
    offset = np.array([x[0] for x in offset])
    offset *= np.cos(forced[1].data['coord_dec'])*180./np.pi*3600.

    # Build the output data frame.
    index = (tuple([x[0] for x in index]),)
    galSimType = instcat['galSimType'].values[index]
    true_mags = instcat[band].values[index]

    flux = np.array(forced[1].data[flux_col].tolist())
    flux_err = np.array(forced[1].data[flux_col + 'Sigma'].tolist())
    df = pd.DataFrame(dict(objectId=forced[1].data['objectId'].tolist(),
                           coord_ra=(forced[1].data['coord_ra']*180./np.pi).tolist(),
                           coord_dec=(forced[1].data['coord_dec']*180./np.pi).tolist(),
                           magnitude=mag_from_adu(flux),
                           magnitude_error=mag_from_adu.error(flux, flux_err),
                           coord_ra_true=instcat['raICRS'].values[index],
                           coord_dec_true=instcat['decICRS'].values[index],
                           magnitude_true=instcat[band].values[index],
                           offset=offset,
                           galSimType=instcat['galSimType'].values[index]))
    return df, visit_name

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    plt.ion()

    repo = 'output'
    visit = 230
    raft = '2,2'
    sensor = '1,1'
    app_mag_file = 'instcats/%(visit)i/app_mags_%(visit)i.pkl' % locals()

    df, visit_name = instcat_comparison(app_mag_file, repo, visit, raft, sensor)
    pointSource = df.query('galSimType=="pointSource"')

    plt.figure()
    plt.hist(pointSource['offset'], bins=50, histtype='step')
    plt.xlabel('offset from true position (arcsec)')
    plt.ylabel('entries / bin')
    plt.title(visit_name)

    plt.figure()
    delta_mag = pointSource['magnitude'] - pointSource['magnitude_true']
    plt.errorbar(pointSource['magnitude_true'], delta_mag,
                 fmt='.')
    plt.xlabel('true magnitude')
    plt.ylabel('measured - true magnitude')
    plt.title(visit_name)

    plt.figure()
    plt.errorbar(pointSource['coord_ra'], pointSource['coord_dec'],
                 fmt='.', label='measured position')
    plt.scatter(pointSource['coord_ra_true'], pointSource['coord_dec_true'],
                s=60, label='true position', facecolors='none',
                edgecolors='r')
    plt.xlabel('RA (deg)')
    plt.ylabel('Dec (deg)')
    plt.legend()
    plt.title(visit_name)
