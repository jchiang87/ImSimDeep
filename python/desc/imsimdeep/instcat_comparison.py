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
import matplotlib.pyplot as plt
import sklearn.neighbors
import lsst.daf.persistence as dp
import desc.imsim

__all__ = ['instcat_comparison', 'plot_instcat_comparison']

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

def plot_instcat_comparison(app_mag_file, repo, visit, raft, sensor,
                            fontsize='x-small', figsize=(8, 8)):
    plt.rcParams['xtick.labelsize'] = fontsize
    plt.rcParams['ytick.labelsize'] = fontsize
    plt.rcParams['figure.figsize'] = figsize
    ccd = 'R%s S%s' % (raft, sensor)
    df, visit_name = desc.imsimdeep.instcat_comparison(app_mag_file, repo,
                                                       visit, raft, sensor)
    pointSource = df.query('galSimType=="pointSource"')
    sersic = df.query('galSimType=="sersic"')
    figure = plt.figure()

    # Plot true and measure positions for point sources only.
    figure.add_subplot(2, 2, 1)
    plt.errorbar(pointSource['coord_ra'], pointSource['coord_dec'],
                 fmt='.', label='measured position')
    plt.scatter(pointSource['coord_ra_true'], pointSource['coord_dec_true'],
                s=60, label='true position', facecolors='none', edgecolors='r')
    field = plt.axis()
    plt.xlabel('RA (deg)', fontsize=fontsize)
    plt.ylabel('Dec (deg)', fontsize=fontsize)
    plt.legend(fontsize=fontsize)
    plt.title('%(visit_name)s, %(ccd)s' % locals(), fontsize=fontsize)

    # Plot (measured-true mag) vs true mag for point sources only.
    figure.add_subplot(2, 2, 2)
    plt.errorbar(pointSource['magnitude_true'],
                 pointSource['magnitude'] - pointSource['magnitude_true'],
                 fmt='.', label='pointSource')
    axis_range = plt.axis()
    plt.plot(axis_range[:2], [0, 0], 'k:')
    plt.xlabel('true magnitude', fontsize=fontsize)
    plt.ylabel('\nmeasured - true magnitude', fontsize=fontsize)
    plt.title('%(visit_name)s, %(ccd)s' % locals(), fontsize=fontsize)

    # Histogram offsets for point sources and galaxies.
    figure.add_subplot(2, 2, 3)
    plt.hist(pointSource['offset'], bins=50, histtype='step',
             label='pointSource')
    plt.hist(sersic['offset'], bins=50, histtype='step', label='sersic')
    plt.yscale('log')
    xmin, xmax, ymin, ymax = plt.axis()
    plt.axis([xmin, xmax, 0.5, ymax])
    plt.xlabel('offset from nearest instcat source (arcsec)', fontsize=fontsize)
    plt.ylabel('\nentries / bin', fontsize=fontsize)
    plt.title('%(visit_name)s, %(ccd)s' % locals(), fontsize=fontsize)
    plt.legend(loc=1, fontsize=fontsize)

    # Quiver plot of positional offsets: measured - true positions
    figure.add_subplot(2, 2, 4)
    X, Y = pointSource['coord_ra_true'], pointSource['coord_dec_true']
    Xp, Yp = pointSource['coord_ra'], pointSource['coord_dec']
    xhat = Xp - X
    yhat = Yp - Y
    # Set arrow lengths in degrees.
    length = pointSource['offset'].values/3600./np.sqrt(xhat**2 + yhat**2)
    U = length*xhat
    V = length*yhat
    q = plt.quiver(X, Y, U, V)
    qk = plt.quiverkey(q, 0.9, 0.9, 1./3600., 'offset (1")', labelpos='N',
                       fontproperties={'size': fontsize})
    plt.axis(field)
    plt.xlabel('RA (deg)', fontsize=fontsize)
    plt.ylabel('Dec (deg)', fontsize=fontsize)
    plt.title('%(visit_name)s, %(ccd)s' % locals(), fontsize=fontsize)

    plt.tight_layout()
