"""
Tools for manipulating instance catalogs read in as pandas DataFrames.
"""
from __future__ import print_function
import os
import sys
import time
import glob
import logging
import subprocess
import numpy as np
import pandas as pd
import lsst.obs.lsstSim as obs_lsstSim
from lsst.sims.coordUtils import chipNameFromRaDec, raDecFromPixelCoords
from lsst.sims.photUtils import LSSTdefaults
from lsst.sims.utils import ObservationMetaData
import desc.imsim

__all__ = ['apply_acceptance_cone', 'select_by_chip_name',
           'obs_metadata', 'LargeInstanceCatalog']

logging.basicConfig(format="%(message)s", level=logging.DEBUG,
                    stream=sys.stdout)
default_logger = logging.getLogger()

def apply_acceptance_cone(objs, ra, dec, radius, logger=default_logger):
    """
    Apply an acceptance cone cut using small angle approximation.

    Parameters
    ----------
    objs : pandas.DataFrame
        DataFrame of phosim objects.
    ra : float
        RA (ICRS) of the center of the acceptance cone in degrees.
    dec : float
        Dec (ICRS) of the center of the acceptance cone in degrees.
    radius : float
        Radius of the acceptance cone.
    logger : logging.Logger, optional
        The logger to use.

    Returns
    -------
    pandas.DataFrame
        The DataFrame containing down-selected objects.
    """
    t0 = time.time()
    my_objs = objs.copy(deep=True)
    my_objs['distance'] = \
        np.sqrt((my_objs['dec'] - dec)**2 +
                (np.cos(my_objs['dec']*np.pi/180.)*(my_objs['ra'] - ra))**2)
    my_objs = my_objs.query('distance<=%.8f' % radius)
    del my_objs['distance']
    logger.debug('apply_acceptance_cone:\n  elapsed time: %f s',
                 time.time()- t0)
    logger.debug('  # objects remaining: %i\n', len(my_objs))
    return my_objs

def select_by_chip_name(objs, chip_name, obs_metadata, camera,
                        logger=default_logger):
    """
    Select only objects that are on the specified chip.
    Parameters
    ----------
    objs : pandas.DataFrame
        DataFrame of phosim objects.
    obs_metadata : lsst.sims.utils.ObservationMetaData
        Obsevation metadata extracted from the phosim commands.
    camera : lsst.afw.cameraGeom.camera.Camera
        The camera instance from lsst.obs.lsstSim.LsstSimMapper().
    logger : logging.Logger, optional
        The logger to use.

    Returns
    -------
    pandas.DataFrame
        The DataFrame containing down-selected objects.
    """
    t0 = time.time()
    my_objs = objs.copy(deep=True)
    my_objs['chip_name'] = chipNameFromRaDec(my_objs['ra'].values,
                                             my_objs['dec'].values,
                                             camera=camera,
                                             obs_metadata=obs_metadata)
    my_objs = my_objs.query('chip_name=="%s"' % chip_name)
    del my_objs['chip_name']
    logger.debug('select_by_chip_name:\n  elapsed time: %f s',
                 time.time()- t0)
    logger.debug('  # objects remaining: %i\n', len(my_objs))
    return my_objs

def obs_metadata(phosim_commands):
    """
    Create an ObservationMetaData instance from phosim commands.

    Parameters
    ----------
    phosim_commands : dict
        Dictionary of phosim instance catalog commands.

    Returns:
    lsst.sims.utils.ObservationMetaData
    """
    return ObservationMetaData(pointingRA=phosim_commands['rightascension'],
                               pointingDec=phosim_commands['declination'],
                               mjd=phosim_commands['mjd'],
                               rotSkyPos=phosim_commands['rotskypos'],
                               bandpassName=phosim_commands['bandpass'],
                               m5=LSSTdefaults().m5(phosim_commands['bandpass']),
                               seeing=phosim_commands['seeing'])

class LargeInstanceCatalog(object):
    """
    Class to manage down-selections of large instance catalogs.
    """
    def __init__(self, instcat_file, num_lines=300000, temp_file_prefix='temp_',
                 logger=default_logger):
        self._temp_file = temp_file_prefix + "_instcat.txt"
        self._sub_files = self.split_catalog(instcat_file, num_lines, logger)
        self._header_file = self.extract_header(instcat_file, temp_file_prefix,
                                                logger)
        self._read_obs_metadata()
        self._camera = obs_lsstSim.LsstSimMapper().camera
        self._logger = logger

    @staticmethod
    def split_catalog(instcat_file, num_lines, logger=default_logger):
        """
        Split the instance catalog into manageably-sized files of just
        the objects.

        Parameters
        ----------
        instcat_file : str
            phosim instance catalog file.
        num_lines : int
            The maximum number of lines in each sub-file.
        logger : logging.Logger, optional
            The logger to use.

        Returns
        -------
        list
            A sorted list of temporary files containing the object lines.
        """
        command = "grep object %(instcat_file)s | split --lines=%(num_lines)i"\
            % locals()
        logger.debug("executing:\n  %s", command)
        t0 = time.time()
        subprocess.call(command, shell=True)
        logger.debug("  elapsed time: %f\n", time.time() - t0)
        return sorted(glob.glob('x??'))

    @staticmethod
    def extract_header(instcat_file, temp_file_prefix, logger=default_logger):
        """
        Get the header to prepend to each sub-file.

        Parameters
        ----------
        instcat_file : str
            phosim instance catalog file.
        temp_file_prefix : str
            The prefix to prepend to the temporary header file.
        logger : logging.Logger, optional
            The logger to use.

        Returns
        -------
        str
            The name of the temporary header file.
        """
        header_file = temp_file_prefix + "_header.txt"
        command = "grep -v object %(instcat_file)s > %(header_file)s" % locals()
        logger.debug("executing:\n  %s", command)
        t0 = time.time()
        subprocess.call(command, shell=True)
        logger.debug("  elapsed time: %f\n", time.time() - t0)
        return header_file

    def clean_up(self):
        """
        Delete all temporary files.
        """
        os.remove(self._header_file)
        for sub_file in self._sub_files:
            os.remove(sub_file)
        os.remove(self._temp_file)

    def cone_select(self, ra, dec, radius):
        """
        Apply acceptance cone selection on the object data.

        Parameters
        ----------
        ra : float
            RA (ICRS) of the center of the acceptance cone in degrees.
        dec : float
            Dec (ICRS) of the center of the acceptance cone in degrees.
        radius : float
            Radius of the acceptance cone.

        Returns
        -------
        tuple(dict, pandas.DataFrame)
            The phosim commands and objects DataFrame.
        """
        return self._apply_selection(ra, dec, radius)

    def chip_select(self, chip_name, radius=0.17):
        """
        Select objects landing on a chip.

        Parameters
        ----------
        chip_name : str
            The chip name, e.g., "R:2,2 S:1,1".
        radius : float, optional
            Radius in degrees of the acceptance cone (centered on the
            middle of the chip) to apply to pre-downselect before
            computing the chip names for each object.

        Returns
        -------
        tuple(dict, pandas.DataFrame)
            The phosim commands and objects DataFrame.
        """
        ra, dec = tuple(raDecFromPixelCoords(2036, 2000, chip_name,
                                             camera=self._camera,
                                             obs_metadata=self._obs_metadata))
        return self._apply_selection(ra, dec, radius, chip_name=chip_name)

    def _read_obs_metadata(self):
        "Read the obsevation metadata for the chip-based acceptance cone."
        commands, objs = desc.imsim.parsePhoSimInstanceFile(self._header_file)
        self._obs_metadata = obs_metadata(commands)

    def _apply_selection(self, ra, dec, radius, chip_name=None):
        "Apply acceptance cone and optional chip-based selection"
        times = [time.time()]
        object_dfs = []
        header_file = self._header_file
        temp_file = self._temp_file
        for sub_file in self._sub_files:
            command = "cat %(header_file)s %(sub_file)s > %(temp_file)s" \
                % locals()
            self._logger.debug("executing:\n  %s\n", command)
            subprocess.call(command, shell=True)
            commands, objs = desc.imsim.parsePhoSimInstanceFile(temp_file)
            times.append(time.time())
            self._logger.debug("parsed instance catalog sub_file:")
            self._logger.debug("  elapsed time: %f\n", times[-1] - times[-2])
            objs = apply_acceptance_cone(objs, ra, dec, radius,
                                         logger=self._logger)
            if chip_name is not None and len(objs) > 0:
                objs = select_by_chip_name(objs, chip_name, self._obs_metadata,
                                           self._camera, logger=self._logger)
            object_dfs.append(objs)
        objs = pd.concat(tuple(object_dfs))
        self._logger.debug('down selection applied:\n  elapsed time: %f s\n',
                           time.time()- times[0])
        return commands, objs
