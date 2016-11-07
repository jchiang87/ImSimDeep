"""
Tools for manipulating instance catalogs read in as pandas DataFrames.
"""
from __future__ import absolute_import, print_function, division
import os
import time
import psutil
from lsst.sims.coordUtils import chipNameFromRaDec
from lsst.sims.photUtils import LSSTdefaults
from lsst.sims.utils import ObservationMetaData
import desc.imsim

__all__ = ['select_by_chip_name', 'obs_metadata']

default_logger = desc.imsim.get_logger("DEBUG")

def mem_use_message(pid=None):
    """
    Return memory usage string.

    Parameters
    ----------
    pid : int, optional
        Process id, e.g., from os.getpid().
    """
    if pid is None:
        pid = os.getpid()
    process = psutil.Process(pid)
    mem_info = process.memory_full_info()
    return "  memory used: %.3f GB\n" % (mem_info.uss/1024.**3)

def select_by_chip_name(objs, chip_name, obs_md, camera, logger=default_logger):
    """
    Select only objects that are on the specified chip.
    Parameters
    ----------
    objs : pandas.DataFrame
        DataFrame of phosim objects.
    obs_md : lsst.sims.utils.ObservationMetaData
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
                                             obs_metadata=obs_md)
    my_objs = my_objs.query('chip_name=="%s"' % chip_name)
    del my_objs['chip_name']
    logger.debug('select_by_chip_name:\n  elapsed time: %f s',
                 time.time()- t0)
    logger.debug('  # objects remaining: %i', len(my_objs))
    logger.debug(mem_use_message())
    return my_objs

def obs_metadata(commands):
    """
    Create an ObservationMetaData instance from phosim commands.

    Parameters
    ----------
    commands : dict
        Dictionary of phosim instance catalog commands.

    Returns:
    lsst.sims.utils.ObservationMetaData
    """
    return ObservationMetaData(pointingRA=commands['rightascension'],
                               pointingDec=commands['declination'],
                               mjd=commands['mjd'],
                               rotSkyPos=commands['rotskypos'],
                               bandpassName=commands['bandpass'],
                               m5=LSSTdefaults().m5(commands['bandpass']),
                               seeing=commands['seeing'])
