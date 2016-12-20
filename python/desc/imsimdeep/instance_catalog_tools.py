"""
Tools for manipulating instance catalogs read in as pandas DataFrames.
"""
from __future__ import absolute_import, print_function, division
import os
import time
import subprocess
import psutil
import lsst.sims.coordUtils as coordUtils
from lsst.sims.photUtils import LSSTdefaults
from lsst.sims.utils import ObservationMetaData
import desc.imsim

__all__ = ['select_by_chip_name', 'obs_metadata', 'instcat_commands',
           'chip_center_coords']

default_logger = desc.imsim.get_logger("DEBUG")

def _cast(value):
    try:
        if '.' in value or 'e' in value:
            return float(value)
        else:
            return int(value)
    except ValueError:
        return value

def instcat_commands(instcat_file, numlines=1000):
    """
    Read the commands from an instance catalog.

    Parameters
    ----------
    instcat_file : str
        The filename of the instance catalog file.
    numlines : int, optional
        The number of lines from the top of the file to read. Default: 1000

    Returns
    -------
    dict
        The PhoSim instance catalog physics commands.
    """
    phosim_commands = dict()
    command = "head -%i %s | grep -v object" % (numlines, instcat_file)
    lines = subprocess.check_output(command, shell=True).split('\n')
    for line in lines:
        if line.startswith('#'):
            continue
        tokens = line.split()
        try:
            phosim_commands[tokens[0]] = _cast(tokens[1])
        except IndexError:
            pass
    phosim_commands['bandpass'] = 'ugrizy'[phosim_commands['filter']]
    return phosim_commands

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

def chip_center_coords(chip_name, obs_md, camera):
    """
    The coordinates of the center of the specified chip.

    Parameters
    ----------
    chip_name : str
        The chip name, e.g., "R:2,2 S:1,1".
    obs_md : lsst.sims.utils.ObservationMetaData
        Obsevation metadata extracted from the phosim commands.
    camera : lsst.afw.cameraGeom.camera.Camera
        The camera instance from lsst.obs.lsstSim.LsstSimMapper().
    """
    corner_pixels = coordUtils.getCornerPixels(chip_name, camera)
    xmid = (corner_pixels[-1][0] - corner_pixels[0][0] + 1)/2
    ymid = (corner_pixels[-1][1] - corner_pixels[0][1] + 1)/2
    return tuple(coordUtils.raDecFromPixelCoords(ymid, xmid, chip_name,
                                                 camera=camera,
                                                 obs_metadata=obs_md))

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
    my_objs['chip_name'] = coordUtils.chipNameFromRaDec(my_objs['ra'].values,
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
