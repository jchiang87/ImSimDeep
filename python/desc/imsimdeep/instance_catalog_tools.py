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
from lsst.obs.lsstSim import LsstSimMapper
from lsst.sims.coordUtils import chipNameFromRaDec
from lsst.sims.photUtils import LSSTdefaults
from lsst.sims.utils import ObservationMetaData
import desc.imsim

__all__ = ['apply_acceptance_cone', 'select_by_chip_name',
           'split_and_prune_instcat']

logging.basicConfig(format="%(message)s", level=logging.DEBUG,
                    stream=sys.stdout)
default_logger = logging.getLogger()

def apply_acceptance_cone(objs, ra, dec, radius, logger=default_logger):
    """
    Apply an acceptance cone cut using small angle approximation.
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

def select_by_chip_name(objs, commands, chip_name, logger=default_logger):
    """
    Select only objects that are on the specified chip.
    """
    t0 = time.time()
    my_objs = objs.copy(deep=True)
    obs_md = ObservationMetaData(pointingRA=commands['rightascension'],
                                 pointingDec=commands['declination'],
                                 mjd=commands['mjd'],
                                 rotSkyPos=commands['rotskypos'],
                                 bandpassName=commands['bandpass'],
                                 m5=LSSTdefaults().m5(commands['bandpass']),
                                 seeing=commands['seeing'])
    my_objs['chip_name'] = chipNameFromRaDec(my_objs['ra'].values,
                                             my_objs['dec'].values,
                                             camera=LsstSimMapper().camera,
                                             obs_metadata=obs_md)
    my_objs = my_objs.query('chip_name=="%s"' % chip_name)
    del my_objs['chip_name']
    logger.debug('select_by_chip_name:\n  elapsed time: %f s',
                 time.time()- t0)
    logger.debug('  # objects remaining: %i\n', len(my_objs))
    return my_objs

def split_and_prune_instcat(instcat_file, ra, dec, radius, chip_name=None,
                            num_lines=300000, temp_file='temp_instcat.txt',
                            clean_up=True, logger=default_logger):
    """
    Split an instance catalog into smaller files, and apply desired
    pruning to each file, concatenating the corresponding DataFrames.
    """
    t0 = time.time()
    # Split the instance catalog into manageably-sized files of just
    # the objects.
    command = \
        "grep object %(instcat_file)s | split --lines=%(num_lines)i" % locals()
    logger.debug("executing:\n  %s", command)
    t0 = time.time()
    subprocess.call(command, shell=True)
    logger.debug("  elapsed time: %f\n", time.time() - t0)

    # Get the header to prepend to each sub-file.
    header_file = "header_file.txt"
    command = "grep -v object %(instcat_file)s > %(header_file)s" % locals()
    logger.debug("executing:\n  %s", command)
    t0 = time.time()
    subprocess.call(command, shell=True)
    logger.debug("  elapsed time: %f\n", time.time() - t0)

    # Loop over split files.
    sub_files = sorted(glob.glob('x??'))

    times = [time.time()]
    object_dfs = []
    for sub_file in sub_files:
        command = "cat %(header_file)s %(sub_file)s > %(temp_file)s" % locals()
        logger.debug("executing:\n  %s\n", command)
        subprocess.call(command, shell=True)
        # Read in the entire instance catalog.
        commands, objs = desc.imsim.parsePhoSimInstanceFile(temp_file)
        times.append(time.time())
        logger.debug("parsed instance catalog sub_file:")
        logger.debug("  elapsed time: %f\n", times[-1] - times[-2])
        objs = apply_acceptance_cone(objs, ra, dec, radius, logger=logger)
        if chip_name is not None and len(objs) > 0:
            objs = select_by_chip_name(objs, commands, chip_name, logger=logger)
        object_dfs.append(objs)

    objs = pd.concat(tuple(object_dfs))

    logger.debug('split_and_prune_instcat:\n  elapsed time: %f s\n',
                 time.time()- t0)

    if clean_up:
        os.remove(header_file)
        for sub_file in sub_files:
            os.remove(sub_file)
        os.remove(temp_file)

    return commands, objs
