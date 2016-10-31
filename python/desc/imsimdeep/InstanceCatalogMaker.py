"""
Code to create instance catalogs via CatSim.
"""
from __future__ import absolute_import, print_function
import sys
import logging
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings('ignore', 'Duplicate object type id', UserWarning)
    warnings.filterwarnings('ignore', 'duplicate object identifie', UserWarning)
    from lsst.sims.catalogs.db import CatalogDBObject
    from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
    from lsst.sims.catUtils.exampleCatalogDefinitions.phoSimCatalogExamples \
        import PhoSimCatalogPoint, PhoSimCatalogSersic2D

__all__ = ['InstanceCatalogMaker']

class InstanceCatalogMaker(object):
    """
    Class for creating instance catalogs.

    Attributes
    ----------
    gen : lsst.sims.catUtils.utils.ObservationMetaDataGenerator

    db_config : dict
        Dictionary of database connection parameters.

    logger : logging.logger
        Logger object.
    """
    star_objs = ['msstars', 'bhbstars', 'wdstars', 'rrlystars', 'cepheidstars']
    gal_objs = ['galaxyBulge', 'galaxyDisk']
    def __init__(self, opsim_db, db_config=None, logger=None):
        """
        Constructor.

        Parameters
        ----------
        opsim_db : str
            sqlite3 db file containing observing plan.

        db_config : dict, optional
            Dictionary of database connection parameters.  Parameters
            for connecting to fatboy.phys.washington.edu from a
            whitelisted machine will be used.

        logger : logging.logger, optional
            Logger object.
        """
        self.gen = ObservationMetaDataGenerator(database=opsim_db,
                                                driver='sqlite')
        if db_config is not None:
            self.db_config = db_config
        else:
            self.db_config = dict(database='LSSTCATSIM',
                                  port=1433,
                                  host='fatboy.phys.washington.edu',
                                  driver='mssql+pymssql')
        if logger is None:
            logging.basicConfig(format="%(message)s", level=logging.INFO,
                                stream=sys.stdout)
            logger = logging.getLogger()
        self.logger = logger

    def make_instance_catalog(self, obsHistID, band, boundLength, outfile=None):
        """
        Method to create instance catalogs.

        Parameters
        ----------
        obsHistID : int
            obsHistID for the desired visit from the opsim db file.

        band : str
            Desired LSST filter to use, ugrizy.

        boundLength : float
            Radius in degrees of sky cone in which to produce objects.

        outfile : str, optional
            File name of the instance catalog to be produced.  If None,
            a default name will be generated, e.g.,
            phosim_input_0000230_r_0.3deg.txt.
        """
        if outfile is None:
            outfile = 'phosim_input_%07i_%s_%.1fdeg.txt' % (obsHistID, band,
                                                            boundLength)
            obs_md = self.gen.getObservationMetaData(obsHistID=obsHistID,
                                                     boundLength=boundLength)[0]
        do_header = True
        for objid in self.star_objs:
            self.logger.info("processing %s", objid)
            db_obj = CatalogDBObject.from_objid(objid, **self.db_config)
            phosim_object = PhoSimCatalogPoint(db_obj, obs_metadata=obs_md)
            if do_header:
                with open(outfile, 'w') as file_obj:
                    phosim_object.write_header(file_obj)
                do_header = False
            phosim_object.write_catalog(outfile, write_mode='a',
                                        write_header=False,
                                        chunk_size=20000)

        for objid in self.gal_objs:
            self.logger.info("processing %s", objid)
            db_obj = CatalogDBObject.from_objid(objid, **self.db_config)
            phosim_object = PhoSimCatalogSersic2D(db_obj, obs_metadata=obs_md)
            phosim_object.write_catalog(outfile, write_mode='a',
                                        write_header=False,
                                        chunk_size=20000)
