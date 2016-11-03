#!/usr/bin/env python
"""
This script will query the CatSim database for stars in the Twinkles
field and produce a reference catalog that can be used to generate
astrometry.net index files.
"""
from __future__ import absolute_import
import argparse
import numpy
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.catalogs.generation.db import CatalogDBObject
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.catUtils.mixins import AstrometryStars, PhotometryStars

class SimulationReference(InstanceCatalog, AstrometryStars, PhotometryStars):
    "Reference stars for simulation astrometry."
    catalog_type = 'simulation_ref_star'
    column_outputs = ['uniqueId', 'raJ2000', 'decJ2000', 'lsst_u', 'lsst_g',
                      'lsst_r', 'lsst_i', 'lsst_z', 'lsst_y', 'isvariable']
    default_columns = [('isresolved', 0, int), ('isvariable', 0, int)]
    default_formats = {'S': '%s', 'f': '%.8f', 'i': '%i'}
    transformations = {'raJ2000': numpy.degrees, 'decJ2000': numpy.degrees}

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate the reference catalog')
    parser.add_argument('opsimDB', help='OpSim database sqlite file')
    parser.add_argument('opshistID', type=int,
                        help='Visit number to use as the center of the field')
    parser.add_argument('-o', '--outfile', type=str,
                        default='simulation_ref.txt',
                        help='Filename of output reference catalog')
    parser.add_argument('--boundLength', type=float, default=1.5,
                        help="Radius of the extraction region in degrees")
    parser.add_argument('--host', type=str,
                        default='fatboy.phys.washington.edu',
                        help='CatSim db host name')
    parser.add_argument('--port', type=int, default=1433, help='db host port')
    args = parser.parse_args()

    # you need to provide ObservationMetaDataGenerator with the connection
    # string to an OpSim output database.  This is the connection string
    # to a test database that comes when you install CatSim.
    generator = ObservationMetaDataGenerator(database=args.opsimDB,
                                             driver='sqlite')
    obsMetaDataResults = generator.getObservationMetaData(fieldRA=(53, 54),
                                                          fieldDec=(-29, -27),
                                                          boundLength=0.3)

    stars = CatalogDBObject.from_objid('allstars', host=args.host,
                                       port=args.port)

    while True:
        try:
            ref_stars = SimulationReference(stars,
                                            obs_metadata=obsMetaDataResults[0])
            break
        except RuntimeError:
            continue
    ref_stars.write_catalog(args.outfile, write_mode='w', write_header=True,
                            chunk_size=20000)
