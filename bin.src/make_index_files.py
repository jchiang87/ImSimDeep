#!/usr/bin/env python
"""
This script will query the CatSim database for stars in the Twinkles
field and produce a reference catalog, and from that, generate
astrometry.net index files.
"""
from __future__ import absolute_import
import os
import argparse
import desc.imsimdeep

parser = argparse.ArgumentParser(description='Generate the astrometry.net index files based on a OpSim visit.')
parser.add_argument('opsim_db', help='OpSim database sqlite file')
parser.add_argument('obsHistID', type=int,
                    help='Visit number for the center of the field')
parser.add_argument('--outfile_root', type=str, default='simulation_ref',
                    help='Root filename of reference catalogs (.txt and .fits)')
parser.add_argument('--boundLength', type=float, default=0.3,
                    help="Radius of the extraction region in degrees")
parser.add_argument('--index_id', type=str, default=None,
                    help="ID string to identify astrometry.net index files.  If None, then construct it from the root of opsim_db + obsHistID")
parser.add_argument('--max_scale', type=int, default=4,
                    help='maximum scale for astrometry.net index files')
parser.add_argument('--host', type=str,
                    default='fatboy.phys.washington.edu',
                    help='CatSim database host name')
parser.add_argument('--port', type=int, default=1433,
                    help='CatSim database host port')
parser.add_argument('--database', type=str, default='LSSTCATSIM',
                    help='CatSim database name')
parser.add_argument('--driver', type=str, default='mssql+pymssql',
                    help='CatSim database driver')
args = parser.parse_args()

db_info = dict(host=args.host, port=args.port, database=args.database,
               driver=args.driver)

refcat_txt = args.outfile_root + '.txt'
refcat_fits = args.outfile_root + '.fits'

if args.index_id is None:
    index_id = (os.path.basename(args.opsim_db).split('.')[0] + '_'
                + str(args.obsHistID) + '_')
else:
    index_id = args.index_id

desc.imsimdeep.make_refcat(args.opsim_db, args.obsHistID, args.boundLength,
                           refcat_txt, catsim_db_info=db_info)

desc.imsimdeep.refcat_to_astrometry_net_input(refcat_txt, outfile=refcat_fits)

desc.imsimdeep.build_index_files(refcat_fits, index_id,
                                 max_scale_number=args.max_scale)
