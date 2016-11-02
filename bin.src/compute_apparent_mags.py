#!/usr/bin/env python
"""
Compute apparent magnitudes for each object in a phosim instance catalog.
"""
from __future__ import absolute_import, print_function
import argparse
import pickle
import numpy as np
import pandas as pd
import desc.imsim
import desc.imsimdeep

parser = argparse.ArgumentParser()
parser.add_argument('instance_catalog', type=str,
                    help='The phosim instance catalog (text or pickled dfs)')
parser.add_argument('outfile', type=str, help='The output filename')
parser.add_argument('--numrows', type=int, default=None,
                    help='Number of rows to read from the instance catalog')
args = parser.parse_args()

try:
    commands, objs = pickle.load(open(args.instance_catalog))
    if args.numrows is not None:
        objs = objs.iloc[:args.numrows-len(commands)]
except IndexError:
    commands, objs = desc.imsim.parsePhoSimInstanceFile(args.instance_catalog,
                                                        numRows=args.numrows)

band = commands['bandpass']
columns = ('objectID', 'ra', 'dec', band)

sed_names = [x[0] for x in objs.groupby('sedName').sedName.unique()]

data_frames = []
with open(args.instance_catalog) as inst_cat:
    for sed_name in sed_names:
        app_mag = desc.imsimdeep.ApparentMagnitude(sed_name)
        my_objs = objs.query("sedName=='%s'" % sed_name)
        nrows = len(my_objs)
        mags = [app_mag(my_objs.iloc[i], band) for i in range(nrows)]
        df = pd.DataFrame(np.zeros((len(my_objs), len(columns))),
                          columns=columns)
        df['objectID'] = pd.to_numeric(my_objs['objectID']).tolist()
        df['ra'] = pd.to_numeric(my_objs['ra']).tolist()
        df['dec'] = pd.to_numeric(my_objs['dec']).tolist()
        df['galSimType'] = my_objs['galSimType'].tolist()
        df[band] = mags
        data_frames.append(df)

my_df = pd.concat(tuple(data_frames), ignore_index=True)
my_df.to_pickle(args.outfile)
