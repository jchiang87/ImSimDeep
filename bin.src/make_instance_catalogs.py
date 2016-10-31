#!/usr/bin/env python
"""
Script to create instance catalogs.
"""
import numpy as np
import desc.imsimdeep

opsim_db = 'kraken_1042_sqlite.db'

ic_maker = desc.imsimdeep.InstanceCatalogMaker(opsim_db)

# Use an r-band visit from Run1.1
obsHistID = 1668469
band = 'r'
#boundLength = 0.2/3600.*2036*sqrt(2.)  # = 0.17; this should cover R22_S11
boundLength = 0.01

ic_maker.make_instance_catalog(obsHistID, band, boundLength)
