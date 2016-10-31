"""
Script to create instance catalogs.
"""
import desc.imsimdeep

opsim_db = 'kraken_1042_sqlite.db'

ic_maker = desc.imsimdeep.InstanceCatalogMaker(opsim_db)

# Selected visits from Twinkles Run 1.1.
obsHistIDs = (1668469, 1648025, 1414156, 1973403, 921297)
band = 'r'
boundLength = 2.5

for obsHistID in obsHistIDs:
    ic_maker.make_instance_catalog(obsHistID, band, boundLength)
