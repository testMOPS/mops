# Observation metadata modules
from lsst.sims.utils import ObservationMetaData
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
# SSM catalog modules
from lsst.sims.catalogs.generation.db import CatalogDBObject
from lsst.sims.catUtils.baseCatalogModels import SolarSystemObj

from lsst.sims.catalogs.measures.instance import InstanceCatalog

# Build sso instance class



# Get opsim data.
opsdb = '/Users/lynnej/opsim/db/enigma_1189_sqlite.db'
generator = ObservationMetaDataGenerator(database=opsdb, driver='sqlite')

obsMetaDataResults = generator.getObservationMetaData(expMJD=(50100.031278, 50100.375172), limit=3)

mySsmDb = CatalogDBObject.from_objid('ssm')

for obsMetaData in obsMetaDataResults:
    print obsMetaData.mjd, obsMetaData.unrefractedRA, obsMetaData.unrefractedDec, obsMetaData.bandpass



