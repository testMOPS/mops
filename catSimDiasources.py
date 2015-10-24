import numpy as np
# Observation metadata modules
from lsst.sims.utils import ObservationMetaData
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
# photometric parameters (exptime lives here for dmag calculation)
from lsst.sims.photUtils import PhotometricParameters
# SSM catalog modules
from lsst.sims.catUtils.baseCatalogModels import SolarSystemObj, CometObj, MBAObj, NEOObj, MiscSolarSystemObj
from lsst.sims.catUtils.mixins import PhotometrySSM, AstrometrySSM
from lsst.sims.catalogs.measures.instance import InstanceCatalog

import time
def dtime(time_prev):
    return (time.time() - time_prev, time.time())

# Build sso instance class


class ssmCat(InstanceCatalog, PhotometrySSM, AstrometrySSM):
    column_outputs = ['objid', 'raJ2000', 'decJ2000', 'sedFilename',
                      'velRa', 'velDec', 'skyVelocity', 'raObserved', 'decObserved',
                      'lsst_u' ,'lsst_g', 'lsst_r', 'lsst_i', 'lsst_z', 'lsst_y',
                      'dmagTrailing', 'dmagDetection']

    transformations = {'raJ2000': np.degrees, 'decJ2000': np.degrees,
                       'velRa': np.degrees, 'velDec': np.degrees}


######

t = time.time()
# Get opsim data.
opsdb = '/Users/lynnej/opsim/db/enigma_1189_sqlite.db'
generator = ObservationMetaDataGenerator(database=opsdb, driver='sqlite')
obsMetaDataResults = generator.getObservationMetaData(expMJD=(50100.031278, 50100.375172), limit=3, boundLength=2.2)
# Need to get exptime above, and move exptime into photParams.
photParams = PhotometricParameters()

dt, t = dtime(t)
print 'To query opsim database: %f seconds' %(dt)

for obs in obsMetaDataResults:
    print obs.mjd, obs.unrefractedRA, obs.unrefractedDec, obs.bandpass, obs.boundType, obs.boundLength
    """
    mySsmDb = NEOObj()
    chunk = mySsmDb.query_columns(obs_metadata=obs)
    for c in chunk:
        print c
        print c.dtype
    """
    mySsmDb = ssmCat(NEOObj(), obs_metadata = obs)
    #for row in mySsmDb.iter_catalog():
    #    print row
    mySsmDb.write_catalog('test')

    dt, t = dtime(t)
    print 'To query solar system objects: %f seconds' %(dt)





