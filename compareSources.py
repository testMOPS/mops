import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def match(maf, catsim):
    times = np.unique(np.concatenate([maf['time'], catsim['obs_mjd']]))
    dra = []
    ddec = []
    for t in times:
        dra_t = []
        ddec_t = []
        tmatch_maf = np.where(np.abs(maf['time'] - t) < 15./60./60./24.0)[0]
        tmatch_catsim = np.where(np.abs(catsim['obs_mjd'] - t) < 15./60./60./24.)[0]
        if len(tmatch_maf) and len(tmatch_catsim) > 0:
            # maf and catsim ids are the same.
            objids = np.unique(np.concatenate([maf['objId'][tmatch_maf], catsim['objId'][tmatch_catsim]]))
            for objid in objids:
                maf_match = np.where(maf['objId'][tmatch_maf] == objid)[0]
                catsim_match = np.where(catsim['objId'][tmatch_catsim] == objid)[0]
                if len(maf_match)>0 and len(catsim_match) >0:
                    dra_t.append(maf['ra'][tmatch_maf][maf_match][0] - catsim['raJ2000'][tmatch_catsim][catsim_match][0])
                    ddec_t.append(maf['dec'][tmatch_maf][maf_match][0] - catsim['decJ2000'][tmatch_catsim][catsim_match][0])
        dra.append(np.array(dra_t)*3600.)
        ddec.append(np.array(ddec_t)*3600.)
    return dra, ddec

# Read my old sources - generated with focal plane included, from catalog propagated to 49353 (in 2008?)
infile = 'neos_s3m_747_allObs.txt'
lj_maf = pd.read_table(infile, delim_whitespace=True)
print lj_maf.columns

# Read my catsim sources - generated without focal plane at the moment, from ?catalog
infile = 'catsim_n747'
lj_catsim = pd.read_csv(infile, skipinitialspace=True)
newcols = lj_catsim.columns.values
newcols[0] = 'objId'
lj_catsim.columns = newcols
print lj_catsim.columns

lj_maf = lj_maf.to_records()
lj_catsim = lj_catsim.to_records()
condition = np.where(np.abs(lj_maf['expMJD'] - 50100.0317) < 10e-5)
condition = np.ones(len(lj_maf['ra']), dtype='bool')

plt.figure()
plt.plot(lj_maf['ra'][condition], lj_maf['dec'][condition], 'ko', label='MAF')
plt.plot(lj_catsim['raJ2000'], lj_catsim['decJ2000'], 'r.', label='CatSim')


dra, ddec = match(lj_maf, lj_catsim)
plt.figure()
for ddra, ddec in zip(dra, ddec):
    plt.plot(ddra, ddec, 'k.')

plt.show()
