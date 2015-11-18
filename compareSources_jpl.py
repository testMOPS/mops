import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from lsst.sims.utils import haversine

match_tol = 0.3 # degrees for max limit for match between jpl/catsim.

def match(catsim, jpl, match_tol=match_tol):
    times = np.unique(np.concatenate([catsim['expMJD'], jpl['epoch_mjd']]))
    dra = []
    ddec = []
    dmag = []
    ra_jpl = np.radians(jpl['ra_deg'])
    dec_jpl = np.radians(jpl['dec_deg'])
    ra_catsim = np.radians(catsim['raJ2000'])
    dec_catsim = np.radians(catsim['decJ2000'])
    for t in times:
        dra_t = []
        ddec_t = []
        dmag_t = []
        tmatch_catsim = np.where(np.abs(catsim['expMJD'] - t) < 15./60./60./24.0)[0]
        tmatch_jpl = np.where(np.abs(jpl['epoch_mjd'] - t) < 15./60./60./24.)[0]
        if len(tmatch_catsim) and len(tmatch_jpl) > 0:
            # catsim and jpl ids are different. Have to match on ra/dec.
            for ra, dec, mag, objid in zip(ra_catsim[tmatch_catsim], dec_catsim[tmatch_catsim],
                                           catsim['magFilter'][tmatch_catsim], catsim['#objid'][tmatch_catsim]):
                sep = haversine(ra_jpl[tmatch_jpl], dec_jpl[tmatch_jpl], ra, dec)
                sep = np.degrees(sep)
                if sep.min() < match_tol:
                    match = np.where(sep == sep.min())[0]
                    dra_t.append(np.degrees(ra-ra_jpl[tmatch_jpl][match])*3600.)
                    ddec_t.append(np.degrees(dec-dec_jpl[tmatch_jpl][match])*3600.)
                    dmag_t.append(mag - jpl['mag'][tmatch_jpl][match])
                    if sep.min() > 0.1:
                        print objid, jpl['object_name'][tmatch_jpl][match]
        if len(dra_t) > 0:
            dra.append(np.array(dra_t).flatten())
            ddec.append(np.array(ddec_t).flatten())
            dmag.append(np.array(dmag_t).flatten())
    return dra, ddec, dmag

# Read my old sources - generated with focal plane included, from catalog propagated to 49353 (in 2008?)
infile = 'n747/catsim_n747'
catsim = pd.read_csv(infile, sep=', ', skipinitialspace=False)
print catsim.columns

# Read JPL sources
infile = 'n747/n747_dets.txt'
jpl = pd.read_table(infile, delim_whitespace=True)
print jpl.columns


catsim = catsim.to_records()
jpl = jpl.to_records()

plt.figure()
plt.plot(catsim['raJ2000'], catsim['decJ2000'], 'r.', label='CatSim')
plt.plot(jpl['ra_deg'], jpl['dec_deg'], 'k.', label='JPL')
plt.legend(loc='lower right', numpoints=1)
plt.title('Catalog contents')


dra, ddec, dmag = match(catsim, jpl)

mag = []
for m in dmag:
    for im in m:
        mag.append(im)
mag = np.array(mag)


plt.figure()
for ddra, dddec in zip(dra, ddec):
    plt.plot(ddra, dddec, 'k.')
plt.xlabel('dRA (arcseconds)')
plt.ylabel('dDec (arcseconds)')
theta = np.radians(np.arange(0, 360, 0.2))
x = np.cos(theta) * match_tol*3600.0
y = np.sin(theta) * match_tol*3600.0
plt.plot(x, y, 'b-')
plt.title('Matches, using max separation tolerance of %.1f Degrees' %(match_tol))

#dra, ddec, dmag = match(catsim, jpl)

plt.figure()
ra_dev = 0
dec_dev = 0
for ddra, dddec in zip(dra, ddec):
    plt.plot(ddra, dddec, 'k.')
    ra_dev += np.abs(np.mean(ddra))
    dec_dev += np.abs(np.mean(dddec))
ra_dev = ra_dev / float(len(dra))
dec_dev = dec_dev / float(len(dra))
plt.xlabel('dRA (arcseconds)')
plt.ylabel('dDec (arcseconds)')
plt.title('Matches, using max separation tolerance of %.1f Degrees' %(match_tol))
plt.figtext(0.6, 0.2, 'ave abs mean: %.2f %.2f' %(ra_dev, dec_dev))


plt.figure()
plt.hist(mag, bins=500)
mag_dev = np.abs(np.mean(mag))
plt.xlabel('dMag')
plt.title('Matches, using max separation tolerance of %.1f deg' %(match_tol))
plt.figtext(0.7, 0.8, 'ave abs mean: %.2f' %(mag_dev))
plt.xlim(-1, 1)

plt.show()
