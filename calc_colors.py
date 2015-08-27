# Requirements: LSST packages photUtils and throughputs
# Get the asteroid spectra and harris_V throughput curve at http://www.astro.washington.edu/users/lynnej/seds/asteroids.tar
# Run this in the directory containing the asteroids SEDs. (note that it assumes there are no extraneous *.dat files in that directory).

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from lsst.sims.photUtils import Bandpass, Sed

filterdir = os.getenv('LSST_THROUGHPUTS_BASELINE')
filterlist = ('u', 'g', 'r', 'i', 'z', 'y')

lsst = {}
for f  in filterlist:
    lsst[f] = Bandpass()
    lsst[f].readThroughput(os.path.join(filterdir, 'total_' + f + '.dat'))
harrisV = Bandpass()
harrisV.readThroughput('harris_V.dat')

seds = glob.glob('*.dat')
seds.remove('harris_V.dat')

print 'Calculating colors for seds: \n', seds

writestring = 'Sed '
for f in filterlist:
    writestring += 'V-%s ' %(f)
print writestring


fig = plt.figure()
for s in seds:
    sed = Sed()
    sed.readSED_flambda(s)
    sname = s.replace('.dat', '')
    vmag = sed.calcMag(harrisV)
    colors = {}
    writestring = '%s ' %(sname)
    for f in filterlist:
        mag = sed.calcMag(lsst[f])
        colors[f] = vmag - mag
        writestring += '%.3f ' %(colors[f])
    print writestring
    plotPointColor = np.random.rand(3)
    plt.plot(colors['g'], colors['r'], linestyle='', marker='o', color=plotPointColor, label='%s' %(sname))
    pylab.annotate(sname, colors['g']+0.005, colors['r']+0.005, color=plotPointColor)

plt.xlabel('V-g')
plt.ylabel('V-r')
#plt.legend(fancybox=True, numpoints=1, fontsize='smaller')
plt.show()


