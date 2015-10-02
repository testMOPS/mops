import sys
import numpy as np
from astropy.time import Time


# Oh my freaking god, there are too many timescales and everyone uses a different one.


if len(sys.argv) < 4:
    print "usage: python convertTimes.py YourTime YourFormat YourTimeScale"
    print " YourFormat ought to be in JD or MJD (but could be ISO as well)."
    print " outputs: your time in many other timescales"
    exit()


yourtime = float(sys.argv[1])
yourformat = sys.argv[2]
yourtimescale = sys.argv[3]

print "Converting %f / %s / %s to other timescales." %(yourtime, yourformat, yourtimescale)

t = Time(yourtime, format=yourformat, scale=yourtimescale)

print 'TT: %f %f' %(t.tt.jd, t.tt.mjd)
print 'TAI: %f %f' %(t.tai.jd, t.tai.mjd)
print 'TDB: %f %f' %(t.tdb.jd, t.tdb.mjd)
print 'UTC: %f %f' %(t.utc.jd, t.utc.mjd)
#print 'UT1: %f %f' %(t.ut1.jd, t.ut1.mjd)


