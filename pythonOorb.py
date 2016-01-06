import os
import argparse
import numpy as np
from itertools import repeat
import pandas as pd
import pyoorb as oo


# Simple script to read an orbit file from oorb and generate an ephemeris at a user-specified (UTC) time.

def pack_oorbArray(orbits):
    """Translate orbital element dictionary (easy for humans) into pyoorb-suitable input orbit array."""
    # Translate orbital elements into array that pyoorb will like.
    # PyOrb wants ::
    # 0: orbitId
    # 1 - 6: orbital elements, using radians for angles
    # 7: element type code, where 2 = cometary - means timescale is TT, too
    # 8: epoch
    # 9: timescale for the epoch; 1= MJD_UTC, 2=UT1, 3=TT, 4=TAI (we use 3)
    # 10: magHv
    # 11: G
    elem_type = np.zeros(len(orbits)) + 2
    epoch_type = np.zeros(len(orbits)) + 3
    gval = np.zeros(len(orbits)) + 0.15
    # Also, the orbitID has to be a float, rather than a string, so substitute if needed.
    if ((isinstance(orbits['objid'][0], float) == True) |
        (isinstance(orbits['objid'][0], int) == True)):
        orbids = orbits['objid']
    else:
        orbids = np.arange(0, len(orbits['objid']), 1)
    # Convert to format for pyoorb, INCLUDING converting inclination, node, argperi to RADIANS
    oorbArray = np.column_stack((orbids, orbits['q'], orbits['e'], np.radians(orbits['i']),
                                    np.radians(orbits['node']), np.radians(orbits['argperi']),
                                    orbits['t_p'], elem_type, orbits['t_0'], epoch_type, orbits['H'], gval))
    return oorbArray


    orbits = pd.read_table('test_orbits.des', sep='\s*', engine='python')
    newcols = orbits.columns.values
    newcols[0] = 'objid'
    orbits.columns = newcols
    dt, t = dtime(t)
    print "Reading %d orbits required %f s" %(len(orbits['q']), dt)
    # convert orbit array to 'packed' form needed in oorb.
    oorbArray = pack_oorbArray(orbits)


def unpackOorbEphs(oorbephems, byObject=True):
    """
    Given oorb ephemeris array (shape = object / times / eph@time),
    Return a numpy array aranged with
     columns = ['delta', 'ra', 'dec', 'mag', 'time', 'timescale', 'dradt', 'ddecdt', 'phase', 'solarelon']
     as the second
    grouped either by object (i.e. length of ra array == length of times) (default)
    or grouped by time (i.e. length of ra array == number of objects) (if byObject not true).
    """
    # Python oorb ephem array format:
    #   [objid][time][ephemeris information @ that time]
    # 0 = distance (geocentric distance)
    # 1 = ra (deg)
    # 2 = dec (deg)
    # 3 = mag
    # 4 = ephem mjd
    # 5 = ephem mjd timescale
    # 6 = dra/dt (deg/day) sky motion
    # 7 = ddec/dt (deg/day) sky motion
    # 8 = phase angle (deg)
    # 9 = solar elongation angle (deg)
    # So usually we want to swap the axes at least, so that instead of all the ephemeris information @ a particular time
    # being the accessible bit of information, we have all the RA values over time for a single object ('byObject')
    # Alternatively, we may want all the RA values for all objects at one time.
    #     This is also an option, by setting 'byObject' to False.
    ephs = np.swapaxes(oorbephems, 2, 0)
    # oorbcols=['delta', 'ra', 'dec', 'magV', 'time', 'timescale', 'dradt', 'ddecdt', 'phase', 'solarelon']
    velocity = np.sqrt(ephs[6]**2 + ephs[7]**2)
    if byObject:
        ephs = np.swapaxes(ephs, 2, 1)
        velocity = np.swapaxes(velocity, 1, 0)
    # Create a numpy recarray. (not using a dAtaframe here, because the numpy recarray is just easier to swap around later).
    ephs = np.rec.fromarrays([ephs[0], ephs[1], ephs[2], ephs[3], ephs[4],
                              ephs[6], ephs[7], ephs[8], ephs[9], velocity],
                              names=['delta', 'ra', 'dec', 'magV', 'time', 'dradt',
                                     'ddecdt', 'phase', 'solarelon','velocity'])
    return ephs

if __name__=="__main__":

    parser = argparse.ArgumentParser(description='Python script to generate ephemerides using oorb')
    parser.add_argument('--orbitFile', type=str, default=None, help='Input orbit file (.des COMETARY format)')
    parser.add_argument('--obsCode', type=int, default=807, help='Observatory code for ephemerides')
    parser.add_argument('--ephTimesFile', type=str, default=None, help='Filename containing MJD-UTC times for ephemerides')
    parser.add_argument('--ephTimes', type=str, default=None, help='List of MJD-UTC times for ephemerides (if >1: in quotes, separated by spaces')
    parser.set_defaults()
    args = parser.parse_args()


    # Read orbits.
    orbits = pd.read_table(args.orbitFile, delim_whitespace=True)
    newcols = orbits.columns.values
    newcols[0] = 'objid'
    orbits.columns = newcols
    # convert orbit array to 'packed' form needed in oorb.
    oorbArray = pack_oorbArray(orbits)

    # Set up oorb.
    ephfile = os.path.join(os.getenv('OORB_DATA'), 'de430.dat')
    oo.pyoorb.oorb_init(ephemeris_fname=ephfile)

    # set observatory code
    obscode = args.obsCode

    # Set up dates to predict ephemerides.
    if args.ephTimes is not None:
        times = np.array(args.ephTimes.split(), float)
        obshistids = np.arange(0, len(times))
        # For pyoorb, we need to tag times with timescales;
        # 1= MJD_UTC, 2=UT1, 3=TT, 4=TAI
        ephTimes = np.array(zip(times, repeat(1, len(times))), dtype='double')

    elif args.ephTimesFile is not None:
        times = pd.read_table(args.ephTimesFile, sep='\s*', engine='python', names=['MJD-UTC'])
        ephTimes = np.array(zip(times['MJD-UTC'].as_matrix(), repeat(1, len(times))), dtype='double')

    else:
        raise Exception('Did not have any ephemeride times')

    print "Using times: ", ephTimes[0]

    # Generate ephemerides.
    oorbephs, err = oo.pyoorb.oorb_ephemeris(in_orbits = oorbArray, in_obscode=obscode, in_date_ephems=ephTimes)

    # Returned ephems contain a 3-D Fortran array of ephemerides, the axes are:
    #   [objid][time][ephemeris information element]
    # the ephemeris information elements are (in order):
    # distance, ra, dec, mag, ephem mjd, ephem mjd timescale, dradt(sky), ddecdt(sky)
    # per object, per date, 8 elements (array shape is OBJ(s)/DATE(s)/VALUES)
    # Note that ra/dec, dradt, etc. are all in DEGREES.
    # First: (to arrange ephems for easier later use)
    # Swap the order of the axes: DATE / Objs / values

    # Unpack ephemerides.
    ephs = unpackOorbEphs(oorbephs, byObject=True)

    # print results?
    print 'DiaID ObshistId ObjID MJD(UTC)    RA    Dec   MJD_UTC MagV FakeSNR'
    diacounter = 0
    for j, o in orbits.iterrows():
        obshistids = np.arange(0, len(ephTimes))
        for i, obshist in enumerate(obshistids):
            print diacounter, obshist, o.objid, ephs['ra'][j][i], ephs['dec'][j][i], ephs['time'][j][i], ephs['magV'][j][i], 5
            diacounter += 1

