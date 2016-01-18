import os
import numpy as np
import matplotlib.pyplot as plt
from itertools import repeat
import pandas as pd
import pyoorb as oo


def setup_oorb():
    # Set up oorb.
    ephfile = os.path.join(os.getenv('OORB_DATA'), 'de405.dat')
    oo.pyoorb.oorb_init(ephemeris_fname=ephfile)

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
    try:
        test_objid = orbits['objid'][0]
    except KeyError:
        test_objid = orbits.objid.head(1)
    if ((isinstance(test_objid, float) == True) |
        (isinstance(test_objid, int) == True)):
        orbids = orbits['objid']
    else:
        orbids = np.arange(0, len(orbits['objid']), 1)
    # Convert to format for pyoorb, INCLUDING converting inclination, node, argperi to RADIANS
    oorbArray = np.column_stack((orbids, orbits['q'], orbits['e'], np.radians(orbits['i']),
                                    np.radians(orbits['node']), np.radians(orbits['argperi']),
                                    orbits['t_p'], elem_type, orbits['t_0'], epoch_type, orbits['H'], gval))
    return oorbArray


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

def read_JPL_dets(filename='n747_dets_new.txt'):
    jpl = pd.read_table(filename, delim_whitespace=True)
    newcols = jpl.columns.values
    newcols[0] = 'objid'
    jpl.columns = newcols
    return jpl

def transform_to_V(filterMags, filters):
    # S colors filter transformation (V to filter) that Peter used.
    # https://github.com/testMOPS/mops/commit/1e50ee9a90124107df61c0f2df695c03f4266882
    transform = {'u':-1.77,
                'g':-0.348,
                'r':0.21332,
                'i':0.37317,
                'z':0.349,
                'y':0.354,
                'w':0.16}
    magV = np.zeros(len(filterMags), float)
    filterlist = np.unique(filters)
    for f in filterlist:
        match = np.where(filters == f)
        magV[match] = filterMags[match] - transform[f]
    return magV

def read_orbits(orbitfile='../S0'):
    # Read NEO orbit file
    orbits = pd.read_table(orbitfile, delim_whitespace=True, skiprows=1)
    newcols = orbits.columns.values
    newcols[0] = 'objid'
    orbits.columns = newcols
    return orbits

def trim_orbits(orbits, jpl):
    uniqueObj = jpl.objid.unique()
    orbMatch = orbits.query('objid in @uniqueObj')
    # Hack to remove objects in JPL list which are not in ours.
    # Not needed anymore really, but can stay in here as a check.
    jplMatch = jpl.query('objid in @orbMatch.objid')
    return orbMatch, jplMatch



if __name__ == '__main__':

    # read jpl data.
    jpl = read_JPL_dets('n747_dets_new.txt')
    # add V mags to JPL data.
    magV = transform_to_V(jpl['mag'].as_matrix(), jpl['filter_id'].as_matrix())
    jpl['magV'] = magV


    # read neo orbits
    orbits = read_orbits('../S0')
    orbits, jpl = trim_orbits(orbits, jpl)

    # predict neo positions at all times in jpl

    times = jpl.epoch_mjd.unique()
    # For pyoorb, we need to tag times with timescales;
    # 1= MJD_UTC, 2=UT1, 3=TT, 4=TAI
    ephTimes = np.array(zip(times, repeat(1, len(times))), dtype='double')
    # Generate ephemerides.
    setup_oorb()
    obscode = 807
    oorbits = pack_oorbArray(orbits)
    oorbephs, err = oo.pyoorb.oorb_ephemeris(in_orbits = oorbits, in_obscode=obscode, in_date_ephems=ephTimes)
    ephs = unpackOorbEphs(oorbephs, byObject=False)
    # For each time, pull out the objects which were actually in the jpl data.
    radiff = []
    decdiff = []
    magdiff = []
    for ephi, time in zip(ephs, times):
        print time
        e = pd.DataFrame(ephi)
        e['objid'] = orbits.objid.as_matrix()
        j = jpl.query('epoch_mjd == @time').sort_values('objid')
        e = e.query('objid in @j.objid').sort_values('objid')
        radiff += list(e.ra.as_matrix() - j.ra_deg.as_matrix())
        decdiff += list(e.dec.as_matrix() - j.dec_deg.as_matrix())
        magdiff += list(e.magV.as_matrix() - j.magV.as_matrix())
    radiff = np.array(radiff).ravel()
    decdiff = np.array(decdiff).ravel()
    magdiff = np.array(magdiff).ravel()
    print radiff.shape
    print decdiff.shape

    f = open('offsets.txt', 'w')
    for dra, ddec, dmag in zip(radiff, decdiff, magdiff):
        print >>f, dra, ddec, dmag
    f.close()

    print np.min(radiff)

    print 'ra difference min/max (deg)', np.min(radiff), np.max(radiff)
    print 'dec difference min/max (deg)', np.min(decdiff), np.max(decdiff)
    print 'Vmag difference min/max', np.min(magdiff), np.max(magdiff)

    plt.figure()
    plt.plot(radiff, decdiff, 'k.')
    plt.xlabel('RA offsets (deg)')
    plt.ylabel('Dec offsets (deg)')

    plt.figure()
    plt.hist(radiff, bins=150)
    plt.xlabel('RA offset (deg)')

    plt.figure()
    plt.hist(decdiff, bins=150)
    plt.xlabel('Dec offset (deg)')

    plt.figure()
    plt.hist(magdiff, bins=150)
    plt.xlabel('Mag offset')

    plt.show()
    
