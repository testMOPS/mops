import sys
from astropy.time import Time

# Horizons returns:
#
#2453645.500000000 = A.D. 2005-Oct-02 00:00:00.0000 (CT)
# EC= 4.921563013560958E-01 QR= 7.847163798437903E-01 IN= 1.906064824393883E+01
# OM= 3.590186205111589E+02 W = 1.065326884835588E+02 Tp=  2453863.572011837270
# N = 5.131331815261347E-01 MA= 2.481000147641282E+02 TA= 2.087004973016340E+02
# A = 1.545192707794188E+00 AD= 2.305669035744587E+00 PR= 7.015722486105969E+02


infile = open(sys.argv[1])
lines = infile.readlines()

if len(sys.argv) > 2:
    objid = sys.argv[2]
else:
    objid = 'test'

if len(sys.argv) > 3:
    H = float(sys.argv[3])
else:
    H = 10.

# line 0
epoch = float(lines[0].split()[0])
oorb_epoch = Time(epoch, format='jd', scale='tt')

# line 1
l = lines[1].split()
e = l[1]
q = l[3]
inc = l[5]

# line 2
l = lines[2].split()
node = l[1]
argperi = l[4]
tperi = float(l[6])
t_p = Time(tperi, format='jd', scale='tt')

print '!!OID FORMAT q e i node argperi t_p H t_0 INDEX N_PAR MOID COMPCODE'
print '%s COM %s %s %s %s %s %.9f %f %.9f 1 6 1 HORIZONS' %(objid, q, e, inc, node, argperi, t_p.tt.mjd, H, oorb_epoch.tt.mjd)
