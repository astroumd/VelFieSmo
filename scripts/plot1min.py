#! /usr/bin/env python
#
#        brute force (3D) grid minimization for beam smearing correction deconvolution
#

import sys, math
import numpy as np
import matplotlib.pyplot as plt

# get the observations, command line arguments are    tablename  Re B R1 m1
# where Re,B,R1 are now in sky units
# tablename still is for a given inclination
table1 = sys.argv[1]
oRe    = float(sys.argv[2])
oB     = float(sys.argv[3])
oR1    = float(sys.argv[4])
om1    = float(sys.argv[5])

print("Table = %s" % table1)
print("Re    = %g" % oRe)
print("B     = %g" % oB)
print("R1    = %g" % oR1)
print("m1    = %g" % om1)

# get the grid data
(m,re,b,vslit,rslit,mslit,vring,rring,mring,vshape,rshape,mshape) = np.loadtxt(table1).T
v1 = vshape;  r1 = rshape;  m1 = mshape

# grid variables
g_rer1 = 60*re/r1       # obs
g_br1  = 60*b/r1        # obs
g_m1   = m1             # obs
g_r1r0 = r1/120         # model
g_m0   = m              # model

# the grid parameters to minimize to
g1 = oRe/oR1
g2 = oB/oR1
g3 = om1
print("Minimizing to: Re/R1= %g B/R1= %g m1= %g" % (g1,g2,g3))

d = np.sqrt( (g_rer1-g1)**2 + (g_br1-g2)**2 + (g_m1-g3)**2 )
idx = np.argsort(d)
nlist = 10
if len(idx) < nlist: nlist=len(idx)

# placeholder for the top "nlist"
g_points = np.random.rand(nlist, 3)

# list the "nlist' smallest distances to (g1,g2,g3)
print("idx, dist, R1/R0;   R0=, m0=")
for j in range(nlist):
    i = idx[j]
    r0 = oR1/g_r1r0[i]
    m0 = g_m0[i]
    g_points[j,0] = g_rer1[i]
    g_points[j,1] = g_br1[i]
    g_points[j,2] = g_m1[i]
    print(i,d[i],"[",g_rer1[i],g_br1[i],g_m1[i],"]",g_r1r0[i],"R0=",r0,"m0=",m0)
