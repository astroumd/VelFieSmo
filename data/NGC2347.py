#! /usr/bin/env python
#
#  compare rotation curves for NGC2347 (see also the generic rotcmp1.py)
#

import sys
import math
import glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, join

#                     # pick a table
table1 = 'run.rcmod'
table2  = 'NGC2347.rotcurtab'
table2a = 'NGC2347.rotcurtab2'
#table3 = 'run.rcslit'
table3 = 'tony.tab'
table4 = '../models/ngc5016/run.rcring'
table5 = 'barolo_2dfit/_2dtrm.txt'
table6 = 'barolo_3dfit/ringlog1.txt'

# edge_pydb/edge_pydb/dat_prof/bbarolo

t10 = 'edge_pydb/edge_pydb/dat_prof/rotcur_levy/HA_natv/NGC5016.Ha.RC.txt'
t11 = 'edge_pydb/edge_pydb/dat_prof/rotcur_levy/CO_natv/NGC5016.CO.RC.txt'
t12 = 'edge_pydb/edge_pydb/dat_prof/rotcur_levy/HA_smo6/NGC5016.Ha.RC.txt'
t13 = 'edge_pydb/edge_pydb/dat_prof/rotcur_levy/CO_smo6/NGC5016.CO.RC.txt'
t14 = 'edge_pydb/edge_pydb/dat_prof/radprof/build/rproftxt/NGC5016.smo7_smob1.rprof.txt'
t15 = 'edge_pydb/edge_pydb/dat_prof/radprof/build/rproftxt/NGC5016.de20_smob1.rprof.txt'
t16 = 'edge_pydb/edge_pydb/dat_prof/radprof/build/rproftxt/NGC5016.de20_smob0.rprof.txt'
t17 = 'edge_pydb/edge_pydb/dat_prof/radprof/build/rproftxt/NGC5016.smo7_smob0.rprof.txt'
t18 = 'edge_pydb/edge_pydb/dat_prof/rotcur_leung/EDGE_CO_vrot/NGC2347_co_vrot.txt'
t19 = 'edge_pydb/edge_pydb/dat_prof/rotcur_leung/EDGE_CO_vrot_bsc/NGC2347_co_vrot.txt'

(r10,v10,dv10,z,z,z,z) = np.loadtxt(t10,comments='%').T
(r11,v11,dv11,z,z,z,z) = np.loadtxt(t11,comments='%').T
(r12,v12,dv12,z,z,z,z) = np.loadtxt(t12,comments='%').T
(r13,v13,dv13,z,z,z,z) = np.loadtxt(t13,comments='%').T
(r18,v18,dv18) = np.loadtxt(t18).T
(r19,v19,dv19) = np.loadtxt(t19).T


def vmodel(r, v0, r0, m0):
    x = r/r0
    v = v0 * x / pow( (1+pow(x,m0)) , 1/m0)
    return v

# observed fit (re/ is needed to deduce r0)
# from tabnllsqfit on tables
# beam=4.5
if False:
    v1 = 290
    r1 =  6.32
    m1 =  2.32
    re = 8.5

    # nearest pt in grid
    v0 = v1
    r0 = 3.0
    m0 = 1.5
else:
    v1 = 292
    r1 =  3.5
    m1 =  2.50
    re = 8.5

    # nearest pt in grid
    # d=0.41
    v0 = v1
    r0 = 1.7
    m0 = 2.0
    # 0.40
    r0 = 1.25
    m0 = 3.0
    # 0.39
    r0 = 1.72
    m0 = 2.5
    # 0.35 - better fit than 0.32 !
    r0 = 1.46
    m0 = 2.75
    # 0.321
    r0 = 1.448
    m0 = 2.5
    # 0.299
    r0 = 1.332
    m0 = 2.6
    # 0.283
    r0 = 1.247
    m0 = 2.7


# this means, for reference
#  Re/R0 = 1.75
#  B/R0  = 0.62

# inc factor correction from 39.9 to 46
ifac = math.sin(39.9/57.266)/math.sin(46/57.266)
ifac = 1.0

v10 = v10 * ifac
v11 = v11 * ifac
v12 = v12 * ifac
v13 = v13 * ifac
v18 = v18 * ifac
v19 = v19 * ifac

rmod  = np.linspace(0.5,21,72)
vmod1 = vmodel(rmod,v1,r1,m1)
vmod0 = vmodel(rmod,v0,r0,m0)


(r1,vsys1,dvsys1,vrot1,dvrot1,pa1,dpa1,inc1,dinc1,xpos1,dxpos1,ypos1,dypos1,npt1,sigma1) = np.loadtxt(table2).T
(r2,vsys2,dvsys2,vrot2,dvrot2,pa2,dpa2,inc1,dinc2,xpos2,dxpos2,ypos2,dypos2,npt2,sigma2) = np.loadtxt(table2a).T
#(rslit,vslit,eslit) = np.loadtxt(table3).T
#(rring,vring,ering) = np.loadtxt(table4).T
#(r2p,r2,vsys2,vrot2,vexp2,pa2,inc2,xpos2,ypos2) = np.loadtxt(table5).T
#(r3p,r3,vrot3,disp3,inc3,pa3,z03p,z03,e3,xpos3,ypos3,vsys3,vrad3) = np.loadtxt(table6).T

# old grep&awk tony table
# (r3,v3,dv3) = np.loadtxt(table3).T
# new official style
bfiles = glob.glob('edge_pydb/edge_pydb/dat_prof/bbarolo/*csv')
r4 = np.array([])
v4 = np.array([])
i4 = np.array([])
for f in bfiles:
    print(f)
    table = Table.read(f, format = 'ascii.ecsv')
    galrows = table[table['bbName']=='NGC2347']
    rt = galrows['bbRad']
    vt = galrows['bbVrot']
    it = galrows['bbInc']
    tmp = galrows['bbPA']
    r4 = np.append(r4,rt.data)
    v4 = np.append(v4,vt.data)
    i4 = np.append(i4,it.data)

# they're all 39.9
v4 = v4 * ifac

plt.figure(1)
plt.subplot(1,1,1)
plt.plot(rmod,vmod1,  '-',label='observed 2D model fit')
plt.plot(rmod,vmod0, '--',label='deconvolved 2D model')
plt.plot(r1,vrot1,   '-o',c='black',label='rotcur-mom1')
plt.plot(r2,vrot2,   '-o',c='black',label='rotcur-cgrad-vcen')
if False:
    plt.scatter(r4,v4   ,c='g',label='BBarolo')
    plt.scatter(r10,v10 ,label='levy-10 inc=39.9')
    plt.scatter(r11,v11 ,label='levy-11')
    plt.scatter(r12,v12 ,label='levy-12')
    plt.scatter(r13,v13 ,label='levy-13')
plt.scatter(r18,v18 ,c='r',label='leung')
plt.scatter(r19,v19 ,c='g',label='leung-BSC')



#plt.plot(rring,vring,'-o',label='rotcur')
#plt.plot(rslit,vslit,'-o',label='slit')
#plt.plot(r2,vrot2,   '-o',label='barolo_2dfit')
#plt.plot(r3,vrot3,   '-o',label='barolo_3dfit')
plt.title('NGC2347 rotation curves')
plt.xlabel('Radius')
plt.ylabel('Velocity')
plt.grid()
#plt.legend(loc='lower right',fontsize = 'small')
plt.legend(loc='lower right',prop={'size':8})
plt.savefig('NGC2347.pdf')
plt.show()
