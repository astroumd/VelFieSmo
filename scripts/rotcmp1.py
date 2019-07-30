#! /usr/bin/env python
#
#  compare rotation curves
#

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

#
do_barolo = True

#                     # pick a table
table0 = 'run.rotmod'
table1 = 'run.rcmod'
table2 = 'run.rotcurtab'
table3 = 'run.rcslit'
table4 = 'run.rcring'
table5 = 'barolo_2dfit/_2dtrm.txt'
table6 = 'barolo_3dfit/ringlog1.txt'
table7 = 'run.rcshape'

(r0, v0) = np.loadtxt(table0).T
(rmod, vmod, zmod) = np.loadtxt(table1).T
(r1,vsys1,dvsys1,vrot1,dvrot1,pa1,dpa1,inc1,dinc1,xpos1,dxpos1,ypos1,dypos1,npt1,sigma1) = np.loadtxt(table2).T
(rslit,vslit,eslit) = np.loadtxt(table3).T
(rring,vring,ering) = np.loadtxt(table4).T
if do_barolo:
    (r2p,r2,vsys2,vrot2,vexp2,pa2,inc2,xpos2,ypos2) = np.loadtxt(table5).T
    (r3p,r3,vrot3,disp3,inc3,pa3,z03p,z03,e3,xpos3,ypos3,vsys3,vrad3) = np.loadtxt(table6).T
(rshape,vshape) = np.loadtxt(table7).T

fmod = interp1d(r0, v0, kind = 'cubic')

#plt.figure(1,figsize=(11, 8.5))    # landscape
plt.figure(1,figsize=(8.5,11))      # portrait

plt.subplot(2,1,1)
#plt.plot(rmod,vmod,  '-o',label='model')
plt.plot(r0,v0,        '-', c='black',label='model')
plt.plot(rshape,vshape,'--',c='black',label='shape')
plt.plot(rring,vring,  '-o',c='red',  label='rotcur')
plt.plot(rslit,vslit,  '-o',c='green',label='slit')
if do_barolo:
    plt.plot(r2,vrot2,   '-o',c='blue',label='barolo_2dfit')
    plt.plot(r3,vrot3,   '-o',c='magenta',label='barolo_3dfit')
plt.xlim(0,800)    
#plt.xlabel('Radius')
plt.ylabel('Velocity')
plt.grid()
plt.legend(loc='lower right',fontsize = 'small')

vmax = 10

plt.subplot(4,1,3)
plt.plot(rmod,vmod-fmod(rmod),      '-', c='black',label='model')
plt.plot(rshape,vshape-fmod(rshape),'--',c='black',label='shape')
plt.plot(rring,vring-fmod(rring),   '-o',c='red',  label='rotcur')
plt.plot(rslit,vslit-fmod(rslit),   '-o',c='green',label='slit')
if do_barolo:
    plt.plot(r2,vrot2-fmod(r2),     '-o',c='blue',   label='barolo_2dfit')
    plt.plot(r3,vrot3-fmod(r3),     '-o',c='magenta',label='barolo_3dfit')
plt.xlim(0,800)
plt.ylim(-vmax,vmax)
#plt.xlabel('Radius')
plt.ylabel('Velocity difference')
plt.grid()
#plt.legend(fontsize = 'x-small')

vmax = 1

plt.subplot(4,1,4)
plt.plot(rmod,vmod-fmod(rmod),      '-', c='black',label='model')
plt.plot(rshape,vshape-fmod(rshape),'--',c='black',label='shape')
plt.plot(rring,vring-fmod(rring),   '-o',c='red',  label='rotcur')
plt.plot(rslit,vslit-fmod(rslit),   '-o',c='green',label='slit')
if do_barolo:
    plt.plot(r2,vrot2-fmod(r2),     '-o',c='blue',   label='barolo_2dfit')
    plt.plot(r3,vrot3-fmod(r3),     '-o',c='magenta',label='barolo_3dfit')
plt.xlim(0,800)
plt.ylim(-vmax,vmax)
plt.xlabel('Radius')
plt.ylabel('Velocity difference')
plt.grid()
#plt.legend(fontsize = 'x-small')



#plt.tight_layout(h_pad=0, w_pad=0, pad=0)

plt.savefig('rotcmp1.pdf')
plt.show()
