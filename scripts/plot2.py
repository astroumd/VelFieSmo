#! /usr/bin/env python
#
#  example of plotting and analysis from mkgalcube by running a series through "make run4"
#
#  grep FITTING run4_*/run.log  | grep -v beam | awk -F: '{print $5}' | awk '{if (NF==12) print $0}' | grep -v "100 120" > plot1.tab
#
#  this produces a table with 12 columns and many rows that are plotted here
#
#  this version makes 12 subplots with various model and observation style relation ships for a graphical interpolation
#  to turn (m1,R1,B) -> (m0,R0) where we ignore V1 (or assume V0=V1)

import sys
import numpy as np
import matplotlib.pyplot as plt

#                     # pick a table via the first command line argument
table1 = sys.argv[1]

# set defaults for remaining arguments
# pick something in case you want to fix those (2,2,1) are the defaults in mkgalcube
mpick = 2
rpick = 1
bpick = 1
if len(sys.argv) > 2: mpick = float(sys.argv[2])
if len(sys.argv) > 3: rpick = float(sys.argv[3])
if len(sys.argv) > 4: bpick = float(sys.argv[4])
    
vpick = 3           # pick which one:   1=slit 2=ring 3=shape

# get data
(m,re,b,vslit,rslit,mslit,vring,rring,mring,vshape,rshape,mshape) = np.loadtxt(table1).T
r0 =  60   # scaling value (see mkgalcube) - we keep them constant
v0 = 100   # scaling value (see mkgalcube) - we keep them constant

# get all uniq input parameters (m, re and beam)
um = np.unique(m)
ur = np.unique(re)
ub = np.unique(b)
print("Unique m    :" , um)
print("Unique re   :" , ur)
print("Unique beam :" , ub)

# pick one of the three fitting methods
if vpick==1:
    v1 = vslit;   r1 = rslit;   m1 = mslit
    title = "Slit_%s" % table1
elif vpick == 2:
    v1 = vring;   r1 = rring;   m1 = mring
    title = "Ring_%s" % table1
elif vpick == 3:    
    v1 = vshape;  r1 = rshape;  m1 = mshape
    title = "Shape_%s" % table1

# remainder are plots    ================================================================

#                                            B/R1 vs. m1  for fixed Re
#plt.figure(1,figsize=(11, 8.5))
plt.figure(1,figsize=(8.5, 11))

plt.subplot(4,3,4)
for um1 in um:
    mask = ((m==um1) & (re==rpick))
    x = np.extract(mask, 60*b/r1)
    y = np.extract(mask, m1)
    i = x.argsort()
    x = x[i]
    y = y[i]
    plt.plot(x,y,'-o',label='m=%g'%um1)
#plt.title("$R_e/R_0=%g$" % (rpick*60/r0))
plt.xlabel('$B/R_1$')
plt.ylabel('$m_1$')
plt.grid()
#plt.legend()

#                                            B/R0 vs. m1  for fixed Re
plt.subplot(4,3,1)
for um1 in um:
    mask = ((m==um1) & (re==rpick))
    x = np.extract(mask, 60*b/r0)
    y = np.extract(mask, m1)
    i = x.argsort()
    x = x[i]
    y = y[i]
    plt.plot(x,y,'-o',label='m=%g'%um1)
plt.title("$R_e/R_0=%g$" % (rpick*60/r0))
plt.xlabel('$B/R_0$')
plt.ylabel('$m_1$')
plt.grid()
#plt.legend()


#                                            R1/R0 vs. m1  for fixed Re
plt.subplot(4,3,3)
for um1 in um:
    mask = ((m==um1) & (re==rpick))
    x = np.extract(mask, r1/r0)
    y = np.extract(mask, m1)
    i = x.argsort()
    x = x[i]
    y = y[i]
    plt.plot(x,y,'-o',label='m=%g'%um1)
plt.title("$R_e/R_0=%g$" % (rpick*60/r0))
plt.xlabel('$R_1/R_0$')
plt.ylabel('$m_1$')
plt.grid()
#plt.legend()

#                                            B/R0 vs. R1/R0 for fixed Re
plt.subplot(4,3,7)
for um1 in um:
    mask = ((m==um1) & (re==rpick))
    x = np.extract(mask, 60*b/r0)
    y = np.extract(mask, r1/r0)
    i = x.argsort()
    x = x[i]
    y = y[i]
    plt.plot(x,y,'-o',label='m=%g'%um1)
#plt.title("$R_e/R_0=%g$" % (rpick*60/r0))
plt.xlabel('$B/R_0$')
plt.ylabel('$R_1/R_0$')
plt.grid()
#plt.legend()



#                                            B/R1 vs. R1/R0 for fixed Re
plt.subplot(4,3,10)
for um1 in um:
    mask = ((m==um1) & (re==rpick))
    x = np.extract(mask, 60*b/r1)
    y = np.extract(mask, r1/r0)
    i = x.argsort()
    x = x[i]
    y = y[i]
    plt.plot(x,y,'-o',label='m=%g'%um1)
#plt.title("$R_e/R_0=%g$" % (rpick*60/r0))
plt.xlabel('$B/R_1$')
plt.ylabel('$R_1/R_0$')
plt.grid()
#plt.legend()



# now fix m

#                                            B/R1 vs. m1  for fixed m
plt.subplot(4,3,5)
for ur1 in ur:
    mask = ((re==ur1) & (m==mpick))
    x = np.extract(mask, 60*b/r1)
    y = np.extract(mask, m1)
    i = x.argsort()
    x = x[i]
    y = y[i]
    plt.plot(x,y,'-o',label='$R_e/R_0=%g$'%(ur1*60/r0))
#plt.title("$m=%g$" % mpick)
plt.xlabel('$B/R_1$')
plt.ylabel('$m_1$')
plt.grid()
#plt.legend()

#                                            B/R0 vs. m1  for fixed m
plt.subplot(4,3,2)
for ur1 in ur:
    mask = ((re==ur1) & (m==mpick))
    x = np.extract(mask, 60*b/r0)
    y = np.extract(mask, m1)
    i = x.argsort()
    x = x[i]
    y = y[i]
    plt.plot(x,y,'-o',label='$R_e/R_0=%g$'%(ur1*60/r0))
plt.title("$m=%g$" % mpick)
plt.xlabel('$B/R_0$')
plt.ylabel('$m_1$')
plt.grid()
#plt.legend()



#                                           B/R1 vs. R1/R0 for fixed m
plt.subplot(4,3,11)
for ur1 in ur:
    mask = ((re==ur1) & (m==mpick))
    x = np.extract(mask, 60*b/r1)
    y = np.extract(mask, r1/r0)
    i = x.argsort()
    x = x[i]
    y = y[i]
    plt.plot(x,y,'-o',label='$R_e/R_0=%g$'%(ur1*60/r0))    
#plt.title("$m=%g$" % mpick)
plt.xlabel('$B/R_1$')
plt.ylabel('$R_1/R_0$')
plt.grid()
#plt.legend()

#                                           B/R0 vs. R1/R0 for fixed m
plt.subplot(4,3,8)
for ur1 in ur:
    mask = ((re==ur1) & (m==mpick))
    x = np.extract(mask, 60*b/r0)
    y = np.extract(mask, r1/r0)
    i = x.argsort()
    x = x[i]
    y = y[i]
    plt.plot(x,y,'-o',label='$R_e/R_0=%g$'%(ur1*60/r0))    
#plt.title("$m=%g$" % mpick)
plt.xlabel('$B/R_0$')
plt.ylabel('$R_1/R_0$')
plt.grid()
#plt.legend()


#                                            R1/R0 vs. m1  for fixed m
plt.subplot(4,3,6)
for ur1 in ur:
    mask = ((re==ur1) & (m==mpick))
    x = np.extract(mask, r1/r0)
    y = np.extract(mask, m1)
    i = x.argsort()
    x = x[i]
    y = y[i]
    plt.plot(x,y,'-o',label='$R_e/R_0=%g$'%(ur1*60/r0))
plt.title("$m=%g$" % mpick)
plt.xlabel('$R_1/R_0$')
plt.ylabel('$m_1$')
plt.grid()
#plt.legend()



#---

# all                         B/R0 vs. R1/R0
plt.subplot(4,3,9)
mask = (re==rpick)
mask = (m > 0)              # example of everything
x  = np.extract(mask, b*60/r0)
y1 = np.extract(mask, rslit/r0)
y2 = np.extract(mask, rring/r0)
y3 = np.extract(mask, rshape/r0)

#plt.scatter(x,y1, label='slit')
#plt.scatter(x,y2, label='ring')
plt.scatter(x,y3, label='shape')
#plt.title('all %s $R_e/R_0$=%g' % (table1,(rpick*60/r0)))
plt.xlabel('$B/R_0$')
plt.ylabel('$R_1/R_0$')
plt.grid()
#plt.legend()


# all                         B/R0 vs. m
plt.subplot(4,3,12)
mask = (re==rpick) 
mask = (m > 0)              # example of everything
x  = np.extract(mask, b*60/r0)
y1 = np.extract(mask, mslit)
y2 = np.extract(mask, mring)
y3 = np.extract(mask, mshape)

#plt.scatter(x,y1, label='slit')
#plt.scatter(x,y2, label='ring')
plt.scatter(x,y3, label='shape')
#plt.title('all %s $R_e/R_0$=%g' % (table1,(rpick*60/r0)))
plt.xlabel('$B/R_0$')
plt.ylabel('$m_1$')
plt.grid()

#plt.legend()

# pad=0.4, w_pad=0.5, h_pad=1.0)
plt.tight_layout(h_pad=0, w_pad=0, pad=0)
plt.savefig('plot2.pdf')
plt.show()
