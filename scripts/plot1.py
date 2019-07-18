#! /usr/bin/env python
#
#  example of plotting and analysis from "make run4"
#
#  grep FITTING run4_*/run.log  | grep -v beam | awk -F: '{print $5}' | awk '{if (NF==12) print $0}' | grep -v "100 120" > plot1.tab
#

import numpy as np
import matplotlib.pyplot as plt

#                     # pick a table
table1 = 'plot1.tab'

vpick = 3             # pick which one:   1=slit 2=ring 3=shape

# get data
(m,re,b,vslit,rslit,mslit,vring,rring,mring,vshape,rshape,mshape) = np.loadtxt(table1).T
r0 = 120
v0 = 100

# get all uniq m's and re's
um = np.unique(m)
ur = np.unique(re)
ub = np.unique(b)
print("Unique m    :" , um)
print("Unique re   :" , ur)
print("Unique beam :" , ub)
# pick something in case you want to fix those (2,2,1) are the defaults in mkgalcube
mpick = 2
rpick = 2
bpick = 1

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

# remainder are plots    

#                                            B/R1 vs. m1  for fixed Re
plt.figure(1)
for um1 in um:
    mask = ((m==um1) & (re==rpick))
    x = np.extract(mask, 60*b/r1)
    y = np.extract(mask, m1)
    i = x.argsort()
    x = x[i]
    y = y[i]
    plt.plot(x,y,'-o',label='m=%g'%um1)
plt.title(title + " $R_e=%g$" % rpick)
plt.xlabel('$B/R_1$')
plt.ylabel('$m_1$')
plt.grid()
#plt.legend()
plt.savefig('fig1a-%s.pdf' % title)
plt.show()

#                                            R1/R0 vs. m1  for fixed Re
plt.figure(2)
for um1 in um:
    mask = ((m==um1) & (re==rpick))
    x = np.extract(mask, r1/r0)
    y = np.extract(mask, m1)
    i = x.argsort()
    x = x[i]
    y = y[i]
    plt.plot(x,y,'-o',label='m=%g'%um1)
plt.title(title + " $R_e=%g$" % rpick)
plt.xlabel('$R_1/R_0$')
plt.ylabel('$m_1$')
plt.grid()
plt.legend()
plt.savefig('fig1b-%s.pdf' % title)
plt.show()

#                                            B vs. R1 for fixed Re
plt.figure(3)
for um1 in um:
    mask = ((m==um1) & (re==rpick))
    x = np.extract(mask, 60*b)
    y = np.extract(mask, r1)
    i = x.argsort()
    x = x[i]
    y = y[i]
    plt.plot(x,y,'-o',label='m=%g'%um1)
plt.title(title + " $R_e=%g$" % rpick)
plt.xlabel('$B$')
plt.ylabel('$R_1$')
plt.grid()
plt.legend()
plt.savefig('fig1c-%s.pdf' % title)
plt.show()



# all
plt.figure(4)
mask = (re==rpick)
# mask = (m > 0)              # example of everything
x  = np.extract(mask, b)
y1 = np.extract(mask, rslit)
y2 = np.extract(mask, rring)
y3 = np.extract(mask, rshape)

plt.scatter(x,y1, label='slit')
plt.scatter(x,y2, label='ring')
plt.scatter(x,y3, label='shape')
plt.title('all %s $R_e$=%g' % (table1,rpick))
plt.xlabel('$B$')
plt.ylabel('$R_1$')
plt.grid()
plt.legend()
plt.savefig('fig1d-All_%s.pdf' % table1)
plt.show()
