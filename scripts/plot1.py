#! /usr/bin/env python
#
#  example of plotting and analysis
#
#  grep FITTING run_*/run.log  | grep -v beam | awk -F: '{print $4}' | awk '{if (NF==11) print $0}' | grep -v "100 120" > plot1.tab
#

import numpy as np
import matplotlib.pyplot as plt

#                     # pick a table
table1 = 'plot1.tab'

pick = 3               # pick which one:   1=slit 2=ring 3=shape

# get data
(m,b,vslit,rslit,mslit,vring,rring,mring,vshape,rshape,mshape) = np.loadtxt(table1).T
r0 = 120
v0 = 100

# get all uniq m's
um = np.unique(m)

# pick one of the three fittings
if pick==1:
    v1 = vslit;   r1 = rslit;   m1 = mslit
    title = "Slit_%s" % table1
elif pick == 2:
    v1 = vring;   r1 = rring;   m1 = mring
    title = "Ring_%s" % table1
elif pick == 3:    
    v1 = vshape;  r1 = rshape;  m1 = mshape
    title = "Shape_%s" % table1

# remainder are plots    

#                                            B/R1 vs. m
plt.figure(1)
for um1 in um:
    x = np.extract(m==um1, 60*b/r1)
    y = np.extract(m==um1, m1)
    i = x.argsort()
    x = x[i]
    y = y[i]
    plt.plot(x,y,'-o',label='m=%g'%um1)
plt.title(title)    
plt.xlabel('$B/R_1$')
plt.ylabel('$m_1$')
plt.grid()
#plt.legend()
plt.savefig('fig1a-%s.pdf' % title)
plt.show()

#                                            R1/R0 vs. m1
plt.figure(2)
for um1 in um:
    x = np.extract(m==um1, r1/r0)
    y = np.extract(m==um1, m1)
    i = x.argsort()
    x = x[i]
    y = y[i]
    plt.plot(x,y,'-o',label='m=%g'%um1)
plt.title(title)
plt.xlabel('$R_1/R_0$')
plt.ylabel('$m_1$')
plt.grid()
plt.legend()
plt.savefig('fig1b-%s.pdf' % title)
plt.show()

#                                            B vs. R1
plt.figure(3)
for um1 in um:
    x = np.extract(m==um1, 60*b)
    y = np.extract(m==um1, r1)
    i = x.argsort()
    x = x[i]
    y = y[i]
    plt.plot(x,y,'-o',label='m=%g'%um1)
plt.title(title)
plt.xlabel('$B$')
plt.ylabel('$R_1$')
plt.grid()
plt.legend()
plt.savefig('fig1c-%s.pdf' % title)
plt.show()



# all
if True:
    plt.figure(4)
    plt.scatter(b, rslit,  label='slit')
    plt.scatter(b, rring,  label='ring')
    plt.scatter(b, rshape, label='shape')
    plt.title('all %s' % table1)
    plt.xlabel('$B$')
    plt.ylabel('$R_1$')
    plt.grid()
    plt.legend()
    plt.savefig('fig1d-All_%s.pdf' % table1)
    plt.show()
