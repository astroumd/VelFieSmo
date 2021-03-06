#! /usr/bin/env python
#
#  example of plotting and analysis from mkgalcube by running a series through "make run4"
#
#  grep FITTING run4_*/run.log  | grep -v beam | awk -F: '{print $5}' | awk '{if (NF==12) print $0}' | grep -v "100 120" > plot1.tab
#
#  this produces a table with 12 columns and many rows that are plotted here

import sys
import numpy as np
import matplotlib.pyplot as plt

#                     # pick a table via the first command line argument
table1 = sys.argv[1]

vpick = 3             # pick which one:   1=slit 2=ring 3=shape

# get data
(m,re,b,vslit,rslit,mslit,vring,rring,mring,vshape,rshape,mshape) = np.loadtxt(table1).T
r0 = 120   # scaling value (see mkgalcube) - we keep them constant
v0 = 100   # scaling value (see mkgalcube) - we keep them constant

# get all uniq input parameters (m, re and beam)
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

# remainder are plots    ================================================================

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
plt.title(title + " $R_e/R_0=%g$" % (rpick*60/r0))
plt.xlabel('$B/R_1$')
plt.ylabel('$m_1$')
plt.grid()
#plt.legend()
plt.savefig('fig1a-%s.pdf' % title)
plt.show()

#                                            B/R0 vs. m1  for fixed Re
plt.figure(1)
for um1 in um:
    mask = ((m==um1) & (re==rpick))
    x = np.extract(mask, 60*b/r0)
    y = np.extract(mask, m1)
    i = x.argsort()
    x = x[i]
    y = y[i]
    plt.plot(x,y,'-o',label='m=%g'%um1)
plt.title(title + " $R_e/R_0=%g$" % (rpick*60/r0))
plt.xlabel('$B/R_0$')
plt.ylabel('$m_1$')
plt.grid()
#plt.legend()
plt.savefig('fig1b-%s.pdf' % title)
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
plt.title(title + " $R_e/R_0=%g$" % (rpick*60/r0))
plt.xlabel('$R_1/R_0$')
plt.ylabel('$m_1$')
plt.grid()
plt.legend()
plt.savefig('fig1c-%s.pdf' % title)
plt.show()

#                                            B/R0 vs. R1/R0 for fixed Re
plt.figure(3)
for um1 in um:
    mask = ((m==um1) & (re==rpick))
    x = np.extract(mask, 60*b/r0)
    y = np.extract(mask, r1/r0)
    i = x.argsort()
    x = x[i]
    y = y[i]
    plt.plot(x,y,'-o',label='m=%g'%um1)
plt.title(title + " $R_e/R_0=%g$" % (rpick*60/r0))
plt.xlabel('$B/R_0$')
plt.ylabel('$R_1/R_0$')
plt.grid()
plt.legend()
plt.savefig('fig1d-%s.pdf' % title)
plt.show()



#                                            B/R1 vs. R1/R0 for fixed Re
plt.figure(3)
for um1 in um:
    mask = ((m==um1) & (re==rpick))
    x = np.extract(mask, 60*b/r1)
    y = np.extract(mask, r1/r0)
    i = x.argsort()
    x = x[i]
    y = y[i]
    plt.plot(x,y,'-o',label='m=%g'%um1)
plt.title(title + " $R_e/R_0=%g$" % (rpick*60/r0))
plt.xlabel('$B/R_1$')
plt.ylabel('$R_1/R_0$')
plt.grid()
plt.legend()
plt.savefig('fig1e-%s.pdf' % title)
plt.show()



# now fix m

#                                            B/R1 vs. m1  for fixed m

plt.figure(5)
for ur1 in ur:
    mask = ((re==ur1) & (m==mpick))
    x = np.extract(mask, 60*b/r1)
    y = np.extract(mask, m1)
    i = x.argsort()
    x = x[i]
    y = y[i]
    plt.plot(x,y,'-o',label='$R_e/R_0=%g$'%(ur1*60/r0))
plt.title(title + " $m=%g$" % mpick)
plt.xlabel('$B/R_1$')
plt.ylabel('$m_1$')
plt.grid()
plt.legend()
plt.savefig('fig1f-%s.pdf' % title)
plt.show()

#                                            B/R0 vs. m1  for fixed m
plt.figure(5)
for ur1 in ur:
    mask = ((re==ur1) & (m==mpick))
    x = np.extract(mask, 60*b/r0)
    y = np.extract(mask, m1)
    i = x.argsort()
    x = x[i]
    y = y[i]
    plt.plot(x,y,'-o',label='$R_e/R_0=%g$'%(ur1*60/r0))
plt.title(title + " $m=%g$" % mpick)
plt.xlabel('$B/R_0$')
plt.ylabel('$m_1$')
plt.grid()
plt.legend()
plt.savefig('fig1g-%s.pdf' % title)
plt.show()



#                                           B/R1 vs. R1/R0 for fixed m
plt.figure(6)
for ur1 in ur:
    mask = ((re==ur1) & (m==mpick))
    x = np.extract(mask, 60*b/r1)
    y = np.extract(mask, r1/r0)
    i = x.argsort()
    x = x[i]
    y = y[i]
    plt.plot(x,y,'-o',label='$R_e/R_0=%g$'%(ur1*60/r0))    
plt.title(title + " $m=%g$" % mpick)
plt.xlabel('$B/R_1$')
plt.ylabel('$R_1/R_0$')
plt.grid()
plt.legend()
plt.savefig('fig1h-%s.pdf' % title)
plt.show()

#                                           B/R0 vs. R1/R0 for fixed m
plt.figure(6)
for ur1 in ur:
    mask = ((re==ur1) & (m==mpick))
    x = np.extract(mask, 60*b/r0)
    y = np.extract(mask, r1/r0)
    i = x.argsort()
    x = x[i]
    y = y[i]
    plt.plot(x,y,'-o',label='$R_e/R_0=%g$'%(ur1*60/r0))    
plt.title(title + " $m=%g$" % mpick)
plt.xlabel('$B/R_0$')
plt.ylabel('$R_1/R_0$')
plt.grid()
plt.legend()
plt.savefig('fig1i-%s.pdf' % title)
plt.show()


#                                            R1/R0 vs. m1  for fixed m
plt.figure(5)
for ur1 in ur:
    mask = ((re==ur1) & (m==mpick))
    x = np.extract(mask, r1/r0)
    y = np.extract(mask, m1)
    i = x.argsort()
    x = x[i]
    y = y[i]
    plt.plot(x,y,'-o',label='$R_e/R_0=%g$'%(ur1*60/r0))
plt.title(title + " $m=%g$" % mpick)
plt.xlabel('$R_1/R_0$')
plt.ylabel('$m_1$')
plt.grid()
plt.legend()
plt.savefig('fig1j-%s.pdf' % title)
plt.show()



#---

# all                         B/R0 vs. R1/R0
plt.figure(4)
mask = (re==rpick)
# mask = (m > 0)              # example of everything
x  = np.extract(mask, b*60/r0)
y1 = np.extract(mask, rslit/r0)
y2 = np.extract(mask, rring/r0)
y3 = np.extract(mask, rshape/r0)

plt.scatter(x,y1, label='slit')
plt.scatter(x,y2, label='ring')
plt.scatter(x,y3, label='shape')
plt.title('all %s $R_e/R_0$=%g' % (table1,(rpick*60/r0)))
plt.xlabel('$B/R_0$')
plt.ylabel('$R_1/R_0$')
plt.grid()
plt.legend()
plt.savefig('fig1k-All_%s.pdf' % table1)
plt.show()


# all                         B/R0 vs. m
plt.figure(4)
mask = (re==rpick)
# mask = (m > 0)              # example of everything
x  = np.extract(mask, b*60/r0)
y1 = np.extract(mask, mslit)
y2 = np.extract(mask, mring)
y3 = np.extract(mask, mshape)

plt.scatter(x,y1, label='slit')
plt.scatter(x,y2, label='ring')
plt.scatter(x,y3, label='shape')
plt.title('all %s $R_e/R_0$=%g' % (table1,(rpick*60/r0)))
plt.xlabel('$B/R_0$')
plt.ylabel('$m_1$')
plt.grid()
plt.legend()
plt.savefig('fig1l-All_%s.pdf' % table1)
plt.show()
