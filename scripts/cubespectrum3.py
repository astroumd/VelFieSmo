#! /usr/bin/env python
#
#  Load a FITS cube , extract the spectrum at a (or reference) pixel
#  and operate and plot some and then more...
#
#  21-apr-2017  Peter Teuben    hackday at "SPECTROSCOPY TOOLS IN PYTHON WORKSHOP" STSCI
#  14-jan-2018  PJT             finally updated to the current astropy.specutils


import os, sys
import numpy as np

from specutils import Spectrum1D
from specutils.io import read_fits
from astropy.units import Quantity

if len(sys.argv) == 1:
    print("Usage: %s fitsfile [xpos ypos]" % sys.argv[0])
    sys.exit(1)

if len(sys.argv) > 3:
    fitsfile = sys.argv[1]
    pos = [int(sys.argv[2]),int(sys.argv[3])]
    s3 = Spectrum1D.read(fitsfile, format="cubetest1", pos=pos)
elif len(sys.argv) == 2:
    fitsfile = sys.argv[1]
    s3 = read_fits.read_fits_spectrum1d(fitsfile, format="cubetest1")



nfreq = len(s3.data)    
restfrq=Quantity(s3.wcs.wcs.restfrq,"Hz")
print("RESTFREQ",restfrq)
s3.rest_value = restfrq

### nooo!!!!   needs to be done in reader since we're immutable
s3.velocity_convention = 'relativistic'
#print("RESTFREQ",s3.rest_value(restfrq))


#  even though native units in Hz, it needs the rest
#s3.to_dispersion("MHz",  rest = restfrq)
f3 = s3.frequency
print(f3[0]," ... ", f3[nfreq-1])

#  km/s
v3 = s3.velocity
print(v3[0]," ... ", v3[nfreq-1])

#  cm
l3 = s3.wavelength
print(l3[0]," ... ", l3[nfreq-1])


# some EW work, make a EW box so it can be overplotted on spectrum

from specutils.analysis import equivalent_width

#s3.to_dispersion("km/s", rest = restfrq)
ew = equivalent_width(s3)

ipeak = s3.flux.argmax()
xpeak = v3[ipeak].value
ypeak = s3.flux[ipeak].value
dx    = ew.value
print("PJT EW",ipeak,xpeak,ypeak,dx)
rect = [(xpeak-0.5*dx, 0.0), dx, ypeak]


# moments around the peak
m = 5
x = v3[ipeak-m:ipeak+m]
y = s3.flux[ipeak-m:ipeak+m]
xmean = (x*y).sum() / y.sum()
xdisp = (x*x*y).sum() / y.sum() - xmean*xmean
print("MOMENTS:",xmean,xdisp)
ymodel = ypeak * np.exp(-0.5*(x-xmean)**2/xdisp)
# print(x,ymodel)

# some plotting
import matplotlib.pyplot as plt
import matplotlib.patches as patches

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
plt.scatter(v3, s3.flux)
plt.plot(v3, s3.flux)
plt.plot(x,ymodel,'r-')
plt.title('Spectrum of %s at [%d,%d]' % (fitsfile, s3.meta['xpos'], s3.meta['ypos']))
plt.xlabel("km/s")
plt.ylabel(s3.unit.name)
ax1.add_patch(patches.Rectangle(rect[0],rect[1],rect[2],hatch='/',fill=False))
plt.show()

#
