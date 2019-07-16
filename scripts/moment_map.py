##Moment map method
#created by Nate
#06/12/2019

##import all the functions we need

import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
from spectral_cube import SpectralCube
import astropy.units as u
from astropy.wcs import WCS

##get data from our hedaerr

filename= '../data/ngc6503.fits'
cube = fits.open(filename)
header = cube[0].header
wcs = WCS(header)
data = cube[0].data

peak_int=np.max(data,axis=0)
rms=np.sqrt(np.mean(data**2,axis=0))
SNR = peak_int/rms

data_masked = data.copy()
data_masked[:,SNR < 4] =  np.nan

scube = SpectralCube(data=data_masked,wcs=wcs)
##convert from frequency to velocity

scube = scube.with_spectral_unit(u.km/u.s,velocity_convention='radio',rest_value=header['RESTFREQ']*u.Hz)

mom0=scube.moment(order=0)
mom1=scube.moment(order=1)
mom2=scube.linewidth_fwhm()

##saving...
name=header['OBJECT']
mom0.write('../data/MomentMaps/'+name+'_mom0.fits',overwrite=True)
mom1.write('../data/MomentMaps/'+name+'_mom1.fits',overwrite=True)
mom2.write('../data/MomentMaps/'+name+'_mom2.fits',overwrite=True)

mom0=mom0.value
mom1=mom1.value
mom2=mom2.value

fig, (ax1,ax2,ax3)=plt.subplots(nrows=1,ncols=3)
fig.suptitle('')

#ax1 = fig.add_subplot(121)
ax1.imshow(mom0,origin='lower',vmin=0.0,vmax=1111110)##??
ax1.set_title('Integrated intensity')

#ax2 = fig.add_subplot(121)
ax2.imshow(mom1,origin='lower',vmin=-250.0,vmax=270)
ax2.set_title('velocity')

#ax3 = fig.add_subplot(121)
ax3.imshow(mom2,origin='lower',vmin=-250.0,vmax=270)
ax3.set_title('velocity dispersion')

plt.show()
