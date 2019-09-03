#! /usr/bin/env python
#
#  example python script to compute gaussfit vs moments in N6503 cube
#  takes about 75" to execute

import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.modeling import models,fitting

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

freq_ax= np.linspace(1,float(header['NAXIS3']),header['NAXIS3'])
freq=(freq_ax-header['CRPIX3'])*header['CDELT3']+header['CRVAL3']

gauss_init = models.Gaussian1D(amplitude=np.nanmedian(data_masked),mean=header['RESTFREQ'],stddev=(np.max(freq)-np.min(freq))/20.)
fitter = fitting.LevMarLSQFitter()

intensity_gauss = np.zeros((data.shape[1],data.shape[2]))
velocity_gauss = np.zeros((data.shape[1],data.shape[2]))
vel_disp_gauss = np.zeros((data.shape[1],data.shape[2]))

for x in range(data.shape[2]):
    for y in range(data.shape[1]):
        if all(np.isnan(data_masked[:,y,x]))==True:
            intensity_gauss[y,x] = np.nan
            velocity_gauss[y,x] = np.nan
            vel_disp_gauss[y,x] = np.nan
        else:
           spec = data_masked[:,y,x] 
           gauss_fit = fitter(gauss_init,freq,spec)
           peak_amp = gauss_fit.amplitude.value
           mean_freq = gauss_fit.mean.value
           disp_freq = gauss_fit.stddev.value

           c = 2.9979E5 #Km/s
           rest_freq=header['RESTFREQ']#Hz
           velocity= c*(1-(mean_freq)/(rest_freq))#Km/s
           vel_disp_FWHM= c*(disp_freq)/(rest_freq)*2.355#Km/s

           intensity_gauss[y,x] = peak_amp #Jy/beam
           velocity_gauss[y,x]  = velocity
           vel_disp_gauss[y,x]  = vel_disp_FWHM

fig, (ax1,ax2,ax3)=plt.subplots(nrows=1,ncols=3)
fig.suptitle('Peak velocity Vs Moment velocity')
print("intensity: %g %g" % (intensity_gauss.min(),intensity_gauss.max()))
print("velocity:  %g %g" % (velocity_gauss.min(), velocity_gauss.max()))
print("vel_disp:  %g %g" % (vel_disp_gauss.min(), vel_disp_gauss.max()))

ax1.imshow(intensity_gauss,origin='lower')
ax2.imshow(velocity_gauss,origin='lower')
ax3.imshow(vel_disp_gauss,origin='lower')

plt.savefig("gaussian_map1.png")
plt.show()
