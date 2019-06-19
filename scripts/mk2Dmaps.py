#script to make 2D maps of intensity, velocity, and velocity dispersion from a 3D data cube
# R. C. Levy
# file created: 2019-06-10
# based on peak_velocity.py, gaussian_maps.py, and moment_maps.py by N. Feleke and R. C. Levy
# Change log:
# 2019-06-10 - file created, RCL
# 2019-06-12 - added more documentation, RCL

#import modules
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS


#open the data cube
filename = '../data/ngc6503.fits' #change this line for each cube
cube = fits.open(filename)
header = cube[0].header
data = cube[0].data

#mask the cube based on the SNR first
signal = np.max(data,axis=0) #use the max of spectrum at each pixel as the signal
noise = np.sqrt(np.mean(data**2,axis=0)) #use the rms to estimate the noise
SNR = signal/noise #compute signal to noise ratio
data_masked = data.copy()
data_masked[:,SNR < 4]=np.nan #mask pixels with SNR < 4

#define a function to derive the peak intensity and velocity maps
def mkpeakmaps(data,header):
	#copied from peak_velocity.py from N. Feleke with minor modifications
	#peak intensity can be found directly from our data by selecting the maximum value from axis=0 
	peak_intensity = np.max(data,axis=0)
	#find channel number where the spectrum a maximum
	idx=np.argmax(data,axis=0)
	#getting the frequency axis
	freq=(idx-header['CRPIX3']+1)*header['CDELT3']+header['CRVAL3']
	#peak velocity corresponds to the velocity at the peak frequency and can be calculated using one of our three core formulas
	c = 2.9979E5 #Km/s
	rest_freq=header['RESTFREQ']#Hz
	velocity= c*(1-(freq)/(rest_freq))#Km/s
	velocity[idx==0]=np.nan
	#get min and max velocities from freq range (for plotting later)
	min_vel = c*(1-np.max(freq_axis/header['RESTFREQ']))
	max_vel = c*(1-np.min(freq_axis/header['RESTFREQ']))
	med_vel = np.mean([np.abs(min_vel),max_vel])
	#get new header for the peak velocity map
	header_2D = header.copy()
	header_2D.remove('HISTORY',remove_all=True)
	header_2D.remove('CTYPE3')
	header_2D.remove('CRVAL3')
	header_2D.remove('CDELT3')
	header_2D.remove('CROTA3')
	header_2D.remove('CRPIX3')
	header_2D.remove('NAXIS3')
	header_2D['NAXIS']=2
	header_2D['BUNIT']='km/s'
	#get the galaxy name from the fits header
	gal_name=header['OBJECT']
	#save the peak velocity map as a fits file
	hdu = fits.PrimaryHDU(data=velocity,header=header_2D)
	hdul = fits.HDUList([hdu])
	hdul.writeto('../data/PeakMaps/'+gal_name+'_peak_velocity.fits',overwrite=True)
	#now save the peak intensity map
	#update the BUNIT in the header, everything else can star the same
	header_2D['BUNIT']='Jy / beam'
	hdu = fits.PrimaryHDU(data=peak_intensity,header=header_2D)
	hdul = fits.HDUList([hdu])
	hdul.writeto('../data/PeakMaps/'+gal_name+'_peak_intensity.fits',overwrite=True)

	#plot all on one plot
	fig=plt.figure(figsize=(8,2))
	#plot peak intensity
	ax=fig.add_subplot(121,projection=wcs,slices=('x','y',1))
	im=ax.imshow(peak_intensity,origin='lower',cmap='magma',vmin=0.,vmax=np.nanmedian(peak_intensity)+3*np.nanstd(peak_intensity))
	cb=fig.colorbar(im,ax=ax,fraction=0.03, pad=0.1)
	cb.set_label('Intensity (Jy/beam)')
	ax.coords[0].set_major_formatter('hh:mm')
	ax.coords[0].set_ticklabel(exclude_overlapping=True)
	ax.set_xlabel('R.A. (J2000)')
	ax.set_ylabel('Dec. (J2000)')
	#plot velocity at peak
	ax=fig.add_subplot(122,projection=wcs,slices=('x','y',1))
	im=ax.imshow(velocity,origin='lower',cmap='RdYlBu_r',vmin=-med_vel,vmax=med_vel)
	cb=fig.colorbar(im,ax=ax,fraction=0.03, pad=0.1)
	cb.set_label('Velocity (km/s)')
	ax.coords[0].set_major_formatter('hh:mm')
	ax.coords[0].set_ticklabel(exclude_overlapping=True)
	ax.set_xlabel('R.A. (J2000)')
	ax.set_ylabel('Dec. (J2000)')
	#add a title
	plt.suptitle(gal_name)
	#adjust subplot spacing
	plt.subplots_adjust(wspace=1.0)
	#save the figure
	plt.savefig('../plots/'+gal_name+'_peakmaps.pdf',bbox_inches='tight')
	plt.close()
	#plt.show()

	#return the arrays containing the maps
	return peak_intensity, velocity


#define a function derive moment maps
def mkmomentmaps(data,header):
	import astropy.units as u
	from spectral_cube import SpectralCube
	#get world coordinate system (wcs) info
	wcs=WCS(header)
	#store data as a SpectralCube
	scube = SpectralCube(data=data,wcs=wcs)
	#convert frequency axis in cube to velocity
	scube = scube.with_spectral_unit(u.km/u.s, velocity_convention='radio', rest_value=header['RESTFREQ']*u.Hz)
	#get min and max velocities covered by frequency range (for plotting later)
	min_vel = np.min(scube.spectral_axis.value)
	max_vel = np.max(scube.spectral_axis.value)
	med_vel = np.mean([np.abs(min_vel),max_vel])
	#use spectral cube to easily get moment maps
	mom0 = scube.moment(order=0) #Jy/beam*km/s
	mom1 = scube.moment(order=1) #km/s
	mom2 = scube.linewidth_fwhm() #km/s
	#grab galaxy name from the fits header
	gal_name = header['OBJECT']
	#save as fits files
	mom0.write('../data/MomentMaps/'+gal_name+'_mom0.fits',overwrite=True)
	mom1.write('../data/MomentMaps/'+gal_name+'_mom1.fits',overwrite=True)
	mom2.write('../data/MomentMaps/'+gal_name+'_mom2.fits',overwrite=True)
	
	#just save the data values
	mom0 = mom0.value
	mom1 = mom1.value
	mom2 = mom2.value

	#plot all on one plot
	fig=plt.figure(figsize=(12,2))
	#plot mom0
	ax=fig.add_subplot(131,projection=wcs,slices=('x','y',1))
	im=ax.imshow(mom0,origin='lower',cmap='magma',vmin=0.,vmax=np.nanmedian(mom0)+3*np.nanstd(mom0))
	cb=fig.colorbar(im,ax=ax,fraction=0.03, pad=0.1)
	cb.set_label('Int. Intensity (Jy/beam km/s)')
	ax.coords[0].set_major_formatter('hh:mm')
	ax.coords[0].set_ticklabel(exclude_overlapping=True)
	ax.set_xlabel('R.A. (J2000)')
	ax.set_ylabel('Dec. (J2000)')
	#plot mom1
	ax=fig.add_subplot(132,projection=wcs,slices=('x','y',1))
	im=ax.imshow(mom1,origin='lower',cmap='RdYlBu_r',vmin=-med_vel,vmax=med_vel)
	cb=fig.colorbar(im,ax=ax,fraction=0.03, pad=0.1)
	cb.set_label('Velocity (km/s)')
	ax.coords[0].set_major_formatter('hh:mm')
	ax.coords[0].set_ticklabel(exclude_overlapping=True)
	ax.set_xlabel('R.A. (J2000)')
	ax.set_ylabel('Dec. (J2000)')
	#plot mom2
	ax=fig.add_subplot(133,projection=wcs,slices=('x','y',1))
	im=ax.imshow(mom2,origin='lower',vmin=0.,vmax=np.nanmedian(mom2))
	cb=fig.colorbar(im,ax=ax,fraction=0.03, pad=0.1)
	cb.set_label('FWHM Velocity Dispersion (km/s)')
	ax.coords[0].set_major_formatter('hh:mm')
	ax.coords[0].set_ticklabel(exclude_overlapping=True)
	ax.set_xlabel('R.A. (J2000)')
	ax.set_ylabel('Dec. (J2000)')
	#add a title
	plt.suptitle(gal_name)
	#adjust subplot spacing
	plt.subplots_adjust(wspace=1.0)
	#save the figure
	plt.savefig('../plots/'+gal_name+'_momentmaps.pdf',bbox_inches='tight')
	plt.close()
	#plt.show()

	#return the arrays containing the maps
	return mom0, mom1, mom2

#define function to fit a gaussian at every pixel
def mkgaussianmaps(data,header):
	from astropy.modeling import models,fitting
	wcs=WCS(header)
	#get frequency axis
	channels = np.linspace(1,header['NAXIS3'],header['NAXIS3'])
	freq_axis = (channels-header['CRPIX3'])*header['CDELT3']+header['CRVAL3'] #Hz
	c=2.9979E5 #km/s
	#get min and max velocities from freq range (for plotting later)
	min_vel = c*(1-np.max(freq_axis/header['RESTFREQ']))
	max_vel = c*(1-np.min(freq_axis/header['RESTFREQ']))
	med_vel = np.mean([np.abs(min_vel),max_vel])
	half_freq_range = (header['NAXIS3']-header['CRPIX3'])*header['CDELT3'] #Hz
	#set up for Gaussian fitting
	gauss_init = models.Gaussian1D(amplitude=np.nanmedian(data),mean=header['RESTFREQ'],stddev=half_freq_range/5)
	fitter = fitting.LevMarLSQFitter()
	#create empty arrays to store maps
	intensity_gauss = np.zeros((data.shape[1],data.shape[2]))
	velocity_gauss = np.zeros((data.shape[1],data.shape[2]))
	dispersion_gauss = np.zeros((data.shape[1],data.shape[2]))
	#loop over all pixels
	#This step will take a while
	for x in range(data.shape[2]):
		for y in range(data.shape[1]):
			if all(np.isnan(data[:,y,x]))==True:
				#don't waste time fitting NaNs
				intensity_gauss[y,x]=np.nan
				velocity_gauss[y,x]=np.nan
				dispersion_gauss[y,x]=np.nan
			else:
				spec = data[:,y,x]
				gauss_fit = fitter(gauss_init,freq_axis,spec)
				#grab the parameters
				peak_intensity = gauss_fit.amplitude.value
				mean_frequency = gauss_fit.mean.value #Hz
				linewidth_sigma = gauss_fit.stddev.value #Hz
				#convert from frequency to velocity
				velocity = c*(1-mean_frequency/header['RESTFREQ']) #km/s
				dispersion_FWHM = c*(linewidth_sigma/header['RESTFREQ'])*2.355 #km/s
				#add these values to the array
				intensity_gauss[y,x]=peak_intensity #Jy/beam
				velocity_gauss[y,x]=velocity #km/s
				dispersion_gauss[y,x]=dispersion_FWHM #km/s
	
	#get new header for the maps
	header_2D = header.copy()
	header_2D.remove('HISTORY',remove_all=True)
	header_2D.remove('CTYPE3')
	header_2D.remove('CRVAL3')
	header_2D.remove('CDELT3')
	header_2D.remove('CROTA3')
	header_2D.remove('CRPIX3')
	header_2D.remove('NAXIS3')
	header_2D['NAXIS']=2
	#get the galaxy name from the fits header
	gal_name=header['OBJECT']
	
	#save the intensity map as a fits file
	header_2D['BUNIT']='Jy / beam'
	hdu = fits.PrimaryHDU(data=intensity_gauss,header=header_2D)
	hdul = fits.HDUList([hdu])
	hdul.writeto('../data/GaussianMaps/'+gal_name+'_intensity_gauss.fits',overwrite=True)
	
	#save the velocity map as a fits file
	header_2D['BUNIT']='km / s'
	hdu = fits.PrimaryHDU(data=velocity_gauss,header=header_2D)
	hdul = fits.HDUList([hdu])
	hdul.writeto('../data/GaussianMaps/'+gal_name+'_velocity_gauss.fits',overwrite=True)
	
	#save the dispersion map as a fits file
	header_2D['BUNIT']='km / s'
	hdu = fits.PrimaryHDU(data=dispersion_gauss,header=header_2D)
	hdul = fits.HDUList([hdu])
	hdul.writeto('../data/GaussianMaps/'+gal_name+'_veldisp_FWHM_gauss.fits',overwrite=True)
	
	#plot all on one plot
	fig=plt.figure(figsize=(12,2))
	#plot intensity from gaussian fit
	#add a subplot with RA and Dec axes
	ax=fig.add_subplot(131,projection=wcs,slices=('x','y',1))
	#plot the intensity map
	im=ax.imshow(intensity_gauss,origin='lower',cmap='magma',vmin=0.,vmax=np.nanmax(data))
	#add a colorbar
	cb=fig.colorbar(im,ax=ax,fraction=0.03, pad=0.1)
	#label the colorbar
	cb.set_label('Intensity (Jy/beam)')
	#display only hh:mm on the x-axis tick labels
	ax.coords[0].set_major_formatter('hh:mm')
	#force the tick labels to not overlap one another
	ax.coords[0].set_ticklabel(exclude_overlapping=True)
	#label the x- and y-axes
	ax.set_xlabel('R.A. (J2000)')
	ax.set_ylabel('Dec. (J2000)')
	#plot velocity from gaussian fit
	ax=fig.add_subplot(132,projection=wcs,slices=('x','y',1))
	im=ax.imshow(velocity_gauss,origin='lower',cmap='RdYlBu_r',vmin=-med_vel,vmax=med_vel)
	cb=fig.colorbar(im,ax=ax,fraction=0.03, pad=0.1)
	cb.set_label('Velocity (km/s)')
	ax.coords[0].set_major_formatter('hh:mm')
	ax.coords[0].set_ticklabel(exclude_overlapping=True)
	ax.set_xlabel('R.A. (J2000)')
	ax.set_ylabel('Dec. (J2000)')
	#plot vel disp from gaussian fit
	ax=fig.add_subplot(133,projection=wcs,slices=('x','y',1))
	im=ax.imshow(dispersion_gauss,origin='lower',vmin=0.,vmax=2*np.nanmedian(dispersion_gauss))
	cb=fig.colorbar(im,ax=ax,fraction=0.03, pad=0.1)
	cb.set_label('FWHM Velocity Dispersion (km/s)')
	ax.coords[0].set_major_formatter('hh:mm')
	ax.coords[0].set_ticklabel(exclude_overlapping=True)
	ax.set_xlabel('R.A. (J2000)')
	ax.set_ylabel('Dec. (J2000)')
	#add a title
	plt.suptitle(gal_name)
	#adjust subplot spacing
	plt.subplots_adjust(wspace=1.0)
	#save the figure
	plt.savefig('../plots/'+gal_name+'_gaussianmaps.pdf',bbox_inches='tight')
	plt.close()
	#plt.show()

	#return the arrays containing the maps
	return intensity_gauss, velocity_gauss, dispersion_gauss



#run the above functions
#make the peak maps
intensity_peak, velocity_peak = mkpeakmaps(data_masked,header)
#make the moment maps
intensity_mom0, velocity_mom1, dispersion_mom2 = mkmomentmaps(data_masked,header)
#make the maps from the gaussian fits
intensity_gauss, velocity_gauss, dispersion_gauss = mkgaussianmaps(data_masked,header)




