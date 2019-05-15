def beam_smearing_correction(path_to_veldispmap,path_to_rotationcurve,spectral_resolution_FWHM,chan_width,PA,inc,RA,Dec,Vsys):
	
	r'''
	Simulate beam smearing and remove from velocity dispersion map using method described in Appendix B of Levy et al. (2018).
	This code: 1) Creates a model model data cube with spectral resolution set by the instrument and an input rotation curve, 
	2) Convolves the data cube to the desired spatial resolution, 
	3) Derives a moment 2 map representing the velocity dispersion due to beam smearing, 
	4) Removes the beam-smeared moment 2 from the input moment 2 in quadrature to derive the beam smearing-corrected velocity dispersion map.
	5) Saves the simulated and beam smearing corrected mom2 maps.
	6) Plots a comparison of the velocity dispersion maps.

	Parameters
	----------
	path_to_veldispmap : str
		Path+filename to the input velocity dispersion fits file, assumes this is a sigma not FWHM, in km/s
	path_to_rotationcurve : str
		Path+filename to the (previously derived) rotation curve, should be a txt file with columns for radius (arcsec), Vrot (km/s), and eVrot (km/s)
	spectral_resolution_FWHM : float
		The FWHM spectral resolution of the data (assumed Gaussian), in km/s
	chan_width : float
		The channel width, in km/s
	PA : float
		The position angle of the approaching side of the major axis measured counter clockwise from North, in degrees
	inc : float
		The inclination angle of the galaxy along the line of sight, in degrees
	RA : float
		The J2000 right ascension of the center of the galaxy, in decimal degrees
	Dec : float
		The J2000 declination of the center of the galaxy, in decimal degrees
	Vsys : float
		The systemic velocity at the center of the galaxy, in km/s

		
	Returns
	-------
	None
	
	Notes
	-----
	Required packages: numpy, astropy, scipy, spectral_cube, pandas, radio_beam, datetime, subprocess, matplotlib
	Author: R. C. Levy (rlevy@astro.umd.edu)
	Based on COBeamSmearingCorr.m, HgBeamSmearingCorr.m, Cdefvar.m, fitpersic.m by R. C. Levy and mkvelfield.m by A. D. Bolatto
	Last updated: 2019-05-14
	Change log:
		2019-05-13 : file created, RCL
		2019-05-14 : initial commit, RCL
	
	Examples
	--------
	>>> from BeamSmearingCorrection_Levy18 import beam_smearing_correction
	>>> path_to_veldispmap = 'mom2.fits'
	>>> path_to_rotationcurve = 'rotcurv.txt'
	>>> spectral_resolution_FWHM = 50. #km/s
	>>> chan_width = 20. #km/s
	>>> PA = 30. #degrees
	>>> inc = 50. #degrees
	>>> RA = 109.0163115 #degrees
	>>> Dec = 64.7107858 #degrees
	>>> Vsys = 4300. #km/s
	>>> beam_smearing_correction(path_to_veldispmap,path_to_rotationcurve,spectral_resolution_FWHM,chan_width,PA,inc,RA,Dec,Vsys)
	
	'''

	#import modules
	import subprocess
	import numpy as np
	import astropy
	from astropy.io import fits
	import warnings
	warnings.simplefilter('ignore',category=fits.verify.VerifyWarning)
	from astropy.wcs import WCS
	import astropy.units as u
	#from astropy.convolution import Gaussian2DKernel, convolve
	import pandas as pd
	from scipy.optimize import curve_fit
	from scipy.interpolate import griddata
	import spectral_cube
	from spectral_cube import SpectralCube
	import radio_beam
	import datetime
	import matplotlib.pyplot as plt
	plt.rcParams['font.family'] = 'serif'

	#store version info
	spectral_cube_version = spectral_cube.__version__
	astropy_version = astropy.__version__
	radio_beam_version = radio_beam.__version__
	ipython_version = subprocess.Popen('ipython --version', shell=True, stdout=subprocess.PIPE).stdout.read().decode().replace('\n','')
	python_version = subprocess.Popen('python --version', shell=True, stdout=subprocess.PIPE).stdout.read().decode().replace('\n','')
	

	#define functions
	def project_coordinates(hdr,RA,Dec,PA,inc):
		#make projected map of coordinates
		#make arrays to hold the x and y pixel coordinates
		x_pix = np.arange(1,hdr['NAXIS1']+1,1)
		y_pix = np.arange(1,hdr['NAXIS2']+1,1)
		#convert from pixels to RA and Dec
		x,y = wcs.all_pix2world(x_pix,y_pix,1,ra_dec_order=True)
		#cast as offsets from the central RA and Dec
		x_cen = (x-RA)*3600. #arcsec
		y_cen = (y-Dec)*3600. #arcsec
		#make a 2D grid
		xx,yy = np.meshgrid(x_cen,y_cen)
		#project
		xr = xx*np.cos(np.radians(PA))+yy*np.sin(np.radians(PA))
		yr = xx*np.sin(np.radians(PA))-yy*np.cos(np.radians(PA))
		yr = yr/np.cos(np.radians(inc))
		#find radial distance of each pixel from the center
		rr = np.sqrt(xr**2+yr**2) #in arcsec
		tt = np.arctan2(yr,xr)
		return rr, tt


	def fit_persic(rotcurve):
		#fit Universal Rotation Curve from Persic et al. 1996 (Eq. 14)
		def URC(x,a,b,c,):
			return a*((-3.856+0.44*np.log10(b))*(1.97*(x/c)**(1.22)/((x/c)**2+0.6084)**(1.43))+1.6*np.exp((-1.5924E-11)*b)*((x/c)**2/((x/c)**2+(1.56E-4)*(b**0.4))))**(0.5)
		R = rotcurve['R_arcsec'].values #arcsec
		Vrot = rotcurve['Vrot_kms'].values #km/s
		fit,cov = curve_fit(URC,R,Vrot,bounds=([0., 1E6, 0.],[1E3, 1E9, 50.]))
		Vrot_model = URC(R,*fit)
		return R, Vrot, Vrot_model

	#load velocity dispersion map
	vel_disp = fits.open(path_to_veldispmap)
	#grab the header
	hdr = vel_disp[0].header
	#store the world coordinate system (wcs) information
	wcs = WCS(hdr)
	#load the original mom2 data
	mom2_data = vel_disp[0].data #km/s, sigma
	mom2_data_FWHM = mom2_data*2.355 #km/s, FWHM

	#load rotation curve
	rotcurve = pd.read_csv(path_to_rotationcurve,sep='\t',comment='%')

	#get smooth model rotation curve by fitting Universal Rotation Curve
	R, Vrot_rc, Vrot_model = fit_persic(rotcurve)

	#make a model velocity field from the model rotation curve
	rr,tt = project_coordinates(hdr,RA,Dec,PA,inc)
	max_rr = np.max(rr)
	max_R = np.max(R)
	min_R = np.min(R)
	#fill the center
	if min_R > 0.:
		R = np.insert(R,0,0.)
		Vrot_rc = np.insert(Vrot_rc,0,0.)
		Vrot_model = np.insert(Vrot_model,0,0.)
	#fill to the edge of the map
	if max_rr > max_R:
		R = np.append(R,max_rr)
		Vrot_rc = np.append(Vrot_rc,Vrot_rc[-1])
		Vrot_model = np.append(Vrot_model,Vrot_model[-1])
	#interpolate into 2D grid
	vvrot = griddata(R,Vrot_model,rr)
	#get model velocity field
	vfield_model = vvrot*np.cos(tt)*np.sin(np.radians(inc))+Vsys

	#create a rotation-only data cube from the model velocity field
	chan_vstart = Vsys-1000. #km/s
	chan_vend = Vsys+1000. #km/s
	chans = np.arange(chan_vstart,chan_vend+chan_width,chan_width)

	#make an empty data cube
	datacube_model = np.zeros((len(chans),hdr['NAXIS2'], hdr['NAXIS1']))

	#convert spectral res from FWHM to sigma
	specres_sigma = spectral_resolution_FWHM/2.355

	#at each spaxel, put a Gaussian centered on the correct velocity channel given by the model vfield
	for i in range(hdr['NAXIS1']):
		for j in range(hdr['NAXIS2']):
			#find index of channel that is closest to the velocity at that pixel
			idx = np.argmin(np.abs(vfield_model[i,j]-chans))
			#put a Gaussian line
			gauss = np.exp(-(chans-vfield_model[i,j])**2/(2*specres_sigma**2))
			datacube_model[:,j,i]=gauss.copy()

	#grab beam info from header
	bmaj = hdr['BMAJ']*3600. #arcsec, beam major axis
	bmin = hdr['BMIN']*3600. #arcsec, baem minor axis
	bpa = hdr['BPA'] #degrees, beam position angle

	#add the spectral axis to the header
	hdr_cube = hdr.copy()
	hdr_cube.remove('HISTORY',remove_all=True)
	hdr_cube['NAXIS']=3
	hdr_cube.insert('NAXIS2',('NAXIS3',len(chans)),useblanks=False,after=True)
	hdr_cube.insert('CTYPE2',('CRPIX3',1E0),useblanks=False,after=True)
	hdr_cube.insert('CRPIX3',('CDELT3',chan_width),useblanks=False,after=True)
	hdr_cube.insert('CDELT3',('CRVAL3',chan_vstart),useblanks=False,after=True)
	hdr_cube.insert('CRVAL3',('CTYPE3','VRAD'),useblanks=False,after=True)
	hdr_cube.insert('CTYPE3',('CUNIT3','km/s'),useblanks=False,after=True)
	wcs_cube = WCS(hdr_cube)

	#convolve 
	#set the original beam half the pixel size
	beam_orig = radio_beam.Beam(major = hdr['CDELT2']/2*3600*u.arcsec, minor=hdr['CDELT2']/2*3600*bmaj/bmin*u.arcsec, pa=bpa*u.deg)
	#final beam is from original mom2 
	beam_final = radio_beam.Beam(major = bmaj*u.arcsec, minor=bmin*u.arcsec, pa=bpa*u.deg)
	#store simulated data cube as spectralcube
	spec_cube = SpectralCube(data=datacube_model,wcs=wcs_cube,beam=beam_orig)
	#convolve to the beam size of the data
	cube_conv = spec_cube.convolve_to(beam_final)

	#derive velocity dispersion (mom2) from convolved data cube
	mom2_beamsmeared_FWHM = cube_conv.linewidth_fwhm().to_value('km/s')

	#remove simulated velocity dispersion due to beam smearing from the data in quadrature
	mom2_corr = np.sqrt(mom2_data_FWHM**2-mom2_beamsmeared_FWHM**2) #km/s FWHM

	#save simulated beam smearing to fits file
	hdr_sim = hdr.copy()
	hdr_sim.remove('HISTORY',remove_all=True)
	hdr_sim['BMAJ'] = hdr['CDELT2']/2 #degrees
	hdr_sim['BMIN'] = hdr['CDELT2']/2*bmin/bmaj #degrees, keep beam axis ratio
	hdr_sim['ORIGIN'] = 'Created with '+python_version+', astropy '+astropy_version+', spectral_cube '+spectral_cube_version+': '+datetime.datetime.utcnow().strftime('%c')+' UTC'
	hdu = fits.PrimaryHDU(mom2_beamsmeared_FWHM)
	hdul  =fits.HDUList([hdu])
	hdul[0].header=hdr_sim
	fname = path_to_veldispmap.rsplit('.fits')[0]+'.SimBeamSmearing.fits'
	hdul.writeto(fname,overwrite=True)


	#save beam smearing-corrected mom2 to fits file
	hdr_corr = hdr.copy()
	hdr_corr.remove('HISTORY',remove_all=True)
	hdr_corr['ORIGIN'] = hdr_corr['ORIGIN']+'; Created with '+python_version+', astropy '+astropy_version+', spectral_cube '+spectral_cube_version+': '+datetime.datetime.utcnow().strftime('%c')+' UTC'
	hdu = fits.PrimaryHDU(mom2_corr)
	hdul  =fits.HDUList([hdu])
	hdul[0].header=hdr_corr
	fname = path_to_veldispmap.rsplit('.fits')[0]+'.BeamSmearingCorr.fits'
	hdul.writeto(fname,overwrite=True)

	#plot the mom2 maps from data, simulated beam smearing, and beam smearing-corrected
	fig = plt.figure(1,figsize=(9,3))
	plt.clf()

	ax0 = plt.subplot(131,projection=wcs)
	im=ax0.imshow(mom2_data_FWHM,origin='lower')
	cb=plt.colorbar(im,fraction=0.046,pad=0.04) #these values ~magically~ keep the colorbar scaled to the plot size regardless of the plot size
	cb.set_label('$\\sigma_{\mathrm{V,FWHM}}$ (km s$^{-1}$)',fontsize=9.)
	cb.ax.tick_params(labelsize=8.)
	ax0.coords[0].set_major_formatter('hh:mm:ss.s')
	ax0.coords[0].set_ticklabel(size=8.,exclude_overlapping=True) 
	ax0.coords[1].set_ticklabel(size=8.,exclude_overlapping=True) 
	ax0.set_title('Data',size=10.)

	ax1 = plt.subplot(132,projection=WCS(hdr_sim))
	im=ax1.imshow(mom2_beamsmeared_FWHM,origin='lower')
	cb=plt.colorbar(im,fraction=0.046,pad=0.04)
	cb.set_label('$\\sigma_{\mathrm{V,FWHM}}$ (km s$^{-1}$)',fontsize=9.)
	cb.ax.tick_params(labelsize=8.)
	ax1.coords[0].set_major_formatter('hh:mm:ss.s')
	ax1.coords[0].set_ticklabel(size=8.,exclude_overlapping=True) 
	ax1.coords[1].set_ticklabel(size=8.,exclude_overlapping=True) 
	ax1.set_title('Simulated Beam Smearing',size=10.)

	ax2 = plt.subplot(133,projection=WCS(hdr_corr))
	im=ax2.imshow(mom2_corr,origin='lower')
	cb=plt.colorbar(im,fraction=0.046,pad=0.04)
	cb.set_label('$\\sigma_{\mathrm{V,FWHM}}$ (km s$^{-1}$)',fontsize=9.)
	cb.ax.tick_params(labelsize=8.)
	ax2.coords[0].set_major_formatter('hh:mm:ss.s')
	ax2.coords[0].set_ticklabel(size=8.,exclude_overlapping=True) 
	ax2.coords[1].set_ticklabel(size=8.,exclude_overlapping=True) 
	ax2.set_title('Beam Smearing Corrected',size=10.)

	plt.subplots_adjust(wspace=1.0)
	fname = path_to_veldispmap.rsplit('.fits')[0]+'.BeamSmearingComparison.pdf'
	plt.savefig(fname,bbox_inches='tight')
	plt.close()

	





