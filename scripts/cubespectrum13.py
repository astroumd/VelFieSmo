#! /usr/bin/env python
#
#  Load a FITS cube , extract the spectrum at a (or reference) pixel
#  and operate and plot some and then more....
#
#
#  22-jun-2017  PJT             summer project - cloned off cubespectrum.py
#  july-2017    Thomas/Peter    various improvements
#
#  @todo
#     - have optional RESTFRQ or RESTFREQ as 3rd argument [done]
#     - output the spectrum in a table, much like testCubeSpectrum.tab [done]
#     - resample the gauss finer (not 5 points but may be 10x more?)
#
#  Note xpos,ypos are 0 based (thus reference pixel default is off by 1 pixel)



import os, sys, math
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import pyspeckit
from scipy.optimize import curve_fit
from astropy.io import fits
from astropy.units import Quantity
c = 299792.458     # [km/s] there should be a way to get 'c' from astropy.units ?


#  set common command line options (unless overriden below)
vlsr = None
    
na = len(sys.argv)
if na == 8:
    fitsfile = sys.argv[1] 
    pos = [int(sys.argv[2]),int(sys.argv[3])]
    restfreq = float(sys.argv[4])* 1e9
    vmin = float(sys.argv[5])
    vmax = float(sys.argv[6])
    vlsr = float(sys.argv[7])
    use_vel = False
elif na == 7:
    # Must be in Km/s
    fitsfile = sys.argv[1] 
    pos = [int(sys.argv[2]),int(sys.argv[3])]
    restfreq = float(sys.argv[4])* 1e9
    vmin = float(sys.argv[5])
    vmax = float(sys.argv[6])
    use_vel = True
elif na == 5: 
    # Must be in GHz
    fitsfile = sys.argv[1]
    pos = [int(sys.argv[2]),int(sys.argv[3])]
    vmin = vmax = None
    restfreq = float(sys.argv[4])* 1e9
    use_vel = True
elif na == 4:
    # Pixel position
    fitsfile = sys.argv[1]
    pos = [int(sys.argv[2]),int(sys.argv[3])]
    restfreq = None
    vmin = vmax = None
    use_vel = False
elif na == 2:
    # Fits file
    fitsfile = sys.argv[1]
    pos = None    
    restfreq = None
    vmin = vmax = None
    use_vel = False 
else:
    print("Usage: %s fitsfile [xpos ypos] [restfreq [vmin vmax] [vlsr]" % sys.argv[0])
    sys.exit(1)


# open the fits file
hdu = fits.open(fitsfile)
print(len(hdu))


# get a reference to the header and data.  Data should be 3dim numpy array now
h = hdu[0].header
d = hdu[0].data.squeeze()
print(d.shape)


#  grab the restfreq, there are at least two ways how this is done
if restfreq == None:
    if 'RESTFRQ' in h:
        restfreq=h['RESTFRQ']
    elif 'RESTFREQ' in h:
        restfreq=h['RESTFREQ']
    else:
        restfreq= h['CRVAL3']
print("RESTFREQ",restfreq/1e9)


if pos == None:
    # the FITS reference pixel is always a good backup
    xpos = int(h['CRPIX1'])
    ypos = int(h['CRPIX2'])
    print("No position given, using reference pixel %g %g" % (xpos,ypos))
else:
    xpos = pos[0]
    ypos = pos[1]


flux     = d[:,ypos,xpos]


nchan    = d.shape[0]       
channeln = np.arange(nchan)
zero     = np.zeros(nchan)


cdelt3 = h['CDELT3']
crval3 = h['CRVAL3']
crpix3 = h['CRPIX3']
# to convert the channel to frequency
channelf = (channeln-crpix3+1)*cdelt3 + crval3
# to convert the Frequency to velocity
#channelv = (1.0-channelf/restfreq) * c
#print (channelf)
#print (channelv)

# what we plot
#channel = channelv 
#channel = channelf
#channel = channeln

# to convert the Frequency to velocity
channelv = (1.0-channelf/restfreq) * c
print (channelv.min(), channelv.max())
print ("min freq", channelf.min()/1e9, "max freq", channelf.max()/1e9)


# to create a spectrum, a table of the flux vs. (rest or sky) frequency 
if vlsr == None:
    xtab = channelf / 1e9                  # sky freqency to GHz
    w = fitsfile + str(pos) + " [sky freq]"
else:  
    xtab = channelf / (1-vlsr/c) / 1e9     # rest freq in GHz
    w = fitsfile + str(pos) + " [rest freq w/ vlsr=%g]" % vlsr
ytab = flux 
np.savetxt('Frequency_Flux.tab',np.c_[xtab,ytab], delimiter='  ', header = (w), comments='#',fmt='%.8f')




def gfit1(xi,yi,m=5):
    """
    moments around a peak
    (also) rely on number of pixels left and right of the peak. Masking optional
    """
    print("GFIT1")
    ipeak = yi.argmax()
    xpeak = xi[ipeak]
    ypeak = yi[ipeak]
    print(ipeak,xpeak,ypeak)
    print(xi[ipeak-m:ipeak+m])
    print(yi[ipeak-m:ipeak+m])
    x = xi[ipeak-m:ipeak+m]
    y = yi[ipeak-m:ipeak+m]
    print(xi.shape)
    print(yi.shape)
    xmean = (x*y).sum() / y.sum()
    xdisp = (x*x*y).sum() / y.sum() - xmean*xmean
    if xdisp > 0:
        xdisp = math.sqrt(xdisp)
    fwhm = 2.355 * xdisp    
    print("MEAN/DISP/FWHM:",xmean,xdisp,fwhm)
    ymodel = ypeak * np.exp(-0.5*(xi-xmean)**2/(xdisp*xdisp))
    return ymodel

def gfit2(x,y):
    """
    rely on masking completely
    moments around a peak    
    """
    print("GFIT2")
    xmean = (x*y).sum() / y.sum()
    xdisp = (x*x*y).sum() / y.sum() - xmean*xmean
    if xdisp > 0:
        xdisp = math.sqrt(xdisp)
    fwhm = 2.355 * xdisp    
    print("MEAN/DISP/FWHM:",xmean,xdisp,fwhm)
    ypeak = y.max()
    print(ypeak)
    ymodel = ypeak * np.exp(-0.5*(x-xmean)**2/(xdisp*xdisp))
    return ymodel

def gfit3(xi,yi):
    """
    relies on masking , use pyspeckit
    """
    # not sure if we need this, or try x = xi
    x = ma.compressed(xi)
    y = ma.compressed(yi)    
    sp = pyspeckit.Spectrum(data=y, xarr=x, error=None, header=h,)
    sp.plotter()
    sp.specfit(fittype='gaussian')
    sp.specfit.plot_fit()
    # sp.baseline()
    print(x)
    print(y)
    # fake a return array 
    return yi

def gfit4(x,y):
    """
        relies on masking , use scipy's curve_fit
    """
    def gauss(x, *p):
        A, mu, sigma, B = p
        return A*np.exp(-(x-mu)**2/(2.*sigma**2)) + B
    #
    # first get some initial estimates
    B = y.min()
    A = y.max() - B
    sigma = (x.max() - x.min() ) /2.0 / 2.5     # this can be done better, moment analysis?
    mu    = (x.max() + x.min() ) /2.0
    p0 = [A, mu, sigma, B]
    # p0 = [0.030544054, 50, 200, 0.029194169]

    print("p0 = ",p0)
    
    coeff, cm = curve_fit(gauss, x, y, p0=p0)
    ymodel = gauss(x, *coeff)
    
    print("Fitted amp            :",coeff[0])
    print("Fitted mean           :",coeff[1])
    print("Fitted sigma and FWHM :",coeff[2], coeff[2]*2.355)
    print("Fitted baseline       :",coeff[3])
    print("Covariance Matrix     :\n",cm)
    # what are now the errors in the fitted values?
    print("error amp     :",math.sqrt(cm[0][0]))
    print("error mean    :",math.sqrt(cm[1][1]))
    print("error sigma   :",math.sqrt(cm[2][2]))
    print("error baseline:",math.sqrt(cm[3][3]))
    
    return ymodel,coeff[1]
    
if use_vel == True:
    plt.figure()  
    if vmin != None:
        # mask
        channelv = ma.masked_outside(channelv,vmin,vmax)
        flux     = ma.masked_array(flux, channelv.mask)
        # make arrays smaller
        channelv = ma.compressed(channelv)
        flux     = ma.compressed(flux)
        # plotting
        # plt.xlim([vmin,vmax])         # technically not needed
        ymodel,vfit = gfit4(channelv,flux)
        plt.plot(channelv,ymodel,label='gfit4')    
    plt.plot(channelv,flux,'o-',markersize=2,label='data')
    # plt.plot(channelv,zero)
    plt.xlabel("Velocity (km/s)")
    plt.ylabel("Flux")
    plt.title(fitsfile +"  @ %g %g" % (xpos,ypos)+ "   %g" % (restfreq/1e9)+ 'Ghz')
    plt.legend()
    plt.show()
else:  
    plt.figure()
    if vlsr != None:
        print ("Gaussian Distribution")
        channelv = ma.masked_outside(channelv,vmin,vmax)
        channelf = ma.masked_array(channelf, channelv.mask)        
        flux     = ma.masked_array(flux,     channelv.mask)
        channelv = ma.compressed(channelv)
        channelf = ma.compressed(channelf)
        flux     = ma.compressed(flux)
        ymodel,ffit = gfit4(channelf/1e9,flux)
        plt.plot(channelf/1e9,ymodel,label='gfit4')
        plt.plot(channelf/1e9,flux,'o-',markersize=2,label='data')
        # compute restfreq
        f0 = ffit / (1-vlsr/c)
        print("Fitted restfreq f0=",f0)  
    else:
        plt.plot(channelf/1e9,flux,'o-',markersize=2,label='data')
        plt.plot(channelf/1e9,zero)
    plt.xlabel("Frequency (GHz)")
    plt.ylabel("Flux")
    plt.title(fitsfile + " @ %g %g" % (xpos,ypos))
    plt.legend()
    plt.show()
    

print("Mean and RMS of %d points: %g %g" % (len(flux),flux.mean(),flux.std()))

















