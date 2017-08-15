import numpy as np
#import pyfits as fits
from astropy.io import fits
from home_grown import read_list
import astro_utilities as au
import shutil
from scipy.signal import fftconvolve
#lis='images_trans.lis'
#tfiles=read_list(lis)
#stack = au.stack_images(tfiles)

#------------- Add point source to deconvolved image -------------

#-- Set file names --
#decon_img = 'acisf15810_repro_evt2_subpix_300-8000.fits.deconvolved'
decon_img = 'hard_2000-10000_deconvolved_stacked.fits_trans'
out_img = decon_img+'.ptsrconly'

#-- Make copy of image to modify and open --
shutil.copy(decon_img,out_img) # preserves header information this way
hdulist = fits.open(out_img,'update')
img = hdulist[0].data
naxis = hdulist[0].header['NAXIS1']
ncounts = img.sum()

#-- Set ptsrc parameters --
sigma = 1. # gaussian width in deconvolved image pixels (pix=0.125")
ptsrc_counts = np.max(img) #set net flux to flux of brightest pixel
print 'ptsrc_counts = '+str(ptsrc_counts)
xsrc = naxis/2#-4.231#naxis/2 # source location taken as center of torus
ysrc = naxis/2#-4.397#naxis/2 #  from lobe+torus model fit

#-- Create ptsrc image array --
minthresh = 1e-10 # force the ptsrc to go to zero
xx,yy=np.mgrid[:naxis,:naxis]

# gaussian (same formulation as in Judith's paper)
r = ((xx-xsrc)**2.0 + (yy-ysrc)**2.0)**0.5
ptsrc = np.exp(-1.0*r**2.0/(2.0*sigma**2.0))
# normalize
ptsrc = ptsrc_counts*ptsrc/ptsrc.sum()
print 'total ptsrc counts = '+str(ptsrc.sum())
# truncate low value pixels
ptsrc[ptsrc<minthresh] = 0.0

#point
#ptsrc = np.zeros((naxis,naxis))
#ptsrc[np.max(xx)/2,np.max(yy)/2] = ptflux

#-- add ptsrc to image --
new_image = img+ptsrc
hdulist[0].data = new_image

#-- write updated image to file --
hdulist.close()

#write image to fits file
#hdu = fits.PrimaryHDU(ptsrc)
#hdulist = fits.HDUList([hdu])
#hdulist.writeto('ptsrc.fits')

#------------- (Re)Convolve New Image with PSF -------------

#-- Make copy of image to modify and open --
conv_img = out_img+'.reconvolved'
shutil.copy(out_img,conv_img) # preserves header information this way
hdulist = fits.open(conv_img,'update')
img = hdulist[0].data
naxis = hdulist[0].header['NAXIS1']
ncounts = img.sum()

#-- Read PSF image --
psf_file = 'psf.fits' #must have same binning as the deconvolved+ptsrc image
psf_hdulist = fits.open(psf_file)
psf_img = psf_hdulist[0].data
psf_naxis = psf_hdulist[0].header['NAXIS1']
psf_hdulist.close()
psf_counts = psf_img.sum()
psf_img_norm = psf_img/psf_counts

#-- Convolve --
conv_img = fftconvolve(img,psf_img,mode='same') #output image is same size as img

#-- Save convolved image --
hdulist[0].data = conv_img
hdulist.close()

#------------- 
# Reconvolve original (no point source) 
#   deconvolved image for comparison 
#-------------

reconv_img = decon_img+'.reconvolved'
shutil.copy(decon_img,reconv_img) # preserves header information this way
hdulist = fits.open(reconv_img,'update')
img = hdulist[0].data
naxis = hdulist[0].header['NAXIS1']
ncounts = img.sum()

#-- Convolve --
reconv_img = fftconvolve(img,psf_img,mode='same') #output image is same size as img

#-- Save convolved image --
hdulist[0].data = reconv_img
hdulist.close()

