#! /usr/bin/env python

#Author: Kari A. Frank 
#Date: November 7, 2014
#Purpose: Align any image with a deconvolved, smoothed, ACIS image.
#Usage: align_radio_image.py targetimg refimg [--threshold threshold] [--clobber clobber]
#
#Input:
#
# targetimg = name of fits file containing the image to be shifted
# refimg = name of fits file containing the reference image       
# 
# threshold = fraction of the maximum image value to suppress.
#             all pixels below threshold will have their values
#             multiplied by 0.4 in the output image.  set to 0
#             to disable. (default = 0.0)
#             
#
# clobber =   specifies whether files should be overwritten
#             if they already exist (same as in ciao tools,
#             default = 'no')
#
#Output:
# - a translated copy of targetimg that should be aligned with refimg
#
#Usage Notes:

#---------------------------------------
#           Import Modules
#---------------------------------------
import argparse,os,sys
#import pyfits as fits
import astropy.io.fits as fits
import home_grown as hg
import numpy as np
import fitting
import astro_utilities as au
#import home_grown as hg #import my own custom functions
from astropy.modeling.models import custom_model_1d
from astropy.modeling.fitting import LevMarLSQFitter
import image_registration as imgreg
import scipy as sp

#from scipy import optimize

#---------------------------------------
#   Parse arguments and set defaults
#---------------------------------------

parser = argparse.ArgumentParser(description='Align the target image with the reference image.')

pwd = os.getcwd()

parser.add_argument('targetimg',help='Target fits image.')
parser.add_argument('refimg',help='Reference fits image.')
parser.add_argument('--threshold',help='Suppression threshold.',default=0.0)
parser.add_argument('--clobber',help='Overwrite existing files.',default='no')

args = parser.parse_args()
args.threshold = float(args.threshold)
#args.burnin = int(args.burnin)
#args.niterations = int(args.niterations)

#---------------------------------------
#        Prep Image Files
#---------------------------------------

## Open and Normalize Reference Image

# Reference image
# read image
refimgfile = fits.open(args.refimg)
refimg = refimgfile[0].data
refimgfile.close()
print refimg.shape

# divide image
refimgtotal = np.sum(refimg)
refimg = refimg/refimgtotal


##---- Begin Loop Over Target Images ----

if os.path.splitext(args.targetimg)[1] == '.lis':
    imgfiles = hg.read_list(args.targetimg)
else:
    imgfiles = [args.targetimg]


for imgfile in imgfiles:

    ## Open Target Image

    # read file
    targetimgfile = fits.open(imgfile)
    targetimg = targetimgfile[0].data
    targetimgfile.close()

    ## Rebin Target Image (if necessary)
    print targetimg.shape
    targetimg = sp.ndimage.zoom(targetimg,refimg.shape[0]/targetimg.shape[0])
    print targetimg.shape

    ## Normalize Target Image

    # divide image
    targetimgtotal = np.sum(targetimg)
    targetimg_norm = targetimg/targetimgtotal


#---------------------------------------
#           Align Image
#---------------------------------------

#    final_dx,final_dy,dxerr,dyerr = imgreg.chi2_shift(targetimg_norm,refimg,upsample_factor=1000)
#    print 'dx, dy, = ',final_dx,final_dy
#    print 'dxerr, dyerr, = ',dxerr,dyerr

#---------------------------------------
#        Save Translated Image
#---------------------------------------

    # Translate Image
#    img_trans = np.roll(targetimg,int(final_dx),axis=1)
#    img_trans = np.roll(img_trans,int(final_dy),axis=0)

    # Truncate lower limit to zero (eliminate negative values)
#    img_trans[img_trans < 0.0] = 0.0

    # Suppress very faint edge emission to zero (soft removal of 
    # faint artificants in bright images)
#    thresh = args.threshold*np.max(img_trans)
#    print 'thresh = ',thresh
#    img_trans[img_trans < thresh] = img_trans[img_trans < thresh]*img_trans[img_trans < thresh]/(args.threshold*np.max(img_trans))

    # Create new fits
#    outfile = imgfile+'_trans'
#    newhdu = fits.PrimaryHDU(img_trans)
#    newhdu.writeto(outfile,clobber=True)

    # Copy header from original target image file
#    stuff = au.transfer_header(imgfile,outfile,outfile)

# end loop over files
