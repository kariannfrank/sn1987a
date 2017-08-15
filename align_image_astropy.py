#! /usr/bin/env python

#Author: Kari A. Frank 
#Date: November 7, 2014
#Purpose: Align a deconvolved, smoothed, ACIS image with a reference image.
#Usage: align_image.py targetimg refimg [--clobber clobber]
#
#Input:
#
# targetimg = name of fits file containing the image to be shifted
# refimg = name of fits file containing the reference image       
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
from home_grown import ls_to_list
import numpy as np
import fitting
import astro_utilities as au
#import home_grown as hg #import my own custom functions
from astropy.modeling.models import custom_model_1d
from astropy.modeling.fitting import LevMarLSQFitter
import image_registration as imgreg

#from scipy import optimize

#---------------------------------------
#   Parse arguments and set defaults
#---------------------------------------

parser = argparse.ArgumentParser(description='Align the target image with the reference image.')

pwd = os.getcwd()

parser.add_argument('targetimg',help='Target fits image.')
parser.add_argument('refimg',help='Reference fits image.')

parser.add_argument('--clobber',help='Overwrite existing files.',default='no')

args = parser.parse_args()

#args.burnin = int(args.burnin)
#args.niterations = int(args.niterations)

#---------------------------------------
#        Prep Image Files
#---------------------------------------

# Create normalized images

# Target image
# read file
targetimgfile = fits.open(args.targetimg)
targetimg = targetimgfile[0].data
targetimgfile.close()
# divide image
targetimgtotal = np.sum(targetimg)
targetimg = targetimg/targetimgtotal

# Reference image
# read image
refimgfile = fits.open(args.refimg)
refimg = refimgfile[0].data
refimgfile.close()
# divide image
refimgtotal = np.sum(refimg)
refimg = refimg/refimgtotal

#---------------------------------------
#           Align Image
#---------------------------------------

# Define model

#@custom_model_1d
def img_align( in_image, dx=10, dy=10):

    # the model is a function which transforms target image into the 
    #  reference image
    
    # Translate Image
    img_trans = np.roll(in_image,int(dx),axis=1)
    img_trans = np.roll(img_trans,int(dy),axis=0)

    return img_trans

align_model = custom_model_1d(img_align)
    


# Fit to reference image
#  (here the 'x' values are targetimg and the 'y' are refimg, 
#   as the function maps x (target) to y (ref))

## using astropy.modeling

# initialize model and fitter
#m_init = img_align()
#m_init = align_model()
#fit  = LevMarLSQFitter()

# fit
#m = fit(m_init,targetimg,refimg)
 
# Print results
#print "best-fit model = ", m

## using image_registration

#offsets_target,offsets_ref = imgreg.tests.compare_methods(targetimg,refimg,noise=0.1)
final_dx,final_dy,dxerr,dyerr = imgreg.chi2_shift(targetimg,refimg,upsample_factor=1000)
print 'dx, dy, = ',final_dx,final_dy
print 'dxerr, dyerr, = ',dxerr,dyerr

#---------------------------------------
#        Save Translated Image
#---------------------------------------

# Translate Image
img_trans = np.roll(targetimg,int(final_dx),axis=1)
img_trans = np.roll(img_trans,int(final_dy),axis=0)

# Create new fits
outfile = args.targetimg+'_trans'
newhdu = fits.PrimaryHDU(img_trans)
newhdu.writeto(outfile,clobber=True)

# Copy header from original target image file
##stuff = au.transfer_header(args.targetimg,outfile,outfile)
