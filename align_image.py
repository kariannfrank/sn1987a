#!/usr/bin/env python

#Author: Kari A. Frank 
#Date: November 7, 2014
#Purpose: Align a deconvolved, smoothed, ACIS image with a reference image.
#Usage: align_image.py targetimg refimg [--threshold threshold] [--rgb rgb --rgbdir rgbdir] [--raw rawdir] [--mincounts mincounts] [--outfile outfile] [--clobber clobber]
#
#Input:
#
# targetimg = name of fits file containing the image to be shifted, 
#             or a .lis text file containing a list of target image
#             file names
#
# refimg = name of fits file containing the reference image       
# 
# threshold = fraction of the maximum image value to suppress.
#             all pixels below threshold will have their values
#             multiplied by 0.4 in the output image.  set to 0
#             to disable. (default = 0.0)(try 0.35 for broadband images)
#             
# rgb =       optionally translate the soft, medium, and hard images 
#             associated with targetimg, by the same deltax,deltay
#             as targetimg (default = 'no')
#
# rgbdir =    directory containing the rgb image files (ignored if rgb='no')
#
# raw =      optionally translate the raw (non-deconvolved) images bye 
#            the same amount (default='no')
# 
# rawdir =   directory containing the raw images (ignored if raw='no')
#
# mincounts = optionally skip images with low counts (default=0)
#
# suffix =   optional suffix for the output (aligned) image file
#
# clobber =   specifies whether files should be overwritten
#             if they already exist (same as in ciao tools,
#             default = 'no')
#
#Output:
# - a translated copy of targetimg that should be aligned with refimg
#
#Usage Notes:
# - cannot be used with the ciao version of python (requires astropy)

#---------------------------------------
#           Import Modules
#---------------------------------------
import argparse,os,re
#import pyfits as fits
import astropy.io.fits as fits
import home_grown as hg
import numpy as np
import fitting
import astro_utilities as au
#import home_grown as hg #import my own custom functions
#from astropy.modeling.models import custom_model_1d
from astropy.modeling.fitting import LevMarLSQFitter
import image_registration as imgreg

#from scipy import optimize

#---------------------------------------
#   Parse arguments and set defaults
#---------------------------------------

parser = argparse.ArgumentParser(description='Align the target image with the reference image.')

pwd = os.getcwd()

parser.add_argument('targetimg',help='Target fits image or .lis file containing list of target image files.')
parser.add_argument('refimg',help='Reference fits image.')
parser.add_argument('--threshold',help='Suppression threshold.',default=0.0,type=float)
parser.add_argument('--rgb',help='Also translate rgb images.',default='no')
parser.add_argument('--rgbdir',help='Directory containing rgb images.',default='')
parser.add_argument('--raw',help='Also translate raw images.',default='no')
parser.add_argument('--rawdir',help='Directory containing raw images.',default='')
parser.add_argument('--mincounts',help='Skip images with less than mincounts total counts.',default=0.0,type=float)
parser.add_argument('--suffix',help='Suffix to append to the end of output file name.',default='')
parser.add_argument('--clobber',help='Overwrite existing files.',default='no')

args = parser.parse_args()
args.threshold = float(args.threshold)
#args.burnin = int(args.burnin)
#args.niterations = int(args.niterations)

#(in future, may make this an argument)
#rgb_bands = ['300-800','800-1200','1200-8000','3000-8000']
rgb_bands = ['300-800','800-1200','1200-8000','2000-10000','3000-8000']

#---------------------------------------
#        Prep Image Files
#---------------------------------------

## Open and Normalize Reference Image

# Reference image
# read image
refimgfile = fits.open(args.refimg)
refimg = refimgfile[0].data
refimgfile.close()
# divide image
refimgtotal = np.sum(refimg)
refimg = refimg/refimgtotal


##---- Begin Loop Over Target Images ----

if os.path.splitext(args.targetimg)[1] == '.lis':
    imgfiles = hg.read_list(args.targetimg,comment='#')
else:
    imgfiles = [args.targetimg]

for imgfile in imgfiles:

    ## Open and Normalize Target Image

    # read file
    targetimgfile = fits.open(imgfile)
    targetimg = targetimgfile[0].data
    targetimgfile.close()
    # divide image
    targetimgtotal = np.sum(targetimg)
    print targetimgtotal
    if targetimgtotal >= args.mincounts:
        targetimg_norm = targetimg/targetimgtotal


#---------------------------------------
#           Align Image
#---------------------------------------

        final_dx,final_dy,dxerr,dyerr = imgreg.chi2_shift(targetimg_norm,refimg,upsample_factor=1000)
        print 'dx, dy, = ',final_dx,final_dy
        print 'dxerr, dyerr, = ',dxerr,dyerr

#---------------------------------------
#        Save Translated Image
#---------------------------------------

    # Translate Image
        img_trans = np.roll(targetimg,int(final_dx),axis=1)
        img_trans = np.roll(img_trans,int(final_dy),axis=0)

    # Truncate lower limit to zero (eliminate negative values)
        img_trans[img_trans < 0.0] = 0.0

    # Suppress very faint edge emission to zero (soft removal of 
    # faint artifacts in bright images)
        thresh = args.threshold*np.max(img_trans)
        print 'thresh = ',thresh
        img_trans[img_trans < thresh] = img_trans[img_trans < thresh]*img_trans[img_trans < thresh]/(args.threshold*np.max(img_trans))

    # Create new fits
        outfile = imgfile+'_trans'+args.suffix
        newhdu = fits.PrimaryHDU(img_trans)
        newhdu.writeto(outfile,clobber=True)

    # Copy header from original target image file
        stuff = au.transfer_header(imgfile,outfile,outfile)

#---------------------------------------
#   Repeat for Other Bands (Same Observation)
#---------------------------------------

        if args.rgb == 'yes':

            # Get original band string
            targetband=re.findall('\d*-\d*',imgfile)[0]
            obsid = re.findall('\d{5,5}_',imgfile)[0]
            print obsid

            for band in rgb_bands:
                print band
                # set up file name(s)
                bandimgname=args.rgbdir+'/'+imgfile.replace(targetband,band)

                # open file
                if os.path.isfile(bandimgname):
                    bandimgfile = fits.open(bandimgname)
                    bandimg = bandimgfile[0].data
                    bandimgfile.close()
                    # divide image
                    bandimgtotal = np.sum(bandimg)
                    print bandimgtotal
                    if bandimgtotal >= args.mincounts:
                        bandimg_norm = bandimg/bandimgtotal

                        # Translate image arrays
                        img_trans = np.roll(bandimg,int(final_dx),axis=1)
                        img_trans = np.roll(img_trans,int(final_dy),axis=0)

                        # Truncate lower limit to zero (eliminate negative values)
                        img_trans[img_trans < 0.0] = 0.0

                        # Suppress very faint edge emission to zero (soft removal of 
                        # faint artifacts in bright images)
                        thresh = args.threshold*np.max(img_trans)
                        print 'thresh = ',thresh
                        img_trans[img_trans < thresh] = img_trans[img_trans < thresh]*img_trans[img_trans < thresh]/(args.threshold*np.max(img_trans))

                        # Create new fits
                        outfile = bandimgname+'_trans'
                        newhdu = fits.PrimaryHDU(img_trans)
                        newhdu.writeto(outfile,clobber=True)

                        # Copy header from original target image file
                        stuff = au.transfer_header(bandimgname,outfile,outfile)

                # Repeat for 'raw' image
                if args.raw=='yes':
                    bandimgname_raw = args.rawdir+'/'+hg.ls_to_list(args.rawdir\
                                      ,'*'+obsid+'*'+band+'*')[0]

                    # open file
                    if os.path.isfile(bandimgname_raw):
                        bandimgfile = fits.open(bandimgname_raw)
                        bandimg = bandimgfile[0].data
                        bandimgfile.close()
                    # divide image
                        bandimgtotal = np.sum(bandimg)
                        if bandimgtotal >= args.mincounts:
                            bandimg_norm = bandimg/bandimgtotal

                            # Convert shifts (raw and smoothed images have different pixel sizes)
                            rawpix_scale = 8.45 # raw pixels 8.45*smoothed pixels (=0.125"/0.0147929")

                            # Translate image arrays
                            img_trans = np.roll(bandimg,int(final_dx/rawpix_scale),axis=1)
                            img_trans = np.roll(img_trans,int(final_dy/rawpix_scale),axis=0)

                            # Truncate lower limit to zero (eliminate negative values)
                            img_trans[img_trans < 0.0] = 0.0

                            # Suppress very faint edge emission to zero (soft 
                            # removal of 
                            # faint artifacts in bright images)
                            thresh = args.threshold*np.max(img_trans)
                            print 'thresh = ',thresh
                            img_trans[img_trans < thresh] = img_trans[img_trans < thresh]*img_trans[img_trans < thresh]/(args.threshold*np.max(img_trans))

                            # Create new fits
                            outfile = bandimgname_raw+'_trans'
                            newhdu = fits.PrimaryHDU(img_trans)
                            newhdu.writeto(outfile,clobber=True)

                            # Copy header from original target image file
                            stuff = au.transfer_header(bandimgname_raw,outfile,outfile)


# end loop over files
