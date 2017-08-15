#! /usr/astro/ciao-4.6/bin/python

#Author: Kari A. Frank 
#Date: November 7, 2014
#Purpose: Align a deconvolved, smoothed, ACIS image with a reference image.
#Usage: align_image.py targetimg refimg targetfpars reffpars [--niterations niteraions] [--burnin burnin] [--clobber clobber]
#
#Input:
#
# targetimg = name of fits file containing the image to be shifted
# refimg = name of fits file containing the reference image       
# 
# targetfpars = name of fpars file for target image
# reffpars = name of fpars file for reference image
#
# niterations = number of iterations to use for the MCMC (default = 1000)
# burnin = number of iterations to disregard (default = 500)
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
import pyfits as fits
import ciao_contrib.runtool as crt #all functions should be prefixed with crt
from home_grown import ls_to_list
import numpy as np
import fitting
import astro_utilities as au
#import home_grown as hg #import my own custom functions

#---------------------------------------
#   Parse arguments and set defaults
#---------------------------------------

parser = argparse.ArgumentParser(description='Align the target image with the reference image.')

pwd = os.getcwd()

parser.add_argument('targetimg',help='Target fits image.')
parser.add_argument('refimg',help='Reference fits image.')

parser.add_argument('targetfpars',help='Target fpars file.')
parser.add_argument('reffpars',help='Reference fpars file.')

parser.add_argument('--niterations',help='Number of MCMC iterations',default=1000)
parser.add_argument('--burnin',help='Number of MCMC iterations to disregard',default=500)

parser.add_argument('--clobber',help='Overwrite existing files.',default='no')

args = parser.parse_args()

args.burnin = int(args.burnin)
args.niterations = int(args.niterations)

#----Set file paths----


#---------------------------------------
#        Prep Image Files
#---------------------------------------

# get obsids
target_obsid = crt.dmkeypar(args.targetimg,"OBS_ID")
ref_obsid = crt.dmkeypar(args.refimg,"OBS_ID")

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

# Set up Priors

# mean of delta priors
dx = 0
dy = 0

# standard deviation of delta priors (size of shifts)
sigma_dx = np.ceil(refimg.shape[0]/5)
sigma_dy = np.ceil(refimg.shape[1]/5)

# Run MCMC
mcmc_pars = fitting.mcmc_align_img(targetimg,refimg,[dx,dy],[sigma_dx,sigma_dy],niterations=args.niterations)

print mcmc_pars

final_dx = np.median(mcmc_pars[args.burnin:,1])
final_dy = np.median(mcmc_pars[args.burnin:,2])
print 'final dx,dy = ',final_dx,final_dy

sigmadx = np.std(mcmc_pars[args.burnin:,1])
sigmady = np.std(mcmc_pars[args.burnin:,2])
print 'sigmadx,sigmady = ',sigmadx,sigmady

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
stuff = au.transfer_header(args.targetimg,outfile,outfile)
