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
#   Read in fpars files and set priors
#---------------------------------------

# read files and headers
target_fpars_file = fits.open(args.targetfpars)
ref_fpars_file = fits.open(args.reffpars)
target_head = fits.getheader(args.targetfpars)
ref_head = fits.getheader(args.reffpars)

# get tables
target_table = target_fpars_file[1].data
ref_table = ref_fpars_file[1].data

# get values
target_x0 = target_table.field('X0')[0]
target_x0err = target_table.field('X0ERR')[0]
target_y0 = target_table.field('Y0')[0]
target_y0err = target_table.field('Y0ERR')[0]

ref_x0 = ref_table.field('X0')[0]
ref_x0err = ref_table.field('X0ERR')[0]
ref_y0 = ref_table.field('Y0')[0]
ref_y0err = ref_table.field('Y0ERR')[0]

# calculate deltaX and deltaY and associated sigmas
dx = target_x0 - ref_x0
#dy = target_y0 - ref_y0
dy = ref_y0 - target_y0
sigma_dx = 5.0*(target_x0err**2.0 + ref_x0err**2.0)**0.5
sigma_dy = 5.0*(target_y0err**2.0 + ref_y0err**2.0)**0.5

print dx,dy,sigma_dx,sigma_dy

#---------------------------------------
#        Prep Image Files
#---------------------------------------

# get obsids
target_obsid = crt.dmkeypar(args.targetimg,"OBS_ID")
ref_obsid = crt.dmkeypar(args.refimg,"OBS_ID")

# Create normalized images

# Target image
# get normalization constant
targetimgfile = fits.open(args.targetimg)
targetimg = targetimgfile[0].data
targetimgfile.close()
targetimgtotal = np.sum(targetimg)
# divide image
crt.dmimgcalc.punlearn()
crt.dmimgcalc.infile = args.targetimg
crt.dmimgcalc.infile2 = 'none'
crt.dmimgcalc.outfile = args.targetimg+'_normalized'
crt.dmimgcalc.op = 'add'
crt.dmimgcalc.weight = 1.0/targetimgtotal
crt.dmimgcalc.clobber = args.clobber
crt.dmimgcalc()

# Reference image
# get normalization constant
refimgfile = fits.open(args.refimg)
refimg = refimgfile[0].data
refimgfile.close()
refimgtotal = np.sum(refimg)
# divide image
crt.dmimgcalc.punlearn()
crt.dmimgcalc.infile = args.refimg
crt.dmimgcalc.infile2 = 'none'
crt.dmimgcalc.outfile = args.refimg+'_normalized'
crt.dmimgcalc.op = 'add'
crt.dmimgcalc.weight = 1.0/refimgtotal
crt.dmimgcalc.clobber = args.clobber
crt.dmimgcalc()

#---------------------------------------
#              Run MCMC
#---------------------------------------

mcmc_pars = fitting.mcmc_align_img(args.targetimg+'_normalized',args.refimg+'_normalized',[dx,dy],[sigma_dx,sigma_dy],niterations=args.niterations)

final_dx = np.median(mcmc_pars[args.burnin:,1])
final_dy = np.median(mcmc_pars[args.burnin:,2])
print 'final dx,dy = ',final_dx,final_dy

sigmadx = np.std(mcmc_pars[args.burnin:,1])
sigmady = np.std(mcmc_pars[args.burnin:,2])
print 'sigmadx,sigmady = ',sigmadx,sigmady

#---------------------------------------
#   Translate Image with Best dx, dy
#---------------------------------------

# copy image
crt.dmcopy.punlearn()
crt.dmcopy.infile = args.targetimg
crt.dmcopy.outfile = args.targetimg+'_trans'
crt.dmcopy.clobber = args.clobber
crt.dmcopy()

# shift image
crt.wcs_update.punlearn()
crt.wcs_update.infile = args.targetimg+'_trans'
crt.wcs_update.wcsfile = args.targetimg+'_trans'
crt.wcs_update.deltax = str(final_dx)
crt.wcs_update.deltay = str(final_dy)
crt.wcs_update.rotang = 0.0
crt.wcs_update.scalefac = 1.0
crt.wcs_update.clobber = 'yes'
crt.wcs_update()
