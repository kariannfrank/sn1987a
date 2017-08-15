#! /usr/bin/env python

#Author: Kari A. Frank
#Date: April 4, 2014
#Purpose: Create copy of a sn1987a model image that has a header
#          obtained from an events file (to provide WCS info)
#
#Usage: add_model_hdr.py evtfile modelfile newfile --source_x source_x --source_y source_y [--clobber CLOBBER]
#
#Input:
#
# evtfile    -- event file from which to create image to get header
# 
# modelfile  -- fits file containing the model image (created with
#               pileupimage.pro
#
# newfile    -- name of file to save the model image with new header
#
# source_x,source_y -- physical coordinates of of the source center in
#                      in the events file (same as used in deconvolve.pro)
#
#
#Output:
# - new fits image file with the model image and header from the 
#   events file image
#
#Usage Notes:
# - only work assuming the default values were used in deconvolve.pro 
#   (40x40 pixel image) and process.pro (rebin_factor = 8.45, final 
#   image 200x200 pixels)
# - event file image is also filtered on energy, 300-8000 eV

#---------------------------------------
#           Import Modules
#---------------------------------------
import argparse
#import ciao_contrib.runtool as crt #all functions should be prefixed with crt
#import astropy.io.fits as fits
import subprocess as sp
from astro_utilities import transfer_header

#---------------------------------------
#   Parse arguments and set defaults
#---------------------------------------

parser = argparse.ArgumentParser(description='Creates copy of model image with a header created from an events file.',epilog='')

#pwd = os.getcwd()

parser.add_argument('evtfile',help='Event file from which to create header.')
parser.add_argument('modelfile',help='Fits file from which to get model image.')
parser.add_argument('newfile',help='Name of new fits file to be created.')

parser.add_argument('--source_x',help='Physical x coordinate of of the source center in the events file',default='4110')
parser.add_argument('--source_y',help='Physical y coordinate of of the source center in the events file',default='4110')

parser.add_argument('--clobber',help='Overwrite existing files.',default='no')

args = parser.parse_args()

#----Set file paths----



#---------------------------------------
#    Create Image from Events file
#---------------------------------------

#----set up image size and binsize----
radius = (40.0/338.0)*100.0/2.0
x_low = str(float(args.source_x) - radius)
x_high = str(float(args.source_x) + radius)
y_low = str(float(args.source_y) - radius)
y_high = str(float(args.source_y) + radius)
binsize = '0.029979674'
binstr = '[bin x='+x_low+':'+x_high+':'+binsize+',y='+y_low+':'+y_high+':'+binsize+'][energy=300:8000]'

#----set up output image----
evtimage = args.evtfile+'_300-8000_img'

#----bin events file----

# old way using ciao python tools
#crt.dmcopy(args.evtfile+binstr,evtimage,clobber=args.clobber)

# new way, using subprocess to access ciao tools
sp.check_call(["dmcopy",args.evtfile+binstr,evtimage,"clobber="+args.clobber])

#---------------------------------------
# Copy header into new model image file
#---------------------------------------


# ciao not compatible with astropy, and no longer has pyfits
print "Transfer Header using anaconda/astropy as follows: \n"
print "from astro_utilities import transfer_header"
print "evtimage = '"+evtimage+"'"
print "modelfile = '"+args.modelfile+"'"
print "newfile = '"+args.newfile+"'"
#print "stuff = transfer_header(evtimage,modelfile,newfile)"

stuff = transfer_header(evtimage,args.modelfile,args.newfile)
