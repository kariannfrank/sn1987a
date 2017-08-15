#! /usr/bin/env python

#Author: Kari A. Frank 
#Date: January 16, 2015
#Purpose: compile total image counts from sn1987a deconvolved images
#Usage: get_image_counts.py infile [--outfile outfile] [--clobber CLOBBER]
#
#Input:
#
# infile: REQUIRED. file name, either a single image/event file 
#         or stack file listing multiple image files (must end .lis if stack)
#          
#
# outfile: output file name, containing the total number of counts
#          in each input image.  format is two columns, 
#          'imagefile  counts'.  default name is the base of infile
#          parameter with _counts.txt appended (e.g.
#          infile='myimages.lis' => outfile='myimages_counts.txt')
#          
# clobber:    specifies whether files should be overwritten
#             if they already exist (same as in ciao tools,
#             default = 'no')
#
#Output:
# - outfile (see above)
#
#Usage Notes:

#---------------------------------------
#           Import Modules
#---------------------------------------
import argparse,os
from home_grown import read_list
import numpy as np
from astropy.io import fits
#import pyfits as fits
#import home_grown as hg #import my own custom functions

#---------------------------------------
#   Parse arguments and set defaults
#---------------------------------------

parser = argparse.ArgumentParser(description='Retrieve total counts from sn1987a image files.')

pwd = os.getcwd()

parser.add_argument('infile',help='Input images file (or list of image files).')

parser.add_argument('--outfile',help='Name of output text file.',default='default')

parser.add_argument('--clobber',help='Overwrite existing files.',default='no')

args = parser.parse_args()

if args.outfile == 'default':
    args.outfile = os.path.splitext(args.infile)[0]+'_counts.txt'

if os.path.splitext(args.infile)[1] == '.lis':
    batchmode = True
else:
    batchmode = False

#----Set file paths----

#--Get Image File List--

if batchmode == True:
    imgfiles = read_list(args.infile,comment='#')
else:
    imgfiles = [args.infile]

#--Open Output Text File--
if os.path.isfile(args.outfile) and args.clobber == 'no':
    print args.outfile+" exists and clobber=no.  Printing to screen only."
    tofile = False
else:
    if os.path.isfile(args.outfile): os.remove(args.outfile) # clobber file
    countfile = open(args.outfile,'w')
    tofile = True

#---------------------------------------
#        Loop Over Observations
#---------------------------------------

for img in imgfiles:

    #--Read in image as array--
    imagefile = fits.open(img)
    image = np.array(imagefile[0].data)

    #--Get total counts--
    counts = np.sum(image)

    #--Close image file--
    imagefile.close()

    #--Write to outfile--
    if tofile is True:
        countfile.write(img+'\t'+str(int(counts)))
    #--Add newline to countfile--
        if img != imgfiles[-1]: countfile.write('\n')

    #--Print to screen--
    print img+'\t'+str(int(counts))

#--close counts file--
if tofile is True:
    countfile.close()

#---------------------------------------
#       Print out final status
#---------------------------------------

print '\nFinished getting counts from '+args.infile+'.\n'
