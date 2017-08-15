#! /usr/bin/env python

#Author: Kari Frank
#Date: June 14, 2013
#Purpose: Run ciao tools to filter out flares.
#
#Usage: chandra_flare.py [--output_dir OUTPUT_DIR]] [--evt2_in EVT2_IN] [--evt2_out EVT2_OUT] [--exclude_region] [--sigma sigma] [--timebin timebin] [--method method] [--clobber CLOBBER]
#
#Input:
#
# output_dir: path to directory where output files will 
#             be saved (default = $PWD)
#
# evt2_in:    optionally provide the evt2 file to be 
#             filtered. if not provided will choose the 
#             first file (in pwd) it finds with 'evt2.fits' 
#             in the file name.
#
# evt2_out:   optionally provide a name for the output
#             evt2 file. default is 'evt2_deflared.fits'
# 
# exclude_region: region file containing the regions to exclude
#                 (e.g. bright sources)(default=unbkg.reg)
#
# sigma:      points falling outside meancountrate+-sigma*stdevcountrate
#             are removed. same as (and given directly to) lc_sigma_clip
#             and lc_clean. default is 3 for lc_sigma_clip and empty for
#             lc_clean.
#
# timebin:    size of the time bin for creating the light curves, in 
#             seconds.  default=50
#
# method:     choose lc_sigma_clip or lc_clean. default is lc_sigma_clip
#
# clobber:       specifies whether files should be overwritten
#                if they already exist (same as in ciao tools,
#                default = 'no')
#
#Output:
# - new filtered evt2 file in output_dir
#
#Usage Notes:
# - CIAO must be initialized prior to running this script (heainit, setciao)
# - This script uses the CIAO installation of python
# - Assumes correct bad pixel file has already been set using acis_set_ardlib
# - Assumes observation was ACIS
# - Input evt2 file must be in the working directory unless full path is 
#   given in evt2_in
# - For default settings, run from within the repro directory
#
# - May get shared libraries errors if running with ciao 4.7 or later.  Try
#   initializing with ciao 4.6 in terminal before running.
#   source /usr/common/defaults/cshrc-modules/ciao-4.6.1.cshrc

#----Import Modules----
import argparse,os
import ciao_contrib.runtool as crt #all functions should be prefixed with crt
import home_grown as hg #import my own custom functions
import lightcurves

#---------------------------------------
#   Parse arguments and set defaults
#---------------------------------------

parser = argparse.ArgumentParser(description='Remove flares from a Chandra evt2 file.',epilog='NOTE: Must initialize CIAO before running this script ("heainit" then "setciao").')

pwd = os.getcwd()+'/'
parser.add_argument('--output_dir',help='Path to directory in which to save filtered evt2 file.',default=pwd)

parser.add_argument('--evt2_in',help='Input evt2 file',default='default')

parser.add_argument('--evt2_out',help='Output evt2 file.',default='evt2_deflared.fits')

parser.add_argument('--exclude_region',help='Region file to exclude.',default='none')

parser.add_argument('--sigma',help='Sigma clipping range for lc_sigma_clip.',default=None)

parser.add_argument('--timebin',help='Size of time bin for light curves',default=50)

parser.add_argument('--method',help='Choose light curve filtering method, either lc_sigma_clip or lc_clean.',default='lc_sigma_clip')

parser.add_argument('--clobber',help='Overwrite existing files.',default='no')

args = parser.parse_args()

#--rename args--
if args.evt2_in == 'default':
    evt2_infile = hg.fetch_file(pwd,pat='repro_evt2.fits',prompt='no')
else:
    evt2_infile = args.evt2_in

if args.sigma != None: args.sigma=float(args.sigma)

evt2_outfile = args.output_dir+args.evt2_out
args.timebin=int(args.timebin)

#--Print out information--
print "\nInput evt2 file: "+evt2_infile
print "Output evt2 file: "+evt2_outfile
print "Exclude region file: "+args.exclude_region
print ""

#---------------------------------------
#             Remove flares
#---------------------------------------

#--Create background light curve--

#timebin = 200 #eventually make this an optional argument
#how to decide on timebin size?
bkg_lc_file = args.output_dir+'bkg_lc.fits'

#- filter on energy -
crt.dmcopy.punlearn()
estr = evt2_infile+'[energy=2400:6000]'
crt.dmcopy(estr,evt2_infile+'_2400-6000')

#- extract light curve -
crt.dmextract.punlearn()
if args.exclude_region != 'none':
    evt2_file_str = evt2_infile+'[exclude sky=region('+args.exclude_region+')][bin time=::'+str(args.timebin)+']'
else:
    evt2_file_str = evt2_infile+'_2400-6000[bin time=::'+str(args.timebin)+']'

crt.dmextract(evt2_file_str,outfile=bkg_lc_file,opt='ltc1',clobber=args.clobber)

#- remove energy filtered event file -
os.remove(evt2_infile+'_2400-6000')

#--Create GTI file--
gti_outfile = args.output_dir+'flare_gti.fits'
if args.method == 'lc_sigma_clip': 
    if args.sigma == None: args.sigma=3
    lightcurves.lc_sigma_clip(bkg_lc_file,outfile=gti_outfile,sigma=args.sigma)
else:
    lightcurves.lc_clean(bkg_lc_file,outfile=gti_outfile,sigma=args.sigma)

#--Apply GTI file--
evt2_file_str = evt2_infile+'[@'+gti_outfile+']'
crt.dmcopy(evt2_file_str,evt2_outfile,opt='all',clobber=args.clobber)

#--Verify flare removal--

#create filtered light curve
#need to check curve visually
crt.dmextract.punlearn()
evt2_file_str = evt2_outfile+'[bin time=::'+str(args.timebin)+']'
crt.dmextract(evt2_file_str,outfile=args.output_dir+'deflared_lc.fits',opt='ltc1',clobber=args.clobber)

raw_input("Press any key to quit ...")
