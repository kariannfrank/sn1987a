#! /usr/astro/ciao-4.5/bin/python

#Author: Kari A. Frank
#Date: January 15, 2015
#Purpose: Create subband event files from the input file, in 
#         the 'standard' 5 sn1987a bands.
#        
#
#Usage: make_subband_evts.py evtfile [--outdir outdir] [--clobber CLOBBER]
#
#Input:
#
# evtfile:   REQUIRED event file, typically the output of 
#            acis_process_events (usually via chandra_repro)
#
# outdir: OPTIONAL provide directory to save all output files. 
#         default is current working directory.
#
# clobber:    specifies whether files should be overwritten
#             if they already exist (same as in ciao tools,
#             default = 'no')
#
#Output:
# - output of multiple calls of dmcopy, filtering on the following bands:
#        300-800, 800-1200, 1200-8000,300-8000, 3000-8000
#   filenames are the same as the input file, with energy range appended
#
# - subband_evtfiles.lis, stack file containing list of all 
#   the created subband event files
#
#Usage Notes:
# - In future, might want to add option to input arbitrary bands
#   

#---------------------------------------
#           Import Modules
#---------------------------------------
import argparse,os
import ciao_contrib.runtool as crt #all functions should be prefixed with crt
#import chandra
#import home_grown as hg #import my own custom functions
#from astropy.io import fits

#---------------------------------------
#   Parse arguments and set defaults
#---------------------------------------

parser = argparse.ArgumentParser(description='Call dmcopy multiple times to filter event file on 5 energy bands.')

pwd = os.getcwd()

parser.add_argument('evtfile',help='Input event file')

parser.add_argument('--outdir',help='Directory to store output files.',default=pwd)

parser.add_argument('--clobber',help='Overwrite existing files.',default='no')

args = parser.parse_args()

verbose = 1

#--- Get file name without extension ---
evtname = os.path.splitext(args.evtfile)[0]

#---------------------------------------
#     Create input stack file
#---------------------------------------

#-- open file --
if os.path.isfile('subbands_in.lis'):
    if args.clobber == 'yes': os.remove('subbands_in.lis')

with open('subbands_in.lis','w') as inlis:

#-- write evtfile name with filters --
    inlis.write(args.evtfile+'[energy=300:800]\n')
    inlis.write(args.evtfile+'[energy=800:1200]\n')
    inlis.write(args.evtfile+'[energy=300:1200]\n')
    inlis.write(args.evtfile+'[energy=1200:8000]\n')
    inlis.write(args.evtfile+'[energy=300:8000]\n')
    inlis.write(args.evtfile+'[energy=500:2000]\n')
    inlis.write(args.evtfile+'[energy=3000:8000]\n')
    inlis.write(args.evtfile+'[energy=2000:10000]')

#---------------------------------------
#     Create output stack file
#---------------------------------------

if os.path.isfile('subbands_out.lis'):
    if args.clobber == 'yes': os.remove('subbands_out.lis')

#-- open file --
with  open('subbands_out.lis','w') as outlis:

#-- write evtfile name with filters --
    outlis.write(evtname+'_300-800.fits\n')
    outlis.write(evtname+'_800-1200.fits\n')
    outlis.write(evtname+'_300-1200.fits\n')
    outlis.write(evtname+'_1200-8000.fits\n')
    outlis.write(evtname+'_300-8000.fits\n')
    outlis.write(evtname+'_500-2000.fits\n')
    outlis.write(evtname+'_3000-8000.fits\n')
    outlis.write(evtname+'_2000-10000.fits')

#---------------------------------------
#          Run dmcopy
#---------------------------------------

crt.dmcopy.punlearn()
crt.dmcopy.infile = '@subbands_in.lis'
crt.dmcopy.outfile = '@subbands_out.lis'
crt.dmcopy(clobber=args.clobber)

#---------------------------------------
#        Delete/rename lis files
#---------------------------------------

os.remove('subbands_in.lis')
#os.remove('subbands_out.lis')
os.rename('subbands_out.lis','subband_evtfiles.lis')
