#! /usr/astro/ciao-4.6/bin/python

#Author: Kari A. Frank 
#Date: June 16, 2014
#Purpose: compile source counts from sn1987a spectra of all observations
#Usage: get_counts.py specdir [--clobber CLOBBER]
#
#Input:
#
# specdir: directory containing spectra
#          
#
# OPTIONAL_ARG1:  description
#                 of optional_arg1
#
# clobber:    specifies whether files should be overwritten
#             if they already exist (same as in ciao tools,
#             default = 'no')
#
#Output:
# - sn1987a_src_counts.txt (list counts for all spectra)
#
#Usage Notes:

#---------------------------------------
#           Import Modules
#---------------------------------------
import argparse,os,sys
import ciao_contrib.runtool as crt #all functions should be prefixed with crt
from home_grown import ls_to_list
#import home_grown as hg #import my own custom functions

#---------------------------------------
#   Parse arguments and set defaults
#---------------------------------------

parser = argparse.ArgumentParser(description='Retrieve counts from sn987a source spectra.')

pwd = os.getcwd()

parser.add_argument('specdir',help='Directory containing the source spectra.',default='../spectra/')

parser.add_argument('--outfile',help='Name of output file.',default='sn1987a_src_counts.txt')

parser.add_argument('--quadrant',help='Which quadrant to gather counts from (se,ne,nw,sw,total)',default='total')

parser.add_argument('--clobber',help='Overwrite existing files.',default='no')

args = parser.parse_args()

#----Set file paths----

#--Get Spectra List--
if args.quadrant == 'total':
    spectra = ls_to_list(args.specdir,ls_args="*src_5grp.pi*")
else:
    spectra = ls_to_list(args.specdir,ls_args="*_"+args.quadrant+"_grp.pi")
print spectra

#--Open Output Text File--
if os.path.isfile(args.outfile) and args.clobber == 'no':
    print args.outfile+" exists and clobber=no.  Quitting."
else:
    outfile = args.outfile
    countfile = open(outfile,'w')

#---------------------------------------
#        Loop Over Observations
#---------------------------------------

    for s in range(len(spectra)):
    
#        counts = crt.dmlist(args.specdir+spectra[s],"TOTCTS")
        counts = crt.dmkeypar(args.specdir+spectra[s],"TOTCTS","echo")
        obsid = crt.dmkeypar(args.specdir+spectra[s],"OBS_ID","echo")

        print obsid, counts
        
        countfile.write(obsid+'\t'+counts+'\n')
        
#--close file--

    countfile.close()

#---------------------------------------
#       Print out final status
#---------------------------------------

