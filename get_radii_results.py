#!/usr/bin/env python

####! /usr/astro/ciao-4.6/bin/python

#Author: Kari A. Frank 
#Date: June 16, 2014
#Purpose: Get radius fit results (radii and errors) and save to file.
#Usage: get_radii_results.py [--results_dir results_dir][--band band][--clobber CLOBBER]
#
#Input:
#
# results_dir: directory containing fpars files
#          
# clobber:    specifies whether files should be overwritten
#             if they already exist (same as in ciao tools,
#             default = 'no')
#
#Output:
# - radii_fits.txt (contains r0,r0err,counts from fpars files and counts from deconvolved images)
#
#Usage Notes:

#---------------------------------------
#           Import Modules
#---------------------------------------
import argparse,os,sys,re
#import ciao_contrib.runtool as crt #all functions should be prefixed with crt
from home_grown import ls_to_list
from astropy.io import fits
#import pyfits as fits
#import home_grown as hg #import my own custom functions
import numpy as np
import sn1987a_time as snt
#from sn1987a_time import get_obs_time
from math import sqrt

#---------------------------------------
#   Parse arguments and set defaults
#---------------------------------------

parser = argparse.ArgumentParser(description='Extract results from new_sn1987a_image_analysis_(parallel).pro for a batch of observations.')

pwd = os.getcwd()

parser.add_argument('--results_dir',help='Directory containing the fpars files.',default='lobes/')
#parser.add_argument('band',help='String containing the energy band.',default='300-8000')
parser.add_argument('--band',help='String containing the energy band (e.g. "300-8000").',default='')

parser.add_argument('--clobber',help='Overwrite existing files.',default='no')

args = parser.parse_args()

#----Handle Files----

#--Get fpars File List--
fparsfiles = ls_to_list(args.results_dir,ls_args="*"+args.band+"*fpars*")
print fparsfiles

#--Check for regions directory--
if not os.path.exists(args.results_dir+'/regions'):
    os.mkdir(args.results_dir+'/regions')

#--Open Output Text File--
outfile = args.results_dir+args.band+'radii_fits.txt'
if os.path.isfile(outfile) and args.clobber == 'no':
    print outfile+" exists and clobber=no.  Quitting."
else:
    radfile = open(outfile,'w')
    radfile.write('obsid\tage\tband\tR0\tR0errl\tR0erru\tSIGR\tSIGRerrl\tSIGRerru\tcounts\tdeconcounts\tSKYNS\tSKYNSerrl\tSKYNSerru\n')

    #--Get counts from spectrum (input)--
#    specfile = '/export/bulk/rice3/kaf33/Chandra_Observations/SN1987A/comparison_dir/imaging/sn1987a_decon_counts.txt'
#    specobsids,speccounts = np.loadtxt(specfile,dtype='str',unpack=True)

#---------------------------------------
#     Loop Over Observations (files)
#---------------------------------------

    for f in range(len(fparsfiles)):
        
        #----Read File----
        thepath = args.results_dir+'/'
        thefile = fparsfiles[f]
        print thepath+thefile
        fitsfile = fits.open(thepath+thefile)
        head = fits.getheader(thepath+thefile)

        #----Get Table----
        table = fitsfile[1].data
        
        #----Get Values----
        R0 = table.field('R0')
        R0err = table.field('R0ERR') #two elements
        counts = table.field('COUNTS')
        SIGR = table.field('SIGR')
        SIGRerr = table.field('SIGRERR') #two elements
        X0 = table.field('X0')
        Y0 = table.field('Y0')

        #----Convert errors----
        R0err_low = sqrt(R0err[0,0]**2.0+R0[0]**2.0/float(counts))
        R0err_high = sqrt(R0err[0,1]**2.0+R0[0]**2.0/float(counts))
        SIGRerr_low = sqrt(SIGRerr[0,0]**2.0+SIGR[0]**2.0/float(counts))
        SIGRerr_high = sqrt(SIGRerr[0,1]**2.0+SIGR[0]**2.0/float(counts))

        #----Convert to arcsec----
        R0_arcsec = float(snt.convert_units(R0,fromunit='r',tounit='arcsec'))
        R0err_arcsec_low = float(snt.convert_units(R0err_low,fromunit='r',tounit='arcsec'))
        R0err_arcsec_high = float(snt.convert_units(R0err_high,fromunit='r',tounit='arcsec'))
        SIGR_arcsec = float(snt.convert_units(SIGR,fromunit='r',tounit='arcsec'))
        SIGRerr_arcsec_low = float(snt.convert_units(SIGRerr_low,fromunit='r',tounit='arcsec'))
        SIGRerr_arcsec_high = float(snt.convert_units(SIGRerr_high,fromunit='r',tounit='arcsec'))
        
        #----Calculate (Projected and Rotated to North-South) Semiminor 
        #         Axis and Errors----
        theta_deg = 43.0 # inclination angle in degrees
        theta_rad = theta_deg*np.pi/180.0 
        phi_deg = 81.2
        phi_rad = phi_deg*np.pi/180.0
        sinratio = np.sin(theta_rad)/np.sin(phi_rad)
        A = R0_arcsec*sinratio
        Aerr_low = A-(R0_arcsec-R0err_arcsec_low)*sinratio
        Aerr_high = (R0_arcsec+R0err_arcsec_high)*sinratio-A

#        #----Get Input Counts (from Spectrum)----
#        incounts = speccounts[np.where(specobsids == obsid)]
#        print obsid

        #----Get Input Counts from header----
        incounts = head['INCOUNTS']

        #----Get source smoothed image and band from header----
        # assumes file name conforms to standard and includes band
        srcimgfile = head['INFILE']
        fband=re.findall(r'[\d]+-[\d]+',srcimgfile)[0]

        #----Get Obsid----
        obsid = head['OBS_ID']

        #----Write to File----
        radfile.write("{zero}\t{half}\t{half1}\t{one:0.3f}\t{two:0.3f}\t{three:0.3f}\t{four:0.3f}\t{five:0.3f}\t{six:0.3f}\t{seven:0.0f}\t{eight}\t{nine:0.3f}\t{ten:0.3f}\t{eleven:0.3f}\n".format(zero=obsid,half=str(snt.get_obs_time(obsid)),half1=fband,one=R0_arcsec,two=R0err_arcsec_low,three=R0err_arcsec_high,four=SIGR_arcsec,five=SIGRerr_arcsec_low,six=SIGRerr_arcsec_high,seven=counts[0],eight=str(incounts),nine=A,ten=Aerr_low,eleven=Aerr_high))
#        radfile.write(obsid+'\t'+str(R0[0])+'\t'+str(R0err[0,0])+'\t'+str(R0err[0,1])+'\t'+str(counts[0])+'\n')

        #----Write Ellipse Region File----
        # in ds9 format, image coordinates
#        regfile = open(args.results_dir+'/regions/'+os.path.splitext(srcimgfile)[0]+'.reg','w')
        regfile = open(args.results_dir+'/regions/'+os.path.splitext(os.path.basename(srcimgfile))[0]+'.reg','w')
        regfile.write('# Filename: '+srcimgfile+'\n')
        regfile.write('global color=white dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 delete=1 include=1 source=1\n')
        regfile.write('image\n')
        X0ds9 = 100.0+X0[0]
        Y0ds9 = 100.0+Y0[0]
        Amajor = R0[0]
        Aminor = R0[0]*np.sin(theta_rad)
        Amajor_low = R0[0]-R0err_low
        Amajor_high = R0[0]+R0err_high
        Aminor_low = Amajor_low*np.sin(theta_rad)
        Aminor_high = Amajor_high*np.sin(theta_rad)
        regfile.write('ellipse('+str(X0ds9)+','+str(Y0ds9)+','+str(Aminor)+','+str(Amajor)+','+str(phi_deg)+') # color=white')
        regfile.write('\nellipse('+str(X0ds9)+','+str(Y0ds9)+','+str(Aminor_low)+','+str(Amajor_low)+','+str(phi_deg)+') # color=white dash=1')
        regfile.write('\nellipse('+str(X0ds9)+','+str(Y0ds9)+','+str(Aminor_high)+','+str(Amajor_high)+','+str(phi_deg)+') # color=white dash=1')

        regfile.close()
        
#--close file--

    radfile.close()

#---------------------------------------
# Create CIAO Elliptical Region Files
#---------------------------------------



#---------------------------------------
#       Print out final status
#---------------------------------------

print "Wrote or updated "+outfile
