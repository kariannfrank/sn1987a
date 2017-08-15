#!/usr/astro/ciao-4.6/bin/python

####!/usr/bin/env python

#Author: Kari A. Frank
#Date: July 9, 2014
#Purpose: Update spectrum file headers with background, rmf, and arf files
#
#Usage: update_spec_headers.py obsid --clobber clobber
#
#Input:
#
# obsid -- obsid (including leading zeros)
# 
#Output:
# - updates piled and unpiled spectrum file headers
# - copied piled spectrum with the associated arf that
#   uses no contamination model
#
#Usage Notes:

#---------------------------------------
#           Import Modules
#---------------------------------------
import argparse
import ciao_contrib.runtool as crt #all functions should be prefixed with crt
import pyfits as fits
#import astropy.io.fits as fits
#import astro_utilities as au

#---------------------------------------
#   Parse arguments and set defaults
#---------------------------------------

parser = argparse.ArgumentParser(description='Update spectrum file header backfile, ancr, and respfile keywords.',epilog='')

#pwd = os.getcwd()

parser.add_argument('obsid',help='Observation ID of current observation (including leading zeros).')

parser.add_argument('--clobber',help='Overwrite existing files.',default='no')

args = parser.parse_args()

#----Set file paths----



#---------------------------------------
#           Set File Names
#---------------------------------------

#----Spectrum Files----
unpiled_spec = 'obs'+args.obsid+'_unpiled_5grp.pi'
piled_spec = 'obs'+args.obsid+'_src_5grp.pi'
nocontam_spec = 'obs'+args.obsid+'_unpiled_5grp_nocontam.pi'

#----rmf, arf, background files----
bkg_file = 'obs'+args.obsid+'_bkg.pi'
arf_file = 'obs'+args.obsid+'.arf'
rmf_file = 'acisf'+args.obsid+'.rmf'
nocontam_arf_file = 'obs'+args.obsid+'_nocontam.arf'

print bkg_file
print arf_file
print rmf_file
print nocontam_arf_file

#---------------------------------------
#  Update piled and unpiled headers
#---------------------------------------

#----Unpiled Spectrum----
#hdr = fits.getheader(unpiled_spec,1)

#hdr['ANCRFILE'] = arf_file
#hdr['RESPFILE'] = rmf_file
#hdr['BACKFILE'] = bkg_file

hdulist = fits.open(unpiled_spec,mode='update')
shdu = hdulist[1]
#print 'Unpiled AREASCAL = ',shdu.header['AREASCAL']
shdu.header['ANCRFILE'] = arf_file
shdu.header['RESPFILE'] = rmf_file
shdu.header['BACKFILE'] = bkg_file
shdu.header['AREASCAL'] = 1.0
hdulist.close(output_verify='ignore') #otherwise yells about bad AREASCAL

#with fits.open(unpiled_spec,mode='update',output_verify='ignore') as hdu1:
#    hdr = hdu1[1].header
#    print hdr['AREASCAL']
#    hdr['ANCRFILE'] = arf_file
#    hdr['RESPFILE'] = rmf_file
#    hdr['BACKFILE'] = bkg_file
#    hdu1[1].header['ANCRFILE'] = arf_file
#    hdu1[1].header['RESPFILE'] = rmf_file
#    hdu1[1].header['BACKFILE'] = bkg_file

#----Piled Spectrum----    
with fits.open(piled_spec,mode='update') as hdu2:
    hdr = hdu2[1].header
    hdr['ANCRFILE'] = arf_file
    hdr['RESPFILE'] = rmf_file
    hdr['BACKFILE'] = bkg_file
    
#fits.writeto(unpiled_spec,)

#---------------------------------------
#      Update no contamination
#---------------------------------------

skip = 1
if skip == 0:

#----Copy Updated UnPiled Spectrum----
    crt.dmcopy(unpiled_spec,nocontam_spec,clobber=args.clobber)

    #with fits.open(nocontam_spec,mode='update') as hdu1:
    #    hdr = hdu1[1].header
    #    print 'NoContam AREASCAL = ',hdr['AREASCAL']
    #    hdr['ANCRFILE'] = nocontam_arf_file
    
    hdulist = fits.open(nocontam_spec,mode='update')
    shdu = hdulist[1]
    #print 'Unpiled AREASCAL = ',shdu.header['AREASCAL']
    shdu.header['ANCRFILE'] = nocontam_arf_file
    hdulist.close(output_verify='ignore') #otherwise yells about bad AREASCAL
