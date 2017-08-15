#!/usr/bin/env python

###!/usr/astro/ciao-4.7/bin/python

#Author: Kari A. Frank
#Date: March 4, 2014
#Purpose: Read in deconvolved, smoothed Chandra images of SNR 1987A and scale them all the same scale based on their 0.3-8.0 keV fluxes.
#Usage: flux_scale_images.py [--clobber CLOBBER]
#
#Input:
#                
# clobber       -- specifies whether files should be overwritten
#                  if they already exist (same as in ciao tools,
#                  default = 'no')
#
#Output:
# - output1
#
#Usage Notes:
# - note1

#---------------------------------------
#           Import Modules
#---------------------------------------
import argparse,os,sys,datetime
#import ciao_contrib.runtool as crt #all functions should be prefixed with crt
import home_grown as hg #import my own custom functions
import numpy as np
import pandas as pd
from astropy.io import fits
#import pyfits as fits

#---------------------------------------
#   Parse arguments and set defaults
#---------------------------------------

parser = argparse.ArgumentParser(description='Scale smoothed images according to 0.5-8.0 keV flux.',epilog='NOTE: Extra note.')

pwd = os.getcwd()

#parser.add_argument('image_file_in',help='File containing list of image files to read in.')

parser.add_argument('--band',help='Image band.',default='300-8000')

parser.add_argument('--clobber',help='Overwrite existing files.',default='no')

args = parser.parse_args()

#----Set file paths----

#sn1987a_dir = '~/Dropbox/Research/SN1987A/comparison_dir/'
#sn1987a_dir = '/export/bulk/rice3/kaf33/Chandra_Observations/SN1987A/comparison_dir/spectra462/'
#sn1987a_dir = '/export/bulk/rice1/kaf33/main/Chandra_Observations/SN1987A/comparison_dir/'
sn1987a_dir = '/data/nasdrive/main/kaf33/sn1987a_monitoring/comparison_dir/'

specdir = sn1987a_dir+'spectra_current/'

imagedir = sn1987a_dir+'images_current/smoothed_images/'

#---------------------------------------
#          Read in Flux Values
#---------------------------------------

#obsid,obsdate,age,grating,caldb,modelname,pcounts,pflux_soft,pflux_soft_low,pflus_soft_high,pflux_hard,pflux_hard_low,pflux_hard_high,pflux_broad,pflux_broad_low,pflux_broad_high,counts,flux_soft,flux_soft_low,flux_soft_high,flux_hard,flux_hard_low,flux_hard_high,flux_broad,flux_broad_low,flux_broad_high,kTcool,kTcool_low,kTcool_high,kThot,kThot_low,kThot_high,taucool,taucool_err,tauhot,tauhot_low,tauhot_high,normcool,normcool_low,normcool_high,normhot,normhot_low,normhot_high,chi2,dof,redchi2,ratiosoft,ratiohard,ratiobroad = np.genfromtxt(sn1987a_dir+'spectra_fits.txt',unpack=True,skip_header=1)

#datadf = pd.read_table(specdir+'current_chandra_spectra_fits.txt',sep=r'\s+',comment='#',na_values=['na','None','-100'],engine='python')
datadf = pd.read_table(specdir+'current_chandra_fluxes.txt',sep=r'\s+',comment='#',na_values=['na','None','-100'],engine='python')

#---------------------------------------
#          Read in Images
#---------------------------------------

#files = hg.read_list(args.image_file_in,comment='#')
#num_images = len(files)

#fname = '_repro_evt2_subpix_'+args.band+'fits_smoothed'
fname = '_evt2_subpix_'+args.band+'fits_smoothed'

for i,obrow in datadf.iterrows():
#for i in range(data['obsid'].shape[0]):

    ob = str(obrow['obsid'])
    print 'obsid = '+ob
#    if len(ob) == 3:
#        a = 'acisf00'+ob+fname
#    if len(ob) == 4:
#        a = 'acisf0'+ob+fname
#    if len(ob) == 5:
#        a = 'acisf'+ob+fname
    a = hg.ls_to_list('./','acisf*'+ob+'_repro*'+fname+'.img')
    a = os.path.splitext(a[0])[0]
    print a

    in_file = a + '.img'
    out_file = imagedir+ a + '_flux' + '.fits'

    #-read image-
    if os.path.isfile(in_file):
        fitshdus = fits.open(in_file)
    #    fitsdata.info()
        img = fitshdus[0].data
        # copy header
        new_hdr = fitshdus[0].header.copy()
        fitshdus.close()

    #-get observation broadband or hardband flux-
        flux = obrow['broad']
#        flux = data['broad'][i]
        if args.band == '3000-8000':
            flux = obrow['hard']
#            flux = data['hard'][i]

    #-calculate scale factor-
        tot_img = np.sum(img)
        scale = flux/tot_img
        if args.band == '3000-8000':
            scale = 5.0*scale
        print 'flux,tot_img,scale = ',flux,tot_img,scale

    #--scale the image--
#    print 'in_file = '+in_file
#    print 'out_file = '+out_file
        if tot_img > 0.0:
#            if os.path.isfile(out_file):
#                print out_file+" exists and clobber=no. Not scaling image."
#            else:
            imgout = img*scale
            fits.writeto(out_file,imgout,new_hdr,clobber=args.clobber)
                
                            
            """
            crt.dmimgcalc.punlearn()
            crt.dmimgcalc.infile=in_file
            crt.dmimgcalc.infile2='none'
            crt.dmimgcalc.outfile=out_file
            crt.dmimgcalc.op = 'imgout=img1*'+str(scale)
            crt.dmimgcalc.clobber = args.clobber
            crt.dmimgcalc()
            """

#---------------------------------------
#       Print out final status
#---------------------------------------

print 'Finished Scaling Images'
