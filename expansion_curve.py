#! /usr/bin/env python

#Author: Kari A. Frank
#Date: June 16, 2014
#Purpose: Plot and fit the radial expansion curve of SN1987A.
#Usage: expansion_curve.py [--infile infile] [--plotmp plotmp] [--plotmcmc plotmcmc] [--mcmcfile mcmcfile] [--clobber clobber]
#
#Input:
#
# infile:     optionally specify file containing the obsids, radii,
#             and errors.  default is infile = 
#             '/export/bulk/rice3/kaf33/Chandra_Observations/SN1987A/comparison_dir/imaging/sn1987a_radii_fits.txt'
#
# plotmp:    optionally plot best fit model from mpfit least-squares fit. 
#               default='yes'
#
# plotmcmc:    optionally plot best fit model from the mcmc fit.  default='no'
#
# mcmcfile     if plotmcmc='yes', mcmc best fit parameters will be read
#              from this file.  default = '/export/bulk/rice3/kaf33/Chandra_Observations/SN1987A/comparison_dir/imaging/expansion_fits/expansion_mcmc_three_2000000-10000000_results.txt'
#
# clobber:    specifies whether files should be overwritten
#             if they already exist (same as in ciao tools,
#             default = 'no')
#
#Output:
# - pdf file of the expansion curve plot
#
#Usage Notes:
# - requires as input a data file with the obsids, radii, and errors,
#   in the format created by get_sn1987a_radii_results.py

#---------------------------------------
#           Import Modules
#---------------------------------------
import argparse,os,sys
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
import fitting as fit

#---------------------------------------
#   Parse arguments and set defaults
#---------------------------------------

parser = argparse.ArgumentParser(description='Plot and fit SN1987A expansion curve.')

pwd = os.getcwd()

#parser.add_argument('required_arg',help='Description and usage for required_arg.',default='default for required_arg')

parser.add_argument('--infile',help='File containing radius data.',default='/export/bulk/rice3/kaf33/Chandra_Observations/SN1987A/comparison_dir/imaging/sn1987a_radii_fits.txt')

parser.add_argument('--plotmp',help='Plot best fit MPFIT model.',default='yes')

parser.add_argument('--plotmcmc',help='Plot best fit MCMC model.',default='no')

parser.add_argument('--mcmcfile',help='File from which to read the MCMC results.',default='/export/bulk/rice3/kaf33/Chandra_Observations/SN1987A/comparison_dir/imaging/expansion_fits/expansion_mcmc_three_2000000-10000000_results.txt')

parser.add_argument('--clobber',help='Overwrite existing files.',default='no')

args = parser.parse_args()

#--constant to convert radius measurements into arcseconds
r_to_arcsec = 0.0147928866

#--constants to convert velocity into km/s--
kpc_to_km = 1000.0*3.0857*10.0**13.0
distance_kpc = 51.4 #Panagia2003
distance_km = distance_kpc*kpc_to_km
arcsec_to_radian = np.pi/(3600.0*180.0)
arcsec_to_km = arcsec_to_radian*distance_km
days_to_s = 24.0*3600.0

arcsecdays_to_kms = arcsec_to_km/days_to_s

#----Set file paths----
#-output plot file-
pfile = '/export/bulk/rice3/kaf33/Chandra_Observations/SN1987A/comparison_dir/imaging/expansion_fits_sept2014/expansion_curve.pdf'
#-file containing fit parameters-
mpfile = '/export/bulk/rice3/kaf33/Chandra_Observations/SN1987A/comparison_dir/imaging/expansion_fits_sept2014/mpfit_results_3.txt'

#---------------------------------------
#        Read and Convert Data
#---------------------------------------

#----Read Raw Data from File----
obsids,ages,r0s,r0err_low_raw,r0err_upp_raw,counts,deconcounts = np.loadtxt(args.infile,unpack=True,skiprows=1)

#----Calculate Errors----
r0err_lows = np.sqrt( r0err_low_raw**2.0 + r0s**2.0/counts  )
r0err_upps = np.sqrt( r0err_upp_raw**2.0 + r0s**2.0/counts  )

#----Convert to Arcsec----
r0_arcsec = r0s*r_to_arcsec
r0err_low_arcsec = r0err_lows*r_to_arcsec
r0err_upp_arcsec = r0err_upps*r_to_arcsec


#---------------------------------------
#        Read MPFIT results
#---------------------------------------

if args.plotmp == 'yes':
    mp_results = np.loadtxt(mpfile,comments='#',skiprows=1,usecols=(1,2,3,4))
    mppars = np.array([0,mp_results[0,0]/arcsecdays_to_kms,mp_results[0,1]/arcsecdays_to_kms,mp_results[0,2],mp_results[0,3]])
    mp_modely = fit.arr_brokenlinear(range(11000),mppars)

#---------------------------------------
#        Read MCMC results
#---------------------------------------

if args.plotmcmc == 'yes':
    mcmc_results = np.loadtxt(args.mcmcfile,comments='#',skiprows=2,usecols=(1,2,3,4))
    mcmcpars = np.array([0,mcmc_results[0,0]/arcsecdays_to_kms,mcmc_results[0,1]/arcsecdays_to_kms,mcmc_results[0,2],mcmc_results[0,3]])
    mcmc_modely = fit.arr_brokenlinear(range(11000),mcmcpars)


#---------------------------------------
#       Plot Expansion Curve
#---------------------------------------

#----Set Up Plot----

pdffile = PdfPages(pfile)
#fig = plt.figure()
fig,ax=plt.subplots()

ax.set_xticks(ages,minor=False)
ax.xaxis.grid(True,which='major')
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(10)
    tick.label.set_rotation('vertical')
plt.subplots_adjust(bottom=0.13)

plt.xlim(4000,11000)
plt.ylim(0.5,0.9)

#--add labels--
plt.xlabel('SN1987A Age [days]')
plt.ylabel('Radius [arcsec]')

#----Plot Meaurements----

plt.errorbar(ages,r0_arcsec,yerr=[r0err_low_arcsec,r0err_upp_arcsec],fmt='.')

plt.scatter(ages,r0_arcsec,marker='o')#,s=deconcounts/100)

#----Plot Best Fit MP Model----
if args.plotmp == 'yes':
    plt.plot(range(11000),mp_modely,linestyle='-',color='green')

#----Plot Best Fit MCMC Model----
if args.plotmcmc == 'yes':
    plt.plot(range(11000),mcmc_modely,linestyle='--',color='green')

#----Close Plot----
pdffile.savefig()
plt.close()
pdffile.close()

#---------------------------------------
#       Print out final status
#---------------------------------------

