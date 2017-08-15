#! /usr/bin/env python

#Author: Kari A. Frank
#Date: June 25, 2014
#Purpose: Plot and fit the light curve and other spectral results of SN1987A.
#Usage: spectra_results.py [--infile_a infile_a] [--clean clean] [--ylog ylog] [--plotfit plotfit] [--dofit dofit] [--helder helder] [--xmm xmm] [--caldb459 caldb459] [--piled piled] [--rosat rosat] [--early early] [--soft soft] [--soft12 soft12] [--hard hard] [--soft13 soft13] [--soft23 soft23] [--broad broad] [--skip462] [--clobber clobber] [--lateonly lateonly] [--agemin agemin] [--agemax agemax] [--instcolor instcolor] [--fluxmax fluxmax] [--ignorechippos ignorechippos] [--official official] [--quadcounts quadcounts]
#
#Input:
#
# infile_b:     optionally specify file with observation and spectra data
#              default is infile_a = 
#             '/export/bulk/rice1/kaf33/main/Chandra_Observations/SN1987A/comparison_dir/spectra459/spectra_fits.txt'
# infile_s:     optionally specify second file with observation and spectra data
#              default is infile_b = 
#             '/export/bulk/rice1/kaf33/main/Chandra_Observations/SN1987A/comparison_dir/spectra462/spectra_fits462.txt'
#
# clean:      optional switch to change the axes labels and ticks to even
#             years and ages (default = 'yes')
#
# ylog:       switch to plot the light curve with log y-axis (default='no')
#
# dofit:      do linear fit of light curve
#
# plotfit:    plot best fit(s) on the contamination absorption plot
#
# helder: plot flux measurements from Helder2012 (default='no')
#
# xmm:    plot flux measurements from XMM (default='yes')]
#
# rosat:  plot flux measurements from ROSAT (default='no')
#
# caldb459:    plot flux measurements from Chandra CALDB 4.5.9 (default='no')
#
# piled:  plot non-pileup corrected flux measurements (default='no')
#
# early:  plot non-pileup corrected flux measurements for observations 
#             1387 and 122 (default='yes')
#
# lateonly: plot only the late-time observations on the light curve
#            (day 8000+)
#
# soft/soft12/soft13/soft23/broad/hard: plot the corresponding band's
#         light curve and/or ratios (default='yes')
#
# official: switch to only include observations in the official list,
#           $sna/comparison_dir/official_obsid_list.txt (default='no')
#
# agemin/agemax: set the minimum and maximum age (in day) to plot
#
# fluxmax: optionally set the maximum of flux axis in light curve 
#            (in 10^-13 erg/cm^2/s)
#
# instcolor: plot different instruments by color and different bands by symbol.
#             only applies to the xmm and chandra (4.6.2) soft and hard bands.
#             (default = 'no')
#
# quadcounts: plot the fraction of counts in each quadrant over time
#
# plotextras: optionally plot information from other wavelengths: optical hot spots,
#             Plait1995 ring size/location, Ng2013 radio break (default=None)
#             options are: 'plait1995','ng2013','larsson2013','hotspots','sugarman2002','all'
#
# clobber:    specifies whether files should be overwritten
#             if they already exist (same as in ciao tools,
#             default = 'no')
#
#Output:
# - pdf file of the expansion curve plot
#
#Usage Notes:


#---------------------------------------
#           Import Modules
#---------------------------------------
import argparse,os,sys
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as dates
import numpy as np
from numpy.lib.recfunctions import append_fields
from sn1987a_time import *
import fitting
#sys.path.append('/astro/research/kaf33/Dropbox/Research/Python_Programs/')
sys.path.append('/Users/kafrank/Dropbox/Research/Python_Programs/')
from fancy_plot import fancy_plot
from scipy.optimize import curve_fit
import analytical_functions as af
import pandas as pd
from sn1987a_plot_helpers import plot_extras

#---------------------------------------
#      Define Helper Functions
#---------------------------------------

def plot_a(ax,age,flux,low=None,high=None,syms='o',colors='b',sizes=4,alphas=1.0,mecs='black',mews=0.5,line='',errorband=False,label='_nolegend_'):        
    if low is not None:
        elow,ehigh = ebars(flux,low,high)
    else:
        elow = None
        ehigh = None

#    fancy_plot(age,flux,yerror=elow,yerror_high=ehigh,syms=syms,colors=colors,sizes=sizes,alphas=alphas,mecs=mecs,errorband=errorband,line=line,mews=mews)
    fancy_plot(age,flux,yerror=np.array(flux-low),yerror_high=np.array(high-flux),syms=syms,colors=colors,sizes=sizes,alphas=alphas,mecs=mecs,errorband=errorband,line=line,mews=mews)

#    ax.errorbar(age,flux,yerr=[elow,ehigh],fmt='.',color=colors[0],markersize=0,markeredgewidth=0)
#    ax.plot(age,flux,syms[0],markersize=4,markeredgewidth=0.5,color=colors[0])

def plot_b(ax,age,flux,flux_low,flux_high,sym='o',color='b',alpha=1):
    ax.errorbar(age,flux,yerr=[flux-flux_low,flux_high-flux],fmt='.',markersize=0,color=color,markeredgewidth=0)
    ax.plot(age,flux,sym,markersize=4,color=color,markeredgewidth=0.5,alpha=alpha)

def ebars(x,x_low,x_high):
    #function to convert error intervals intor error bars
    errlow = x-x_low
    errhigh = x_high - x
    return errlow,errhigh
    
#-function to derive vectors of errors on ratios-
def ratio_err(top,bottom,top_low,top_high,bottom_low,bottom_high):
    #uses simple propagation of errors (partial derivatives)
    #note it returns error interval

    #-make sure input is numpy arrays-
    top = np.array(top).astype(float)
    top_low = np.array(top_low)
    top_high = np.array(top_high)
    bottom = np.array(bottom).astype(float)
    bottom_low = np.array(bottom_low)
    bottom_high = np.array(bottom_high)
    
    #-calculate errorbars-
    top_errlow = np.subtract(top,top_low)
    top_errhigh = np.subtract(top_high,top)
    bottom_errlow = np.subtract(bottom,bottom_low)
    bottom_errhigh = np.subtract(bottom_high,bottom)

    #-calculate ratio_low-
    ratio_low = top/bottom - ( (top_errlow/bottom)**2.0 + (top/(bottom)**2.0*bottom_errlow)**2.0  )**0.5
#    ratio_low = np.divide(top,bottom) - ( np.divide(top_errlow,bottom)**2.0 + (np.divide(top,bottom**2.0)*bottom_errlow)**2.0  )**0.5
    #-calculate ratio_high-
    ratio_high = ( (top_errhigh/bottom)**2.0 + (top/(bottom)**2.0*bottom_errhigh)**2.0  )**0.5 + top/bottom


    # return two vectors, err_low and err_high
    return ratio_low,ratio_high

def set_axis(ax,x=None,twin=False,title=None,xlab=None,grid=False,
             clean='yes'):

#    font_size = 10
    font_size = 13
#    fontangle = 50
    fontangle = 40

    if (twin == False) and (clean == 'yes'):
        fontangle = 0

    if clean == 'no':
        font_size = 9

    if twin == False:
        ax0 = ax
        ax0.xaxis.grid(grid,which='major')
        if x is not None:
            ax0.set_xticks(x,minor=False)
    else:
        ax0 = ax.twiny()
#        x_min,x_max = plt.xlim()
        ax0.set_xlim(ax.get_xlim())
        if x is not None:
            ax0.set_xticks(x,minor=False)
        if xlab is not None:
            ax0.set_xticklabels(xlab)
        
    if title is not None:
        ax0.set_xlabel(title)

    for tick in ax0.xaxis.get_major_ticks():
        if twin == False:
            tick.label.set_fontsize(font_size)
            tick.label.set_rotation(fontangle)
        if twin == True:
            tick.label2.set_fontsize(font_size)
            tick.label2.set_rotation(fontangle)

    #ax.set_yscale('log')

def start_plot(xtitle='',ytitle='',xmin=None,xmax=None,ymin=None,ymax=None,
               ylog=False,clean='yes'):

    #-initialize main plot-
    fig,ax1=plt.subplots()

    #-set axis titles-
    plt.xlabel(xtitle)
    plt.ylabel(ytitle)

    #-set bottom axis-
    if ylog == True: plt.yscale('log')
    set_axis(ax1,x=botx,grid=mgrid,clean=clean)
    if xmin is not None or xmax is not None:
        ax1.set_xlim(xmin,xmax)
    if ymin is not None or ymax is not None:
        ax1.set_ylim(ymin,ymax)

    #-add extra space on bottom-
    plt.subplots_adjust(bottom=bott,top=topp)

    return ax1

#---------------------------------------
#   Parse arguments and set defaults
#---------------------------------------

parser = argparse.ArgumentParser(description='Plot and fit SN1987A light curve and other spectral results.')

pwd = os.getcwd()

#parser.add_argument('required_arg',help='Description and usage for required_arg.',default='default for required_arg')

#rice
sn1987a_dir = '/data/nasdrive/main/kaf33/sn1987a_monitoring/'
#brevis
#sn1987a_dir = '/Users/kafrank/Research/SN1987A/'

parser.add_argument('--infile_b',help='File containing observation and spectra data.',default='')

parser.add_argument('--infile_a',help='File containing observation and spectra data.',default=sn1987a_dir+'comparison_dir/spectra_current/spectra465_frank2015a_fits_official.txt')

parser.add_argument('--clean',help='Determines x-axis ticks and labels.',default='yes')

parser.add_argument('--ylog',help='Plot light curve with log scale y-axis.',default='no')

parser.add_argument('--plotfit',help='Plot best fit contamination.',default='no')

parser.add_argument('--dofit',help='Run mcmc on contamination absorption.',default='no')

parser.add_argument('--helder',help='Plot fluxes from Helder2012',default='no')

parser.add_argument('--xmm',help='Plot fluxes from XMM EPIC-pn',default='yes')

parser.add_argument('--rosat',help='Plot fluxes from ROSAT',default='no')

parser.add_argument('--caldb459',help='Plot fluxes from Chandra CALDB 4.5.9',default='no')

parser.add_argument('--piled',help='Plot non-pileup corrected fluxes.',default='no')

parser.add_argument('--early',help='Plot (non-pileup corrected) fluxes of observations 1387 and 122.',default='yes')

parser.add_argument('--lateonly',help='Plot only day 8000+ on light curve.',default='no')

parser.add_argument('--soft',help='Plot 0.5-2.0 keV band light curve.',default='yes')
parser.add_argument('--soft12',help='Plot 1.0-2.0 keV band light curve.',default='yes')
parser.add_argument('--soft13',help='Plot 1.0-3.0 keV band light curve.',default='no')
parser.add_argument('--hard',help='Plot 3.0-8.0 keV band light curve.',default='yes')
parser.add_argument('--broad',help='Plot 0.5-8.0 keV band light curve.',default='yes')
parser.add_argument('--soft23',help='Plot 2.0-3.0 keV band light curve.',default='no')

parser.add_argument('--skip462',help='Do not plot the caldb462 light curve.',default='no')

parser.add_argument('--official',help='Include only obsids in the official list.',default='no')

parser.add_argument('--agemin',help='Minimum age for x-axes.',default=0)
parser.add_argument('--agemax',help='Maximum age for x-axes.',default=0)

parser.add_argument('--fluxmax',help='Maximum flux for light curve y-axis.',default=0)

parser.add_argument('--instcolor',help='Switch meaning of symbols and colors.',default='no')

parser.add_argument('--quadcounts',help='Plot fraction of counts per quadrant.',default='yes')

parser.add_argument('--ignorechippos',help='Plot normal borders for observations with offset chip positions.',default='no')

parser.add_argument('--plotextras',help='Overplot information from other wavelengths.',default=None)

#not currently implemented!
parser.add_argument('--clobber',help='Overwrite existing files.',default='no')

args = parser.parse_args()


#----Set dependent arguments----

args.agemin = int(args.agemin)
args.agemax = int(args.agemax)
args.fluxmax = float(args.fluxmax)

if args.ylog == 'yes': 
    args.ylog = True
else:
    args.ylog = False

#-check if all chandra is skipped-
skipchandra = 'no'
if args.skip462 == 'yes' and args.caldb459 == 'no' and args.helder == 'no': skipchandra = 'yes'

#-if rosat and/or xmm only, only plot soft band-
#(will prevent adding extraneous bands to legends)
if skipchandra == 'yes':
    args.broad = 'no'
    args.soft12 = 'no'
    args.soft13 = 'no'
    args.soft23 = 'no'
    if args.xmm == 'no': args.hard = 'no'

#----Set file paths----
helderfile = '/astro/research/kaf33/Dropbox/Research/SN1987A/Helder_table1.txt'
#helderfile = '/Users/kafrank/Dropbox/Research/SN1987A/Helder_table1.txt'

#plotfile = sn1987a_dir+'comparison_dir/spectra_results.pdf'
plotfile = './spectra465_frank2015a_results.pdf'

posteriorfile = sn1987a_dir+'comparison_dir/posteriors.pdf'

mcmcfile = sn1987a_dir+'comparison_dir/mcmc_contamfit.txt'

countratefile = sn1987a_dir+'comparison_dir/countrate_ratios.txt'

xmmfile = sn1987a_dir+'comparison_dir/xmm_results.txt'

rosatfile = sn1987a_dir+'comparison_dir/rosat_results.txt'

earlyfile = sn1987a_dir+'comparison_dir/chandra_early_fluxes.txt'

obsfile = sn1987a_dir+'comparison_dir/chandra_observations.txt'

radiusfile =sn1987a_dir+'comparison_dir/imaging/sn1987a_radii_fits.txt'

radiusfitfile = sn1987a_dir+'comparison_dir/imaging/expansion_fits/sn1987a_mpfit_radii_fits.txt'

expansionfitfile = sn1987a_dir+'comparison_dir/imaging/expansion_fits_sept2014/mpfit_results_3.txt'

#expansionfitfile = sn1987a_dir+'/comparison_dir/imaging/expansion_fits/expansion_mcmc_three_2000000-10000000_results.txt'

officialfile = sn1987a_dir+'comparison_dir/official_obsid_list.txt'

quadcountfile_front = sn1987a_dir+'comparison_dir/spectra_current/sn1987a_'
quadcountfile_back = '_counts.txt'
nwcountsfile = quadcountfile_front+'nw'+quadcountfile_back
necountsfile = quadcountfile_front+'ne'+quadcountfile_back
swcountsfile = quadcountfile_front+'sw'+quadcountfile_back
secountsfile = quadcountfile_front+'se'+quadcountfile_back
totalcountsfile = quadcountfile_front+'total'+quadcountfile_back

#---------------------------------------
#        Read and Convert Data
#---------------------------------------

#----Read CALDB4.6.2 Results from File----

a_data = np.genfromtxt(args.infile_a,skip_header=0,names=True,comments='#',dtype=None)
if args.infile_b != '':
    b_data = np.genfromtxt(args.infile_b,skip_header=0,names=True,comments='#',dtype=None)
else:
    b_data = np.genfromtxt(args.infile_a,skip_header=0,names=True,comments='#',dtype=None)
#columns names:
#obsid	date		age	grating	caldb	model 	pcounts	psoft psoftlow psofthigh	phard phardlow phardhigh	pbroad pbroadlow pbroadhigh	psoft12 psoft12low psoft12high	Counts	soft softlow softhigh	hard hardlow hardhigh	broad broadlow broadhigh	soft12 soft12low soft12high	soft23 soft23low soft23high	soft13 soft13low soft13high	kTcool	kTcoollow kTcoolhigh	kThot	kThotlow kThothigh	Taucool Taucoolerr	tauhot	tauhotlow tauhothigh	normcool normcoollow normcoolhigh	normhot normhotlow normhothigh	chi2	dof	redchi2		csoft12 csoft12low csoft12high	csoft csoftlow csofthigh	csoft23 csoft23low csoft23high	csoft13 csoft13low csoft13high	chard chardlow chardhigh

#----Read Radius Measurements----
if False:
#    r_data = np.genfromtxt(radiusfile,skip_header=0,names=True,comments='#',dtype=None)
#columsn:
#obsid	age	radius	radiuslow	radiushigh	counts	deconcounts

#--calculate errors--
    r_data['radiuslow'] = np.sqrt(r_data['radiuslow']**2.0 + r_data['radius']**2.0/r_data['counts'])
    r_data['radiushigh'] = np.sqrt(r_data['radiushigh']**2.0 + r_data['radius']**2.0/r_data['counts'])

#--convert errorbars to error intervals--
    r_data['radiuslow'] = r_data['radius'] - r_data['radiuslow']
    r_data['radiushigh'] = r_data['radius'] + r_data['radiushigh']

#--convert values to arcsec--
    r_data['radius'] = convert_units(r_data['radius'],'r','arcsec')
    r_data['radiuslow'] = convert_units(r_data['radiuslow'],'r','arcsec')
    r_data['radiushigh'] = convert_units(r_data['radiushigh'],'r','arcsec')

#-read radius fit results-
    rfit_data = np.genfromtxt(expansionfitfile,skip_header=0,names=True,dtype=None)
#columns:
#stat  v_early  v_late  intercept  changepoint
#rows:
#median, mean, standard deviation

#-read fitted radius data points-
    rmodel_data = np.genfromtxt(radiusfitfile,skip_header=0,names=True,comments='#',dtype=None)

#--convert errorbars to error intervals--
    rmodel_data['radiuslow'] = rmodel_data['radius'] - rmodel_data['radiuslow']
    rmodel_data['radiushigh'] = rmodel_data['radius'] + rmodel_data['radiushigh']

#----Get Observation Info----

#--Get Observation Years (and ages if necessary)--
a_obsyears = [str(get_obs_time(int(obs),get='year')) for obs in a_data['obsid']]
a_data = append_fields(a_data,names='year',data=a_obsyears)
if 'age' not in a_data.dtype.names:
    a_age = [get_obs_time(obs) for obs in a_data['obsid']]
    a_data = append_fields(a_data,names='age',data=a_age)

b_obsyears = [str(get_obs_time(int(obs),get='year')) for obs in b_data['obsid']]
b_data = append_fields(b_data,names='year',data=b_obsyears)
if 'age' not in b_data.dtype.names:
    b_age = [get_obs_time(obs) for obs in b_data['obsid']]
    b_data = append_fields(b_data,names='age',data=b_age)

#r_obsyears = [str(get_obs_time(int(obs),get='year')) for obs in r_data['obsid']]
#r_data = append_fields(r_data,names='year',data=r_obsyears)
#if 'age' not in r_data.dtype.names:
#    r_age = [get_obs_time(obs) for obs in r_data['obsid']]
#    r_data = append_fields(r_data,names='age',data=r_age)

#--Read Main Info File--
obs_info = np.genfromtxt(obsfile,skip_header=0,names=True,comments='#',dtype=None)
#column names:
#obsid	date	age	pi	configuration	grating	exposure	frametime	simoffset	yoffset	zoffset

#-get years-
obsyears = [str(get_obs_time(int(obs),get='year')) for obs in obs_info['obsid']]
obs_info = append_fields(obs_info,names='year',data=obsyears)

#--Get FrameTimes--
a_frames = [obs_info['frametime'][np.where(obs_info['obsid'] == obs)] for obs in a_data['obsid']]
a_data = append_fields(a_data,names='frametime',data=np.array(a_frames)[:,0])

b_frames = [obs_info['frametime'][np.where(obs_info['obsid'] == obs)] for obs in b_data['obsid']]
b_data = append_fields(b_data,names='frametime',data=np.array(b_frames)[:,0])

#r_frames = [obs_info['frametime'][np.where(obs_info['obsid'] == obs)] for obs in r_data['obsid']]
#r_data = append_fields(r_data,names='frametime',data=np.array(r_frames)[:,0])

#--Get SimZ--

#-convert 'default' to '0.0'-
obs_info['simoffset'][np.where(obs_info['simoffset'] == 'default')] = '0.0'
#print i_sim,i_obsid

#-associate with a,b,r data-
a_sims = [obs_info['simoffset'][np.where(obs_info['obsid'] == obs)] for obs in a_data['obsid']]
a_data = append_fields(a_data,names='sim',data=np.array(a_sims)[:,0])

b_sims = [obs_info['simoffset'][np.where(obs_info['obsid'] == obs)] for obs in b_data['obsid']]
b_data = append_fields(b_data,names='sim',data=np.array(b_sims)[:,0])

#r_sims = [obs_info['simoffset'][np.where(obs_info['obsid'] == obs)] for obs in r_data['obsid']]
#r_data = append_fields(r_data,names='sim',data=np.array(r_sims)[:,0])

#--get gratings for radius observations--
#r_gratings = [obs_info['grating'][np.where(obs_info['obsid'] == obs)] for obs in r_data['obsid']]
#r_data = append_fields(r_data,names='grating',data=np.array(r_gratings)[:,0])

#----Calculate Contamination Abssorption and Append to Data Structure----
#unabsorbed/absorbed flux 
#(without contamination models / with contamination models)

#--save light curve data to text file for easy reference--

lightcurve = pd.DataFrame(np.ma.filled(a_data,fill_value=-9999))
lightcurve.to_csv(args.infile_a+'_fluxes.txt',sep='\t',columns=['obsid','age','grating','counts','soft','softlow','softhigh','soft12','soft12low','soft12high','hard','hardlow','hardhigh','broad','broadlow','broadhigh'],index=False)

if args.infile_a == sn1987a_dir+'comparison_dir/spectra462/spectra_fits462.txt':
    abs_soft = np.divide(a_data['csoft'],a_data['soft'])
    abs_soft12 = np.divide(a_data['csoft12'],a_data['soft12'])
    abs_soft23 = np.divide(a_data['csoft23'],a_data['soft23'])
    abs_soft13 = np.divide(a_data['csoft13'],a_data['soft13'])
    abs_hard = np.divide(a_data['chard'],a_data['hard'])
    abs_soft_low,abs_soft_high = ratio_err(a_data['csoft'],a_data['soft'],a_data['csoftlow'],a_data['csofthigh'],a_data['softlow'],a_data['softhigh'])
    abs_soft12_low,abs_soft12_high = ratio_err(a_data['csoft12'],a_data['soft12'],a_data['csoft12low'],a_data['csoft12high'],a_data['soft12low'],a_data['soft12high'])
    abs_soft23_low,abs_soft23_high = ratio_err(a_data['csoft23'],a_data['soft23'],a_data['csoft23low'],a_data['csoft23high'],a_data['soft23low'],a_data['soft23high'])
    abs_soft13_low,abs_soft13_high = ratio_err(a_data['csoft13'],a_data['soft13'],a_data['csoft13low'],a_data['csoft13high'],a_data['soft13low'],a_data['soft13high'])
    abs_hard_low,abs_hard_high = ratio_err(a_data['chard'],a_data['hard'],a_data['chardlow'],a_data['chardhigh'],a_data['hardlow'],a_data['hardhigh'])

    a_data = append_fields(a_data,names=['abssoft','abssoft12','abssoft23','abssoft13','abshard'],data=[abs_soft,abs_soft12,abs_soft23,abs_soft13,abs_hard])
    a_data = append_fields(a_data,names=['abssoftlow','abssoft12low','abssoft23low','abssoft13low','abshardlow'],data=[abs_soft_low,abs_soft12_low,abs_soft23_low,abs_soft13_low,abs_hard_low])
    a_data = append_fields(a_data,names=['abssofthigh','abssoft12high','abssoft23high','abssoft13high','abshardhigh'],data=[abs_soft_high,abs_soft12_high,abs_soft23_high,abs_soft13_high,abs_hard_high])

    #----Calculate Optical Depth and Append to Data Structure----
    #errors are error intervals
    depth_soft = -1.0*np.log(abs_soft)
    depth_soft_low = -1.0*np.log(abs_soft_low)
    depth_soft_high = -1.0*np.log(abs_soft_high)
    depth_soft12 = -1.0*np.log(abs_soft12)
    depth_soft12_low = -1.0*np.log(abs_soft12_low)
    depth_soft12_high = -1.0*np.log(abs_soft12_high)
    depth_soft13 = -1.0*np.log(abs_soft13)
    depth_soft13_low = -1.0*np.log(abs_soft13_low)
    depth_soft13_high = -1.0*np.log(abs_soft13_high)
    depth_soft23 = -1.0*np.log(abs_soft23)
    depth_soft23_low = -1.0*np.log(abs_soft23_low)
    depth_soft23_high = -1.0*np.log(abs_soft23_high)
    depth_hard = -1.0*np.log(abs_hard)
    depth_hard_low = -1.0*np.log(abs_hard_low)
    depth_hard_high = -1.0*np.log(abs_hard_high)

    a_data = append_fields(a_data,names=['depthsoft','depthsoftlow','depthsofthigh','depthsoft12','depthsoft12low','depthsoft12high','depthsoft13','depthsoft13low','depthsoft13high','depthsoft23','depthsoft23low','depthsoft23high','depthhard','depthhardlow','depthhardhigh'],data=[depth_soft,depth_soft_low,depth_soft_high,depth_soft12,depth_soft12_low,depth_soft12_high,depth_soft13,depth_soft13_low,depth_soft13_high,depth_soft23,depth_soft23_low,depth_soft23_high,depth_hard,depth_hard_low,depth_hard_high])

#----Calculate Pileup Fractions and Append to Data----
if 'psoft' in a_data.dtype.names:
    ratiosoft = a_data['soft']/a_data['psoft']
    ratiohard = a_data['hard']/a_data['phard']
    ratiobroad = a_data['broad']/a_data['pbroad']
    a_data = append_fields(a_data,names=['softpileup','hardpileup','broadpileup'],data=[ratiosoft,ratiohard,ratiobroad])

    soft_err_low,soft_err_high = ratio_err(a_data['soft'],a_data['psoft'],a_data['softlow'],a_data['softhigh'],a_data['psoftlow'],a_data['psofthigh'])
    hard_err_low,hard_err_high = ratio_err(a_data['hard'],a_data['phard'],a_data['hardlow'],a_data['hardhigh'],a_data['phardlow'],a_data['phardhigh'])
    broad_err_low,broad_err_high = ratio_err(a_data['broad'],a_data['pbroad'],a_data['broadlow'],a_data['broadhigh'],a_data['pbroadlow'],a_data['pbroadhigh'])
    a_data = append_fields(a_data,names=['softpileuplow','softpileuphigh','hardpileuplow','hardpileuphigh','broadpileuplow','broadpileuphigh'],data=[soft_err_low,soft_err_high,hard_err_low,hard_err_high,broad_err_low,broad_err_high])

#----Read Helder2012 Results----
if args.helder == 'yes':
    h_data = np.genfromtxt(helderfile,names=True,dtype=None,comments='#')
#columnn names:
#obsid		configuration	soft	softlow	softhigh	frametimes	softpileup	hard    hardlow    hardhigh   hardpileup

#-convert errorbars to error intervals for consistency with other data-
    h_data['softlow']=h_data['soft']-h_data['softlow']
    h_data['softhigh']=h_data['softhigh']+h_data['soft']
    h_data['hardlow']=h_data['hard']-h_data['hardlow']
    h_data['hardhigh']=h_data['hardhigh']+h_data['hard']

#--get associated ages--
    h_age = [get_obs_time(obs) for obs in h_data['obsid']]
    h_data = append_fields(h_data,names='age',data=h_age)

#--get chip position--
    h_sims = [obs_info['simoffset'][np.where(obs_info['obsid'] == obs)] for obs in h_data['obsid']]
    h_data = append_fields(h_data,names='sim',data=np.array(h_sims)[:,0])

#----Read CountRates and Calculate Ratios----
#(estimate of contamination absorption)
count_data = np.genfromtxt(countratefile,names=True,dtype=None,comments='#')
#column names: (em=emitted flux, meas=measured (contamination absorbed) flux)
#obsid	date	age	grating	caldb	model	emsoft	emsoft12	emsoft13	emhard	meassoft	meassoft12	meassoft13	meashard

count_ratio_soft = count_data['meassoft']/count_data['emsoft']
count_ratio_soft12 = count_data['meassoft12']/count_data['emsoft12']
count_ratio_soft13 = count_data['meassoft13']/count_data['emsoft13']
count_ratio_hard = count_data['meashard']/count_data['emhard']

a_data = append_fields(a_data,names=['ratiosoft','ratiosoft12','ratiosoft13','ratiohard'],data=[count_ratio_soft,count_ratio_soft12,count_ratio_soft13,count_ratio_hard])

#----Read XMM Results----
#notes:
# the hard flux is 3-10 keV
# obsids are read as numbers, not strings, so leading zeros are stripped
xmm_data = np.genfromtxt(xmmfile,names=True,comments='#',dtype=None)

#column names:
#obsid		date		age	filter	totalexp filtexp	soft	softlow		softhigh	hard	        hardlow 	        hardhigh

#fluxes are 10^-13 ergs/s/cm^2
#errors are lower and upper errorbars

#-convert errorbars to error intervals for consistency with other data-
xmm_data['softlow']=xmm_data['soft']-xmm_data['softlow']
xmm_data['softhigh']=xmm_data['softhigh']+xmm_data['soft']
xmm_data['hardlow']=xmm_data['hard']-xmm_data['hardlow']
xmm_data['hardhigh']=xmm_data['hardhigh']+xmm_data['hard']

#----Read ROSAT Results----
#notes:
rosat_data = np.genfromtxt(rosatfile,names=True,comments='#',dtype=None)

#column names:
#age		configuration	soft		softlow		softhigh

#fluxes are 10^-13 ergs/s/cm^2
#errors are lower and upper errorbars

#-convert errorbars to error intervals for consistency with other data-
rosat_data['softlow']=rosat_data['soft']-rosat_data['softlow']
rosat_data['softhigh']=rosat_data['softhigh']+rosat_data['soft']

#fluxes are all in units of 10^-13 ergs/s/cm^2
#errors are 68%
#norms in units of 10^-3
#pileup flux ratios are unpiled/piled

#----Read Early Chandra Fluxes----
#notes:
early_data = np.genfromtxt(earlyfile,names=True,comments='#',dtype=None)

#column names:
#obsid data age grating caldb model counts soft softlow softhigh hard hardlow hardhigh broad broadlow broadhigh soft12 soft12low soft12high

#fluxes are 10^-13 ergs/s/cm^2
#errors are error interval limits

#----Read Quadrant Counts Files----
if args.quadcounts == 'yes':
    nwcounts = pd.read_table('sn1987a_nw_counts.txt',index_col=0,comment='#',header=None,names=['obsid','nwcounts'])
    necounts = pd.read_table('sn1987a_ne_counts.txt',index_col=0,comment='#',header=None,names=['obsid','necounts'])
    swcounts = pd.read_table('sn1987a_sw_counts.txt',index_col=0,comment='#',header=None,names=['obsid','swcounts'])
    secounts = pd.read_table('sn1987a_se_counts.txt',index_col=0,comment='#',header=None,names=['obsid','secounts'])
    
    # combine into one table
    ctable = pd.merge(nwcounts,necounts,left_index=True,right_index=True)
    ctable = pd.merge(ctable,secounts,left_index=True,right_index=True)
    ctable = pd.merge(ctable,swcounts,left_index=True,right_index=True)

    # add total counts columns
    ctable['counts'] = ctable['nwcounts']+ctable['necounts']+ctable['swcounts']+ctable['secounts']

    # remove HETG observations for epochs with both HETG and bare 
    # ACIS (bare ACIS had more counts, and we don't care about pileup here)
    #    ctable.drop([14698,9144])
#    ctable = ctable[ctable.index != 10855]
#    ctable = ctable[ctable.index != 9144]
#    print ctable

    # calculate count fractions
    ctable['nwfraction']=ctable['nwcounts']/ctable['counts']
    ctable['nefraction']=ctable['necounts']/ctable['counts']
    ctable['sefraction']=ctable['secounts']/ctable['counts']
    ctable['swfraction']=ctable['swcounts']/ctable['counts']

    # sum count fractions
    ctable['fractionsum'] = ctable['nwfraction']+ctable['nefraction']+ctable['swfraction']+ctable['sefraction']

    # calculate fraction errors (poisson)
    ctable['nwfraction_errlow'],ctable['nwfraction_errhigh'] = ratio_err(ctable['nwcounts'],ctable['counts'],(ctable['nwcounts']-np.sqrt(ctable['nwcounts'])),ctable['nwcounts']+np.sqrt(ctable['nwcounts']),ctable['counts']-np.sqrt(ctable['counts']),ctable['counts']+np.sqrt(ctable['counts'])) 
    ctable['nefraction_errlow'],ctable['nefraction_errhigh'] = ratio_err(ctable['necounts'],ctable['counts'],ctable['necounts']-np.sqrt(ctable['necounts']),ctable['necounts']+np.sqrt(ctable['necounts']),ctable['counts']-np.sqrt(ctable['counts']),ctable['counts']+np.sqrt(ctable['counts'])) 
    ctable['swfraction_errlow'],ctable['swfraction_errhigh'] = ratio_err(ctable['swcounts'],ctable['counts'],ctable['swcounts']-np.sqrt(ctable['swcounts']),ctable['swcounts']+np.sqrt(ctable['swcounts']),ctable['counts']-np.sqrt(ctable['counts']),ctable['counts']+np.sqrt(ctable['counts'])) 
    ctable['sefraction_errlow'],ctable['sefraction_errhigh'] = ratio_err(ctable['secounts'],ctable['counts'],ctable['secounts']-np.sqrt(ctable['secounts']),ctable['secounts']+np.sqrt(ctable['secounts']),ctable['counts']-np.sqrt(ctable['counts']),ctable['counts']+np.sqrt(ctable['counts'])) 

    # add observation age
    obsyears = [str(get_obs_time(int(obs),get='year')) for obs in ctable.index]
    age = [get_obs_time(obs) for obs in ctable.index]
    ctable['age'] = age

    #--get gratings for radius observations--
    cgratings = [obs_info['grating'][np.where(obs_info['obsid'] == obs)] for 
                 obs in ctable.index]
    ctable['grating'] = cgratings

    #--sort by observation age--
    ctable.sort_values('age',inplace=True)

#---------------------------------------
#---------------------------------------
#                PLOTS
#---------------------------------------
#---------------------------------------

#---------------------------------------
#        Global Plot Properties
#---------------------------------------

#--Set whitespace on bottom, top--
bott = 0.13
topp = 0.87

#--Initialize Plot and Set Axes--
if args.agemin == 0:
    if args.rosat == 'yes':
        args.agemin = 1000
    else:
        args.agemin = 5000
    if args.lateonly == 'yes':
        args.agemin = 7800
    if args.early == 'yes' and args.rosat == 'no':
        args.agemin = 4500
if skipchandra == 'yes' and args.xmm == 'no' and args.rosat == 'yes': 
    if args.agemax == 0: args.agemax = 4500
    fluxmax = 1.5
else:
    if args.agemax == 0: args.agemax = 11000
    if args.fluxmax == 0:
        if min(a_data['age']) <= args.agemax: 
            if args.broad == 'yes':
                #fluxmax = 120
                fluxmax = 1.1*max(a_data['broad'][np.where(a_data['age']<=args.agemax)])       
            else:
    #        fluxmax = 100
                fluxmax = 1.1*max(a_data['soft'][np.where(a_data['age']<=args.agemax)])
        else: #before main chandra, use xmm or rosat
            if min(xmm_data['age']) <= args.agemax: 
                fluxmax = 1.1*max(xmm_data['soft'][np.where(xmm_data['age']<=args.agemax)])
            else:
                fluxmax = 1.1*max(rosat_data['soft'][np.where(rosat_data['age']<=args.agemax)])     
    else: fluxmax = args.fluxmax

if args.ylog == True:
    if args.rosat == 'yes':
        fluxmin = 0.07
    else:
        fluxmin = 0.3
else:
    fluxmin = 0.0

#----Set Axis Labels and Ticks----

#-calculate required top axis labels automatically based on agemin/agemax-
#print 'agemin = ',args.agemin
years = range(convert_time(args.agemin,get='year')+1,convert_time(args.agemax,get='year')+1)
year_ticks = [str(yr) for yr in years]
regular_ages = [convert_time(yr+'-01-01',get='age',informat='date') for yr in year_ticks]
#print 'yeartick = ',year_ticks
#print 'regular ages = ',regular_ages

#radyears = range(1999,2016)
#radius_year_ticks = [str(yr) for yr in radyears]
#radius_regular_ages = [convert_time(yr+'-01-01',get='age') for yr in radius_year_ticks]

spectrayears = ['2001','2002','2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014']
spectra_year_ticks = [str(yr) for yr in spectrayears]
spectra_regular_ages = [convert_time(yr+'-01-01',get='age',informat='date') for yr in spectra_year_ticks]

if args.clean == 'yes':
    lightcurve_topx = regular_ages
    lightcurve_toplab = year_ticks
    topx = spectra_regular_ages
    toplab = spectra_year_ticks
#    ext_topx = radius_regular_ages
#    ext_toplab = radius_year_ticks
    botx = None
    toptitle = 'Year'
    mgrid = False
#    frame = False
    frame = True #but sets edge and facecolor to white
else:
    topx = obs_info['age']#a_data['age']
    toplab = obs_info['year']#a_data['year']
    lightcurve_topx = topx
    lightcurve_toplab = toplab
    ext_topx = topx
    ext_toplab = toplab
    toptitle = 'Observation Year'
    botx = obs_info['age']#a_data['age']
    mgrid = True
    frame = True

#-set fonts-
font = {'family':'serif','weight':'normal','size':15}
matplotlib.rc('font',**font)
matplotlib.rcParams['text.latex.preamble'] = [r"\boldmath"]
plt.rc('text',usetex=True)

#----Set Symbols and Colors----

#--find indices--

#-detector position-
a_offz_i = np.where(a_data['sim'] == '-8.42')
b_offz_i = np.where(b_data['sim'] == '-8.42')
#r_offz_i = np.where(r_data['sim'] == '-8.42')
if args.helder == 'yes': h_offz_i = np.where(h_data['sim'] == '-8.42')

#-ACIS grating-
a_hetg_i = np.where(a_data['grating'] == 'HETG')
a_bare_i = np.where(a_data['grating'] == 'NONE')
b_hetg_i = np.where(b_data['grating'] == 'HETG')
b_bare_i = np.where(b_data['grating'] == 'NONE')
#r_hetg_i = np.where(r_data['grating'] == 'HETG')
#r_bare_i = np.where(r_data['grating'] == 'NONE')
if args.helder == 'yes':
    h_hetg_i = np.where(h_data['configuration'] == 'HETG')
    h_bare_i = np.where(h_data['configuration'] == 'S3')
a_letg_i = np.where(a_data['grating'] == 'LETG')
b_letg_i = np.where(b_data['grating'] == 'LETG')

#--make lists for each--

#-colors = energy band-
softcolor = 'green'#'red'
soft12color = 'purple'#'green'
soft13color = 'blue'
soft23color = 'magenta'
hardcolor = 'cyan'
broadcolor = 'black'#'purple'

#--colors = instrument,symbols = band--
#(used only for chandra and xmm, soft and hard bands, if instcolor='yes')
#(supersedes the bandcolor and instrument symbols)
xmmsoftcolor = 'red'
chandrasoftcolor = 'royalblue'
xmmsoftsym = '^'
chandrasoftsym = 'o'
xmmhardcolor = 'red'
chandrahardcolor = 'royalblue'
xmmhardsym = 'd'
chandrahardsym = 's'

#-number of observations-
a_nobs = a_data.shape[0]
b_nobs = b_data.shape[0]
#r_nobs = r_data.shape[0]
xmm_nobs = xmm_data.shape[0]
rosat_nobs = rosat_data.shape[0]
if args.helder == 'yes': h_nobs = h_data.shape[0]
print 'nobs = ',a_nobs
#print 'xmm_nobs = ',xmm_nobs

#-symbol = instrument (ACIS-HETG,ACIS-NONE,XMM-EPIC-pn)-
hetg_sym = 'o' #HETG = circle
bare_sym = 'd' #bare acis = thin diamond
letg_sym = 'p' #LETG = pentagon
a_syms = [hetg_sym]*a_nobs
for i in list(a_bare_i[:][0]):
    a_syms[i] = bare_sym 
for i in list(a_letg_i[:][0]):
    a_syms[i] = letg_sym
b_syms = [hetg_sym]*b_nobs 
for i in list(b_bare_i[:][0]):
    b_syms[i] = bare_sym 
for i in list(b_letg_i[:][0]):
    b_syms[i] = letg_sym 
#r_syms = [hetg_sym]*r_nobs 
#for i in list(r_bare_i[:][0]):
#    r_syms[i] = bare_sym 
if args.helder == 'yes':
    h_syms = [hetg_sym]*h_nobs 
    for i in list(h_bare_i[:][0]):
        h_syms[i] = bare_sym 
xmm_syms = '^' #xmm = triangle
rosat_syms = 'x' #rosat = x
early_syms = [hetg_sym,bare_sym]

chandrahardsyms = ['s']*a_nobs
chandrasoftsyms = ['o']*a_nobs
#for i in list(a_bare_i[:][0]):
#    chandrahardsyms[i] = 'H'
#    chandrasoftsyms[i] = 'p'

chandrahardcolors = ['royalblue']*a_nobs
chandrasoftcolors = ['royalblue']*a_nobs
for i in list(a_bare_i[:][0]):
    chandrahardcolors[i] = 'navy'
    chandrasoftcolors[i] = 'navy'

#-markeredgecolor = simz (detector z)-
if args.instcolor == 'no':
    offcolor = 'gray'
else: 
    offcolor = 'black'
a_mecs = ['black']*a_nobs
b_mecs = ['black']*b_nobs
#r_mecs = ['black']*r_nobs
if args.helder == 'yes': h_mecs = ['black']*h_nobs
if args.ignorechippos == 'no':
    for i in list(a_offz_i[:][0]):
        a_mecs[i] = offcolor #offset observations have offcolor outline
    for i in list(b_offz_i[:][0]):
        b_mecs[i] = offcolor #offset observations have offcolor outline
#    for i in list(r_offz_i[:][0]):
#        r_mecs[i] = offcolor #offset observations have offcolor outline
    if args.helder == 'yes':
        for i in list(h_offz_i[:][0]):
            h_mecs[i] = offcolor #offset observations have offcolor outline
xmm_mecs = 'black'
rosat_mecs = softcolor
early_mecs = 'black'

#-markeredgewidth = simz (detector z)-
offmew = 1.5
a_mews = [0.5]*a_nobs
b_mews = [0.5]*b_nobs
#r_mews = [0.5]*r_nobs
if args.helder == 'yes': h_mews = [0.5]*h_nobs
if args.ignorechippos == 'no':
    for i in list(a_offz_i[:][0]):
        a_mews[i] = offmew #offset observations have different outline
    for i in list(b_offz_i[:][0]):
        b_mews[i] = offmew #offset observations have different outline
#    for i in list(r_offz_i[:][0]):
#        r_mews[i] = offmew #offset observations have different outline
    if args.helder == 'yes': 
        for i in list(h_offz_i[:][0]):
            h_mews[i] = offmew #offset observations have different outline
#h_mews = 0.5
xmm_mews = 1.0
rosat_mews = 0.7
early_mews = 0.5

#-markersize - # = frametime-
xmm_sizes = 8
#a_sizes = [f*3 for f in a_data['frametime']]
#b_sizes = [f*3 for f in b_data['frametime']]
a_sizes = 8
b_sizes = 4
#r_sizes = 4
rosat_sizes = 8
early_sizes = 8

#-alpha (transparency) = caldb-
a_alpha = 1.0
if args.skip462 == 'no':
    b_alpha = 0.6
else:
    b_alpha = 1.0
#r_alpha = 1.0
xmm_alpha = [1.0]*xmm_nobs
xmm_alpha[-1] = 0.7 #set last point to different alpha (since reduced myself)
if args.skip462 == 'no':
    h_alpha = 0.6
else:
    if args.caldb459 == 'yes':
        h_alpha = 0.6
    else:
        h_alpha = 1.0
rosat_alpha = 1.0
early_alpha = 1.0

piled_alpha = 0.1

#---------------------------------------
#        Light Curve(s)
#---------------------------------------

#----Set Up Plot----

pdffile = PdfPages(plotfile)
#fig = plt.figure()

#-initialize main plot-
#fig,ax1=plt.subplots()
#ax1.set_xlim(xmin,xmax)
#ax1.set_ylim(ymin,ymax)

#-set axis titles-
#plt.xlabel('SN1987A Age [days]')
#plt.ylabel('Flux [10$^{-13}$ erg cm$^{-2}$ s$^{-1}$]')

#-set bottom axis-
#set_axis(ax1,x=botx,grid=mgrid)

#-add extra space on bottom-
#plt.subplots_adjust(bottom=bott,top=topp)

#--Initialize Plot--
ax1 = start_plot(xtitle='SN1987A Age [days]',ytitle='Flux [10$^{-13}$ erg cm$^{-2}$ s$^{-1}$]',xmin=args.agemin,xmax=args.agemax,ymin=fluxmin,ymax=fluxmax,ylog=args.ylog,clean=args.clean)


if args.skip462 == 'no':
    if args.soft == 'yes':
#--Soft Flux--
        if args.instcolor == 'yes': #color=instrument,symbol=band (effects xmm and chandra, soft and hard only)
            plot_a(ax1,a_data['age'],a_data['soft'],low=a_data['softlow'],high=a_data['softhigh'],syms=chandrasoftsyms,colors=chandrasoftcolors,alphas=a_alpha,mecs=a_mecs,mews=a_mews,sizes=a_sizes)
        else:
            plot_a(ax1,a_data['age'],a_data['soft'],low=a_data['softlow'],high=a_data['softhigh'],syms=a_syms,colors=softcolor,alphas=a_alpha,mecs=a_mecs,mews=a_mews,sizes=a_sizes)
#    if args.dofit == 'yes': #do linear fit to light curve
#        lc_fit = np.polyfit(a_data['age'],a_data['soft'],1)
#        fit_fn = np.poly1d(lc_fit) #make function out of best fit
#        plt.plot(a_data['age'],fit_fn(a_data['age']),'-')#plot best fit
        if args.dofit == 'yes': #do powerlaw fit to light curve
            fit_i = np.array([-11,-10,-9,-8,-7,-6,-4,-3,-2,-1])
#            fit_i = range(3,a_nobs-1)
            fit_params, fit_covars = curve_fit(af.powerlaw,a_data['age'][fit_i],a_data['soft'][fit_i],p0=[1.0,0.5,-8000.0])
            sigmas = [np.sqrt(fit_covars[0,0]),np.sqrt(fit_covars[1,1]),np.sqrt(fit_covars[2,2])]
            print 'fit_params = ',fit_params
            print 'fit_covars = ',fit_covars
            fitx=np.linspace(min(a_data['age'][fit_i]),11000.0)
            plt.plot(fitx,af.powerlaw(fitx,fit_params[0],fit_params[1],fit_params[2]),'--',color=softcolor)
            plt.plot(fitx,af.powerlaw(fitx,fit_params[0]+sigmas[0],fit_params[1]+sigmas[1],fit_params[2]+sigmas[1]),'--',color=softcolor,alpha=0.5)            
            plt.plot(fitx,af.powerlaw(fitx,fit_params[0]-sigmas[0],fit_params[1]-sigmas[1],fit_params[2]-sigmas[1]),'--',color=softcolor,alpha=0.5)      

            future_ages = np.array([10257,10441,10623,10807])
            print 'fit(future_ages) = ',af.powerlaw(future_ages,fit_params[0],fit_params[1],fit_params[2])
            plt.plot(future_ages,af.powerlaw(future_ages,fit_params[0],fit_params[1],fit_params[2]),'*',color=softcolor)            
#            plt.plot(future_ages,af.powerlaw(future_ages,fit_params[0]+sigmas[0],fit_params[1]+sigmas[1],fit_params[2]+sigmas[1]),'*',color=softcolor,alpha=0.6)            
#            plt.plot(future_ages,af.powerlaw(future_ages,fit_params[0]-sigmas[0],fit_params[1]-sigmas[1],fit_params[2]-sigmas[1]),'*',color=softcolor,alpha=0.6)            
            

    if args.soft12 == 'yes':
#--Soft12 Flux (1.0-2.0keV)--
        plot_a(ax1,a_data['age'],a_data['soft12'],low=a_data['soft12low'],high=a_data['soft12high'],syms=a_syms,colors=soft12color,alphas=a_alpha,mecs=a_mecs,mews=a_mews,sizes=a_sizes)

    if args.soft13 == 'yes':
#--Soft13 Flux (1.0-3.0keV)--
        plot_a(ax1,a_data['age'],a_data['soft13'],low=a_data['soft13low'],high=a_data['soft13high'],syms=a_syms,colors=soft13color,alphas=a_alpha,mecs=a_mecs,mews=a_mews,sizes=a_sizes)

    if args.hard == 'yes':
#--Hard Flux--
        if args.instcolor == 'yes':
            plot_a(ax1,a_data['age'],a_data['hard'],low=a_data['hardlow'],high=a_data['hardhigh'],syms=chandrahardsyms,colors=chandrahardcolors,alphas=a_alpha,mecs=a_mecs,mews=a_mews,sizes=a_sizes)
        else:
            plot_a(ax1,a_data['age'],a_data['hard'],low=a_data['hardlow'],high=a_data['hardhigh'],syms=a_syms,colors=hardcolor,alphas=a_alpha,mecs=a_mecs,mews=a_mews,sizes=a_sizes)
        #print 'hard fluxes = ',a_data['hard']
    if args.broad == 'yes':
#--Broad Flux--
        plot_a(ax1,a_data['age'],a_data['broad'],low=a_data['broadlow'],high=a_data['broadhigh'],syms=a_syms,colors=broadcolor,alphas=a_alpha,mecs=a_mecs,mews=a_mews,sizes=a_sizes)


#--CALDB4.5.9--
if args.caldb459 == 'yes':
    #--Soft Flux--
    if args.soft == 'yes': plot_a(ax1,b_data['age'],b_data['soft'],low=b_data['softlow'],high=b_data['softhigh'],syms=b_syms,colors=softcolor,alphas=b_alpha,mecs=b_mecs,mews=b_mews,sizes=b_sizes)

    #--Soft12 Flux (1.0-2.0keV)--
    if args.soft12 == 'yes': plot_a(ax1,b_data['age'],b_data['soft12'],low=b_data['soft12low'],high=b_data['soft12high'],syms=b_syms,colors=soft12color,alphas=b_alpha,mecs=b_mecs,mews=b_mews,sizes=b_sizes)

    #--Soft13 Flux (1.0-3.0keV)--
    if args.soft13 == 'yes': plot_a(ax1,b_data['age'],b_data['soft13'],low=b_data['soft13low'],high=b_data['soft13high'],syms=b_syms,colors=soft13color,alphas=b_alpha,mecs=b_mecs,mews=b_mews,sizes=b_sizes)

    #--Hard Flux--
    if args.hard == 'yes': plot_a(ax1,b_data['age'],b_data['hard'],low=b_data['hardlow'],high=b_data['hardhigh'],syms=b_syms,colors=hardcolor,alphas=b_alpha,mecs=b_mecs,mews=b_mews,sizes=b_sizes)

    #--Broad Flux--
    if args.broad == 'yes': plot_a(ax1,b_data['age'],b_data['broad'],low=b_data['broadlow'],high=b_data['broadhigh'],syms=b_syms,colors=broadcolor,alphas=b_alpha,mecs=b_mecs,mews=b_mews,sizes=b_sizes)

#--Early Observations (obs 1387 and 122)--
if args.early == 'yes':

    #--Soft Flux--
    if args.soft == 'yes': plot_a(ax1,early_data['age'],early_data['soft'],low=early_data['softlow'],high=early_data['softhigh'],syms=early_syms,colors=softcolor,alphas=early_alpha,mecs=early_mecs,mews=early_mews,sizes=early_sizes)

    #--Soft12 Flux (1.0-2.0keV)--
    if args.soft12 == 'yes': plot_a(ax1,early_data['age'],early_data['soft12'],low=early_data['soft12low'],high=early_data['soft12high'],syms=early_syms,colors=soft12color,alphas=early_alpha,mecs=early_mecs,mews=early_mews,sizes=early_sizes)

    #--Hard Flux--
    if args.hard == 'yes': plot_a(ax1,early_data['age'],early_data['hard'],low=early_data['hardlow'],high=early_data['hardhigh'],syms=early_syms,colors=hardcolor,alphas=early_alpha,mecs=early_mecs,mews=early_mews,sizes=early_sizes)

    #--Broad Flux--
    if args.broad == 'yes': plot_a(ax1,early_data['age'],early_data['broad'],low=early_data['broadlow'],high=early_data['broadhigh'],syms=early_syms,colors=broadcolor,alphas=early_alpha,mecs=early_mecs,mews=early_mews,sizes=early_sizes)


#--Piled Data--
if args.piled == 'yes':
    #-Soft Flux-
    if args.soft == 'yes': plot_a(ax1,a_data['age'],a_data['psoft'],low=a_data['psoftlow'],high=a_data['psofthigh'],syms=a_syms,colors=softcolor,alphas=piled_alpha,mecs=a_mecs,mews=a_mews,sizes=a_sizes)

    #-Soft12 Flux (1.0-2.0keV)-
    if args.soft12 == 'yes': plot_a(ax1,a_data['age'],a_data['psoft12'],low=a_data['psoft12low'],high=a_data['psoft12high'],syms=a_syms,colors=soft12color,alphas=piled_alpha,mecs=a_mecs,mews=a_mews,sizes=a_sizes)

    #-Soft13 Flux (1.0-3.0keV)-
#    if args.soft13 == 'yes': plot_a(ax1,a_data['age'],a_data['psoft13'],low=a_data['psoft13low'],high=a_data['psoft13high'],syms=a_syms,colors=soft13color,alphas=piled_alpha,mecs=a_mecs,mews=a_mews,sizes=a_sizes)

    #-Hard Flux-
    if args.hard == 'yes': plot_a(ax1,a_data['age'],a_data['phard'],low=a_data['phardlow'],high=a_data['phardhigh'],syms=a_syms,colors=hardcolor,alphas=piled_alpha,mecs=a_mecs,mews=a_mews,sizes=a_sizes)

    #-Broad Flux-
    if args.broad == 'yes': plot_a(ax1,a_data['age'],a_data['pbroad'],low=a_data['pbroadlow'],high=a_data['pbroadhigh'],syms=a_syms,colors=broadcolor,alphas=piled_alpha,mecs=a_mecs,mews=a_mews,sizes=a_sizes)

#--Plot XMM Fluxes--
if args.xmm == 'yes':
    if args.soft == 'yes': 
        if args.instcolor == 'yes':
            plot_a(ax1,xmm_data['age'],xmm_data['soft'],low=xmm_data['softlow'],high=xmm_data['softhigh'],syms=xmmsoftsym,colors=xmmsoftcolor,alphas=xmm_alpha,mecs=xmm_mecs,mews=xmm_mews,sizes=xmm_sizes)
        else:
            plot_a(ax1,xmm_data['age'],xmm_data['soft'],low=xmm_data['softlow'],high=xmm_data['softhigh'],syms=xmm_syms,colors=softcolor,alphas=xmm_alpha,mecs=xmm_mecs,mews=xmm_mews,sizes=xmm_sizes)
    if args.dofit == 'yes': #do linear fit to light curve
        lc_fit = np.polyfit(xmm_data['age'][-4:],xmm_data['soft'][-4:],1)
        fit_fn = np.poly1d(lc_fit) #make function out of best fit
        plt.plot(xmm_data['age'][-4:],fit_fn(xmm_data['age'][-4:]),'--',color='blue')#plot best fit
        print 'fit = ',lc_fit

    if args.hard == 'yes': 
        if args.instcolor == 'yes':
            plot_a(ax1,xmm_data['age'],xmm_data['hard'],low=xmm_data['hardlow'],high=xmm_data['hardhigh'],syms=xmmhardsym,colors=xmmhardcolor,alphas=xmm_alpha,mecs=xmm_mecs,mews=xmm_mews,sizes=xmm_sizes)
        else: 
            plot_a(ax1,xmm_data['age'],xmm_data['hard'],low=xmm_data['hardlow'],high=xmm_data['hardhigh'],syms=xmm_syms,colors=hardcolor,alphas=xmm_alpha,mecs=xmm_mecs,mews=xmm_mews,sizes=xmm_sizes)

#--Plot ROSAT Fluxes--
if args.rosat == 'yes':
    if args.soft == 'yes': plot_a(ax1,rosat_data['age'],rosat_data['soft'],low=rosat_data['softlow'],high=rosat_data['softhigh'],syms=rosat_syms,colors=softcolor,alphas=rosat_alpha,mecs=rosat_mecs,mews=rosat_mews,sizes=rosat_sizes)

#--Plot Helder2012 Measurements--

if args.helder == 'yes':
    if args.soft == 'yes': plot_a(ax1,h_data['age'],h_data['soft'],low=h_data['softlow'],high=h_data['softhigh'],syms=h_syms,colors=softcolor,alphas=h_alpha,mecs=h_mecs,mews=h_mews)
    if args.hard == 'yes': plot_a(ax1,h_data['age'],h_data['hard'],low=h_data['hardlow'],high=h_data['hardhigh'],syms=h_syms,colors=hardcolor,alphas=h_alpha,mecs=h_mecs,mews=h_mews)
        
#--Plot Extras--
if args.plotextras != None:
    plot_extras(names=args.plotextras)

        
#--Add Legend--
emptyx = [0]
emptyy = [0]
#plt.scatter(emptyx,emptyy,color='black',marker='^',label='CALDB 4.5.9')
if args.instcolor == 'yes':
    plt.scatter(emptyx,emptyy,color=xmmsoftcolor,marker=xmmsoftsym,label='XMM     0.5 - 2.0 keV',linewidths=xmm_mews)
    plt.scatter(emptyx,emptyy,color=chandrasoftcolor,marker=chandrasoftsym,label='Chandra 0.5 - 2.0 keV',linewidth=1.0)
    plt.scatter(emptyx,emptyy,color=xmmhardcolor,marker=xmmhardsym,label='XMM     3.0 - 10.0 keV',linewidths=xmm_mews)
    plt.scatter(emptyx,emptyy,color=chandrahardcolor,marker=chandrahardsym,label='Chandra 3.0 - 8.0 keV',linewidths=1.0)
else: 
    if skipchandra == 'no':
        if args.agemax >= 4500: #only include if chandra observations plotted
            plt.scatter(emptyx,emptyy,facecolor='white',marker=hetg_sym,label='ACIS (w/HETG)',alpha=a_alpha,linewidth=0.5)
            plt.scatter(emptyx,emptyy,facecolor='white',marker=bare_sym,label='ACIS (no grating)',alpha=a_alpha,linewidth=0.5)
            plt.scatter(emptyx,emptyy,facecolor='white',marker=letg_sym,label='LETG',alpha=a_alpha,linewidth=0.5)
    #    plt.scatter(emptyx,emptyy,color='black',marker=bare_sym,label='ACIS',alpha=a_alpha)
    if args.xmm == 'yes':
        plt.scatter(emptyx,emptyy,facecolor='white',marker=xmm_syms,label='EPIC-pn',alpha=xmm_alpha[0],linewidth=0.5)
    if args.rosat == 'yes':
        plt.scatter(emptyx,emptyy,color='black',marker=rosat_syms,label='ROSAT',alpha=rosat_alpha)
    if args.caldb459 == 'yes' and args.skip462 == 'no':
        plt.scatter(emptyx,emptyy,color='black',marker='s',label='CALDB 4.6.2',alpha=a_alpha)
        if args.helder == 'yes': plt.scatter(emptyx,emptyy,color='gray',marker='s',label='CALDB 4.4.10/4.5.9',alpha=b_alpha)
        else: plt.scatter(emptyx,emptyy,color='gray',marker='s',label='CALDB 4.5.9',alpha=b_alpha)
    if args.piled == 'yes':
        plt.scatter(emptyx,emptyy,color='black',marker='s',label='non-pileup corrected',alpha=piled_alpha)
    if args.broad == 'yes': plt.scatter(emptyx,emptyy,color=broadcolor,marker='s',label='0.5-8.0 keV')
    if args.soft == 'yes': plt.scatter(emptyx,emptyy,color=softcolor,marker='s',label='0.5-2.0 keV')
    if args.soft12 == 'yes': plt.scatter(emptyx,emptyy,color=soft12color,marker='s',label='1.0-2.0 keV')
    if args.soft13 == 'yes': plt.scatter(emptyx,emptyy,color=soft13color,marker='s',label='1.0-3.0 keV')
    if args.hard == 'yes': plt.scatter(emptyx,emptyy,color=hardcolor,marker='s',label='3.0-8.0 keV')
    #if skipchandra == 'no': plt.scatter(emptyx,emptyy,facecolor='none',edgecolor=offcolor,linewidth=offmew,marker='s',label='Low Chip Y',alpha=1)

leg = plt.legend(loc='upper left',fontsize='small',frameon=frame,scatterpoints=1,numpoints=1,markerscale=1.5)
lframe = leg.get_frame()
lframe.set_facecolor('white')
lframe.set_edgecolor('white')

#lab1,lab2,lab3,lab4,lab5,lab6 = leg.get_texts()
#lab4.set_color(softcolor)
#lab5.set_color(soft12color)
#lab6.set_color(hardcolor)

if args.ylog == True: ax1.set_yscale('log')

#--Add Top (year) x-axis--
#note this must come after plotting things on the first axis!

set_axis(ax1,x=lightcurve_topx,twin=True,title=toptitle,xlab=lightcurve_toplab,clean=args.clean)#,xlab=b_obsyears)

#--Close Plot--
pdffile.savefig()
plt.close()

#---------------------------------------
#          Pileup Fractions
#---------------------------------------

#----Initialize Plot----
ax1 = start_plot(xtitle='SN1987A Age [days]',ytitle='Unpiled Flux/Piled Flux',xmin=args.agemin,xmax=args.agemax,ymin=0.5,ymax=1.8,clean=args.clean)

#fig,ax1=plt.subplots()

#plt.xlabel('SN1987A Age [days]')
#plt.ylabel('Unpiled/Piled Flux')

#set_axis(ax1,x=botx,grid=mgrid)

#plt.subplots_adjust(bottom=bott,top=topp)

#-set limits and axis titles-
#xmin=5000
#xmax=11000
#plt.xlim(xmin,xmax)
#ax1.set_ylim(0.5,1.8)

#--plot new fractions--
if 'psoft' in a_data.dtype.names:
    if args.soft == 'yes': plot_a(ax1,a_data['age'],a_data['softpileup'],low=a_data['softpileuplow'],high=a_data['softpileuphigh'],syms=a_syms,colors=softcolor,alphas=a_alpha,mecs=a_mecs,line='-',mews=a_mews,errorband=True)
    if args.hard == 'yes': plot_a(ax1,a_data['age'],a_data['hardpileup'],low=a_data['hardpileuplow'],high=a_data['hardpileuphigh'],syms=a_syms,colors=hardcolor,alphas=a_alpha,mecs=a_mecs,line='-',mews=a_mews,errorband=True)
    if args.broad == 'yes': plot_a(ax1,a_data['age'],a_data['broadpileup'],low=a_data['broadpileuplow'],high=a_data['broadpileuphigh'],syms=a_syms,colors=broadcolor,alphas=a_alpha,mecs=a_mecs,line='-',mews=a_mews,errorband=True)

#--plot Helder2012 fractions--
if args.helder == 'yes':
    if args.soft == 'yes': plot_a(ax1,h_data['age'],h_data['softpileup'],syms=h_syms,colors=softcolor,alphas=h_alpha,mecs=h_mecs,line='-')
    if args.hard == 'yes': plot_a(ax1,h_data['age'],h_data['hardpileup'],syms=h_syms,colors=hardcolor,alphas=h_alpha,mecs=h_mecs,line='-')

    plt.plot([0],[0.4],markersize=0,label='Helder2012',alpha=h_alpha,marker=h_syms[0])

#--Plot Extras--
if args.plotextras != None:
    plot_extras(names=args.plotextras)
    
#--plot legend--
if skipchandra == 'no':
    plt.plot([0],[0.4],hetg_sym,label='with HETG',alpha=a_alpha,color='white')
    plt.plot([0],[0.4],bare_sym,label='bare ACIS',alpha=a_alpha,color='white')
if args.soft == 'yes': plt.plot([0],[0.4],'s',label='0.5-2.0 keV',color=softcolor)
if args.broad == 'yes':plt.plot([0],[0.4],'s',label='0.5-8.0 keV',color=broadcolor)
if args.hard == 'yes':plt.plot([0],[0.4],'s',label='3.0-8.0 keV',color=hardcolor)
plt.scatter(emptyx,emptyy,facecolor='none',edgecolor=offcolor,linewidth=offmew,marker='s',label='Low Chip Y',alpha=1)

leg = plt.legend(numpoints=1,frameon=frame,fontsize='xx-small',scatterpoints=1)
lframe = leg.get_frame()
lframe.set_facecolor('white')
lframe.set_edgecolor('white')

#--Plot Extras--
if args.plotextras != None:
    plot_extras(names=args.plotextras)

#set_axis(ax1,b_age,twin=True,title='Observation Year',xlab=b_obsyears)
set_axis(ax1,x=topx,twin=True,title=toptitle,xlab=toplab,clean=args.clean)

#--close plot--
pdffile.savefig()
plt.close()

#---------------------------------------
#          Hardness Ratios
#---------------------------------------
#(hard/soft flux)

ax1 = start_plot(xtitle='SN1987A Age [days]',ytitle='Hard Flux/Soft Flux',xmin=args.agemin,xmax=args.agemax,ymin=0.0,clean=args.clean)

#--plot soft ratios--
rlow,rhigh = ratio_err(a_data['hard'],a_data['soft'],a_data['hardlow'],a_data['hardhigh'],a_data['softlow'],a_data['softhigh'])
plot_a(ax1,a_data['age'],np.divide(a_data['hard'],a_data['soft']),low=rlow,high=rhigh,syms=a_syms,alphas=a_alpha,mecs=a_mecs,mews=a_mews,colors=softcolor,sizes=a_sizes)

#--plot soft12 ratios--
rlow,rhigh = ratio_err(a_data['hard'],a_data['soft12'],a_data['hardlow'],a_data['hardhigh'],a_data['soft12low'],a_data['soft12high'])
plot_a(ax1,a_data['age'],np.divide(a_data['hard'],a_data['soft12']),low=rlow,high=rhigh,syms=a_syms,alphas=a_alpha,mecs=a_mecs,mews=a_mews,colors=soft12color,sizes=a_sizes)

#--plot Helder2012 fractions--
#plt.plot(hage,hflux_hard/hflux_soft,marker='o',color='green',label='Soft=0.5-2.0 keV, Helder2012')

#--plot legend--
plt.plot([0],[0.4],hetg_sym,label='with HETG',alpha=a_alpha,color='white')
plt.plot([0],[0.4],bare_sym,label='bare ACIS',alpha=a_alpha,color='white')
plt.plot([0],[0.4],'s',label='F$_{3.0-8.0keV}$ / F$_{0.5-2.0keV}$',color=softcolor)
plt.plot([0],[0.4],'s',label='F$_{3.0-8.0keV}$ / F$_{1.0-2.0keV}$',color=soft12color)
plt.scatter(emptyx,emptyy,facecolor='none',edgecolor=offcolor,linewidth=offmew,marker='s',label='Low Chip Y',alpha=1)

#--Plot Extras--
if args.plotextras != None:
    plot_extras(names=args.plotextras)

#--plot legend--
#plt.scatter([0],[0.05],color='black',marker=None,label='Hard=3.0-8.0 keV')
#plt.legend(numpoints=1,frameon=frame,fontsize='xx-small')

leg = plt.legend(numpoints=1,frameon=frame,fontsize='xx-small',scatterpoints=1)
lframe = leg.get_frame()
lframe.set_facecolor('white')
lframe.set_edgecolor('white')

set_axis(ax1,x=topx,twin=True,title=toptitle,xlab=toplab,clean=args.clean)

#--close plot--
pdffile.savefig()
plt.close()

#---------------------------------------
#       Contamination Absorption
#---------------------------------------

if args.infile_a == sn1987a_dir+'comparison_dir/spectra462/spectra_fits462.txt':

#----Set Up Plot----
    ax1 = start_plot(xtitle='SN1987A Age [days]',ytitle='Absorbed Flux / Unabsorbed Flux',xmin=args.agemin,xmax=args.agemax,clean=args.clean)

#--plot with errorbands--
    emptyy = [min(a_data['abssoft'])]
    if args.soft == 'yes': 
        plot_a(ax1,a_data['age'],a_data['abssoft'],low=a_data['abssoftlow'],high=a_data['abssofthigh'],syms=a_syms,colors=softcolor,alphas=a_alpha,mecs=a_mecs,errorband=True,mews=a_mews,line='-',sizes=a_sizes)
        plt.scatter(emptyx,emptyy,color=softcolor,marker='s',label='0.5-2.0 keV')
    if args.soft23 == 'yes': 
        plot_a(ax1,a_data['age'],a_data['abssoft23'],low=a_data
    ['abssoft23low'],high=a_data['abssoft23high'],syms=a_syms,colors=soft23color,alphas=a_alpha,mecs=a_mecs,errorband=True,mews=a_mews,line='-',sizes=a_sizes)
        plt.scatter(emptyx,emptyy,color=soft23color,marker='s',label='2.0-3.0 keV')
    if args.soft13 == 'yes': 
        plot_a(ax1,a_data['age'],a_data['abssoft13'],low=a_data['abssoft13low'],high=a_data['abssoft13high'],syms=a_syms,colors=soft13color,alphas=a_alpha,mecs=a_mecs,errorband=True,mews=a_mews,line='-',sizes=a_sizes)
        plt.scatter(emptyx,emptyy,color=soft13color,marker='s',label='1.0-3.0 keV')
    if args.soft12 == 'yes': 
        plot_a(ax1,a_data['age'],a_data['abssoft12'],low=a_data['abssoft12low'],high=a_data['abssoft12high'],syms=a_syms,colors=soft12color,alphas=a_alpha,mecs=a_mecs,errorband=True,mews=a_mews,line='-',sizes=a_sizes)
        plt.scatter(emptyx,emptyy,color=soft12color,marker='s',label='1.0-2.0 keV')
    if args.hard == 'yes': 
        plot_a(ax1,a_data['age'],a_data['abshard'],low=a_data['abshardlow'],high=a_data['abshardhigh'],syms=a_syms,colors=hardcolor,alphas=a_alpha,mecs=a_mecs,errorband=True,mews=a_mews,line='-',sizes=a_sizes)
        plt.scatter(emptyx,emptyy,color=hardcolor,marker='s',label='3.0-8.0 keV')

    plt.scatter(emptyx,emptyy,facecolor='none',edgecolor=offcolor,linewidth=offmew,marker='s',label='Low Chip Y',alpha=1)

    #--plot ratios from count rates--
    #y_min,y_max = plt.ylim()
    #plt.plot(count_age,count_ratio_soft,'^',color='blue')
    #plt.plot(count_age,count_ratio_soft12,'^',color='green')
    #plt.plot(count_age,count_ratio_soft13,'^',color='cyan')
    #plt.plot(count_age,count_ratio_hard,'^',color='red')
    #plt.scatter(emptyx,[y_min],color='black',marker='^',label='Count Rate Ratios')

    #--plot legend--
    leg = plt.legend(numpoints=1,frameon=frame,fontsize='xx-small',loc='upper left',scatterpoints=1,markerscale=1.0)
    lframe = leg.get_frame()
    lframe.set_facecolor('white')
    lframe.set_edgecolor('white')


    #--plot top axis--
    set_axis(ax1,x=topx,twin=True,title=toptitle,xlab=toplab,clean=args.clean)

    #--close plot--
    pdffile.savefig()
    plt.close()

    #---------------------------------------
    #       Optical Depth
    #---------------------------------------

    #----Set Up Plot----
    ax1 = start_plot(xtitle='SN1987A Age [days]',ytitle='Optical Depth',xmin=args.agemin,xmax=args.agemax,clean=args.clean)

    #--plot with errorbands--
    emptyy = [min(a_data['depthhard'])]
    if args.soft == 'yes': 
        plot_a(ax1,a_data['age'],a_data['depthsoft'],low=a_data['depthsoftlow'],high=a_data['depthsofthigh'],syms=a_syms,colors=softcolor,alphas=a_alpha,mecs=a_mecs,errorband=True,mews=a_mews,line='-')
        plt.scatter(emptyx,emptyy,color=softcolor,marker='s',label='0.5-2.0 keV')
    if args.soft12 == 'yes': 
        plot_a(ax1,a_data['age'],a_data['depthsoft12'],low=a_data['depthsoft12low'],high=a_data['depthsoft12high'],syms=a_syms,colors=soft12color,alphas=a_alpha,mecs=a_mecs,errorband=True,mews=a_mews,line='-')
        plt.scatter(emptyx,emptyy,color=soft12color,marker='s',label='1.0-2.0 keV')
    if args.soft13 == 'yes': 
        plot_a(ax1,a_data['age'],a_data['depthsoft13'],low=a_data['depthsoft13low'],high=a_data['depthsoft13high'],syms=a_syms,colors=soft13color,alphas=a_alpha,mecs=a_mecs,errorband=True,mews=a_mews,line='-')
        plt.scatter(emptyx,emptyy,color=soft13color,marker='s',label='1.0-3.0 keV')
    if args.soft23 == 'yes': 
        plot_a(ax1,a_data['age'],a_data['depthsoft23'],low=a_data['depthsoft23low'],high=a_data['depthsoft23high'],syms=a_syms,colors=soft23color,alphas=a_alpha,mecs=a_mecs,errorband=True,mews=a_mews,line='-')
        plt.scatter(emptyx,emptyy,color=soft23color,marker='s',label='2.0-3.0 keV')
    if args.hard == 'yes': 
        plot_a(ax1,a_data['age'],a_data['depthhard'],low=a_data['depthhardlow'],high=a_data['depthhardhigh'],syms=a_syms,colors=hardcolor,alphas=a_alpha,mecs=a_mecs,errorband=True,mews=a_mews,line='-')
        plt.scatter(emptyx,emptyy,color=hardcolor,marker='s',label='3.0-8.0 keV')

    plt.scatter(emptyx,emptyy,facecolor='none',edgecolor=offcolor,linewidth=offmew,marker='s',label='Low Chip Y',alpha=1)

    #--plot ratios from count rates--
    #y_min,y_max = plt.ylim()
    #plt.plot(count_age,count_ratio_soft,'^',color='blue')
    #plt.plot(count_age,count_ratio_soft12,'^',color='green')
    #plt.plot(count_age,count_ratio_soft13,'^',color='cyan')
    #plt.plot(count_age,count_ratio_hard,'^',color='red')
    #plt.scatter(emptyx,[y_min],color='black',marker='^',label='Count Rate Ratios')

    #--plot legend--
    leg = plt.legend(numpoints=1,frameon=frame,fontsize='small',loc='upper left',scatterpoints=1)
    lframe = leg.get_frame()
    lframe.set_facecolor('white')
    lframe.set_edgecolor('white')


    #plt.plot(a_age,a_abs_hard+1.0,marker='o',color='red',markersize=psize)
    #plt.plot(a_age,a_abs_soft23+1.0,marker='o',color='magenta',markersize=psize)
    #plt.plot(a_age,a_abs_soft13+1.0,marker='o',color='cyan',markersize=psize)
    #plt.plot(a_age,a_abs_soft12+1.0,marker='o',color='green',markersize=psize)
    #plt.plot(a_age,a_abs_soft+1.0,marker='o',color='blue',markersize=psize)

    #--plot top axis--
    set_axis(ax1,x=topx,twin=True,title=toptitle,xlab=toplab,
             clean=args.clean)

    #--close plot--
    pdffile.savefig()
    plt.close()


#---------------------------------------
#     Hot/Cool Normalization Ratios
#---------------------------------------

#--Initialize Plot--
ax1 = start_plot(xtitle='SN1987A Age [days]',ytitle='Norm$_{Hot}$/Norm$_{Cool}$',xmin=args.agemin,xmax=args.agemax,ymax=0.5,ymin=0.0,clean=args.clean)

#--calculate ratio--
#normratio = b_normhot/b_normcool
normratio_low,normratio_high = ratio_err(a_data['normhot'],a_data['normcool'],a_data['normhotlow'],a_data['normhothigh'],a_data['normcoollow'],a_data['normcoolhigh'])

#--plot ratios--
plot_a(ax1,a_data['age'],a_data['normhot']/a_data['normcool'],low=normratio_low,high=normratio_high,syms=a_syms,colors='black',alphas=a_alpha,mecs=a_mecs,mews=a_mews,sizes=a_sizes)
#ax1.errorbar(b_age,normratio,yerr=[normratio_low,normratio_high],marker='o')

#--Plot Extras--
if args.plotextras != None:
    plot_extras(names=args.plotextras)

#--plot legend--
#plt.legend(numpoints=1,frameon=frame,fontsize='xx-small',loc='upper left')

#--set top axis--
set_axis(ax1,x=topx,twin=True,title=toptitle,xlab=toplab,clean=args.clean)

#--save and close plot--
pdffile.savefig()
plt.close()

#---------------------------------------
#       kT (hot and cold) vs Age
#---------------------------------------

#--Initialize Plot--
ax1 = start_plot(xtitle='SN1987A Age [days]',ytitle='kT [keV]',xmin=args.agemin,xmax=args.agemax,ymin=0.0,ymax=3.5,clean=args.clean)

#--plot temperatures--
#legend
#blue = cool
#red = hot
plot_a(ax1,a_data['age'],a_data['kThot'],low=a_data['kThotlow'],high=a_data['kThothigh'],syms=a_syms,colors='red',alphas=a_alpha,mecs=a_mecs,mews=a_mews,sizes=a_sizes)
plot_a(ax1,a_data['age'],a_data['kTcool'],low=a_data['kTcoollow'],high=a_data['kTcoolhigh'],syms=a_syms,colors='blue',alphas=a_alpha,mecs=a_mecs,mews=a_mews,sizes=a_sizes)

#--Plot Extras--
if args.plotextras != None:
    plot_extras(names=args.plotextras)

#--plot legend--
#plt.legend(numpoints=1,frameon=frame,fontsize='xx-small',loc='upper right')

#--set top axis--
set_axis(ax1,x=topx,twin=True,title=toptitle,xlab=toplab,clean=args.clean)

#--save and close plot--
pdffile.savefig()
plt.close()

#---------------------------------------
#     Hot/Cool kT Ratios
#---------------------------------------

#--Initialize Plot--
ax1 = start_plot(xtitle='SN1987A Age [days]',ytitle='kT$_{Hot}$/kT$_{Cool}$',xmin=args.agemin,xmax=args.agemax,ymax=17.0,clean=args.clean)

#--calculate ratio--
#kTratio = b_kThot/b_kTcool
kTratio_low,kTratio_high = ratio_err(a_data['kThot'],a_data['kTcool'],a_data['kThotlow'],a_data['kThothigh'],a_data['kTcoollow'],a_data['kTcoolhigh'])

#--plot ratios--
plot_a(ax1,a_data['age'],a_data['kThot']/a_data['kTcool'],low=kTratio_low,high=kTratio_high,syms=a_syms,colors='black',alphas=a_alpha,mecs=a_mecs,mews=a_mews,sizes=a_sizes)

#--Plot Extras--
if args.plotextras != None:
    plot_extras(names=args.plotextras)

#--plot legend--
#plt.legend(numpoints=1,frameon=frame,fontsize='xx-small',loc='upper left')

#--set top axis--
set_axis(ax1,x=topx,twin=True,title=toptitle,xlab=toplab,clean=args.clean)

#--save and close plot--
pdffile.savefig()
plt.close()


#---------------------------------------
#     Ionization Age vs Age
#---------------------------------------

#--Initialize Plot--
ax1 = start_plot(xtitle='SN1987A Age [days]',ytitle='Ionization Age [10$^{11}$ s cm$^{-3}$]',xmin=args.agemin,xmax=args.agemax,ymin=0.0,ymax=3.5,clean=args.clean)

#--plot tau--
plot_a(ax1,a_data['age'],a_data['tauhot'],low=a_data['tauhotlow'],high=a_data['tauhothigh'],syms=a_syms,colors='black',alphas=a_alpha,mecs=a_mecs,mews=a_mews,sizes=a_sizes)

#--Plot Extras--
if args.plotextras != None:
    plot_extras(names=args.plotextras)

#--set top axis--
set_axis(ax1,x=topx,twin=True,title=toptitle,xlab=toplab,clean=args.clean)

#--plot legend--
#plt.legend(numpoints=1,frameon=frame,fontsize='xx-small',loc='upper left')

#--save and close plot--
pdffile.savefig()
plt.close()

#---------------------------------------
#     Quadrant Count Fractions
#---------------------------------------

if args.quadcounts == 'yes':

    #--Set up colors and symbols--
    nwcolor = 'purple'
    swcolor = 'blue'
    secolor = 'green'
    necolor = 'orangered'
    
    count_nobs = ctable.shape[0]
    countsyms = [hetg_sym]*count_nobs
    count_bare_i = np.where(ctable['grating'] == 'NONE')
    for i in list(count_bare_i[:][0]):
        countsyms[i] = bare_sym 

    c_alphas = 0.6
#    c_sizes = a_sizes
    c_sizes = 6
    band_alpha = 0.1
    
    #--Initialize Plot--
    ax1 = start_plot(xtitle='SN1987A Age [days]',ytitle='Fraction of Counts in Quadrant',xmin=args.agemin,xmax=args.agemax,ymin=0.1,ymax=0.42,clean=args.clean)

    #--plot fractions--

    fancy_plot(ctable['age'],ctable['nwfraction'],yerror=
               ctable['nwfraction']-ctable['nwfraction_errlow'],syms=
               countsyms,colors=nwcolor,sizes=c_sizes,errorband=True,
               alphas=c_alphas,mecs = nwcolor,bandalpha=band_alpha)
    fancy_plot(ctable['age'],ctable['nefraction'],yerror=
               ctable['nefraction']-ctable['nefraction_errlow'],syms=
               countsyms,colors=necolor,sizes=c_sizes,errorband=True,
               alphas=c_alphas,mecs=necolor,bandalpha=band_alpha)
    fancy_plot(ctable['age'],ctable['swfraction'],yerror=
               ctable['swfraction']-ctable['swfraction_errlow'],syms=
               countsyms,colors=swcolor,sizes=c_sizes,errorband=True,
               alphas=c_alphas,mecs=swcolor,bandalpha=band_alpha)
    fancy_plot(ctable['age'],ctable['sefraction'],yerror=
               ctable['sefraction']-ctable['sefraction_errlow'],syms=
               countsyms,colors=secolor,sizes=c_sizes,errorband=True,
               alphas=c_alphas,mecs=secolor,bandalpha=band_alpha)

    #--Plot Extras--
    if args.plotextras != None:
        plot_extras(names=args.plotextras)
    
    #--plot legend--
#    plt.scatter(emptyx,emptyy,color=nwcolor,marker='',label='NW',
#                linewidths=1.0)
#    plt.scatter(emptyx,emptyy,color=necolor,marker='',label='NE',
#                linewidths=1.0)
#    plt.scatter(emptyx,emptyy,color=swcolor,marker='',label='SW',
#                linewidths=1.0)
#    plt.scatter(emptyx,emptyy,color=secolor,marker='',label='SE',
#                linewidths=1.0)

 #   leg = plt.legend(loc='upper right',fontsize='small',frameon=frame,
 #                    scatterpoints=1,numpoints=1,markerscale=1.5)

#    lab1,lab2,lab3,lab4 = leg.get_texts()
#    lab1.set_color(nwcolor)
#    lab2.set_color(necolor)
#    lab3.set_color(swcolor)
#    lab4.set_color(secolor)

    #--set top axis--
    set_axis(ax1,x=lightcurve_topx,twin=True,title=toptitle,
             xlab=lightcurve_toplab,clean=args.clean)
    #set_axis(ax1,x=topx,twin=True,title=toptitle,xlab=toplab)
    
    #--plot pie chart as legend--

    #-create inset axes-
    pax = plt.axes([0.7,0.67,0.2,0.2])
    #-plot pie chart figure-
    wedges,texts=pax.pie([25,25,25,25],labels=['NW','NE','SE','SW'],
            colors=[nwcolor,necolor,secolor,swcolor],radius=1.0)#,
#            labeldistance=0.4)
    pax.set_aspect('equal')
    #-set label properties-
    yoffset = 0.35
    xoffset = 0.43
    labely = [yoffset-0.05,yoffset-0.05,-yoffset,-yoffset]
    labelx = [xoffset,-xoffset,-xoffset,xoffset]
    q = 0
    for t in texts:
        t.set_horizontalalignment('center')
        t.set_verticalalignment('center')
        t.set_y(labely[q])
        t.set_x(labelx[q])
        t.set_color('white')
        q = q+1
    for w in wedges:
        w.set_linewidth(0)

    #--save and close plot--
    pdffile.savefig()
    plt.close()

#---------------------------------------
#          Expansion Curve
#---------------------------------------
plotradius = 'no'
if plotradius == 'yes':

    #--Initialize Plot--
    ax1 = start_plot(xtitle='SN1987A Age [days]',ytitle='Radius [arcsec]',xmin=4000,xmax=args.agemax,ymin=0.5,ymax=0.85,clean=args.clean)

    #--plot radii--
    #plot_a(ax1,r_data['age'],r_data['radius'],low=r_data['radiuslow'],high=r_data['radiushigh'],syms=r_syms,colors='black',alphas=r_alpha,mecs=r_mecs,mews=r_mews,sizes=5.0)
    plot_a(ax1,r_data['age'],r_data['radius'],low=r_data['radiuslow'],high=r_data['radiushigh'],syms='o',colors='black',alphas=r_alpha,mecs='black',mews=0.5,sizes=5.0)

    #--plot best fit--
    #columns:
    #stat  v_early(km/s)  v_late(km/s)  intercept(arcsec)  changepoint(day)
    #rows:
    #median, mean, standard deviation
    fitparams = [0.0,convert_units(rfit_data['v_early'][0],'km/s','arcsec/day'),convert_units(rfit_data['v_late'][0],'km/s','arcsec/day'),rfit_data['intercept'][0],rfit_data['changepoint'][0]]
    expx = [0,rfit_data['changepoint'][0],11000]
    expy = fitting.arr_brokenlinear(expx,fitparams)

    plot_a(ax1,expx,expy,syms='',line='-',colors='black')

    #--plot fitted data points with error band--
    #plot_a(ax1,rmodel_data['age'],rmodel_data['radius'],low=rmodel_data['radiuslow'],high=rmodel_data['radiushigh'],syms='',colors='blue',alphas=r_alpha,mecs='blue',mews=0.5,errorband=True)

    verbose = 'yes'
    if verbose == 'yes':
        change_radius = fitting.brokenlinear(fitparams[4],fitparams)
        print 'change_radius [arcsec,pc,cm] = ',change_radius,convert_units(change_radius,'arcsec','pc'),convert_units(change_radius,'arcsec','cm')
        opt_radius = fitting.brokenlinear(4999,fitparams)
        print 'day 4999 radius [arcsec,pc,cm]',opt_radius,convert_units(opt_radius,'arcsec','pc'),convert_units(opt_radius,'arcsec','cm')

        lastrad = r_data['radius'][-1]
        print 'last radius [arcsec,pc,cm] = ',lastrad,convert_units(lastrad,'arcsec','pc'),convert_units(lastrad,'arcsec','cm')
        lastrad_errlow = r_data['radiuslow'][-1]
        print 'last radius err low = ',lastrad_errlow,convert_units(lastrad_errlow,'arcsec','pc'),convert_units(lastrad_errlow,'arcsec','cm')
        lastrad_errhigh = r_data['radiushigh'][-1]
        print 'last radius err high = ',lastrad_errhigh,convert_units(lastrad_errhigh,'arcsec','pc'),convert_units(lastrad_errhigh,'arcsec','cm')
        print 'changepoint = ',fitparams[4],' days'

        day = 5789
        changerad = r_data['radius'][np.where(r_data['age'] == day)]
        changerad_low = r_data['radiuslow'][np.where(r_data['age'] == day)]
        changerad_high = r_data['radiushigh'][np.where(r_data['age'] == day)]
        print 'changerad [arcsec,pc,cm] = ',changerad,convert_units(changerad,'arcsec','pc'),convert_units(changerad,'arcsec','cm')
        print 'changerad_low = ',changerad_low,convert_units(changerad_low,'arcsec','pc'),convert_units(changerad_low,'arcsec','cm')
        print 'changerad_high = ',changerad_high,convert_units(changerad_high,'arcsec','pc'),convert_units(changerad_high,'arcsec','cm')

        day = 5036
        optrad = r_data['radius'][np.where(r_data['age'] == day)]
        optrad_low = r_data['radiuslow'][np.where(r_data['age'] == day)]
        optrad_high = r_data['radiushigh'][np.where(r_data['age'] == day)]


    #--plot ring width as bands--
    print expx
    print opt_radius
    #plt.fill_between(expx,changerad[0],lastrad,color='red',alpha=0.3)
    #plt.fill_between(expx,optrad[0],lastrad,color='blue',alpha=0.2)
    #plt.plot(expx,[changerad[0]]*3,line='-',color='red')
    #plt.plot(expx,[lastrad]*3,line='-',color='blue')
    #plt.plot(expx,[optrad[0]]*3,line='-',color='orange')

    #--plot legend--
    #plt.legend(numpoints=1,frameon=frame,fontsize='xx-small',loc='upper left')
    emptyx = [3000]
    emptyy = [0.6]
    plt.scatter(emptyx,emptyy,color='black',marker=hetg_sym,label='with HETG',alpha=a_alpha)
    plt.scatter(emptyx,emptyy,color='black',marker=bare_sym,label='bare ACIS',alpha=a_alpha)
    plt.plot(emptyx,emptyy,markerfacecolor='none',markeredgecolor=offcolor,markeredgewidth=offmew,marker='s',label='Low Chip Y',alpha=1,linestyle='None')

    #plt.legend(loc='lower right',fontsize='xx-small',frameon=frame,scatterpoints=1,numpoints=1,markerscale=1)

    #--set top axis--
    set_axis(ax1,x=ext_topx,twin=True,title=toptitle,xlab=ext_toplab,clean=args.clean)

    #--save and close plot--
    pdffile.savefig()
    plt.close()




#---------------------------------------
#             Wrap Up
#---------------------------------------

#----Close File----
pdffile.close()


