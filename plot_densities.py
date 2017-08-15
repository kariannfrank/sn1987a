#! /usr/bin/env python

#Author: Kari A. Frank
#Date: June 25, 2014
#Purpose: Calculate and plot the densities of the SN1987A equatorial ring based on the ACIS spectral fits and radius fitting results.
#Usage: plot_densities.py [--infile infile] [--radiusfile radiusfile] [--clean clean] [--ylog ylog] [--plotfit plotfit] [--dofit dofit] [--clobber clobber] [--lateonly lateonly] [--agemin agemin] [--agemax agemax] [--fluxmax fluxmax] [--ignorechippos ignorechippos] [--official yes]
#
#Input:
#
# infile:     optionally specify file with observation and spectra data
#              default is infile =  'spectra465_frank2015a_fits.txt'
#
# radiusfile: text file containing radius fit results, of same format created
#             by get_radii_results.py (default='radii_fits.txt')
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
# lateonly: plot only the late-time observations on the light curve
#            (day 8000+)
#
# official: switch to only include observations in the official list,
#           $sna/comparison_dir/official_obsid_list.txt (default='no')
#
# agemin/agemax: set the minimum and maximum age (in day) to plot
#
# fluxmax: optionally set the maximum of flux axis in light curve 
#            (in 10^-13 erg/cm^2/s)
#
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
from astro_utilities import norm_to_em
import fitting
#sys.path.append('/astro/research/kaf33/Dropbox/Research/Python_Programs/')
sys.path.append('/Users/kafrank/Dropbox/Research/Python_Programs/')
from fancy_plot import fancy_plot
from scipy.optimize import curve_fit
import analytical_functions as af

#---------------------------------------
#      Define Helper Functions
#---------------------------------------

def plot_a(ax,age,flux,low=None,high=None,syms='o',colors='b',sizes=4,alphas=1.0,mecs='black',mews=0.5,line='',errorband=False,label='_nolegend_'):        
    if low is not None:
        elow,ehigh = ebars(flux,low,high)
    else:
        elow = None
        ehigh = None
    fancy_plot(age,flux,yerror=elow,yerror_high=ehigh,syms=syms,colors=colors,sizes=sizes,alphas=alphas,mecs=mecs,errorband=errorband,line=line,mews=mews)
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
    top = np.array(top)
    top_low = np.array(top_low)
    top_high = np.array(top_high)
    bottom = np.array(bottom)
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

def set_axis(ax,x=None,twin=False,title=None,xlab=None,grid=False):

#    font_size = 10
    font_size = 13
#    fontangle = 50
    fontangle = 40

    if twin == False:
        fontangle = 0

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
        else:
            tick.label2.set_fontsize(font_size)
            tick.label2.set_rotation(fontangle)
    #ax.set_yscale('log')

def start_plot(xtitle='',ytitle='',xmin=None,xmax=None,ymin=None,ymax=None,ylog=False):

    #-initialize main plot-
    fig,ax1=plt.subplots()

    #-set axis titles-
    plt.xlabel(xtitle)
    plt.ylabel(ytitle)

    #-set bottom axis-
    if ylog == True: plt.yscale('log')
    set_axis(ax1,x=botx,grid=mgrid)
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
#sn1987a_dir = '/export/bulk/rice1/kaf33/main/Chandra_Observations/SN1987A/'
#brevis
sn1987a_dir = '/Users/kafrank/Research/SN1987A/'

parser.add_argument('--infile',help='File containing observation and spectra data.',default=sn1987a_dir+'comparison_dir/images_current/radii_norms.txt')
parser.add_argument('--radiusfile',help='Text file with image fitting results.',default=sn1987a_dir+'comparison_dir/images_current/radii_fits_300-8000.txt')

parser.add_argument('--clean',help='Determines x-axis ticks and labels.',default='yes')

parser.add_argument('--ylog',help='Plot light curve with log scale y-axis.',default='no')

parser.add_argument('--plotfit',help='Plot best fit contamination.',default='no')

parser.add_argument('--dofit',help='Run mcmc on contamination absorption.',default='no')

parser.add_argument('--early',help='Plot (non-pileup corrected) fluxes of observations 1387 and 122.',default='yes')

parser.add_argument('--lateonly',help='Plot only day 8000+ on light curve.',default='no')

parser.add_argument('--official',help='Include only obsids in the official list.',default='no')

parser.add_argument('--agemin',help='Minimum age for x-axes.',default=0)
parser.add_argument('--agemax',help='Maximum age for x-axes.',default=0)

parser.add_argument('--fluxmax',help='Maximum flux for light curve y-axis.',default=0)

parser.add_argument('--instcolor',help='Switch meaning of symbols and colors.',default='no')

parser.add_argument('--ignorechippos',help='Plot normal borders for observations with offset chip positions.',default='no')

#not currently implemented!
parser.add_argument('--clobber',help='Overwrite existing files.',default='no')

args = parser.parse_args()

# set constants
z = 0.00095
dist_kpc = 51.4
dist_km = 1.586*10.0**18.
dist_cm = dist_km*10.**5.


#----Set dependent arguments----

args.agemin = int(args.agemin)
args.agemax = int(args.agemax)
args.fluxmax = float(args.fluxmax)

if args.ylog == 'yes': 
    args.ylog = True
else:
    args.ylog = False

#----Set file paths----

plotfile = './spectra465_frank2015a_densities.pdf'

obsfile = sn1987a_dir+'comparison_dir/chandra_observations.txt'

#radiusfile =sn1987a_dir+'comparison_dir/imaging/sn1987a_radii_fits.txt'

radiusfitfile = sn1987a_dir+'comparison_dir/imaging/expansion_fits/sn1987a_mpfit_radii_fits.txt'

#---------------------------------------
#        Read and Convert Data
#---------------------------------------

#----Read spectral results from file----

#a_data = np.genfromtxt(args.infile,skip_header=0,names=True,comments='#',dtype=None)

#columns names:
#obsid	date		age	grating	caldb	model 	pcounts	psoft psoftlow psofthigh	phard phardlow phardhigh	pbroad pbroadlow pbroadhigh	psoft12 psoft12low psoft12high	Counts	soft softlow softhigh	hard hardlow hardhigh	broad broadlow broadhigh	soft12 soft12low soft12high	soft23 soft23low soft23high	soft13 soft13low soft13high	kTcool	kTcoollow kTcoolhigh	kThot	kThotlow kThothigh	Taucool Taucoolerr	tauhot	tauhotlow tauhothigh	normcool normcoollow normcoolhigh	normhot normhotlow normhothigh	chi2	dof	redchi2		csoft12 csoft12low csoft12high	csoft csoftlow csofthigh	csoft23 csoft23low csoft23high	csoft13 csoft13low csoft13high	chard chardlow chardhigh


#----Read Radius Measurements----
r_data = np.genfromtxt(args.infile,skip_header=0,names=True,comments='#',dtype=None)
#columns:
#obsid   age     band    R0      R0errl  R0erru  SIGR    SIGRerrl        SIGRerru        counts  deconcounts
#     SKYNS   SKYNSerrl       SKYNSerru
# (units in arcsec)

#--convert errorbars to error intervals--
#r_data['R0errl'] = r_data['R0'] - r_data['R0errl']
#r_data['R0erru'] = r_data['R0'] + r_data['R0erru']


#----Get Observation Info----

#--Get Observation Years (and ages if necessary)--
r_obsyears = [str(get_obs_time(int(obs),get='year')) for obs in r_data['obsid']]
r_data = append_fields(r_data,names='year',data=r_obsyears)
if 'age' not in r_data.dtype.names:
    r_age = [get_obs_time(obs) for obs in r_data['obsid']]
    r_data = append_fields(r_data,names='age',data=r_age)

#--Read Main Info File--
obs_info = np.genfromtxt(obsfile,skip_header=0,names=True,comments='#',dtype=None)
#column names:
#obsid	date	age	pi	configuration	grating	exposure	frametime	simoffset	yoffset	zoffset

#-get years-
obsyears = [str(get_obs_time(int(obs),get='year')) for obs in obs_info['obsid']]
obs_info = append_fields(obs_info,names='year',data=obsyears)

#--Get FrameTimes--
r_frames = [obs_info['frametime'][np.where(obs_info['obsid'] == obs)] for obs in r_data['obsid']]
r_data = append_fields(r_data,names='frametime',data=np.array(r_frames)[:,0])

#--Get SimZ--

#-convert 'default' to '0.0'-
obs_info['simoffset'][np.where(obs_info['simoffset'] == 'default')] = '0.0'
#print i_sim,i_obsid

r_sims = [obs_info['simoffset'][np.where(obs_info['obsid'] == obs)] for obs in r_data['obsid']]
r_data = append_fields(r_data,names='sim',data=np.array(r_sims)[:,0])

#--get gratings for radius observations--
#r_gratings = [obs_info['grating'][np.where(obs_info['obsid'] == obs)] for obs in r_data['obsid']]
#r_data = append_fields(r_data,names='grating',data=np.array(r_gratings)[:,0])


#---------------------------------------
#        Estimate Densities
#---------------------------------------

#--Reduce to only include 300-8000 band radii--
#r_data = r_data[np.where(r_data['band']=='300-8000')]
#r0_data = r_data['R0']['band'=='300-8000']
#r_data = np.ma.compress_rows(r_data)


#--Calculate Emission Measures--
coolem = norm_to_em(r_data['norm1']/1000.,dist_cm,z)
hotem = norm_to_em(r_data['norm2']/1000.,dist_cm,z)

#coolem_low = norm_to_em(r_data['normcoollow'],dist_cm,z)
#coolem_high = norm_to_em(r_data['normcoolhigh'],dist_cm,z)

#hotem_low = norm_to_em(a_data['normhotlow'],dist_cm,z)
#hotem_high = norm_to_em(a_data['normhothigh'],dist_cm,z)

#--Calculate volumes (in cm^3)--
volumes = np.pi*convert_units(r_data['SigR'],'arcsec','cm')**2.0*2.0*np.pi*convert_units(r_data['Radius'],'arcsec','cm')
#volumes = np.pi*convert_units(0.2,'arcsec','cm')**2.0*2.0*np.pi*convert_units(r_data['Radius'],'arcsec','cm')

#--Calculate densities, electron--
#cooldensities = 1.21*(coolem/(volumes*0.05)/1.21)**0.5
cooldensities = 1.21*(coolem/volumes/1.21)**0.5
hotdensities = 1.21*(hotem/volumes/1.21)**0.5

#--Append density columns
r_data = append_fields(r_data,names='cooldensity',data=cooldensities)
#a_data = append_fields(a_data,names='cooldensitylow',data=cooldensities_low)
#a_data = append_fields(a_data,names='cooldensityhigh',data=cooldensities_high)
r_data = append_fields(r_data,names='hotdensity',data=hotdensities)
#a_data = append_fields(a_data,names='hotdensitylow',data=hotdensities_low)
#a_data = append_fields(a_data,names='hotdensityhigh',data=hotdensities_high)

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
    args.agemin = 5000
        
if args.agemax == 0: args.agemax = 11000
        
#if args.fluxmax == 0:
#    fluxmax = 1.1*max(r_data['soft'][np.where(a_data['age']<=args.agemax)])       
#else: fluxmax = args.fluxmax

if args.ylog == True:
    fluxmin = 0.3
else:
    fluxmin = 0.0

#----Set Axis Labels and Ticks----

#-calculate required top axis labels automatically based on agemin/agemax-
print 'agemin = ',args.agemin
years = range(convert_time(args.agemin,get='year')+1,convert_time(args.agemax,get='year')+1)
year_ticks = [str(yr) for yr in years]
regular_ages = [convert_time(yr+'-01-01',get='age',informat='date') for yr in year_ticks]
print 'yeartick = ',year_ticks
print 'regular ages = ',regular_ages

radyears = range(1999,2016)
radius_year_ticks = [str(yr) for yr in radyears]
radius_regular_ages = [convert_time(yr+'-01-01',get='age',informat='date') for yr in radius_year_ticks]

spectrayears = ['2001','2002','2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014']
spectra_year_ticks = [str(yr) for yr in spectrayears]
spectra_regular_ages = [convert_time(yr+'-01-01',get='age',informat='date') for yr in spectra_year_ticks]

if args.clean == 'yes':
    lightcurve_topx = regular_ages
    lightcurve_toplab = year_ticks
    topx = spectra_regular_ages
    toplab = spectra_year_ticks
    ext_topx = radius_regular_ages
    ext_toplab = radius_year_ticks
    botx = None
    toptitle = 'Year'
    mgrid = False
    frame = False
else:
    topx = obs_info['age']#r_data['age']
    toplab = obs_info['year']#r_data['year']
    lightcurve_topx = topx
    lightcurve_toplab = toplab
    ext_topx = topx
    ext_toplab = toplab
    toptitle = 'Observation Year'
    botx = obs_info['age']#r_data['age']
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
r_offz_i = np.where(r_data['sim'] == '-8.42')

#-ACIS grating-
r_hetg_i = np.where(r_data['grating'] == 'HETG')
r_bare_i = np.where(r_data['grating'] == 'NONE')

#--make lists for each--

#-colors = energy band-
hotcolor = 'red'
coolcolor = 'blue'

#-number of observations-
r_nobs = r_data.shape[0]
print 'nobs = ',r_nobs

#-symbol = instrument (ACIS-HETG,ACIS-NONE,XMM-EPIC-pn)-
hetg_sym = 'o' #HETG = circle
bare_sym = 'd' #bare acis = thin diamond
letg_sym = 'p' #LETG = hexagon
r_syms = [hetg_sym]*r_nobs 
for i in list(r_bare_i[:][0]):
    r_syms[i] = bare_sym 

chandrahardsyms = ['s']*r_nobs
chandrasoftsyms = ['o']*r_nobs
#for i in list(a_bare_i[:][0]):
#    chandrahardsyms[i] = 'H'
#    chandrasoftsyms[i] = 'p'

#-markeredgecolor = simz (detector z)-
if args.instcolor == 'no':
    offcolor = 'gray'
else: 
    offcolor = 'black'
r_mecs = ['black']*r_nobs
if args.ignorechippos == 'no':
    for i in list(r_offz_i[:][0]):
        r_mecs[i] = offcolor #offset observations have offcolor outline

#-markeredgewidth = simz (detector z)-
offmew = 1.5
r_mews = [0.5]*r_nobs
if args.ignorechippos == 'no':
    for i in list(r_offz_i[:][0]):
        r_mews[i] = offmew #offset observations have different outline

#-markersize - # = frametime-
#a_sizes = [f*3 for f in a_data['frametime']]
r_sizes = 6

#-alpha (transparency) = caldb-
r_alpha = 1.0

#---------------------------------------
#        Plots
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
ax1 = start_plot(xtitle='SN1987A Age [days]',ytitle='Density [cm$^{-3}$]',xmin=args.agemin,xmax=args.agemax,ylog=args.ylog)


#--Plot age vs cool densities--
plot_a(ax1,r_data['age'],r_data['cooldensity'],syms=r_syms,colors=coolcolor,alphas=r_alpha,mecs=r_mecs,mews=r_mews,sizes=r_sizes)

#--Plot age vs hot densities--
plot_a(ax1,r_data['age'],r_data['hotdensity'],syms=r_syms,colors=hotcolor,alphas=r_alpha,mecs=r_mecs,mews=r_mews,sizes=r_sizes)

if args.ylog == True: ax1.set_yscale('log')

#--Add Top (year) x-axis--
#note this must come after plotting things on the first axis!

set_axis(ax1,x=lightcurve_topx,twin=True,title=toptitle,xlab=lightcurve_toplab)#,xlab=b_obsyears)

#--Close Plot--
pdffile.savefig()
plt.close()

#----Plot emission measures----

#--Initialize Plot--
ax1 = start_plot(xtitle='SN1987A Age [days]',ytitle='Emission Measures',xmin=args.agemin,xmax=args.agemax,ylog=args.ylog)

#--Plot age vs cool densities--
plot_a(ax1,r_data['age'],coolem,syms=r_syms,colors=coolcolor,alphas=r_alpha,mecs=r_mecs,mews=r_mews,sizes=r_sizes)

#--Plot age vs hot densities--
plot_a(ax1,r_data['age'],hotem,syms=r_syms,colors=hotcolor,alphas=r_alpha,mecs=r_mecs,mews=r_mews,sizes=r_sizes)

if args.ylog == True: ax1.set_yscale('log')

#--Add Top (year) x-axis--
#note this must come after plotting things on the first axis!

set_axis(ax1,x=lightcurve_topx,twin=True,title=toptitle,xlab=lightcurve_toplab)#,xlab=b_obsyears)

#--Close Plot--
pdffile.savefig()
plt.close()

#----Plot emission measures----

#--Initialize Plot--
ax1 = start_plot(xtitle='SN1987A Age [days]',ytitle='Torus Volume [cm$^3$]',xmin=args.agemin,xmax=args.agemax,ylog=args.ylog)

#--Plot age vs cool densities--
plot_a(ax1,r_data['age'],volumes,syms=r_syms,colors='black',alphas=r_alpha,mecs=r_mecs,mews=r_mews,sizes=r_sizes)

if args.ylog == True: ax1.set_yscale('log')

#--Add Top (year) x-axis--
#note this must come after plotting things on the first axis!

set_axis(ax1,x=lightcurve_topx,twin=True,title=toptitle,xlab=lightcurve_toplab)#,xlab=b_obsyears)

#--Close Plot--
pdffile.savefig()
plt.close()

#---------------------------------------
#             Wrap Up
#---------------------------------------

#----Close File----
pdffile.close()


