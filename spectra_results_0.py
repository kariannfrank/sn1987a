#! /usr/bin/env python

#Author: Kari A. Frank
#Date: June 25, 2014
#Purpose: Plot and fit the radial expansion curve of SN1987A.
#Usage: spectra_results.py [--infile infile] [--clean clean] [--plotfit plotfit] [--dofit dofit] [--plothelder plothelder] [--clobber clobber]
#
#Input:
#
# infile_a:     optionally specify file with observation and spectra data
#              default is infile_a = 
#             '/export/bulk/rice3/kaf33/Chandra_Observations/SN1987A/comparison_dir/spectra459/spectra_fits.txt'
# infile_b:     optionally specify second file with observation and spectra data
#              default is infile_b = 
#             '/export/bulk/rice3/kaf33/Chandra_Observations/SN1987A/comparison_dir/spectra462/spectra_fits462.txt'
#
# clean:      optional switch to change the axes labels and ticks to even
#             years and ages (default = 'no')
#
# dofit:      run mcmc fitting on the contamination absorption data
#
# plotfit:    plot best fit(s) on the contamination absorption plot
#
# plothelder: plot flux measurements from Helder2012 (default='no')
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
from sn1987a_time import *
import fitting
sys.path.append('/astro/research/kaf33/Dropbox/Research/Python_Programs/')
from fancy_plot import fancy_plot


#---------------------------------------
#      Define Helper Functions
#---------------------------------------

def plot_a(ax,age,flux,flux_low,flux_high,syms='o',colors='b',sizes=4,alphas=1,mecs='black'):        
    elow,ehigh = ebars(flux,flux_low,flux_high)
    fancy_plot(age,flux,yerror=elow,yerror_high=ehigh,syms=syms,colors=colors,sizes=sizes,alphas=alphas,mecs=mecs)
#    ax.errorbar(age,flux,yerr=[elow,ehigh],fmt='.',color=color,markersize=0,markeredgewidth=0)
#    ax.plot(age,flux,sym,markersize=4,markeredgewidth=0.5,color=color)

def plot_b(ax,age,flux,flux_low,flux_high,sym='o',color='b',alpha=1):
    ax.errorbar(age,flux,yerr=[flux-flux_low,flux_high-flux],fmt='.',markersize=0,color=color,markeredgewidth=0)
    ax.plot(age,flux,sym,markersize=4,color=color,markeredgewidth=0.5,alpha=alpha)

def ebars(x,x_low,x_high):
    errlow = x-x_low
    errhigh = x_high - x
    return errlow,errhigh
    
#-function to derive vectors of errors on ratios-
def ratio_err(top,bottom,top_low,top_high,bottom_low,bottom_high):
    #uses simple propagation of errors (partial derivatives)
    #note it returns errorbars, not interval

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
    ratio_low  = np.sqrt( np.square(np.divide(top_errlow,bottom)) + np.square( np.multiply(np.divide(top,np.square(bottom)),bottom_errlow)) )
    #-calculate ratio_high-
    ratio_high = np.sqrt( np.square(np.divide(top_errhigh,bottom)) + np.square( np.multiply(np.divide(top,np.square(bottom)),bottom_errhigh)) )
#    ratio_high = ((top_errhigh/bottom)**2.0 + (top/(bottom**2.0))*bottom_errhigh)**2.0)**0.5

    # return two vectors, err_low and err_high
    return ratio_low,ratio_high

def set_axis(ax,x=None,twin=False,title=None,xlab=None,grid=False):

    if twin == False:
        ax0 = ax
        ax0.xaxis.grid(grid,which='major')
        if x != None:
            ax0.set_xticks(x,minor=False)
    else:
        ax0 = ax.twiny()
#        x_min,x_max = plt.xlim()
        ax0.set_xlim(ax.get_xlim())
        if x != None:
            ax0.set_xticks(x,minor=False)
        if xlab != None:
            ax0.set_xticklabels(xlab)
        
    if title != None:
        ax0.set_xlabel(title)

    for tick in ax0.xaxis.get_major_ticks():
        if twin == False:
            tick.label.set_fontsize(10)
            tick.label.set_rotation(50)
        else:
            tick.label2.set_fontsize(10)
            tick.label2.set_rotation(50)
    #ax.set_yscale('log')

#---------------------------------------
#   Parse arguments and set defaults
#---------------------------------------

parser = argparse.ArgumentParser(description='Plot and fit SN1987A light curve and other spectral results.')

pwd = os.getcwd()

#parser.add_argument('required_arg',help='Description and usage for required_arg.',default='default for required_arg')

parser.add_argument('--infile_a',help='File containing observation and spectra data.',default='/export/bulk/rice3/kaf33/Chandra_Observations/SN1987A/comparison_dir/spectra459/spectra_fits.txt')

parser.add_argument('--infile_b',help='File containing observation and spectra data.',default='/export/bulk/rice3/kaf33/Chandra_Observations/SN1987A/comparison_dir/spectra462/spectra_fits462.txt')

parser.add_argument('--clean',help='Determines x-axis ticks and labels.',default='no')

parser.add_argument('--plotfit',help='Plot best fit contamination.',default='no')

parser.add_argument('--dofit',help='Run mcmc on contamination absorption.',default='no')

parser.add_argument('--plothelder',help='Plot fluxes from Helder2012',default='no')

#not currently implemented!
parser.add_argument('--clobber',help='Overwrite existing files.',default='no')


args = parser.parse_args()

#----Set file paths----
helderfile = '/astro/research/kaf33/Dropbox/Research/SN1987A/Helder_table1.txt'

plotfile = '/export/bulk/rice3/kaf33/Chandra_Observations/SN1987A/comparison_dir/spectra_results.pdf'

posteriorfile = '/export/bulk/rice3/kaf33/Chandra_Observations/SN1987A/comparison_dir/posteriors.pdf'

mcmcfile = '/export/bulk/rice3/kaf33/Chandra_Observations/SN1987A/comparison_dir/mcmc_contamfit.txt'

countratefile = '/export/bulk/rice3/kaf33/Chandra_Observations/SN1987A/comparison_dir/countrate_ratios.txt'

xmmfile = '/export/bulk/rice3/kaf33/Chandra_Observations/SN1987A/comparison_dir/xmm_results.txt'

obsfile = '/export/bulk/rice3/kaf33/Chandra_Observations/SN1987A/comparison_dir/chandra_observations.txt'

#---------------------------------------
#        Read and Convert Data
#---------------------------------------

#----Read CALDB4.5.9 Results from File----

#'p' prefix indicates piled spectral results
a_obsid,a_obsdate,a_age,a_grating,a_caldb,a_modelname,a_pcounts,a_pflux_soft,a_pflux_soft_low,a_pflux_soft_high,a_pflux_hard,a_pflux_hard_low,a_pflux_hard_high,a_pflux_broad,a_pflux_broad_low,a_pflux_broad_high,a_pflux_soft12,a_pflux_soft12_low,a_pflux_soft12_high,a_counts,a_flux_soft,a_flux_soft_low,a_flux_soft_high,a_flux_hard,a_flux_hard_low,a_flux_hard_high,a_flux_broad,a_flux_broad_low,a_flux_broad_high,a_flux_soft12,a_flux_soft12_low,a_flux_soft12_high,a_flux_soft23,a_flux_soft23_low,a_flux_soft23_high,a_flux_soft13,a_flux_soft13_low,a_flux_soft13_high,a_kTcool,a_kTcool_low,a_kTcool_high,a_kThot,a_kThot_low,a_kThot_high,a_taucool,a_taucool_err,a_tauhot,a_tauhot_low,a_tauhot_high,a_normcool,a_normcool_low,a_normcool_high,a_normhot,a_normhot_low,a_normhot_high,a_chi2,a_dof,a_redchi2,a_cflux_soft12,a_cflux_soft12_low,a_cflux_soft12_high,a_cflux_soft,a_cflux_soft_low,a_cflux_soft_high,a_cflux_soft23,a_cflux_soft23_low,a_cflux_soft23_high,a_cflux_soft13,a_cflux_soft13_low,a_cflux_soft13_high,a_cflux_hard,a_cflux_hard_low,a_cflux_hard_high = np.genfromtxt(args.infile_a,unpack=True,skip_header=1,comments='#',dtype=str)#'str','str','int','str','str','str','int','float','float','float','float','float','float','float','float','float','float','float','float','int','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float'))

a_obsid,a_obsdate,a_age,a_grating,a_caldb,a_modelname,a_pcounts,a_pflux_soft,a_pflux_soft_low,a_pflux_soft_high,a_pflux_hard,a_pflux_hard_low,a_pflux_hard_high,a_pflux_broad,a_pflux_broad_low,a_pflux_broad_high,a_pflux_soft12,a_pflux_soft12_low,a_pflux_soft12_high,a_counts,a_flux_soft,a_flux_soft_low,a_flux_soft_high,a_flux_hard,a_flux_hard_low,a_flux_hard_high,a_flux_broad,a_flux_broad_low,a_flux_broad_high,a_flux_soft12,a_flux_soft12_low,a_flux_soft12_high,a_flux_soft23,a_flux_soft23_low,a_flux_soft23_high,a_flux_soft13,a_flux_soft13_low,a_flux_soft13_high,a_kTcool,a_kTcool_low,a_kTcool_high,a_kThot,a_kThot_low,a_kThot_high,a_taucool,a_taucool_err,a_tauhot,a_tauhot_low,a_tauhot_high,a_normcool,a_normcool_low,a_normcool_high,a_normhot,a_normhot_low,a_normhot_high,a_chi2,a_dof,a_redchi2,a_cflux_soft12,a_cflux_soft12_low,a_cflux_soft12_high,a_cflux_soft,a_cflux_soft_low,a_cflux_soft_high,a_cflux_soft23,a_cflux_soft23_low,a_cflux_soft23_high,a_cflux_soft13,a_cflux_soft13_low,a_cflux_soft13_high,a_cflux_hard,a_cflux_hard_low,a_cflux_hard_high = np.genfromtxt(args.infile_a,skip_header=1,names=True,comments='#',dtype=None)#'str','str','int','str','str','str','int','float','float','float','float','float','float','float','float','float','float','float','float','int','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float','float'))

#----Read CALDB4.6.2 Results from File----

#'p' prefix indicates piled spectral results
b_obsid,b_obsdate,b_age,b_grating,b_caldb,b_modelname,b_pcounts,b_pflux_soft,b_pflux_soft_low,b_pflus_soft_high,b_pflux_hard,b_pflux_hard_low,b_pflux_hard_high,b_pflux_broad,b_pflux_broad_low,b_pflux_broad_high,b_pflux_soft12,b_pflux_soft12_low,b_pflux_soft12_high,b_counts,b_flux_soft,b_flux_soft_low,b_flux_soft_high,b_flux_hard,b_flux_hard_low,b_flux_hard_high,b_flux_broad,b_flux_broad_low,b_flux_broad_high,b_flux_soft12,b_flux_soft12_low,b_flux_soft12_high,b_flux_soft23,b_flux_soft23_low,b_flux_soft23_high,b_flux_soft13,b_flux_soft13_low,b_flux_soft13_high,b_kTcool,b_kTcool_low,b_kTcool_high,b_kThot,b_kThot_low,b_kThot_high,b_taucool,b_taucool_err,b_tauhot,b_tauhot_low,b_tauhot_high,b_normcool,b_normcool_low,b_normcool_high,b_normhot,b_normhot_low,b_normhot_high,b_chi2,b_dof,b_redchi2,b_cflux_soft12,b_cflux_soft12_low,b_cflux_soft12_high,b_cflux_soft,b_cflux_soft_low,b_cflux_soft_high,b_cflux_soft23,b_cflux_soft23_low,b_cflux_soft23_high,b_cflux_soft13,b_cflux_soft13_low,b_cflux_soft13_high,b_cflux_hard,b_cflux_hard_low,b_cflux_hard_high = np.genfromtxt(args.infile_b,unpack=True,skip_header=1,comments='#',dtype=str)#(str,str,int,str,str,str,int,float,float,float,float,float,float,float,float,float,float,float,float,int,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float))#,dtype=None)


#----Get Observation Info----

#--Get Observation Years--
#print obsid.shape
#obsyears = [obs.split('-')[0] for obs in obsdate]

a_obsyears = [str(get_obs_time(int(obs),get='year')) for obs in a_obsid]
b_obsyears = [str(get_obs_time(int(obs),get='year')) for obs in b_obsid]

#--Read Main Info File--
i_obsid,i_date,i_age,i_pi,i_config,i_grating,i_exp,i_frame,i_sim,i_y,i_z = np.genfromtxt(obsfile,unpack=True,skip_header=1,comments='#',dtype=str)

#--Get FrameTimes--
a_frames = [i_frame[np.where(i_obsid == obs)] for obs in a_obsid]
b_frames = [i_frame[np.where(i_obsid == obs)] for obs in b_obsid]

#--Get SimZ--

#-convert 'default' to '0.0'-
i_sim[np.where(i_sim == 'default')] = '0.0'
#print i_sim,i_obsid

#-associate with a,b data-
a_sim = [i_sim[np.where(i_obsid == obs)] for obs in a_obsid]
b_sim = [i_sim[np.where(i_obsid == obs)] for obs in b_obsid]

#----Calculate Contamination Absorption----
#unabsorbed/absorbed flux 
#(without contamination models / with contamination models)

#a_cflux_soft = a_cflux_soft.astype(float)
#a_cflux_soft_low = a_cflux_soft_low.astype(float)
#a_cflux_soft_high = a_cflux_soft.astype(float)

a_abs_soft = np.divide(a_cflux_soft,a_flux_soft)
a_abs_soft12 = np.divide(a_cflux_soft12,a_flux_soft12)
a_abs_soft23 = np.divide(a_cflux_soft23,a_flux_soft23)
a_abs_soft13 = np.divide(a_cflux_soft13,a_flux_soft13)
a_abs_hard = np.divide(a_cflux_hard,a_flux_hard)
print type(a_abs_soft)
b_abs_soft = np.divide(b_cflux_soft,b_flux_soft)
b_abs_soft12 = np.divide(b_cflux_soft12,b_flux_soft12)
b_abs_soft23 = np.divide(b_cflux_soft23,b_flux_soft23)
b_abs_soft13 = np.divide(b_cflux_soft13,b_flux_soft13)
b_abs_hard = np.divide(b_cflux_hard,b_flux_hard)

b_abs_soft_low,b_abs_soft_high = ratio_err(b_cflux_soft,b_flux_soft,b_cflux_soft_low,b_cflux_soft_high,b_flux_soft_low,b_flux_soft_high)
b_abs_soft12_low,b_abs_soft12_high = ratio_err(b_cflux_soft12,b_flux_soft12,b_cflux_soft12_low,b_cflux_soft12_high,b_flux_soft12_low,b_flux_soft12_high)
b_abs_soft23_low,b_abs_soft23_high = ratio_err(b_cflux_soft23,b_flux_soft23,b_cflux_soft23_low,b_cflux_soft23_high,b_flux_soft23_low,b_flux_soft23_high)
b_abs_soft13_low,b_abs_soft13_high = ratio_err(b_cflux_soft13,b_flux_soft13,b_cflux_soft13_low,b_cflux_soft13_high,b_flux_soft13_low,b_flux_soft13_high)
b_abs_hard_low,b_abs_hard_high = ratio_err(b_cflux_hard,b_flux_hard,b_cflux_hard_low,b_cflux_hard_high,b_flux_hard_low,b_flux_hard_high)

#----Calculate Optical Depth----
depth_soft = -1.0*np.log(b_abs_soft)
depth_soft_low = depth_soft+1.0*np.log(b_abs_soft_low)
depth_soft_high = -1.0*np.log(b_abs_soft_high)-depth_soft

depth_soft12 = -1.0*np.log(b_abs_soft12)
depth_soft12_low = depth_soft12+1.0*np.log(b_abs_soft12_low)
depth_soft12_high = -1.0*np.log(b_abs_soft12_high)-depth_soft12

depth_soft13 = -1.0*np.log(b_abs_soft13)
depth_soft13_low = depth_soft13+1.0*np.log(b_abs_soft13_low)
depth_soft13_high = -1.0*np.log(b_abs_soft13_high)-depth_soft13

depth_soft23 = -1.0*np.log(b_abs_soft23)
depth_soft23_low = depth_soft23+1.0*np.log(b_abs_soft23_low)
depth_soft23_high = -1.0*np.log(b_abs_soft23_high)-depth_soft23

depth_hard = -1.0*np.log(b_abs_hard)
depth_hard_low = depth_hard+1.0*np.log(b_abs_hard_low)
depth_hard_high = -1.0*np.log(b_abs_hard_high)-depth_hard

#----Calculate Pileup Fractions from 462 results----
ratiosoft = b_flux_soft/b_pflux_soft
ratiohard = b_flux_hard/b_pflux_hard
ratiobroad = b_flux_broad/b_pflux_broad

#----Read Helder2012 Results----
hobsid,inst,hflux_soft,hflux_soft_low,hflux_soft_high,frame,hratiosoft,hflux_hard,hflux_hard_low,hflux_hard_high,hratiohard = np.genfromtxt(helderfile,unpack=True,comments='#')

#--get associated ages--
hage = [get_obs_time(obs) for obs in hobsid]


#----Read CountRates and Calculate Ratios----
#(estimate of contamination absorption)
count_obsid,count_date,count_age,count_grating,count_caldb,count_model,count_em_soft,count_em_soft12,count_em_soft13,count_em_hard,count_meas_soft,count_meas_soft12,count_meas_soft13,count_meas_hard = np.genfromtxt(countratefile,unpack=True,comments='#')

count_ratio_soft = count_meas_soft/count_em_soft
count_ratio_soft12 = count_meas_soft12/count_em_soft12
count_ratio_soft13 = count_meas_soft13/count_em_soft13
count_ratio_hard = count_meas_hard/count_em_hard

#----Read XMM Results----
#note the hard flux is 3-10 keV
xmm_obsid,xmm_date,xmm_age,xmm_filter,xmm_totexp,xmm_filtexp,xmm_flux_soft,xmm_flux_soft_low,xmm_flux_soft_high,xmm_flux_hard,xmm_flux_hard_low,xmm_flux_hard_high = np.genfromtxt(xmmfile,unpack=True,comments='#')


#np.loadtxt(args.infile,unpack=True,skiprows=1)

#fluxes are all in units of 10^-13 ergs/s/cm^2
#errors are 68%
#norms in units of 10^-3
#flux ratios are unpiled/piled

#---------------------------------------
#---------------------------------------
# Extrapolate Contamination Absorption
#---------------------------------------
#---------------------------------------

#paste from excerpt file if needed

#---------------------------------------
#---------------------------------------
#                PLOTS
#---------------------------------------
#---------------------------------------

#---------------------------------------
#        Global Plot Properties
#---------------------------------------

#x-axis (age) range
xmin = 5000
xmax = 10000
ymin = 0
ymax = 120

bott = 0.13
topp = 0.87

#-get jan 1 ages for top axes-
year_ticks = ['2001','2002','2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014']
#use convert_time() to convert jan1,year into age, use as tick locations
regular_ages = [convert_time(yr+'-01-01',get='age') for yr in year_ticks]

if args.clean == 'yes':
    topx = regular_ages
    toplab = year_ticks
    botx = None
    toptitle = 'Year'
    mgrid = False
    frame = False
else:
    topx = b_age
    toplab = b_obsyears
    toptitle = 'Observation Year'
    botx = b_age
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
b_offz_i = np.where(b_sim == '-8.42')
a_offz_i = np.where(a_sim == '-8.42')

#-ACIS grating-
a_hetg_i = np.where(a_grating == 'HETG')
b_hetg_i = np.where(b_grating == 'HETG')
a_bare_i = np.where(a_grating == 'NONE')
b_bare_i = np.where(b_grating == 'NONE')

#--make lists for each--

#-colors = energy band-
softcolor = 'blue'
soft12color = 'green'
soft13color = 'cyan'
soft23color = 'yellow'
hardcolor = 'red'
broadcolor = 'purple'

#-symbol = instrument (ACIS-HETG,ACIS-NONE,XMM-EPIC-pn)-
a_syms = ['o']*a_obsid.shape[0] #HETG = circle
b_syms = ['o']*b_obsid.shape[0]
xmm_syms = 'diamond' #xmm = diamond
a_syms[a_bare_i] = '^' #bare ACIS = triangle
b_syms[b_bare_i] = '^'

#-markeredgecolor = simz (detector z)-
a_mec = ['black']*a_obsid.shape[0]
b_mec = ['black']*b_obsid.shape[0]
xmm_mec = 'black'
a_mec[a_offz_i] = 'red' #offset observations have red outline
b_med[b_offz_i] = 'red'

#-markersize = frametime-
xmm_sizes = 4
a_sizes = [f*2.0 for f in a_frames]
b_sizes = [f*2.0 for f in b_frames]

#-alpha (transparency) = caldb-
b_alpha = 1.0
a_alpha = 0.5
xmm_alpha = [1.0]*xmm_obsid.shape[0]
xmm_alpha[-1] = 0.5 #set last point to different alpha (since reduced myself)

#---------------------------------------
#        Light Curve(s)
#---------------------------------------

#----Set Up Plot----

pdffile = PdfPages(plotfile)
#fig = plt.figure()

#-initialize main plot-
fig,ax1=plt.subplots()
ax1.set_xlim(xmin,xmax)
ax1.set_ylim(ymin,ymax)

#-set axis titles-
plt.xlabel('SN1987A Age [days]')
plt.ylabel('Flux [10$^{-13}$ erg cm$^{-2}$ s$^{-1}$]')

set_axis(ax1,x=botx,grid=mgrid)

#-add extra space on bottom-
plt.subplots_adjust(bottom=bott,top=topp)

#--Soft Flux--
#plot_a(ax1,a_age,a_flux_soft,a_flux_soft_low,a_flux_soft_high,sym='^')
#plot_b(ax1,b_age,b_flux_soft,b_flux_soft_low,b_flux_soft_high,sym='o')
plot_a(ax1,a_age,a_flux_soft,a_flux_soft_low,a_flux_soft_high,syms=a_syms,colors=softcolor,alphas=a_alpha,mecs=a_mecs)
plot_a(ax1,b_age,b_flux_soft,b_flux_soft_low,b_flux_soft_high,syms=b_syms,colors=softcolor,alphas=b_alpha,mecs=b_mecs)

#-plot extrapolated fluxes-
#plt.plot(ext_ages,ext_flux_soft,'*',color='blue',markersize=4,markeredgewidth=0.5)

#-plot xmm fluxes-
plt.errorbar(xmm_age,xmm_flux_soft,yerr=[xmm_flux_soft_low,xmm_flux_soft_high],fmt='o',alpha=0.3,markersize=3)

#--Soft12 Flux (1.0-2.0keV)--
plot_a(ax1,a_age,a_flux_soft12,a_flux_soft12_low,a_flux_soft12_high,sym='^',color='g')
plot_b(ax1,b_age,b_flux_soft12,b_flux_soft12_low,b_flux_soft12_high,sym='o',color='g')

#plt.plot(ext_ages,ext_flux_soft12,'*',color='g',markersize=4,markeredgewidth=0.5)

#--Soft13 Flux (1.0-3.0keV)--
plot_a(ax1,a_age,a_flux_soft13,a_flux_soft13_low,a_flux_soft13_high,sym='^',color='cyan')
plot_b(ax1,b_age,b_flux_soft13,b_flux_soft13_low,b_flux_soft13_high,sym='o',color='cyan')

#plt.plot(ext_ages,ext_flux_soft13,'*',color='cyan',markersize=4,markeredgewidth=0.5)

#--Hard Flux--
plot_a(ax1,a_age,a_flux_hard,a_flux_hard_low,a_flux_hard_high,sym='^',color='r')
plot_b(ax1,b_age,b_flux_hard,b_flux_hard_low,b_flux_hard_high,sym='o',color='r')

#-plot xmm fluxes-
plt.errorbar(xmm_age,xmm_flux_hard,yerr=[xmm_flux_hard_low,xmm_flux_hard_high],fmt='o',alpha=0.3,markersize=3,color='r')

#--Broad Flux--
plot_a(ax1,a_age,a_flux_broad,a_flux_broad_low,a_flux_broad_high,sym='^',color='purple')
plot_b(ax1,b_age,b_flux_broad,b_flux_broad_low,b_flux_broad_high,sym='o',color='purple')

#plt.scatter(ages,r0_arcsec,marker='o')#,s=deconcounts/100)
#--Plot Helder2012 Measurements--

if args.plothelder == 'yes':
    plt.errorbar(hage,hflux_soft,yerr=[hflux_soft_low,hflux_soft_high],fmt='o',markersize=2.0,color='green',label='Helder2012')

    plt.errorbar(hage,hflux_hard,yerr=[hflux_hard_low,hflux_hard_high],fmt='d',markersize=2.0,color='green',label='_nolegend')

    
#plt.legend(numpoints=1,loc='upper left',fontsize='small',frameon=True)

#--Add Legend--
emptyx = [0]
emptyy = [0]
plt.scatter(emptyx,emptyy,color='black',marker='^',label='CALDB 4.5.9')
plt.scatter(emptyx,emptyy,color='black',marker='o',label='CALDB 4.6.2')
plt.scatter(emptyx,emptyy,color='purple',marker='s',label='0.5-8.0 keV')
plt.scatter(emptyx,emptyy,color='blue',marker='s',label='0.5-2.0 keV')
plt.scatter(emptyx,emptyy,color='green',marker='s',label='1.0-2.0 keV')
plt.scatter(emptyx,emptyy,color='cyan',marker='s',label='1.0-3.0 keV')
plt.scatter(emptyx,emptyy,color='red',marker='s',label='3.0-8.0 keV')
plt.scatter(emptyx,emptyy,color='gray',marker='o',label='EPIC-pn',alpha=0.3)

plt.legend(loc='upper left',fontsize='xx-small',frameon=frame,scatterpoints=1,numpoints=1)

#--Add Top (year) x-axis--
#note this must come after plotting things on the first axis!

set_axis(ax1,x=topx,twin=True,title=toptitle,xlab=toplab)#,xlab=b_obsyears)

#ax2 = ax1.twiny()
#ax2.set_xlim(xmin,xmax)
#ax2.set_xticks(b_age,minor=False)
#ax2.set_xticklabels(b_obsyears)
#ax2.set_xlabel("Observation Year")
#for tick in ax2.xaxis.get_major_ticks():
#    tick.label2.set_fontsize(10)
#    tick.label2.set_rotation('vertical')

#--Close Plot--
pdffile.savefig()
plt.close()

#---------------------------------------
#          Pileup Fractions
#---------------------------------------

#----Initialize Plot----
fig,ax1=plt.subplots()

plt.xlabel('SN1987A Age [days]')
plt.ylabel('Unpiled/Piled Flux')

set_axis(ax1,x=botx,grid=mgrid)

plt.subplots_adjust(bottom=bott,top=topp)

#-set limits and axis titles-
plt.xlim(xmin,xmax)
#plt.ylim(0.2,0.6)

#--plot new fractions--
plt.plot(b_age,ratiosoft,marker='o',color='blue')
plt.plot(b_age,ratiohard,marker='d',color='blue')
plt.plot(b_age,ratiobroad,marker='s',color='blue')

#--plot Helder2012 fractions--
plt.plot(hage,hratiosoft,marker='o',color='green',label='_nolegend_')
plt.plot(hage,hratiohard,marker='d',color='green',label='_nolegend_')
#plt.plot(hage,hratiobroad,marker='s',color='green')

plt.plot([0],[0.4],markersize=0,label='Helder2012',color='green')
plt.plot([0],[0.4],markersize=0,label='CALDB 4.6.2',color='blue')
plt.plot([0],[0.4],label='0.5-2.0 keV',color='black',marker='o')
plt.plot([0],[0.4],label='0.5-8.0 keV',color='black',marker='s')
plt.plot([0],[0.4],label='3.0-8.0 keV',color='black',marker='d')

plt.legend(numpoints=1,frameon=frame,fontsize='xx-small')

#set_axis(ax1,b_age,twin=True,title='Observation Year',xlab=b_obsyears)
set_axis(ax1,x=topx,twin=True,title=toptitle,xlab=toplab)

#--close plot--
pdffile.savefig()
plt.close()

#---------------------------------------
#          Hardness Ratios
#---------------------------------------
#(hard/soft flux)

fig,ax1=plt.subplots()

plt.xlabel('SN1987A Age [days]')
plt.ylabel('Hard/Soft Flux')

set_axis(ax1,x=botx,grid=mgrid)

#--set tick labels and grid lines--
#ax1.set_xticks(b_age,minor=False)
#ax1.xaxis.grid(True,which='major')
#for tick in ax1.xaxis.get_major_ticks():
#    tick.label.set_fontsize(10)
#    tick.label.set_rotation('vertical')

plt.subplots_adjust(bottom=bott,top=topp)

#--set limits and axis titles--
plt.xlim(xmin,xmax)

#--plot new ratios--
plt.plot(b_age,b_flux_hard/b_flux_soft,marker='o',color='blue',label='Soft=0.5-2.0 keV')

#--plot hard/soft2--
plt.plot(b_age,b_flux_hard/b_flux_soft12,marker='^',color='blue',label='Soft=1.0-2.0 keV')

#--plot Helder2012 fractions--
plt.plot(hage,hflux_hard/hflux_soft,marker='o',color='green',label='Soft=0.5-2.0 keV, Helder2012')

#--plot legend--
plt.scatter([0],[0.05],color='black',marker=None,label='Hard=3.0-8.0 keV')
plt.legend(numpoints=1,frameon=frame,fontsize='xx-small')

set_axis(ax1,x=topx,twin=True,title=toptitle,xlab=toplab)

#--close plot--
pdffile.savefig()
plt.close()

#---------------------------------------
#       Contamination Absorption
#---------------------------------------

#----Set Up Plot----

#pdffile = PdfPages(plotfile)
#fig = plt.figure()

#-initialize main plot-
fig,ax1=plt.subplots()
ax1.set_xlim(xmin,xmax)
#ax1.set_ylim(0.99,1.01)

#-set axis titles-
plt.xlabel('SN1987A Age [days]')
plt.ylabel('Absorbed Flux / Unabsorbed Flux')

set_axis(ax1,x=botx,grid=mgrid)

#-add extra space on bottom-
plt.subplots_adjust(bottom=bott,top=topp)

#--legend--
#4.5.9 = triangle ('^')
#4.6.2 = circle ('o')
#blue = 0.5-2.0 keV
#magenta = 2.0-3.0 keV
#cyan = 1.0 - 3.0 keV
#green = 1.0-2.0 keV
#red = 3.0-8.0 keV

psize = 5
al = 1
er_al = 0.1

#-plot with regular errobars-
#plt.errorbar(b_age,b_abs_hard,yerr=[b_abs_hard_low,b_abs_hard_high],marker='o',color='red',label='3.0-8.0 keV',markersize=psize,alpha=al)
#plt.errorbar(b_age,b_abs_soft23,yerr=[b_abs_soft23_low,b_abs_soft23_high],marker='o',color='magenta',label='2.0-3.0 keV',markersize=psize,alpha=al)
#plt.errorbar(b_age,b_abs_soft13,yerr=[b_abs_soft13_low,b_abs_soft13_high],marker='o',color='cyan',label='1.0-3.0 keV',markersize=psize,alpha=al)
#plt.errorbar(b_age,b_abs_soft12,yerr=[b_abs_soft12_low,b_abs_soft12_high],marker='o',color='green',label='1.0-2.0 keV',markersize=psize,alpha=al)
#plt.errorbar(b_age,b_abs_soft,yerr=[b_abs_soft_low,b_abs_soft_high],marker='o',color='blue',label='0.5-2.0 keV',markersize=psize,alpha=al)

#-plot with errors as shaded region
plt.plot(b_age,b_abs_hard,marker='o',color='red',label='3.0-8.0 keV',markersize=psize,alpha=al)
plt.plot(b_age,b_abs_soft23,marker='o',color='magenta',label='2.0-3.0 keV',markersize=psize,alpha=al)
plt.plot(b_age,b_abs_soft13,marker='o',color='cyan',label='1.0-3.0 keV',markersize=psize,alpha=al)
plt.plot(b_age,b_abs_soft12,marker='o',color='green',label='1.0-2.0 keV',markersize=psize,alpha=al)
plt.plot(b_age,b_abs_soft,marker='o',color='blue',label='0.5-2.0 keV',markersize=psize,alpha=al)

plt.fill_between(b_age,b_abs_soft-b_abs_soft_low,b_abs_soft+b_abs_soft_high,alpha=er_al,color='blue')
plt.fill_between(b_age,b_abs_soft23-b_abs_soft23_low,b_abs_soft23+b_abs_soft23_high,alpha=er_al,color='magenta')
plt.fill_between(b_age,b_abs_soft13-b_abs_soft13_low,b_abs_soft13+b_abs_soft13_high,alpha=er_al,color='cyan')
plt.fill_between(b_age,b_abs_soft12-b_abs_soft12_low,b_abs_soft12+b_abs_soft12_high,alpha=er_al,color='green')
plt.fill_between(b_age,b_abs_hard-b_abs_hard_low,b_abs_hard+b_abs_hard_high,alpha=er_al,color='red')

#--plot ratios from count rates--
y_min,y_max = plt.ylim()
plt.plot(count_age,count_ratio_soft,'^',color='blue')
plt.plot(count_age,count_ratio_soft12,'^',color='green')
plt.plot(count_age,count_ratio_soft13,'^',color='cyan')
plt.plot(count_age,count_ratio_hard,'^',color='red')
plt.scatter(emptyx,[y_min],color='black',marker='^',label='Count Rate Ratios')

#--plot legend--
plt.legend(numpoints=1,frameon=frame,fontsize='xx-small',loc='upper left',scatterpoints=1)

#plt.plot(a_age,a_abs_hard+1.0,marker='o',color='red',markersize=psize)
#plt.plot(a_age,a_abs_soft23+1.0,marker='o',color='magenta',markersize=psize)
#plt.plot(a_age,a_abs_soft13+1.0,marker='o',color='cyan',markersize=psize)
#plt.plot(a_age,a_abs_soft12+1.0,marker='o',color='green',markersize=psize)
#plt.plot(a_age,a_abs_soft+1.0,marker='o',color='blue',markersize=psize)

#--plot fit results--
if args.plotfit == 'yes':
    if args.dofit == 'no':
        #read from file
        slope,intercept = np.loadtxt(mcmcfile,unpack=True,skiprows=1)
        a = slope[0]
        b = intercept[0]
        mcmc_x = range(10000)
        mcmc_y = [a*x+b for x in mcmc_x]
    else: #-generate fitted data to overplot on contamination curve---
        medpars = np.array([0,soft_stats[0,0],soft_stats[0,1],soft_stats[0,2],soft_stats[0,3]])        
        mcmc_x = np.array(range(10000))
        mcmc_y = fitting.arr_linear(mcmc_x,medpars)

#    plt.plot(mcmc_x,mcmc_y,color='blue')
    plt.plot(poly_x,soft_poly_y,linestyle='--',color='blue')
    plt.plot(poly_x,soft12_poly_y,linestyle='--',color='green')
    plt.plot(poly_x,soft13_poly_y,linestyle='--',color='cyan')

#--plot top axis--
set_axis(ax1,x=topx,twin=True,title=toptitle,xlab=toplab)

#--close plot--
pdffile.savefig()
plt.close()

#---------------------------------------
#       Optical Depth
#---------------------------------------

#----Set Up Plot----

#pdffile = PdfPages(plotfile)
#fig = plt.figure()

#-initialize main plot-
fig,ax1=plt.subplots()
ax1.set_xlim(xmin,xmax)
#ax1.set_ylim(0.99,1.01)

#-set axis titles-
plt.xlabel('SN1987A Age [days]')
plt.ylabel('Optical Depth')

set_axis(ax1,x=botx,grid=mgrid)

#-add extra space on bottom-
plt.subplots_adjust(bottom=bott,top=topp)

#--legend--
#4.5.9 = triangle ('^')
#4.6.2 = circle ('o')
#blue = 0.5-2.0 keV
#magenta = 2.0-3.0 keV
#cyan = 1.0 - 3.0 keV
#green = 1.0-2.0 keV
#red = 3.0-8.0 keV

psize = 5
al = 1
er_al = 0.1

#-plot with regular errobars-
#plt.errorbar(b_age,b_abs_hard,yerr=[b_abs_hard_low,b_abs_hard_high],marker='o',color='red',label='3.0-8.0 keV',markersize=psize,alpha=al)
#plt.errorbar(b_age,b_abs_soft23,yerr=[b_abs_soft23_low,b_abs_soft23_high],marker='o',color='magenta',label='2.0-3.0 keV',markersize=psize,alpha=al)
#plt.errorbar(b_age,b_abs_soft13,yerr=[b_abs_soft13_low,b_abs_soft13_high],marker='o',color='cyan',label='1.0-3.0 keV',markersize=psize,alpha=al)
#plt.errorbar(b_age,b_abs_soft12,yerr=[b_abs_soft12_low,b_abs_soft12_high],marker='o',color='green',label='1.0-2.0 keV',markersize=psize,alpha=al)
#plt.errorbar(b_age,b_abs_soft,yerr=[b_abs_soft_low,b_abs_soft_high],marker='o',color='blue',label='0.5-2.0 keV',markersize=psize,alpha=al)

#-plot with errors as shaded region
plt.plot(b_age,depth_hard,marker='o',color='red',label='3.0-8.0 keV',markersize=psize,alpha=al)
plt.plot(b_age,depth_soft23,marker='o',color='magenta',label='2.0-3.0 keV',markersize=psize,alpha=al)
plt.plot(b_age,depth_soft13,marker='o',color='cyan',label='1.0-3.0 keV',markersize=psize,alpha=al)
plt.plot(b_age,depth_soft12,marker='o',color='green',label='1.0-2.0 keV',markersize=psize,alpha=al)
plt.plot(b_age,depth_soft,marker='o',color='blue',label='0.5-2.0 keV',markersize=psize,alpha=al)

plt.fill_between(b_age,depth_soft-depth_soft_low,depth_soft+depth_soft_high,alpha=er_al,color='blue')
plt.fill_between(b_age,depth_soft23-depth_soft23_low,depth_soft23+depth_soft23_high,alpha=er_al,color='magenta')
plt.fill_between(b_age,depth_soft13-depth_soft13_low,depth_soft13+depth_soft13_high,alpha=er_al,color='cyan')
plt.fill_between(b_age,depth_soft12-depth_soft12_low,depth_soft12+depth_soft12_high,alpha=er_al,color='green')
plt.fill_between(b_age,depth_hard-depth_hard_low,depth_hard+depth_hard_high,alpha=er_al,color='red')

#--plot ratios from count rates--
y_min,y_max = plt.ylim()
#plt.plot(count_age,count_ratio_soft,'^',color='blue')
#plt.plot(count_age,count_ratio_soft12,'^',color='green')
#plt.plot(count_age,count_ratio_soft13,'^',color='cyan')
#plt.plot(count_age,count_ratio_hard,'^',color='red')
#plt.scatter(emptyx,[y_min],color='black',marker='^',label='Count Rate Ratios')

#--plot legend--
plt.legend(numpoints=1,frameon=frame,fontsize='xx-small',loc='upper left',scatterpoints=1)

#plt.plot(a_age,a_abs_hard+1.0,marker='o',color='red',markersize=psize)
#plt.plot(a_age,a_abs_soft23+1.0,marker='o',color='magenta',markersize=psize)
#plt.plot(a_age,a_abs_soft13+1.0,marker='o',color='cyan',markersize=psize)
#plt.plot(a_age,a_abs_soft12+1.0,marker='o',color='green',markersize=psize)
#plt.plot(a_age,a_abs_soft+1.0,marker='o',color='blue',markersize=psize)

#--plot top axis--
set_axis(ax1,x=topx,twin=True,title=toptitle,xlab=toplab)

#--close plot--
pdffile.savefig()
plt.close()


#---------------------------------------
#     Hot/Cool Normalization Ratios
#---------------------------------------

#--initialize main plot--
fig,ax1=plt.subplots()
ax1.set_xlim(xmin,xmax)
#ax1.set_ylim(0.99,1.01)

#-set axis titles-
plt.xlabel('SN1987A Age [days]')
plt.ylabel('Norm$_{Hot}$/Norm$_{Cool}$')

#-set up axis-
set_axis(ax1,x=botx,grid=mgrid)

#-add extra space on bottom-
plt.subplots_adjust(bottom=bott,top=topp)

#--calculate ratio--
normratio = b_normhot/b_normcool
normratio_low,normratio_high = ratio_err(b_normhot,b_normcool,b_normhot_low,b_normhot_high,b_normcool_low,b_normcool_high)

#--plot ratios--
ax1.errorbar(b_age,normratio,yerr=[normratio_low,normratio_high],marker='o')

#--plot legend--
#plt.legend(numpoints=1,frameon=frame,fontsize='xx-small',loc='upper left')

#--set top axis--
set_axis(ax1,x=topx,twin=True,title=toptitle,xlab=toplab)

#--save and close plot--
pdffile.savefig()
plt.close()

#---------------------------------------
#       kT (hot and cold) vs Age
#---------------------------------------

#--initialize main plot--
fig,ax1=plt.subplots()
ax1.set_xlim(xmin,xmax)
#ax1.set_ylim(0.99,1.01)

#-set axis titles-
plt.xlabel('SN1987A Age [days]')
plt.ylabel('kT [keV]')

#-set up axis-
set_axis(ax1,x=botx,grid=mgrid)

#-add extra space on bottom-
plt.subplots_adjust(bottom=bott,top=topp)

#--plot temperatures--
#legend
#blue = cool
#red = hot
ax1.errorbar(b_age,b_kTcool,yerr=[b_kTcool-b_kTcool_low,b_kTcool_high-b_kTcool],marker='o',color='blue',label='Cool Component')
ax1.errorbar(b_age,b_kThot,yerr=[b_kThot-b_kThot_low,b_kThot_high-b_kThot],marker='o',color='red',label='Hot Component')

#--plot legend--
plt.legend(numpoints=1,frameon=frame,fontsize='xx-small',loc='upper right')

#--set top axis--
set_axis(ax1,x=topx,twin=True,title=toptitle,xlab=toplab)

#--save and close plot--
pdffile.savefig()
plt.close()

#---------------------------------------
#     Hot/Cool kT Ratios
#---------------------------------------

#--initialize main plot--
fig,ax1=plt.subplots()
ax1.set_xlim(xmin,xmax)
#ax1.set_ylim(0.99,1.01)

#-set axis titles-
plt.xlabel('SN1987A Age [days]')
plt.ylabel('kT$_{Hot}$/kT$_{Cool}$')

#-set up axis-
set_axis(ax1,x=botx,grid=mgrid)

#-add extra space on bottom-
plt.subplots_adjust(bottom=bott,top=topp)

#--calculate ratio--
kTratio = b_kThot/b_kTcool
kTratio_low,kTratio_high = ratio_err(b_kThot,b_kTcool,b_kThot_low,b_kThot_high,b_kTcool_low,b_kTcool_high)

#--plot ratios--
ax1.errorbar(b_age,kTratio,yerr=[kTratio_low,kTratio_high],marker='o')

#--plot legend--
#plt.legend(numpoints=1,frameon=frame,fontsize='xx-small',loc='upper left')

#--set top axis--
set_axis(ax1,x=topx,twin=True,title=toptitle,xlab=toplab)

#--save and close plot--
pdffile.savefig()
plt.close()


#---------------------------------------
#     Ionization Age vs Age
#---------------------------------------

#--initialize main plot--
fig,ax1=plt.subplots()
ax1.set_xlim(xmin,xmax)
#ax1.set_ylim(0.99,1.01)

#-set axis titles-
plt.xlabel('SN1987A Age [days]')
plt.ylabel('Ionization Age [10$^{11}$ s cm$^{-3}$]')

#-set up axis-
set_axis(ax1,x=botx,grid=mgrid)

#-add extra space on bottom-
plt.subplots_adjust(bottom=bott,top=topp)

#--plot tau--
ax1.errorbar(b_age,b_tauhot,yerr=[b_tauhot-b_tauhot_low,b_tauhot_high-b_tauhot],marker='o',color='red')

#--set top axis--
set_axis(ax1,x=topx,twin=True,title=toptitle,xlab=toplab)

#--plot legend--
#plt.legend(numpoints=1,frameon=frame,fontsize='xx-small',loc='upper left')

#--save and close plot--
pdffile.savefig()
plt.close()

#---------------------------------------
#             Wrap Up
#---------------------------------------

#----Close File----
pdffile.close()


