#! /usr/bin/env python

###! /usr/astro/ciao-4.6/bin/python

#Author: Kari A. Frank
#Date: 2015-04-03
#Purpose: Plot results from the radius fitting analysis.
#Usage: $PYDIR/sn1987a/radius_results.py [--radiusfile radiusfile] [--agemin agemin] [--agemax agemax] [--bands bands] [--mincounts mincounts] [--ignorechipy ignorechipy] [--ignoregrating ignoregrating] [--clean clean] [--bandcolors bandcolors] [--plotfile plotfile] [--plotfit plotfit] [--clobber clobber]
#
#Input:
#
# radiusfile: text file containing radius fit results, of same format created
#             by get_radii_results.py (default='radii_fits.txt')
#
# plotfile:   name of output plot file (default=radiusfile+'.pdf', e.g.
#             'radii_fits.pdf')
#
# agemin/agemax: minimum/maximum ages to include (default is unrestricted)
#
# bands:      specify which bands to include in eV, as string separated by
#             by ',', e.g. '300-800,800-1200,1200-8000'.  Also accepts
#             'soft','broad','hard', and 'all', defined as follows. 
#             Default is 'all'.
#                 'soft' = '500-2000'
#                 'hard' = '3000-8000'
#                 'broad'= '300-8000'
#                 'all'  = '300-800,800-1200,1200-8000,300-8000,3000-8000,
#                              2000-10000'
#
# mincounts:  optional minimum required counts to include (default=0)
#
# clean:      optional switch to change the axes labels and ticks to even
#             years and ages (default = 'yes')
#
# ignorechipy: optional switch to not plot gray (or thick) borders
#                for observations with offset chip positions. default='yes' 
#
# ignoregrating: do not use different symbols for observations with gratings
#                (default='no')
#
# plotfit:    optionally overplot the best fit broken-linear function, as output 
#             by expansion_fit.pro (default='no'). must have already run expansion_fit.pro
#             for the relevant band(s), and output files must be in same directory.
#
# plotextras: optionally plot information from other wavelengths: optical hot spots,
#             Plait1995 ring size/location, Ng2013 radio break (default=None)
#             options are: 'plait1995','ng2013','larsson2013','hotspots','sugarman2002','all','ng2013radii'
#
# clobber:    specifies whether files should be overwritten
#             if they already exist (same as in ciao tools,
#             default = 'no')
#
#Output:
# - output1
#
#Usage Notes:
# - note1

#---------------------------------------
#           Import Modules
#---------------------------------------
import argparse,os,sys
#import ciao_contrib.runtool as crt #all functions should be prefixed with crt
import home_grown as hg #import my own custom functions
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as dates
import numpy as np
from numpy.lib.recfunctions import append_fields
from sn1987a_time import *
from fancy_plot import fancy_plot
sys.path.insert(0,'/astro/research/kaf33/Dropbox/Research/Python_Programs/xmc_analysis/')
#sys.path.insert(0,'/Users/kafrank/Dropbox/Research/Python_Programs/xmc_analysis/')
from xmcplot import plot_param_hist
import pandas as pd
#from sn1987a_plot_helpers import plot_extras

#---------------------------------------
#---------------------------------------
#      DEFINE HELPER FUNCTIONS
#---------------------------------------
#---------------------------------------

#---------------------------------------
#         Set Plot Axis
#---------------------------------------
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

#---------------------------------------
#          Initialize Plot
#---------------------------------------
def start_plot(xtitle='',ytitle='',xmin=None,xmax=None,ymin=None,ymax=None,ylog=False):

    #-initialize main plot-
    fig,ax1=plt.subplots()

    #-set axis titles-
    plt.xlabel(xtitle)
    plt.ylabel(ytitle)

    #-set bottom axis-
    if ylog == True: plt.yscale('log')
    set_axis(ax1,x=botx,grid=mgrid)
    if xmin is not None and xmax is not None:
        ax1.set_xlim(xmin,xmax)
    if ymin is not None and ymax is not None:
        ax1.set_ylim(ymin,ymax)

    #-add extra space on bottom-
    plt.subplots_adjust(bottom=bott,top=topp)

    return ax1

#---------------------------------------
#          Plot Legend
#---------------------------------------
def plot_legend(labels=None,colors=['black'],markers=['o'],location='upper left',
                fontsize='medium',frameon=False,scatterpoints=1,numpoints=1,
                markerscale=1.5,mews=None,mecs=None):

    #-convert defaults to lists-
    if mews is None: mews = [None]*len(labels)
    if mecs is None: mecs = [None]*len(labels)

    #-dummy position arrays-
    emptyx=[0]
    emptyy=emptyx

    #-plot dummy points to define labels-
    for p in range(len(labels)):
        plt.scatter(emptyx,emptyy,color=colors[p],marker=markers[p],label=labels[p],edgecolor=mecs[p],linewidth=mews[p])

    #-plot legend-
    leg = plt.legend(loc=location,fontsize=fontsize,frameon=frameon,
                     scatterpoints=scatterpoints,numpoints=numpoints,markerscale=markerscale)


#---------------------------------------
#---------------------------------------
#             MAIN PROGRAM
#---------------------------------------
#---------------------------------------

#---------------------------------------
#   Parse arguments and set defaults
#---------------------------------------

parser = argparse.ArgumentParser(description='Make expansion and other plots from results of SN1987A image fitting.',epilog='')

pwd = os.getcwd()

#parser.add_argument('required_arg',help='Description and usage for required_arg.',default='default for required_arg')

parser.add_argument('--radiusfile',help='Text file with image fitting results.',default='radii_fits.txt')

parser.add_argument('--plotfile',help='Output plot file.',default='default')

parser.add_argument('--bands',help='Bands to include in plotting and fitting.',default='all')

parser.add_argument('--agemin',help='Minimum age to include (in days since SN).',default=0)
parser.add_argument('--agemax',help='Maximum age to include (in days since SN).',default=None)

parser.add_argument('--mincounts',help='Include only measurements from images iwth mincounts counts.',default=0)

parser.add_argument('--clean',help='Determines x-axis ticks and labels.',default='yes')

parser.add_argument('--bandcolors',help='List of color names to override default band color choices.',default=None)

parser.add_argument('--ignorechipy',help='Plot normal borders for observations with offset chip positions.',default='yes')
parser.add_argument('--ignoregrating',help='Do not use different symbols for observations with gratings.',default='no')

parser.add_argument('--plotfit',help='Overplot the best-fit broken-linear functions for each band.',default='no')
parser.add_argument('--plotextras',help='Overplot information from other wavelengths.',default=None)

parser.add_argument('--clobber',help='Overwrite existing files.',default='no')

args = parser.parse_args()

#----Convert default types----
args.mincounts = int(args.mincounts)
args.agemin = int(args.agemin)
if args.agemax != None: args.agemax = int(args.agemax)
    
#----Parse bands----
if args.bands=='all': args.bands='300-800,800-1200,1200-8000,300-8000,3000-8000,2000-10000,500-2000,300-1200'
if args.bands=='broad': args.bands='300-8000'
if args.bands=='hard': args.bands='3000-8000'
if args.bands=='soft': args.bands='500-2000'
bands=args.bands.split(',')
print "Bands = ",bands

bandproperties = pd.DataFrame([['red','0.3-0.8 keV'],['purple','0.8-1.2 keV'],['magenta','0.3-1.2 keV'],['blue','0.5-2 keV'],['cyan','1.2-8 keV'],['green','3-8 keV'],['yellow','2-10 keV'],['black','0.3-8 keV']],index=['300-800','800-1200','300-1200','500-2000','1200-8000','3000-8000','2000-10000','300-8000'],columns=['color','label'])

# remove any duplicate bands
bands = list(set(bands))

# check for invalid bands
for b in bands:
    if b not in bandproperties.index:
        print "Warning: "+b+" is not a valid band. Skipping "+b+"."
        bands.remove(b)

# for later conversion using pandas
#bandcolors = bandproperties['color'][bands]
#bandlabels = bandproperties['label'][bands]
        
softcolor = 'red' #300-800
mediumcolor = 'green' #800-1200
softmediumcolor = 'yellow' #300-1200
hardcolor = 'blue' #1200-8000
veryhardcolor = 'cyan' #3000-8000
broadcolor = 'black' #300-8000
widehardcolor = 'blue' #'purple' #2000-10000
widesoftcolor = 'green'#'magenta' #500-2000
bandcolors = ['']*len(bands)
bandlabels = ['']*len(bands)
for bi in range(len(bands)):
    if bands[bi] == '300-800':
        bandlabels[bi] = '0.3-0.8 keV'
    elif bands[bi] == '800-1200':
        bandlabels[bi] = '0.8-1.2 keV'
    elif bands[bi] == '1200-8000':
        bandlabels[bi] = '1.2-8.0 keV'
    elif bands[bi] == '300-8000':
        bandlabels[bi] = '0.3-8.0 keV' 
    elif bands[bi] == '3000-8000':
        bandlabels[bi] = '3.0-8.0 keV'
    elif bands[bi] == '2000-10000':
        bandlabels[bi] = '2.0-10.0 keV'
    elif bands[bi] == '300-1200':
        bandlabels[bi] = '0.3-1.2 keV'
    elif bands[bi] == '500-2000':
        bandlabels[bi] = '0.5-2.0 keV'
    else:
        print "WARNING: unrecognized band."

if args.bandcolors == None:
    for bi in range(len(bands)):
        if bands[bi] == '300-800':
            bandcolors[bi] = softcolor
        elif bands[bi] == '800-1200':
            bandcolors[bi] = mediumcolor
        elif bands[bi] == '1200-8000':
            bandcolors[bi] = hardcolor
        elif bands[bi] == '300-8000':
            bandcolors[bi] = broadcolor
        elif bands[bi] == '3000-8000':
            bandcolors[bi] = veryhardcolor
        elif bands[bi] == '2000-10000':
            bandcolors[bi] = widehardcolor
        elif bands[bi] == '300-1200':
            bandcolors[bi] = softmediumcolor
        elif bands[bi] == '500-2000':
            bandcolors[bi] = widesoftcolor
        else:
            print "WARNING: unrecognized band."
    

#----Set file paths----

sn1987a_dir = '/data/nasdrive/main/kaf33/sn1987a_monitoring/'
#sn1987a_dir = '/Users/kafrank/Research/SN1987A/'
#sn1987a_dir = '/media/backupdisk/kaf33/main/Chandra_Observations/SN1987A/'
obsfile = sn1987a_dir+'comparison_dir/chandra_observations.txt'
radiusfilename = os.path.splitext(args.radiusfile)[0]
if args.plotfile == 'default':
    args.plotfile = radiusfilename+'.pdf'

#---------------------------------------
#         Read In Data
#---------------------------------------

#----Read Radius Fit Results----

#-into array with named columns----

fitdata = np.genfromtxt(args.radiusfile,skip_header=0,names=True,comments='#',dtype=None)

#----Get Observation Information----

#--Get Observation Years--
obsyears = [str(get_obs_time(int(obs),get='year')) for obs in fitdata['obsid']]
fitdata = append_fields(fitdata,names='year',data=obsyears)

#--Read Main Info File--
obsinfo = np.genfromtxt(obsfile,skip_header=0,names=True,comments='#',dtype=None)
#column names:
#obsid	date	age	pi	configuration	grating	exposure	frametime	simoffset	yoffset	zoffset

#-get years-
obsyears = [str(get_obs_time(int(obs),get='year')) for obs in obsinfo['obsid']]
obsinfo = append_fields(obsinfo,names='year',data=obsyears)

#--Get FrameTimes--
frames = [obsinfo['frametime'][np.where(obsinfo['obsid'] == obs)] for obs in fitdata['obsid']]
fitdata = append_fields(fitdata,names='frametime',data=np.array(frames)[:,0])

#--Get SimZ--

#-convert 'default' to '0.0'-
obsinfo['simoffset'][np.where(obsinfo['simoffset'] == 'default')] = '0.0'
#print i_sim,i_obsid

#-associate with a,b,r data-
sims = [obsinfo['simoffset'][np.where(obsinfo['obsid'] == obs)] for obs in fitdata['obsid']]
fitdata = append_fields(fitdata,names='sim',data=np.array(sims)[:,0])

#--get gratings--
gratings = [obsinfo['grating'][np.where(obsinfo['obsid'] == obs)] for obs in fitdata['obsid']]
fitdata = append_fields(fitdata,names='grating',data=np.array(gratings)[:,0])

#---------------------------------------
#     Filter Out Unwanted Rows
#---------------------------------------

#----Filter on bands----
fitdata=fitdata[np.logical_or.reduce([fitdata['band']==b for b in bands])]
#print fitdata

#----Filter on counts----
fitdata=fitdata[np.where(fitdata['deconcounts'] >= args.mincounts)]

#----Filter on age----
fitdata=fitdata[np.where(fitdata['age'] >= args.agemin)]
if args.agemax != None: fitdata=fitdata[np.where(fitdata['age'] <= args.agemax)]

#---------------------------------------
# MAKE FUNCTION Fit Radius Expansion to Broken Powerlaw
#---------------------------------------

# Runs expansion_pro.fit for each band:
# IDL> expansion_fit,infile=args.radiusfile,band='300-8000'
# (change band string as appropriate)

## currently this runs, but the arguments are not passed (expansion_fit.pro just uses defaults) 

#if args.plotfit == 'yes':
#    for b in bands:

        # Build shell command string
        # command line: "idl -e expansion_fit -args infile=<>,band=<>,mincounts=<>"
#        cmd = "idl -quiet -e expansion_fit -args infile='"+args.radiusfile+"',band='"+b+"',mincounts="+str(args.mincounts)

        # Execute commad line
#        print cmd
#        os.system(cmd)
    
    
#---------------------------------------
#       Plot Data vs Age
#---------------------------------------

#--Set whitespace on bottom, top-- 
bott = 0.13
topp = 0.87

#--Initialize Plot and Set Axes--
if args.agemax == None: args.agemax = 200+np.max(fitdata['age'])
if args.agemin == 0: args.agemin = np.min(fitdata['age'])-200
rmin=0.54
rmax=0.87
#rmin=np.min(fitdata['R0'])-np.max(fitdata['R0errl'])
#rmax=np.max(fitdata['R0erru'])+np.max(fitdata['R0'])
sigrmin=np.min(fitdata['SIGR'])-np.max(fitdata['SIGRerrl'])
sigrmax=np.max(fitdata['SIGRerru'])+np.max(fitdata['SIGR'])
nsmin=np.min(fitdata['SKYNS'])-np.max(fitdata['SKYNSerrl'])
nsmax=np.max(fitdata['SKYNSerru'])+np.max(fitdata['SKYNS'])

#if args.ylog == True:

#----Set Axis Labels and Ticks----

#-calculate required top axis labels automatically based on agemin/agemax-
print 'agemin = ',args.agemin
years = range(convert_time(args.agemin,get='year')+1,convert_time(args.agemax,get='year')+1)
year_ticks = [str(yr) for yr in years]
regular_ages = [convert_time(yr+'-01-01',get='age',informat='date') for yr in year_ticks]
#print 'yeartick = ',year_ticks
#print 'regular ages = ',regular_ages

#radyears = range(1999,2016)
#radius_year_ticks = [str(yr) for yr in radyears]
#radius_regular_ages = [convert_time(yr+'-01-01',get='age') for yr in radius_year_ticks]

if args.clean == 'yes':
    topx = regular_ages
    toplab = year_ticks
    botx = None
    toptitle = 'Year'
    mgrid = False
    frame = False
else:
    topx = obsinfo['age']#a_data['age']
    toplab = obsinfo['year']#a_data['year']
    toptitle = 'Observation Year'
    botx = obsinfo['age']#a_data['age']
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
offz_i = np.where(fitdata['sim'] == '-8.42')

#-ACIS grating-
hetg_i = np.where(fitdata['grating'] == 'HETG')
bare_i = np.where(fitdata['grating'] == 'NONE')
letg_i = np.where(fitdata['grating'] == 'LETG')

#-band-
soft_i   = np.where(fitdata['band'] == '300-800')
medium_i = np.where(fitdata['band'] == '800-1200')
hard_i   = np.where(fitdata['band'] == '1200-8000')
broad_i  = np.where(fitdata['band'] == '300-8000')
veryhard_i = np.where(fitdata['band'] == '3000-8000')
widehard_i = np.where(fitdata['band'] == '2000-10000')
softmedium_i = np.where(fitdata['band'] == '300-1200')
widesoft_i = np.where(fitdata['band'] == '500-2000')

#-initialize legend lists-
legcolors=[]
leglabels=[]
legmews=[]
legmecs=[]
legsyms=[]

#--make lists for each--

#-number of observations-
nobs = fitdata.shape[0]
print 'nobs = ',nobs

#-colors = energy band-
colors = [broadcolor]*nobs
for i in list(soft_i[:][0]):
    colors[i] = softcolor
for i in list(medium_i[:][0]):
    colors[i] = mediumcolor
for i in list(hard_i[:][0]):
    colors[i] = hardcolor
for i in list(veryhard_i[:][0]):
    colors[i] = veryhardcolor
for i in list(widehard_i[:][0]):
    colors[i] = widehardcolor
for i in list(widesoft_i[:][0]):
    colors[i] = widesoftcolor
for i in list(softmedium_i[:][0]):
    colors[i] = softmediumcolor
legcolors = bandcolors #defined at beginning of main
leglabels = bandlabels
legmews = [0.5]*len(bandcolors)
legmecs = ['black']*len(bandcolors)
legsyms = ['s']*len(bandcolors)


#-symbol = instrument (ACIS-HETG,ACIS-NONE,ACIS-LETG)-
hetg_sym = 'o' #HETG = circle
bare_sym = 'd' #bare acis = thin diamond
letg_sym = 'p' #LETG = hexagon
syms = [hetg_sym]*nobs
if args.ignoregrating == 'no':
    for i in list(bare_i[:][0]):
        syms[i] = bare_sym 
    for i in list(letg_i[:][0]):
        syms[i] = letg_sym
    if hetg_sym in syms: 
        legcolors+=['white'] #append to legend
        leglabels+=['ACIS+HETG']
        legmecs+=['black']
        legmews+=[0.5]
        legsyms+=[hetg_sym]
    if letg_sym in syms:
        legcolors+=['white'] #append to legend
        leglabels+=['ACIS+LETG']
        legmecs+=['black']
        legmews+=[0.5]
        legsyms+=[letg_sym]
    if bare_sym in syms:
        legcolors+=['white'] #append to legend
        leglabels+=['ACIS w/o grating']
        legmecs+=['black']
        legmews+=[0.5]
        legsyms+=[bare_sym]
else:
    legsyms = [hetg_sym]*len(bandcolors) #make color symbols match plot symbols

#-markeredgecolor = simz (detector z)-
offcolor = 'gray'
mecs = ['black']*nobs
if args.ignorechipy == 'no':
    for i in list(offz_i[:][0]):
        mecs[i] = offcolor #offset observations have offcolor outline

#-markeredgewidth = simz (detector z)-
offmew = 3.0
mews = [0.5]*nobs
if args.ignorechipy == 'no':
    for i in list(offz_i[:][0]):
        mews[i] = offmew #offset observations have different outline
if offcolor in mecs:
    legcolors+=['white'] #marker will look like border only
    leglabels+=['Low Chip Y']
    legmecs+=[offcolor]
    legmews+=[offmew]
    legsyms+=['o']

#-markersize - # = frametime-
#sizes = [f*3 for f in fitdata['frametime']]
sizes = 5

#-alpha (transparency)- 
if len(bands) > 1:
    alphas = [0.5]*nobs
else:
    alphas = [1.0]*nobs

print 'Markers Defined'

#----Set Up Plot File----
#pdffile = PdfPages(args.plotfile)

with PdfPages(args.plotfile) as pdffile:

    #----R0 vs age----

    #--Initialize Plot--
    ax1 = start_plot(xtitle='SN1987A Age [days]',ytitle='Radius [arcsec]',xmin=args.agemin,xmax=args.agemax,ymin=rmin,ymax=rmax)

    #----Plot----
    fancy_plot(ax1,fitdata['age'],fitdata['R0'],yerror=fitdata['R0errl'],yerror_high=fitdata['R0erru'],syms=syms,colors=colors,sizes=sizes,alphas=alphas,mecs=mecs,mews=mews,errorband=False)

    #----Plot broken-linear fits----
    if args.plotfit == 'yes':
        for bi in range(len(bands)):
            b = bands[bi]
            #--read best-fit model file--
            modeldata = np.genfromtxt(radiusfilename+'_'+b+'_mpfit_radii.txt',skip_header=0,names=True,comments='#',dtype=None)
            #--overplot model--
            fancy_plot(ax1,modeldata['age'],modeldata['radius'],yerror=modeldata['radiuslow'],yerror_high=modeldata['radiushigh'],syms='None',colors=bandcolors[bi],line='-',errorband=True) 

    #----Other noteable data points----
    # [age,radius(arcsec from center)]

    if args.plotextras != None:

        plot_extras(names=args.plotextras)
        
        ## # hot spot appearances in HST, Sugarman2002 table 1
        ## hotspot_ages = [2933,4283,4337,4337,4440,4440,4440,4725,4816,4816,4999,4999]
        ## hotspot_radii = [0.56,0.673,0.607,0.702,0.565,0.697,0.707,0.607,0.535,0.629,0.531,0.713]

        ## fancy_plot(hotspot_ages,hotspot_radii,syms='*',colors='white')
        ## legcolors+=['white'] #marker will look like border only
        ## leglabels+=['HST Hot Spots']
        ## legmecs+=['black']
        ## legmews+=[0.5]
        ## legsyms+=['*']
        
        ## #Larsson2013 -- inner debris (HST) morphology transitioned from core to edge-brightened due to
        ## #  energy source shifting from radioactive decay to X-ray illumination
        ## debris_transition_age = [5500,5500] 
        ## debris_transition_radius = [0.0,2.0]#plot vertical line

        ## fancy_plot(debris_transition_age,debris_transition_radius,syms=None,colors='black',line='--')
        ## #legcolors+=[None] #marker will look like border only
        ## #leglabels+=['Debris Transition']
        ## #legmecs+=[None]
        ## #legmews+=[0.]
        ## #legsyms+=[None]
        
        ## # Width and location of optical ring from UV flash, Plait1995
        ## # plot horizontal band
        ## plaitwidth = np.array([0.121/2.0,0.121/2.0])
        ## plait_radii = [0.86,0.86]
        ## plait_ages = [args.agemin,args.agemax]

        ## fancy_plot(plait_ages,plait_radii,colors='gray',syms=None,errorband=True,yerror=plaitwidth)

        ## # Location of ring from Sugarman2002 (table 3)
        ## sugarman_radii = [0.829,0.829]
        ## sugarman_ages = [args.agemin,args.agemax]
        
        ## #Ng2013 Radio expansion measurements: using the torus model found the following all around day 7600:
        ## # - break in expansion (decrease in velocity)
        ## # - break in torus opening angle, from constant ~40deg to decreasing
        ## # - break in asymmetry; the morphology was ~40% asymmetric, but became steadily more symmetric at the break
        ## # - break in the light curve; the slope flattened, deviating from previous exponential growth
        ## radiobreak_radii = [0.0,2.0] #plot vertical line
        ## radiobreak_age = [7600,7600]

        ## fancy_plot(radiobreak_age,radiobreak_radii,syms=None,colors='black',line=':')
                
    #----Plot Legend----
    plot_legend(labels=leglabels,colors=legcolors,markers=legsyms,mecs=legmecs,mews=legmews,fontsize='small')

    #--set top axis--
    set_axis(ax1,x=topx,twin=True,title=toptitle,xlab=toplab)

    #--save and close plot--
    pdffile.savefig()
    plt.close()

    print 'Radius vs Age plotted'

    #----North-South radius vs age----

    #--Initialize Plot--
    ax1 = start_plot(xtitle='SN1987A Age [days]',ytitle='North-South Radius [arcsec]',xmin=args.agemin,xmax=args.agemax,ymin=nsmin,ymax=nsmax)

    #----Plot----
    fancy_plot(ax1,fitdata['age'],fitdata['SKYNS'],yerror=fitdata['SKYNSerrl'],yerror_high=fitdata['SKYNSerru'],syms=syms,colors=colors,sizes=sizes,alphas=alphas,mecs=mecs,mews=mews,errorband=False)

    # --plot extras--
    if args.plotextras != None:

        plot_extras(names=args.plotextras)
            
    #----Plot Legend----
    plot_legend(labels=leglabels,colors=legcolors,markers=legsyms,mecs=legmecs,mews=legmews,fontsize='small')

    #--set top axis--
    set_axis(ax1,x=topx,twin=True,title=toptitle,xlab=toplab)

    #--save and close plot--
    pdffile.savefig()
    plt.close()

    print 'N-S Radius vs Age plotted'

    #----SigR vs age----

    #--Initialize Plot--
    ax1 = start_plot(xtitle='SN1987A Age [days]',ytitle='Ring Width [arcsec]',xmin=args.agemin,xmax=args.agemax,ymin=sigrmin,ymax=sigrmax)

    #----Plot----
    fancy_plot(ax1,fitdata['age'],fitdata['SIGR'],yerror=fitdata['SIGRerrl'],yerror_high=fitdata['SIGRerru'],syms=syms,colors=colors,sizes=sizes,alphas=alphas,mecs=mecs,mews=mews,errorband=False)

    #--set top axis--
    set_axis(ax1,x=topx,twin=True,title=toptitle,xlab=toplab)

    #--save and close plot--
    pdffile.savefig()
    plt.close()


    #---------------------------------------
    #       Plot Radius Histograms
    #---------------------------------------

    #----Create list of band indices----
    band_i = [list(soft_i[0]),list(medium_i[0]),list(hard_i[0]),list(broad_i[0]),list(veryhard_i[0]),list(widehard_i[0])] #np.where returns single element tuple (the element is an array)
    band_strings = ['300-800','800-1200','1200-8000','300-8000','3000-8000','2000-10000']
    bandcolors = [softcolor,mediumcolor,hardcolor,broadcolor,veryhardcolor,widehardcolor]
    #print band_i
    bi = 0

    #----Loop over bands----
    for b in band_i:
        if len(b) > 0:
            plot_param_hist(fitdata['R0'][b],band_strings[bi],pdffile,color=bandcolors[bi])
        bi = bi+1

    #---------------------------------------
    #       Plot Radius vs SigR
    #---------------------------------------

    #--Initialize Plot--
    ax1 = start_plot(xtitle='Radius [arcsec]',ytitle='Ring Width [arcsec]',xmin=rmin,xmax=rmax,ymin=sigrmin,ymax=sigrmax)

    #----Plot----
    fancy_plot(ax1,fitdata['R0'],fitdata['SIGR'],yerror=fitdata['SIGRerrl'],yerror_high=fitdata['SIGRerru'],xerror=fitdata['R0errl'],xerror_high=fitdata['R0erru'],syms=syms,colors=colors,sizes=sizes,alphas=alphas,mecs=mecs,mews=mews,errorband=False)

    #--save and close plot--
    pdffile.savefig()
    plt.close()

    print 'Radius vs SigR plotted'

#---------------------------------------
#             Wrap Up
#---------------------------------------

#----Close File----
#pdffile.close()
