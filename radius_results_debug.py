#! /usr/bin/env python

###! /usr/astro/ciao-4.6/bin/python

#Author: Kari A. Frank
#Date: 2015-04-03
#Purpose: 
#Usage: $PYDIR/sn1987a/radius_results.py [--radiusfile radiusfile] [--agemin agemin] [--agemax agemax] [--bands bands] [--mincounts mincounts] [--ignorechippos ignorechippos] [--ignoregrating ignoregrating] [--clean clean] [--clobber clobber]
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
#                 'soft' = '300-800,800-1200'
#                 'hard' = '3000-8000'
#                 'broad'= '300-8000'
#                 'all'  = '300-800,800-1200,1200-8000,300-8000,3000-8000,
#                              2000-10000'
#
# mincounts:  optional minimum required counts to include
#
# clean:      optional switch to change the axes labels and ticks to even
#             years and ages (default = 'yes')
#
# ignorechipy: optional switch to not plot gray (or thick) borders
#                for observations with offset chip positions. default='no' 
#
# ignoregrating: do not use different symbols for observations with gratings
#                (default='no')
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
from xmcplot import plot_param_hist

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
def plot_legend(labels=None,colors=['black'],markers=['o'],location='upper left',fontsize='medium',frameon=False,scatterpoints=1,numpoints=1,markerscale=1.5,mews=None,mecs=None):

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
    leg = plt.legend(loc=location,fontsize=fontsize,frameon=frameon,scatterpoints=scatterpoints,numpoints=numpoints,markerscale=markerscale)


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

parser.add_argument('--mincounts',help='Include only measurements from images iwth mincounts counts.',default=300)

parser.add_argument('--clean',help='Determines x-axis ticks and labels.',default='yes')

parser.add_argument('--ignorechipy',help='Plot normal borders for observations with offset chip positions.',default='no')
parser.add_argument('--ignoregrating',help='Do not use different symbols for observations with gratings.',default='no')

parser.add_argument('--clobber',help='Overwrite existing files.',default='no')

args = parser.parse_args()

#----Convert default types----
args.mincounts = int(args.mincounts)
args.agemin = int(args.agemin)
if args.agemax != None: args.agemax = int(args.agemax)

#----Parse bands----
if args.bands=='all': args.bands='300-800,800-1200,1200-8000,300-8000,3000-8000,2000-10000'
if args.bands=='broad': args.bands='300-8000'
if args.bands=='hard': args.bands='3000-8000'
if args.bands=='soft': args.bands='300-800,800-1200'
bands=args.bands.split(',')
print "Bands = ",bands
softcolor = 'red' #300-800
mediumcolor = 'green' #800-1200
hardcolor = 'blue' #1200-8000
veryhardcolor = 'cyan' #3000-8000
broadcolor = 'purple' #300-8000
widehardcolor = 'magenta' #2000-10000
bandcolors = ['']*len(bands)
bandlabels = ['']*len(bands)
for bi in range(len(bands)):
    if bands[bi] == '300-800':
        bandcolors[bi] = softcolor
        bandlabels[bi] = '0.3-0.8 keV'
    elif bands[bi] == '800-1200':
        bandcolors[bi] = mediumcolor
        bandlabels[bi] = '0.8-1.2 keV'
    elif bands[bi] == '1200-8000':
        bandcolors[bi] = hardcolor
        bandlabels[bi] = '1.2-8.0 keV'
    elif bands[bi] == '300-8000':
        bandcolors[bi] = broadcolor
        bandlabels[bi] = '0.3-8.0 keV' 
    elif bands[bi] == '3000-8000':
        bandcolors[bi] = veryhardcolor
        bandlabels[bi] = '3.0-8.0 keV'
    elif bands[bi] == '2000-10000':
        bandcolors[bi] = widehardcolor
        bandlabels[bi] = '2.0-10.0 keV'
    else:
        print "WARNING: unrecognized band."
        
#----Set file paths----

sn1987a_dir = '/export/bulk/rice1/kaf33/main/Chandra_Observations/SN1987A/'
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

#---------------------------------------
#       Plot Data vs Age
#---------------------------------------

#--Set whitespace on bottom, top--
bott = 0.13
topp = 0.87

#--Initialize Plot and Set Axes--
if args.agemax == None: args.agemax = 200+np.max(fitdata['age'])
if args.agemin == 0: args.agemin = np.min(fitdata['age'])-200
#rmin=0.5
#rmax=0.9
rmin=np.min(fitdata['R0'])-np.max(fitdata['R0errl'])
rmax=np.max(fitdata['R0erru'])+np.max(fitdata['R0'])
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
regular_ages = [convert_time(yr+'-01-01',get='age') for yr in year_ticks]
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
        legcolors+=['black'] #append to legend
        leglabels+=['ACIS+HETG']
        legmecs+=['black']
        legmews+=[0.5]
        legsyms+=[hetg_sym]
    if letg_sym in syms:
        legcolors+=['black'] #append to legend
        leglabels+=['ACIS+LETG']
        legmecs+=['black']
        legmews+=[0.5]
        legsyms+=[letg_sym]
    if bare_sym in syms:
        legcolors+=['black'] #append to legend
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
#sizes = 5
sizes = 10

#-alpha (transparency)- 
if len(bands) > 1:
    alphas = [0.5]*nobs
else:
    alphas = [1.0]*nobs

print 'Markers Defined'

#----Set Up Plot File----
#pdffile = PdfPages(args.plotfile)

#with PdfPages(args.plotfile) as pdffile:

#----R0 vs age----

#--Initialize Plot--
ax1 = start_plot(xtitle='SN1987A Age [days]',ytitle='Radius [arcsec]',xmin=args.agemin,xmax=args.agemax,ymin=rmin,ymax=rmax)

#----Plot----
fancy_plot(fitdata['age'],fitdata['R0'],yerror=fitdata['R0errl'],yerror_high=fitdata['R0erru'],syms=syms,colors=colors,sizes=sizes,alphas=alphas,mecs=mecs,mews=mews,errorband=False)

#----Plot Legend----
plot_legend(labels=leglabels,colors=legcolors,markers=legsyms,mecs=legmecs,mews=legmews,fontsize='small')

#--set top axis--
set_axis(ax1,x=topx,twin=True,title=toptitle,xlab=toplab)

#--save and close plot--
plt.show()
#pdffile.savefig()
#plt.close()

print 'Radius vs Age plotted'

#---------------------------------------
#             Wrap Up
#---------------------------------------

#----Close File----
#pdffile.close()
