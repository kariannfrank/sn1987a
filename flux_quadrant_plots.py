# coding: utf-8
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sn1987a.sn1987a_time import *
import sn1987a.sn1987a_plots as snplt
from sn1987a.broken_line import broken_line
from sn1987a.flux_function import flux_function
from fancy_plot import fancy_plot

# switches
plot_quadrants = False
plot_halves = True

sn1987a_dir = '/data/nasdrive/main/kaf33/sn1987a_monitoring/comparison_dir/'

#-load chandra observation results and information into a single dataframe-
obsinfo = pd.read_table('../chandra_observations.txt',sep=r'\s*',engine='python')
#obsinfo.drop('grating',inplace=True,axis=1)
#obsinfo.drop('age',inplace=True,axis=1)
df = pd.read_table(sn1987a_dir+'spectra_current/spectra465_frank2015a_fits_official.txt_fluxes.txt',
                   sep=r'\s*',comment='#',na_values='None',engine='python')
df.drop('grating',inplace=True,axis=1)
df.drop('age',inplace=True,axis=1)
obsinfo.set_index('obsid',inplace=True)
df.set_index('obsid',inplace=True)
df = pd.concat([df,obsinfo],join='inner',axis=1)

#-read in the quadrant fluxes-
qdf = pd.read_table('spectra_franka2015a_quad_fits_official.txt',sep=r'\s*',comment='#',engine='python')
qdf.set_index('obsid',inplace=True)
qdf = pd.concat([qdf,obsinfo],join='inner',axis=1) #add obsinfo

#-add extra errors-
# 20% = 10% astrometry + 10% pileup
# - in future, estimate pileup correction to be same as in the full 
#   ring spectra, and correct the fluxes according to this ratio
astroerr = 0.3
syserr = astroerr

qdf['nw_softlow'] = qdf['nw_softlow']-qdf['nw_soft']*syserr
qdf['nw_softhigh'] = qdf['nw_softhigh']+qdf['nw_soft']*syserr
qdf['ne_softlow'] = qdf['ne_softlow']-qdf['ne_soft']*syserr
qdf['ne_softhigh'] = qdf['ne_softhigh']+qdf['ne_soft']*syserr
qdf['se_softlow'] = qdf['se_softlow']-qdf['se_soft']*syserr
qdf['se_softhigh'] = qdf['se_softhigh']+qdf['se_soft']*syserr
qdf['sw_softlow'] = qdf['sw_softlow']-qdf['sw_soft']*syserr
qdf['sw_softhigh'] = qdf['sw_softhigh']+qdf['sw_soft']*syserr

qdf['nw_hardlow'] = qdf['nw_hardlow']-qdf['nw_hard']*syserr
qdf['nw_hardhigh'] = qdf['nw_hardhigh']+qdf['nw_hard']*syserr
qdf['ne_hardlow'] = qdf['ne_hardlow']-qdf['ne_hard']*syserr
qdf['ne_hardhigh'] = qdf['ne_hardhigh']+qdf['ne_hard']*syserr
qdf['se_hardlow'] = qdf['se_hardlow']-qdf['se_hard']*syserr
qdf['se_hardhigh'] = qdf['se_hardhigh']+qdf['se_hard']*syserr
qdf['sw_hardlow'] = qdf['sw_hardlow']-qdf['sw_hard']*syserr
qdf['sw_hardhigh'] = qdf['sw_hardhigh']+qdf['sw_hard']*syserr

qdf['nw_broadlow'] = qdf['nw_broadlow']-qdf['nw_broad']*syserr
qdf['nw_broadhigh'] = qdf['nw_broadhigh']+qdf['nw_broad']*syserr
qdf['ne_broadlow'] = qdf['ne_broadlow']-qdf['ne_broad']*syserr
qdf['ne_broadhigh'] = qdf['ne_broadhigh']+qdf['ne_broad']*syserr
qdf['se_broadlow'] = qdf['se_broadlow']-qdf['se_broad']*syserr
qdf['se_broadhigh'] = qdf['se_broadhigh']+qdf['se_broad']*syserr
qdf['sw_broadlow'] = qdf['sw_broadlow']-qdf['sw_broad']*syserr
qdf['sw_broadhigh'] = qdf['sw_broadhigh']+qdf['sw_broad']*syserr

# calculate east and west fluxes and errors
qdf['west_broad'] = qdf['nw_broad']+qdf['sw_broad']
qdf['west_broadlow'] = qdf['west_broad']-((qdf['nw_broad']-qdf['nw_broadlow'])**2.0+(qdf['sw_broad']-qdf['sw_broadlow'])**2.0)**0.5
qdf['west_broadhigh'] = qdf['west_broad']+((qdf['nw_broadhigh']-qdf['nw_broad'])**2.0+(qdf['sw_broadhigh']-qdf['sw_broad'])**2.0)**0.5
qdf['east_broad'] = qdf['ne_broad']+qdf['se_broad']
qdf['east_broadlow'] = qdf['east_broad']-((qdf['ne_broad']-qdf['ne_broadlow'])**2.0+(qdf['se_broad']-qdf['se_broadlow'])**2.0)**0.5
qdf['east_broadhigh'] = qdf['east_broad']+((qdf['ne_broadhigh']-qdf['ne_broad'])**2.0+(qdf['se_broadhigh']-qdf['se_broad'])**2.0)**0.5
qdf['west_soft'] = qdf['nw_soft']+qdf['sw_soft']
qdf['west_softlow'] = qdf['west_soft']-((qdf['nw_soft']-qdf['nw_softlow'])**2.0+(qdf['sw_soft']-qdf['sw_softlow'])**2.0)**0.5
qdf['west_softhigh'] = qdf['west_soft']+((qdf['nw_softhigh']-qdf['nw_soft'])**2.0+(qdf['sw_softhigh']-qdf['sw_soft'])**2.0)**0.5
qdf['east_soft'] = qdf['ne_soft']+qdf['se_soft']
qdf['east_softlow'] = qdf['east_soft']-((qdf['ne_soft']-qdf['ne_softlow'])**2.0+(qdf['se_soft']-qdf['se_softlow'])**2.0)**0.5
qdf['east_softhigh'] = qdf['east_soft']+((qdf['ne_softhigh']-qdf['ne_soft'])**2.0+(qdf['se_softhigh']-qdf['se_soft'])**2.0)**0.5
qdf['west_hard'] = qdf['nw_hard']+qdf['sw_hard']
qdf['west_hardlow'] = qdf['west_hard']-((qdf['nw_hard']-qdf['nw_hardlow'])**2.0+(qdf['sw_hard']-qdf['sw_hardlow'])**2.0)**0.5
qdf['west_hardhigh'] = qdf['west_hard']+((qdf['nw_hardhigh']-qdf['nw_hard'])**2.0+(qdf['sw_hardhigh']-qdf['sw_hard'])**2.0)**0.5
qdf['east_hard'] = qdf['ne_hard']+qdf['se_hard']
qdf['east_hardlow'] = qdf['east_hard']-((qdf['ne_hard']-qdf['ne_hardlow'])**2.0+(qdf['se_hard']-qdf['se_hardlow'])**2.0)**0.5
qdf['east_hardhigh'] = qdf['east_hard']+((qdf['ne_hardhigh']-qdf['ne_hard'])**2.0+(qdf['se_hardhigh']-qdf['se_hard'])**2.0)**0.5



#-calculate hardness ratio and errors, add as columns-

# full ring
df['hratio'] = df['hard']/df['soft']
df['hratiolow'],df['hratiohigh'] = snplt.ratio_err(df['hard'],df['soft'],df['hardlow'],df['hardhigh'],df['softlow'],df['softhigh'])

# quadrants
qdf['nw_hratio'] = qdf['nw_hard']/qdf['nw_soft']
qdf['nw_hratiolow'],qdf['nw_hratiohigh'] = snplt.ratio_err(qdf['nw_hard'],qdf['nw_soft'],qdf['nw_hardlow'],qdf['nw_hardhigh'],qdf['nw_softlow'],qdf['nw_softhigh'])
qdf['ne_hratio'] = qdf['ne_hard']/qdf['ne_soft']
qdf['ne_hratiolow'],qdf['ne_hratiohigh'] = snplt.ratio_err(qdf['ne_hard'],qdf['ne_soft'],qdf['ne_hardlow'],qdf['ne_hardhigh'],qdf['ne_softlow'],qdf['ne_softhigh'])
qdf['se_hratio'] = qdf['se_hard']/qdf['se_soft']
qdf['se_hratiolow'],qdf['se_hratiohigh'] = snplt.ratio_err(qdf['se_hard'],qdf['se_soft'],qdf['se_hardlow'],qdf['se_hardhigh'],qdf['se_softlow'],qdf['se_softhigh'])
qdf['sw_hratio'] = qdf['sw_hard']/qdf['sw_soft']
qdf['sw_hratiolow'],qdf['sw_hratiohigh'] = snplt.ratio_err(qdf['sw_hard'],qdf['sw_soft'],qdf['sw_hardlow'],qdf['sw_hardhigh'],qdf['sw_softlow'],qdf['sw_softhigh'])

#-load xmm observation results and information into a single dataframe-
xmmdf = pd.read_table(sn1987a_dir+'xmm_fluxes_nopileup.txt',sep='\t',comment='#',na_values='None')#,dtype={'obsid':str,'date':str,'age':int,'filter':str,'totalexp':float,'filtexp':float,'soft':float,'softlow':float,'softhigh':float,'hard':float,'hardlow':float,'hardhigh':float})
xmmdf.set_index('obsid',inplace=True)

#---plot lightcurves---
pdffile=PdfPages('quadrant_lightcurves.pdf')

fullcolor = 'black'
nwcolor = 'purple'#'mediumpurple'
necolor = 'darkgreen'#'mediumseagreen'
swcolor = 'orangered'
secolor = 'royalblue'
eastcolor = 'darkcyan'
westcolor = 'mediumvioletred'

qsize = 4

qleglabels = ['NW','NE','SE','SW']
qlegcolors = [nwcolor,necolor,secolor,swcolor]
hleglabels = ['W','E','','']
hlegcolors = [westcolor,eastcolor,eastcolor,westcolor]

legsyms = ['d','o','p','^'] #shape for each grating and for PN
leglabels = ['ACIS (no grating)','ACIS (w/HETG)','LETG','EPIC-pn']
legcolors = ['white','white','white','white']
legmecs = ['black','black','black','black']

#-broad band-

# full ring chandra
fig,ax=snplt.standard_age_plot(df,'broad',ytitle='0.3-8 keV Flux [10$^{-13}$ erg cm$^{-2}$ s$^{-1}$]',yerrlow='broadlow',yerrhigh='broadhigh',ymax=200.0,color=fullcolor,gratings='grating',agemin=4400,agemax=10600,ylog=True,ymin=0.4)

# quadrants chandra
if plot_quadrants == True:
    snplt.standard_age_plot(qdf,'nw_broad',yerrlow='nw_broadlow',yerrhigh='nw_broadhigh',color=nwcolor,overplot=True,ax=ax,gratings='grating',symsize=qsize,line='-')
    snplt.standard_age_plot(qdf,'ne_broad',yerrlow='ne_broadlow',yerrhigh='ne_broadhigh',color=necolor,overplot=True,ax=ax,gratings='grating',symsize=qsize,line='-')
    snplt.standard_age_plot(qdf,'se_broad',yerrlow='se_broadlow',yerrhigh='se_broadhigh',color=secolor,overplot=True,ax=ax,gratings='grating',symsize=qsize,line='-')
    snplt.standard_age_plot(qdf,'sw_broad',yerrlow='sw_broadlow',yerrhigh='sw_broadhigh',color=swcolor,overplot=True,ax=ax,gratings='grating',symsize=qsize,line='-')

# halves chandra
if plot_halves == True:
    snplt.standard_age_plot(qdf,'east_broad',yerrlow='east_broadlow',yerrhigh='east_broadhigh',overplot=True,ax=ax,gratings='grating',line='-',color=eastcolor)
    snplt.standard_age_plot(qdf,'west_broad',yerrlow='west_broadlow',yerrhigh='west_broadhigh',overplot=True,ax=ax,gratings='grating',line='-',color=westcolor)

if plot_quadrants == True:
    legax = snplt.plot_pie_legend(ax,labels=qleglabels,colors=qlegcolors,radius=0.7,pos='bottomright')
if plot_halves == True and plot_quadrants == False:
    legax = snplt.plot_pie_legend(ax,labels=hleglabels,colors=hlegcolors,radius=0.7,pos='bottomright')    
snplt.plot_legend(ax,labels=leglabels[:-1],colors=legcolors[:-1],markers=legsyms[:-1],mecs=legmecs[:-1],fontsize='small')

pdffile.savefig()
plt.close()

#-soft band-
# full ring chandra
fig,ax=snplt.standard_age_plot(df,'soft',ytitle='0.5-2 keV Flux [10$^{-13}$ erg cm$^{-2}$ s$^{-1}$]',yerrlow='softlow',yerrhigh='softhigh',ymax=200.0,color=fullcolor,gratings='grating',agemin=4400,agemax=10600,ylog=True,ymin=0.2)

# quadrants chandra
if plot_quadrants == True:
    snplt.standard_age_plot(qdf,'nw_soft',yerrlow='nw_softlow',yerrhigh='nw_softhigh',color=nwcolor,overplot=True,ax=ax,gratings='grating',symsize=qsize,line='-')
    snplt.standard_age_plot(qdf,'ne_soft',yerrlow='ne_softlow',yerrhigh='ne_softhigh',color=necolor,overplot=True,ax=ax,gratings='grating',symsize=qsize,line='-')
    snplt.standard_age_plot(qdf,'se_soft',yerrlow='se_softlow',yerrhigh='se_softhigh',color=secolor,overplot=True,ax=ax,gratings='grating',symsize=qsize,line='-')
    snplt.standard_age_plot(qdf,'sw_soft',yerrlow='sw_softlow',yerrhigh='sw_softhigh',color=swcolor,overplot=True,ax=ax,gratings='grating',symsize=qsize,line='-')

# full ring xmm
snplt.standard_age_plot(xmmdf,'soft',yerrlow='softlow',yerrhigh='softhigh',color=fullcolor,overplot=True,ax=ax,sym='^')

# halves chandra
if plot_halves == True:
    snplt.standard_age_plot(qdf,'east_soft',yerrlow='east_softlow',yerrhigh='east_softhigh',overplot=True,ax=ax,gratings='grating',line='-',color=eastcolor)
    snplt.standard_age_plot(qdf,'west_soft',yerrlow='west_softlow',yerrhigh='west_softhigh',overplot=True,ax=ax,gratings='grating',line='-',color=westcolor)

# legend
if plot_quadrants == True:
    legax = snplt.plot_pie_legend(ax,labels=qleglabels,colors=qlegcolors,radius=0.7,pos='bottomright')
if plot_halves == True and plot_quadrants == False:
    legax = snplt.plot_pie_legend(ax,labels=hleglabels,colors=hlegcolors,radius=0.7,pos='bottomright')
snplt.plot_legend(ax,labels=leglabels,colors=legcolors,markers=legsyms,mecs=legmecs,fontsize='small')

pdffile.savefig()
plt.close()

#-hard band-
# full ring chandra
fig,ax=snplt.standard_age_plot(df,'hard',ytitle='3-8 keV Flux [10$^{-13}$ erg cm$^{-2}$ s$^{-1}$]',yerrlow='hardlow',yerrhigh='hardhigh',ymax=30.0,color=fullcolor,gratings='grating',agemin=4400,agemax=10600,ylog=True,ymin=0.1)

# quadrants chandra
if plot_quadrants == True:
    snplt.standard_age_plot(qdf,'nw_hard',yerrlow='nw_hardlow',yerrhigh='nw_hardhigh',color=nwcolor,overplot=True,ax=ax,gratings='grating',symsize=qsize,line='-')
    snplt.standard_age_plot(qdf,'ne_hard',yerrlow='ne_hardlow',yerrhigh='ne_hardhigh',color=necolor,overplot=True,ax=ax,gratings='grating',symsize=qsize,line='-')
    snplt.standard_age_plot(qdf,'se_hard',yerrlow='se_hardlow',yerrhigh='se_hardhigh',color=secolor,overplot=True,ax=ax,gratings='grating',symsize=qsize,line='-')
    snplt.standard_age_plot(qdf,'sw_hard',yerrlow='sw_hardlow',yerrhigh='sw_hardhigh',color=swcolor,overplot=True,ax=ax,gratings='grating',symsize=qsize,line='-')

# full ring xmm
snplt.standard_age_plot(xmmdf,'hard',yerrlow='hardlow',yerrhigh='hardhigh',color=fullcolor,overplot=True,ax=ax,sym='^')

# halves chandra
if plot_halves == True:
    snplt.standard_age_plot(qdf,'east_hard',yerrlow='east_hardlow',yerrhigh='east_hardhigh',overplot=True,ax=ax,gratings='grating',line='-',color=eastcolor)
    snplt.standard_age_plot(qdf,'west_hard',yerrlow='west_hardlow',yerrhigh='west_hardhigh',overplot=True,ax=ax,gratings='grating',line='-',color=westcolor)

# legend
if plot_quadrants == True:
    legax = snplt.plot_pie_legend(ax,labels=qleglabels,colors=qlegcolors,radius=0.7,pos='bottomright')
if plot_halves == True and plot_quadrants == False:
    legax = snplt.plot_pie_legend(ax,labels=hleglabels,colors=hlegcolors,radius=0.7,pos='bottomright')
snplt.plot_legend(ax,labels=leglabels,colors=legcolors,markers=legsyms,mecs=legmecs,fontsize='small')

pdffile.savefig()
plt.close()
#pdffile.close()

#-----plot hardness ratio-----

#pdffile=PdfPages('hardness_ratio.pdf')

# full ring chandra
fig,ax1=snplt.standard_age_plot(df,'hratio',ytitle='F$_{3.0-8.0 keV}$ / F$_{0.5-2.0keV}$',yerrlow='hratiolow',yerrhigh='hratiohigh',color=fullcolor,ymin=0.0,gratings='grating',agemin=4400,agemax=10600,errinterval=True,ymax=0.5)

# quadrants chandra
snplt.standard_age_plot(qdf,'nw_hratio',yerrlow='nw_hratiolow',yerrhigh='nw_hratiohigh',color=nwcolor,overplot=True,ax=ax1,gratings='grating')
snplt.standard_age_plot(qdf,'ne_hratio',yerrlow='ne_hratiolow',yerrhigh='ne_hratiohigh',color=necolor,overplot=True,ax=ax1,gratings='grating')
snplt.standard_age_plot(qdf,'se_hratio',yerrlow='se_hratiolow',yerrhigh='se_hratiohigh',color=secolor,overplot=True,ax=ax1,gratings='grating')
snplt.standard_age_plot(qdf,'sw_hratio',yerrlow='sw_hratiolow',yerrhigh='sw_hratiohigh',color=swcolor,overplot=True,ax=ax1,gratings='grating')

# legend
legax = snplt.plot_pie_legend(ax,labels=qleglabels,colors=qlegcolors,radius=0.7)
snplt.plot_legend(ax,labels=leglabels[:-1],colors=legcolors[:-1],markers=legsyms[:-1],mecs=legmecs[:-1],fontsize='small')


pdffile.savefig()
plt.close()

#-----plot quadrant fluxes as fraction of total-----

## # quadrants chandra
## snplt.standard_age_plot(qdf,'nw_hard',yerrlow='nw_hardlow',yerrhigh='nw_hardhigh',color=nwcolor,overplot=True,ax=ax,gratings='grating')
## snplt.standard_age_plot(qdf,'ne_hard',yerrlow='ne_hardlow',yerrhigh='ne_hardhigh',color=necolor,overplot=True,ax=ax,gratings='grating')
## snplt.standard_age_plot(qdf,'se_hard',yerrlow='se_hardlow',yerrhigh='se_hardhigh',color=secolor,overplot=True,ax=ax,gratings='grating')
## snplt.standard_age_plot(qdf,'sw_hard',yerrlow='sw_hardlow',yerrhigh='sw_hardhigh',color=swcolor,overplot=True,ax=ax,gratings='grating')

## # legend
## legax = snplt.plot_pie_legend(ax,labels=qleglabels,colors=qlegcolors,radius=0.7)
## snplt.plot_legend(ax,labels=leglabels,colors=legcolors,markers=legsyms,mecs=legmecs,fontsize='small')


## pdffile.savefig()
## plt.close()

pdffile.close()

