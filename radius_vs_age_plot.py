# coding: utf-8
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib
import matplotlib.pyplot as plt
from sn1987a.sn1987a_time import *
import sn1987a.sn1987a_plots as snplt
from fancy_plot import fancy_plot

radio = True

#-load chandra observation results and information into a single dataframe-
obsinfo = pd.read_table('../chandra_observations.txt',sep='\t')
obsinfo.drop('age',inplace=True,axis=1)
df = pd.read_table('radii_fits_official.txt',sep=r'\s+',comment='#',na_values='None',engine='python')
df.set_index('obsid',inplace=True)

#-load expansion fit parameters-
broadparamdf = pd.read_table('radii_fits_official_300-8000_mpfit_params.txt',sep=r'\s*',comment='#',na_values=None,engine='python')

# plot parameters
agemin = 4500
agemax = 10900
ymin = 0.54
ymax = 0.9
alpha = 1.0

### plot broad band expansion curve ###
broaddf = df[df.band == '300-8000']
obsinfo.set_index('obsid',inplace=True)
broaddf = pd.concat([broaddf,obsinfo],join='inner',axis=1)

#-load data from the broken linear fit-
fitdf = pd.read_table('radii_fits_official_300-8000_mpfit_line.txt',sep='\t',comment='#',na_values=None)
#fit2df = pd.read_table('radii_fits_300-8000_300-8000_mpfit_line.txt',sep='\t',comment='#',na_values=None)
#fitdf = pd.read_table('radii_fits_official_300-8000_mpfit_line.txt',sep='\t',comment='#',na_values=None)

#-load data from Svet's density profile fit-
svetdf = pd.read_table('sn1987a_chandra_2016_density_expansion_fit.txt',header=6,sep=r'\s+',engine='python')

#-plot stuff-
pdffile=PdfPages('radius_vs_age_broad.pdf')

fig,ax=snplt.standard_age_plot(broaddf,'R0',ytitle='Radius [arcsec]',yerrlow='R0errl',yerrhigh='R0erru',errinterval=False,gratings='grating',ymin=ymin,ymax=ymax,agemin=agemin,agemax=agemax)

fancy_plot(ax,fitdf['age'],fitdf['radius'],colors='black',line='-',errorband=False,sizes=0.0000)
#fancy_plot(ax,fit2df['age'],fit2df['radius'],colors='black',line='--',errorband=False,sizes=0.0000)

#fancy_plot(ax,svetdf['age'],svetdf['fitradius'],colors='black',line='--',errorband=False,sizes=0.0000)

# legend
legsyms = ['s','d','o'] #square for each band, shape for each grating and for PN
leglabels = ['0.3-8 keV','ACIS (no grating)','ACIS (w/HETG)']
legcolors = ['black','white','white']
legmecs = ['black','black','black']
leg1 = snplt.plot_legend(ax,labels=leglabels,colors=legcolors,markers=legsyms,mecs=legmecs,fontsize='small',location='upper left')

# velocities annotation
v1 = 'v$_{early}='+str(int(round(broadparamdf['v_early'][0])))+'\pm'+str(int(round(broadparamdf['v_early'][1])))+'$ km s$^{-1}$'
v2 = 'v$_{late}~='+str(int(round(broadparamdf['v_late'][0])))+'\pm'+str(int(round(broadparamdf['v_late'][1])))+'$ km s$^{-1}$'
ax.annotate(v1+'\n'+v2,xycoords='figure fraction',textcoords='figure fraction',xy=(0.9,0.1),xytext=(0.5,0.22))

pdffile.savefig()
plt.close()
pdffile.close()

### plot rgb expansion curve ###

#-load chandra observation results and information into a single dataframe-
usoftdf = df[df.band=='300-800']
usoftdf = pd.concat([usoftdf,obsinfo],join='inner',axis=1)

softdf = df[df.band=='500-2000']
softdf = pd.concat([softdf,obsinfo],join='inner',axis=1)

harddf = df[df.band=='2000-10000']
harddf = pd.concat([harddf,obsinfo],join='inner',axis=1)

#-load fit parameters-
usoftparamdf = pd.read_table('radii_fits_official_300-800_mpfit_params.txt',sep=r'\s*',comment='#',na_values=None,engine='python')
softparamdf = pd.read_table('radii_fits_official_500-2000_mpfit_params.txt',sep=r'\s*',comment='#',na_values=None,engine='python')
hardparamdf = pd.read_table('radii_fits_official_2000-10000_mpfit_params.txt',sep=r'\s*',comment='#',na_values=None,engine='python')

#-load data from the broken linear fit-
usoftfitdf = pd.read_table('radii_fits_official_300-800_mpfit_line.txt',sep='\t',comment='#',na_values=None)
softfitdf = pd.read_table('radii_fits_official_500-2000_mpfit_line.txt',sep='\t',comment='#',na_values=None)
hardfitdf = pd.read_table('radii_fits_official_2000-10000_mpfit_line.txt',sep='\t',comment='#',na_values=None)

#-plot stuff-
pdffile=PdfPages('radius_vs_age_rgb.pdf')

fig,ax1=snplt.standard_age_plot(usoftdf,'R0',ytitle='Radius [arcsec]',yerrlow='R0errl',yerrhigh='R0erru',errinterval=False,gratings='grating',ymin=ymin,ymax=ymax,agemin=agemin,agemax=agemax,color='red',alphas=alpha)
snplt.standard_age_plot(softdf,'R0',yerrlow='R0errl',yerrhigh='R0erru',errinterval=False,gratings='grating',overplot=True,ax=ax1,color='green',alphas=alpha)
snplt.standard_age_plot(harddf,'R0',yerrlow='R0errl',yerrhigh='R0erru',errinterval=False,gratings='grating',overplot=True,ax=ax1,color='blue',alphas=alpha)
#snplt.standard_age_plot(broaddf,'R0',yerrlow='R0errl',yerrhigh='R0erru',errinterval=False,gratings='grating',overplot=True,ax=ax1,color=None,mecs='black')

fancy_plot(ax1,usoftfitdf['age'],usoftfitdf['radius'],colors='red',line='-',errorband=False,sizes=0.0000,mecs='red')
fancy_plot(ax1,softfitdf['age'],softfitdf['radius'],colors='green',line='-',errorband=False,sizes=0.0000,mecs='green')
fancy_plot(ax1,hardfitdf['age'],hardfitdf['radius'],colors='blue',line='-',errorband=False,sizes=0.0000,mecs='blue')

if radio==True:
# ng2013 ring radii
    radiodf = pd.read_table('../ng2013_ring_radii.txt',sep=r'\s*',comment='#',na_values='None',engine='python')
    snplt.standard_age_plot(radiodf,'majorradius',yerrlow='majorradiuserr',yerrhigh='majorradiuserr',errinterval=False,overplot=True,ax=ax1,color='black',sym='x',mews=1.5)
#    snplt.standard_age_plot(radiodf,'minorradius',yerrlow='minorradiuserr',yerrhigh='minorradiuserr',errinterval=False,overplot=True,ax=ax1,color='black',sym='x')

# legend
legsyms = ['s','s','s','d','o'] #square for each band, shape for each grating and for PN
leglabels = ['0.3-0.8 keV','0.5-2 keV','2-10 keV','ACIS (no grating)','ACIS (w/HETG)']
legcolors = ['red','green','blue','white','white']
legmecs = ['red','green','blue','black','black','black','black']
if radio == True:
    legsyms = ['s','s','s','d','o','x'] #square for each band, shape for each grating and for PN
    leglabels = ['0.3-0.8 keV','0.5-2 keV','2-10 keV','ACIS (no grating)','ACIS (w/HETG)','ATCA 9 GHz']
    legcolors = ['red','green','blue','white','white','black']
    legmecs = ['red','green','blue','black','black','black','black']
    vradio = 'v$_{ATCA}=3890\pm50$ km s$^{-1}$'
    
snplt.plot_legend(ax1,labels=leglabels,colors=legcolors,markers=legsyms,mecs=legmecs,fontsize='small')

# velocities annotation
sv1 = 'v$_{early}='+str(int(round(softparamdf['v_early'][0])))+'\pm'+str(int(round(softparamdf['v_early'][1])))+'$ km s$^{-1}$'
sv2 = 'v$_{late}~='+str(int(round(softparamdf['v_late'][0])))+'\pm'+str(int(round(softparamdf['v_late'][1])))+'$ km s$^{-1}$'
usv1 = 'v$_{early}='+str(int(round(usoftparamdf['v_early'][0])))+'\pm'+str(int(round(usoftparamdf['v_early'][1])))+'$ km s$^{-1}$'
usv2 = 'v$_{late}~='+str(int(round(usoftparamdf['v_late'][0])))+'\pm'+str(int(round(usoftparamdf['v_late'][1])))+'$ km s$^{-1}$'
hv1 = 'v$_{early}='+str(int(round(hardparamdf['v_early'][0])))+'\pm'+str(int(round(hardparamdf['v_early'][1])))+'$ km s$^{-1}$'
hv2 = 'v$_{late}~='+str(int(round(hardparamdf['v_late'][0])))+'\pm'+str(int(round(hardparamdf['v_late'][1])))+'$ km s$^{-1}$'

if radio == True:
    ax1.text((agemax-agemin)/2.0+agemin,(ymax-ymin)/4.0+ymin+0.04,vradio,color='black')
ax1.text((agemax-agemin)/2.0+agemin,(ymax-ymin)/4.0+ymin,usv1+'\n'+usv2,color='red')
ax1.text((agemax-agemin)/2.0+agemin,(ymax-ymin)/4.0+ymin-0.04,sv1+'\n'+sv2,color='green')
ax1.text((agemax-agemin)/2.0+agemin,(ymax-ymin)/4.0+ymin-0.08,hv1+'\n'+hv2,color='blue')

pdffile.savefig()
plt.close()
pdffile.close()

