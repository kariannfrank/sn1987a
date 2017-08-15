"""
 Misc. Functions for plotting sn1987a data


 Contains the following functions:

 set_plot_fonts
 set_axis
 start_plot
 top_year_axis
 plot_extras
 
""" 
#----------------------------------------------------------
# Import Common Modules
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from sn1987a_time import convert_time
from fancy_plot import fancy_plot

#----------------------------------------------------------
def ratio_err(top,bottom,top_low,top_high,bottom_low,bottom_high):
    """
    function to derive vectors of errors on ratios

    uses simple propagation of errors (partial derivatives)
    returns lower and upper errorbars
    """
    
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

#----------------------------------------------------------
def set_plot_fonts(fontsize=15,fontfamily='serif',fontweight='normal'):
   """Set the font properties for a plot"""
   
   font = {'family':fontfamily,'weight':fontweight,'size':fontsize}
   matplotlib.rc('font',**font)
   matplotlib.rcParams['text.latex.preamble'] = [r"\boldmath"]
   plt.rc('text',usetex=True)

#----------------------------------------------------------
def set_axis(ax,x=None,twin=False,title=None,xlab=None,grid=False,
	     clean=True):	
    """Set up a bottom or top time axis for SN1987A plot"""

#    font_size = 10
    font_size = 13
#    fontangle = 50
    fontangle = 40

    if (twin == False) and (clean == True):
        fontangle = 0

    if clean == False:
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

    return ax0

#----------------------------------------------------------
def start_plot(xtitle='',ytitle='',xmin=None,xmax=None,ymin=None,ymax=None,
               ylog=False,clean=True,figsize=None):
    """Initialize a SN1987A plot with a time x-axis on bottom"""
    
    #-initialize main plot-
    fig,ax1=plt.subplots()

    #-set size of plot-
    if figsize != None: fig.set_size_inches(figsize[0],figsize[1])
    
    #-set axis titles-
    plt.xlabel(xtitle)
    plt.ylabel(ytitle)

    #-set bottom axis-
    if ylog == True: plt.yscale('log')
    if clean is True:
        grid = False
    else:
        grid = True
    set_axis(ax1,grid=grid,clean=clean)
    if xmin is not None or xmax is not None:
        ax1.set_xlim(xmin,xmax)
        ax1.set_autoscalex_on(False)
    if ymin is not None or ymax is not None:
        ax1.set_ylim(ymin,ymax)
        ax1.set_autoscaley_on(False)
                
    #-add extra space on bottom-
    fig.subplots_adjust(bottom=0.13,top=0.87)
    
    return fig,ax1

#----------------------------------------------------------
def top_year_axis(ax,agemin=4400,agemax=10600):
     """Add a top x-axis with year ticks to a plot with a bottom age axis"""
     
     #--get years in range--
     years = [str(yr) for yr in range(convert_time(agemin,get='year')+1,
				      convert_time(agemax,get='year')+1)]


     #--get sn1987a ages associated with the years--
     ages = [convert_time(yr+'-01-01',get='age',informat='date') 
             for yr in years]

     ax2 = set_axis(ax,x=ages,twin=True,title='Year',xlab=years)

     return ax2
 
#----------------------------------------------------------
def plot_legend(ax,labels=None,colors=['black'],markers=['o'],
                location='upper left',
                fontsize='medium',frameon=False,scatterpoints=1,numpoints=1,
                markerscale=1.5,mews=None,mecs=None,ncol=1,alphas=None,
                mealphas=None,
                textcolors=None,textshifts=None):
    import numbers

    #-convert defaults to lists-
    if mews is None or isinstance(mews,numbers.Number): 
        mews = [mews]*len(labels)
    if mecs is None: mecs = [None]*len(labels)
    if alphas is None: 
        alphas = [1.0]*len(labels)
    if mealphas is None:
        mealphas = alphas
    if isinstance(textcolors,str):
        textcolors = [textcolors]*len(labels)
    if isinstance(textshifts,numbers.Number):
        textshifts = [textshifts]*len(labels)

    mecs = [matplotlib.colors.colorConverter.to_rgba(mecs[i],alpha = mealphas[i]) for i in xrange(len(labels))]
 
    #-dummy position arrays-
    emptyx=[0]
    emptyy=emptyx

    #-plot dummy points to define labels-
    for p in range(len(labels)):
#        ax.scatter(emptyx,emptyy,color=colors[p],marker=markers[p],label=labels[p],edgecolor=mecs[p],linewidth=mews[p],alpha=alphas[p])
        ax.plot(emptyx,emptyy,linestyle='None',color=colors[p],
                marker=markers[p],label=labels[p],markeredgecolor=mecs[p],
                markeredgewidth=mews[p],alpha=alphas[p])

    #-plot legend-
    leg = ax.legend(loc=location,fontsize=fontsize,frameon=frameon,
                     scatterpoints=scatterpoints,numpoints=numpoints,
                     markerscale=markerscale,ncol=ncol,handletextpad=0)

    #-change text colors-
    if textcolors is not None or textshifts is not None:
        labs = leg.get_texts()
        lc = 0
        for l in labs:
            if textcolors is not None: l.set_color(textcolors[lc])
            if textshifts is not None: l.set_position((textshifts[lc],0))
            lc = lc + 1

    return leg
    
#----------------------------------------------------------
def plot_pie_legend(ax,labels=['NW','NE','SE','SW'],
                    colors=['red','orange','blue','green'],
                    radius=1.0,sizes=[25,25,25,25],fontsizes=None,
                    yoffset=0.35,xoffset=0.43,textcolor='white',
                    pos='topright',rec=None,textcolors=None,linewidth=0):
    """Plot a legend in the form of a pie chart"""

    if textcolors is None:
        textcolors = [textcolor]*len(labels)
    if isinstance(fontsizes,float): fontsizes = int(fontsizes)
    if isinstance(fontsizes,int):
        fontsizes = [fontsizes]*len(labels)

    if len(labels) == 2:
        startangle=90
    else:
        startangle=0

    #-create inset axes-
    #  rec = [x,y,xsize,ysize]
    if rec is None: # if rec is provided, then it overrides the pos argument
        if pos == 'topright':
            rec = [0.7,0.67,0.2,0.2]
        elif pos == 'topleft':
            rec = [0.1,0.67,0.2,0.2] 
        elif pos == 'bottomright':
            rec = [0.7,0.15,0.2,0.2]
        elif pos == 'bottomleft':
            rec = [0.1,0.2,0.2,0.2] 
        else:
            rec = [0.7,0.67,0.2,0.2] # top right
    
    pax = plt.axes(rec)

    #-plot pie chart figure-
    wedges,texts=pax.pie(sizes,labels=labels,
                         colors=colors,radius=radius,startangle=startangle)
    pax.set_aspect('equal') # force a circle

    #-set label properties-
    yoffset=yoffset
    xoffset=xoffset
    labely = [yoffset-0.05,yoffset-0.05,-yoffset,-yoffset]
    labelx = [xoffset,-xoffset,-xoffset,xoffset]
    q = 0
    for t in texts:
        t.set_horizontalalignment('center')
        t.set_verticalalignment('center')
        t.set_y(labely[q])
        t.set_x(labelx[q])
        t.set_color(textcolors[q])
        if fontsizes is not None: t.set_fontsize(fontsizes[q])
        q = q+1
    for w in wedges:
        w.set_linewidth(linewidth)

    return pax

#----------------------------------------------------------
def plot_extras(names='all',agemin=4000,agemax=11000):
	"""
	Author: Kari A. Frank
	Date: September 22, 2015
	Purpose: Plot lines on already open age vs Y plots (where Y is 
	      typically either radius or flux

	Input:
	 names : string
	     comma separated list of names of which extras to plot. 
	     default = 'all'. options are:
	     'plait1995'- plots position and width of the ER from 
	          Plait1995 (in this case Y must be radius)
	     'ng2013' - plots the radio break age day=7600 from the 
	          Ng2013 torus fits to radio images, as vertical line
	     'sugarman2002' - plots location of the ring from Sugarman 2002 
	          (in this case Y must be radius)
	     'hotspots' - plots the radius and day of appearance for the 
	          HST hot spots in Sugarman 2002
	     'larsson2014' - plots the age of transition from radioactive 
	         decay to X-ray heating as inner debris energy source
	      'chandrabroadbreak' - plot as vertical line the best fit 
	           age of the break point in expansion curve
	      'plaitwidth' - plot distance from chandrabroadbreak radius 
	           to width of plait ring

	Output:
	 - plots to an already open plot

	Usage Notes:
	 - will not add anything to the plot legend
	"""

        dummy_y = [0.00000001,1000.] #dummy y values for vertical lines
        dummy_sym = '.'
        
        if names == 'all' or 'hotspots' in names:
                # hot spot appearances in HST, Sugarman2002 table 1
                hotspot_ages = [2933,4283,4337,4337,4440,4440,4440,4725,
				4816,4816,4999,4999]
                hotspot_radii = [0.56,0.673,0.607,0.702,0.565,0.697,0.707,
				 0.607,0.535,0.629,0.531,0.713]

                fancy_plot(hotspot_ages,hotspot_radii,syms='*',
			   colors='white')
                #        legcolors+=['white'] #marker will look like border only
                #        leglabels+=['HST Hot Spots']
                #        legmecs+=['black']
                #        legmews+=[0.5]
                #        legsyms+=['*']

        if names == 'all' or 'larsson2014' in names:
                #Larsson2013 -- inner debris (HST) morphology transitioned 
		#  from core to edge-brightened due to
                #  energy source shifting from radioactive decay to X-ray 
		#  illumination
                debris_transition_age = [5500,5500] 
                debris_transition_radius = dummy_y#plot vertical line

                fancy_plot(debris_transition_age,debris_transition_radius,
			   colors='black',line='--',syms=dummy_sym)
                #legcolors+=[None] #marker will look like border only
                #leglabels+=['Debris Transition']
                #legmecs+=[None]
                #legmews+=[0.]
                #legsyms+=[None]

        if names == 'all' or 'plait1995' in names:
                # Width and location of optical ring from UV flash, 
		#  Plait1995
                # plot horizontal band
                plaitwidth = np.array([0.121/2.0,0.121/2.0])
                plait_radii = [0.86,0.86]
                plait_ages = [agemin,agemax]
                fancy_plot(plait_ages,plait_radii,colors='gray',
			   syms=dummy_sym,errorband=True,yerror=plaitwidth)

        if names == 'all' or 'plaitwidth' in names:
                # width of Plait ring, from chandra break point
                plaitwidth = np.array([0.121/2.0,0.121/2.0])
                plait_ages = [-1.0,20000.0]
                chandrabreakradii = [0.73+0.121/2.0,0.73+0.121/2.0]
                fancy_plot(plait_ages,chandrabreakradii,colors='purple',
			   syms=dummy_sym,errorband=True,yerror=plaitwidth)
                
        if names == 'all' or 'sugarman2002' in names:
                # Location of ring from Sugarman2002 (table 3)
                sugarman_radii = [0.829,0.829]
                sugarman_ages = [agemin,agemax]
                fancy_plot(plait_ages,plait_radii,colors='gray',
			   syms=dummy_sym,errorband=True,yerror=plaitwidth)
                
        if names == 'all' or 'ng2013' in names:
                #Ng2013 Radio expansion measurements: using the torus 
		#  model found the following all around day 7600:
                # - break in expansion (decrease in velocity)
                # - break in torus opening angle, from constant ~40deg 
		#   to decreasing
                # - break in asymmetry; the morphology was ~40% asymmetric,
		#   but became steadily more symmetric at the break
                # - break in the light curve; the slope flattened, 
		#   deviating from previous exponential growth
                radiobreak_radii = dummy_y #plot vertical line
                radiobreak_age = [7600,7600]

                fancy_plot(radiobreak_age,radiobreak_radii,syms=dummy_sym,
			   colors='black',line=':')

        if names == 'all' or 'chandrabroadbreak' in names:
                # plot best fit age of the break in expansion from the 
		#  300-8000 chandra images
                chandrabreak_radii = dummy_y #plot vertical line
                chandrabreak_age = [5962.77,5962.77]

                fancy_plot(chandrabreak_age,chandrabreak_radii,
			   syms=dummy_sym,colors='black',line='--')

        if (names == 'all' or 'ng2013radii' in names or 'ng2013fluxes' 
	    in names):
            ngradiifile = ('/astro/research/kaf33/Dropbox/Research/SN1987A'
			   '/ng2013_ring_radii.txt')
            ng2013_ringfits = pd.read_table(ngradiifile,comment='#',
					    header=0,names=
					    ['age','flux','fluxerr','majorradius',
					     'majorradiuserr','minorradius',
					     'minorradiuserr','asymmetry',
					     'asymmetryerr','phi','phierr','chi2',
					     'dof'],index_col=False)
             # columns: age(index)	     	flux	fluxerr		majorradius	majorradiuserr	minorradius	minorradiuserr	asymmetry	asymmetryerr	phi	phierr	chi2 dof                        
                
        if names == 'all' or 'ng2013radii' in names:

#             plt.plot(ng2013_ringfits['age'],ng2013_ringfits[])
            fancy_plot(np.array(ng2013_ringfits['age']),
		       np.array(ng2013_ringfits['majorradius']),
		       yerror=np.array(ng2013_ringfits['majorradiuserr']),
		       colors='black',syms='x')
            fancy_plot(np.array(ng2013_ringfits['age']),
		       np.array(ng2013_ringfits['minorradius']),
		       yerror=np.array(ng2013_ringfits['minorradiuserr']),
		       colors='black',syms='x')
    
        if names == 'all' or 'ng2013fluxes' in names:
            ngradiifile = '~/Dropbox/Research/SN1987A/ng2013_ring_radii.txt'
            ng2013_ringfits = pd.read_table(ngradiifile,comment='#',
					    header=0,names=
					    ['age','flux','fluxerr','majorradius',
					     'majorradiuserr','minorradius',
					     'minorradiuserr','asymmetry',
					     'asymmetryerr','phi','phierr','chi2',
					     'dof'],index_col=False)
             # columns: age(index)	     	flux	fluxerr		majorradius	majorradiuserr	minorradius	minorradiuserr	asymmetry	asymmetryerr	phi	phierr	chi2 dof

#             plt.plot(ng2013_ringfits['age'],ng2013_ringfits[])
            fancy_plot(np.array(ng2013_ringfits['age']),
		       np.array(ng2013_ringfits['flux']/3.0),
		       yerror=np.array(ng2013_ringfits['fluxerr']),
		       colors='black',syms='x')

#----------------------------------------------------------
def standard_age_plot(df,y,ages='age',gratings=None,simoffsets=None,
                      agemin=4400,agemax=10600,ymin=0.0,ymax=None,
                      plotextras=None,clean=True,color='black',
                      overplot=False,ax=None,symsize=8,ytitle=None,
                      sym='o',yerrlow=None,yerrhigh=None,errinterval=True,
                      ylog=False,errorband=False,line='',linewidth=1,
                      figsize=None,mews=None,mecs=None,alphas=1.0,
                      mealphas=None,linecolor=None,
                      **kwargs):
     """
     Create a standard plot of age vs Y for SN1987A data
     
     Author: Kari A. Frank
     Date: 2015-12-18

     Input

     df : pd.DataFrame
       dataframe containing columns for, at a minimum, observation ages and
       y values
     ages : string
       name of the column in df containing the observation ages
     y : string
       name of the column in df containing the y values to be plotted
     yerrlow, yerrhigh : string
       name of the columns containing the upper and lower errors on y
     errinterval : bool
       if True, then assume errors are the lower and upper endpoints of the 
       error interval and convertthem to error bars before plotting.
     gratings : string
       name of the column in df specifying the gratings used for each 
       observation represented in x (and y). if provided, each grating
       type will be plotted with a different symbol
     simoffsets : string
       name of column in df specifying detector sim offsets for each 
       observation. if provided, the symbol outlines for these 
       observations will be gray instead of the usual black
     agemin,agemax : numerical
       minimum and maximum ages, in days for the x-axis
     ymin, ymax : numerical
       optionally specify the minimum and/or maximum values for the y-axis
     plotextras : string
       comma separated list of names for the extras to plot (passed directly
       to plot_extras)
     clean : bool
       if False, will plot x ticks and a vertical line across the plot at 
       every data point to more easily associate specific observations with 
       a specific data point
     color : string
       specify a color for all the points. to plot with more than one color,
       call this function again with different y and overplot=True.
     overplot : bool
       if True, will skip creating a new figure and just plot the 
       provided x and y on the current figure. must also provide 
       the ax argument.
     ax : axis object
       used only if overplot=True, to make sure it overplotted on 
       the correct axis
     sym : string
       specify the symbol to use (will be ignored if gratings is given)
     symsize : numerical
       specify the size of the symbols (default 8 is good for 
       presentations, use 4 to more clearly see error bars)
     ylog : bool
       plot the y axix on log scale if True (default=False)
     errorband : bool
       if True will plot the errors as a shaded band rather than bars 
       (passed directly to fancy_plot). default=False
     line : string
       format of the line to connect data points. same as in pyplot.plot 
       and fancy_plot. default = '' (no line)
     linewidth : numeric
       width of the line to connect data points. same as in pyplot.plot 
       and fancy_plot. default = 1
     linecolor : string
       color of the line to connect data points. same as in pyplot.plot 
       and fancy_plot. default = None
     figsize : 2-element numerical tuple
       tuple of the form (5,3) to specify the size of the plot in 
       inches (x,y)
     mews : numerical
       passed to fancy_plot to set the markeredge width
     mecs : string
       passed to fancy_plot to set the markeredge color
     alphas : numerical
       passed directly to fancy_plot to set the symbol alpha
     mealphas : numerical
       passed directly to fancy_plot to set the markeredge alpha. default
       is to be the same as the marker (face) alpha
     Notes
     - In future, may add ability to also specify plotting a second 
       column as y values, optionally with a different color and symbols
     - To save the plot to a pdf file, initialize the pdffile before calling
       this function, e.g.:

       from matplotlib.backends.backend_pdf import PdfPages
       pdffile = PdfPages(plotfile)
       
       and close it afterwards:
       pdffile.savefig()
       plt.close()
       pdffile.close()

     """
     
     #----Set up symbols----
     if gratings is not None:
        df['symbols'] = 'd'
        # change according to grating
        df.loc[(df[gratings].values=='HETG'),'symbols']='o'
        df.loc[(df[gratings].values=='LETG'),'symbols']='*'
            
     #----Set up marker edge colors and widths----
     if mews is None:
         df['mews'] = 0.5
     else:
         df['mews'] = mews
     if mecs is None:
         df['mecs'] = color
     else:
         df['mecs'] = mecs
     if mealphas is None:
         mealphas = alphas
#     if simoffsets is not None:
#          # make sure column is strings
#          df[simoffsets] = pd.Series([str(off) for off in df[simoffsets]])
#          df['mecs'][df[simoffsets] == '-8.42'] = 'gray'
#     mecs = df['mecs'].values

     #----Set plot fonts----
     set_plot_fonts()

     #----Initialize Plot----
     if overplot is False:
        if ytitle is None:
           ytitle = y # label it with the column name
        fig,ax1 = start_plot(xtitle = 'SN1987A Age [days]',ytitle=ytitle,
                         xmin=agemin,xmax=agemax,ymin=ymin,ymax=ymax,
                         ylog=ylog,clean=clean,figsize=figsize)
     else:
        ax1 = ax

     #----Drop rows with missing y data----
     df = df.dropna(subset=[y])

     if gratings is not None:
         sym = df['symbols'].values

     #----Convert Error Intervals to Error Bars----
     if errinterval is True:
        if yerrlow is not None:
            yerrlow = df[y] - df[yerrlow]
        if yerrhigh is not None: 
            yerrhigh = df[yerrhigh] - df[y]
     else:
        if yerrlow is not None:
            yerrlow = df[yerrlow]
        if yerrhigh is not None:
            yerrhigh = df[yerrhigh]

     #----Plot the data----
     fancy_plot(ax1,df[ages],df[y],yerror=yerrlow,yerror_high=yerrhigh,
                syms=sym,mecs=df['mecs'].values,mews=df['mews'].values,
                line=line,linewidth=linewidth,errorband=errorband,
                colors=color,sizes=symsize,linealpha=1.0,bandalpha=0.1,
                alphas=alphas,ealphas=mealphas,mealphas=mealphas)


     #----Plot Extras----
     if plotextras is not None:
          plot_extras(names=plotextras)

     #----Plot Legend----
     # if overplotting, then need to plot the legend separately after all
     # calls to standard_age_plot
#     plot_legend(labels=leglabels,colors=legcolors,markers=legsyms,mecs=legmecs,
#                mews=legmews,fontsize='small')
     
     #----Plot Top Year Axis---
     if overplot is False:
          top_year_axis(ax1,agemin=agemin,agemax=agemax)

     #----Return----
     if overplot is False:
        return fig,ax1
     else:
        return None

