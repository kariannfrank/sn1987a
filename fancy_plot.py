def fancy_plot(ax,x,y,xerror=None,yerror=None,xerror_high=None,
               yerror_high=None,line='',linewidth=1,
               syms='o',colors='blue',sizes=4,alphas=1.0,mecs='black',
               mews=0.5,errorband=False,bandalpha=0.1,linealpha=1.0,
               ealphas=None,mealphas=None):
    """
    Function to make plots with different symbols/colors/sizes for each point

    Author: Kari A. Frank
    Date: July 19, 2014
    Purpose: Function to make plots with different symbols/colors/sizes 
             for each point.

    Input:

     x,y       -- 1D lists or numpy arrays of the x and y data to be 
                  plotted (floats or ints) [required]

     xerror,yerror -- 1D lists or numpy arrays of the errorbars on the x and 
                  y values. if the corresponding xerror_high/yerror_high are 
                  also provided, xerror/yerror will be treated as the lower 
                  errorbars [optional]

     xerror_high,yerror_high -- 1D lists or numpy arrays of the upper 
                  errorbars on the x and y values. if given, the 
                  corresponding x/yerror must also be provided. [optional]

     syms      -- list or single string of symbols (e.g. 'o','^').  if 
                  a single string is provided, all points will use the 
                  same symbol.  if a list, it must have the same number 
                  of elements as x and y, and will be used in the same 
                  order. [optional, default='o']

     colors    -- list or single string of colors (e.g. 'g' or 'green').  
                  works the same as syms. [optional, default='black']

     sizes     -- list or single int of symbols sizes.  works the same as 
                  syms. [optional,default=4]

     alphas    -- scalar or list of scalars specifying alpha (opacity) of the
                  symbols and errorbars. [optional, default=1] 

     mecs      -- scalar or list of scalars specifying markeredgecolor of the
                  symbols. [optional, default='black'] 

     mews      -- scalar or list of scalars specifying markeredgewidth of the
                  symbols. [optional, default=0.5] 

     errorband -- switch to plot the errors as a shaded band rather
                  than errorbars. [optional, default=False]

     bandalpha -- optionally specify alpha for the errorband

     line      -- string of format for connecting line 
                  [optional, default='' (no line)]
    
     linewidth -- float of the linewidth [optional, default=1]

     bandalpha -- optionally specify alpha for the line

     ealphas    -- optionally specify a different alpha for the errorbars;
                  this allows plotting of solid errorbars with translucent
                  markers. [ optional, default=mealphas]

     mealphas   -- optionally specify a different alpha for the marker edges.
                   [optional, default = alphas]

    Output:
     - plot to an already initialized figure (to allow overplotting)

    Usage Notes:
     - can add further functionality of the plot and errorbar functions as
       needed (e.g. allowing errorbars to be different color than the point,
       changing edgecolor, allowing upper/lower limits)
     - errorbar color and alpha will be the same as the markeredgecolor 
       (allows for unfilled 'white' points to have visible errorbars.)
    """

#---------------------------------------
#           Import Modules
#---------------------------------------

    import matplotlib
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

#---------------------------------------
#   Parse arguments and set defaults
#---------------------------------------

#---Set to Default Values---
    if syms is None:
        syms = ['o' for i in x]
    if colors is None:
        colors = ['black' for i in x]
    if sizes is None:
        sizes = [4 for i in x]
    if mealphas is None:
        mealphas = alphas
    if ealphas == None:
        ealphas = mealphas

#---Convert Scalars to Lists---

    if isinstance(syms,str):
        syms = [syms]*len(x)
    if isinstance(colors,str):
        colors = [colors]*len(x)
    if isinstance(sizes,list) == False:
        sizes = [sizes]*len(x)
    if isinstance(alphas,int):
        alphas = float(alphas)
    if isinstance(alphas,float):
        alphas = [alphas]*len(x)
    if isinstance(ealphas,int):
        ealphas = float(ealphas)
    if isinstance(ealphas,float):
        ealphas = [ealphas]*len(x)
    if isinstance(mecs,str):
        mecs = [mecs]*len(x)
    if isinstance(mews,float):
        mews = [mews]*len(x)
    if isinstance(mealphas,float):
        mealphas = [mealphas]*len(x)

    mecs = [matplotlib.colors.colorConverter.to_rgba(mecs[i],alpha = mealphas[i]) for i in range(len(mecs))]

#---Convert to numpy arrays---
    if yerror is not None:
        yerror = np.array(yerror)
    if yerror_high is not None:
        yerror_high = np.array(yerror_high)
    if xerror is not None:
        xerror = np.array(xerror)
    if xerror_high is not None:
        xerror_high = np.array(xerror_high)

    x = np.array(x)
    y = np.array(y)

#---------------------------------------
#            Plot Values
#---------------------------------------

#---Set up errors---

    if xerror_high is None:
        xerror_high = xerror
    if yerror_high is None:
        yerror_high = yerror

    if xerror is not None:
        xerr = [np.array([xerror[i],xerror_high[i]]).reshape(2,1) for i in
                range(len(x))]
        xerr = np.array(xerr)
    else:
        xerr = np.array([None]*len(x))

    if yerror is not None:
        yerr = [np.array([yerror[i],yerror_high[i]]).reshape(2,1) for i in 
                range(len(x))]
        yerr = np.array(yerr)
    else:
        yerr = np.array([None]*len(y))

#---If connecting line, plot data with no symbols to create line---

    if line != '':

        #-set line color-
        if all(col == mecs[0] for col in mecs) == True:
            lcolor = mecs[0] #default = markeredgecolor if all same
        elif all(col == colors[0] for col in colors) == True:
            lcolor = colors[0] #second choice = marker color is all same
        else:
            lcolor = 'black' #else all black

        #-set line alpha-
        if linealpha is None:
            if all(al == alphas[0] for al in alphas) == True:
                linealpha = alphas[0]
            else:
                linealpha = max(alphas)
    
        #-plot data with no symbols (for line only)-
        ax.plot(x,y,line,color=lcolor,alpha=linealpha,linewidth=linewidth)
        
#---Iterate over each point----
    for i in range(len(x)):
        if errorband == False:            
            if ealphas[i] == alphas[i]:
                ax.errorbar(x[i],y[i],yerr=yerr[i],fmt=syms[i],
                            color=colors[i],
                            alpha=alphas[i],markersize=sizes[i],
                            markeredgewidth
                            =mews[i],markeredgecolor=mecs[i],ecolor=mecs[i])
            else:
                # plot errorbars only (set marker color to white)
                ax.errorbar(x[i],y[i],yerr=yerr[i],fmt=syms[i],
                            color='white',
                            alpha=ealphas[i],markersize=sizes[i],
                            markeredgewidth
                            =mews[i],markeredgecolor=mecs[i],ecolor=mecs[i])
                # plot markers only
                ax.plot(x[i],y[i],marker=syms[i],color=colors[i],
                        alpha=alphas[i],markersize=sizes[i],
                        markeredgewidth
                        =mews[i],markeredgecolor=mecs[i])
        else:
            ax.plot(x[i],y[i],marker=syms[i],color=colors[i],alpha=
                    alphas[i],markersize=sizes[i],markeredgewidth=mews[i],
                    markeredgecolor=mecs[i])

#---Plot Errorband if Set---
    if errorband == True:
        if yerror is None:
            print 'ERROR: No errors provided, cannot plot errorband.'
        if all(col == colors[0] for col in colors) == False:
            print ("Warning: points have different colors. Errorband "
                   "will be plotted gray.")
            ercolor = 'black'
        else:
            ercolor = colors[0]

#        plt.fill_between(x,y-yerr[:,0,0],y+yerr[:,0,1],color=ercolor,
#                                     alpha=bandalpha)

        ax.fill_between(x,y-yerr[0,0,:],y+yerr[0,1,:],color=ercolor,
                                     alpha=bandalpha)


    return None
########## END ##########

