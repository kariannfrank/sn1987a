#Module of functions for fitting sn1987a data
#
#Contains the following functions:
#
# likelihood
# arr_brokenlinear
# brokenlinear
# arr_linear
# linear
# mcmc
#
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#  Import Common Modules

import numpy as np

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#
# likelihood()
#
#Author: Kari A. Frank
#Date: June 19, 2016
#Purpose: Calculate the likelihood of the data given theh parameters.
#Usage: result = likelihood(x,y,sigma,params)
#
#Input:
#  
# x = x values of measured data (numpy float array)
# y = y values of measured data (numpy float array)
# sigma = error values of measured data (numpy float array)
# params = broken lineaer model parameters, in format [junk,a1,a2,b1,tb]
#          (numpy float array)
#
#Usage Notes
# - x,y,sigma must have same number of elements
# - assumes that the errors are gaussian and symmetric.  
#   switching to asymmetric errors, e.g. by using piecewise 
#   truncated gaussians, would be done by changing the likelihood 
#   function, but has not been implemented.

def likelihood(x,y,sigma,params,func='brokenlinear'):

    #-get predicted radius values for each observation-
    if func == 'brokenlinear':
        F = arr_brokenlinear(x,params)
    if func == 'linear':
        F = arr_linear(x,params)

    #-calculate chi2-
    chi2 = np.sum( ((F-y)/sigma)**2.0 )
    
    #-calculate likelihood-
    N = y.shape[0] #number of observations
    Li = (2.0*np.pi)**(-N/2.0)*np.exp(-1.0*chi2/2.0)
    
    return Li
    
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# 
# arr_brokenlinear()
# brokenlinear()
#
# Model (F_K) Function

# takes as input 1D numpy arrays
# params must be in form [junk,a1,a2,b1,tb]

def arr_brokenlinear(xarr,params):

    Farr = np.array([brokenlinear(i,params) for i in xarr])
    
    return Farr 


def brokenlinear(x,params):

#    breakpoint = (params[4]-params[3])/(params[1]-params[2])

    #-for a2=change in slope-
    #a2 = params[1]+params[2]
    #-for a2=second slope-
    a2 = params[2]

    #-determine second intercept-
    b2 = params[3]+(params[1]-a2)*params[4]

   
    if x < params[4]:
        value = params[1]*x+params[3]
    else:
        value = a2*x+b2

    #-test with non-broken line (ignores a2 and tb)-
    #value = params[1]*x+params[3]
    
    return value

#--------------------------------------------------------------------------
# 
# arr_linear()
# linear()
#
# Model (F_K) Function

# takes as input 1D numpy arrays
# params must be in form [junk,a1,b1] (a1=slope,b1=intercept)

def arr_linear(xarr,params):

    Farr = np.array([linear(i,params) for i in xarr])
    
    return Farr 


def linear(x,params):

    #-slope-
    a1 = params[1]

    #-intercept-
    b1 = params[2]

    value = a1*x+b1
  
    return value


#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#
# mcmc()
#
#Author: Kari A. Frank
#Date: June 19, 2016
#Purpose: Fit the given x,y data to a broken linear model using MCMC
#Usage: result = mcmc(x,y,sigma,limits,initial)
#
#Input:
#  
# x = x values of measured data (numpy float array)
# y = y values of measured data (numpy float array)
# sigma = error values of measured data (numpy float array)
# limits = upper and lower limits for model 
#          parameters, in format [[0,0],[a1_low,a1_high],[a2_low,a2_high],...]
#          (numpy float array). can also use different model, array must still 
#          have the same dimensions (extra elements will be ignored). for 
#          a given parameter, setting low=high will cause it to be treated as
#          a fixed parameter.
#          
# initial = array of initial values for each model parameter, used as first row 
#           in the output table,e.g. [0.0,a1_0,a2_0,b1_0,tb_0] 
#           (numpy float array)
# func = string containing name of the model function to use.  must be one of 
#        the functions defined below. default='brokenlinear' 
#
#Output:
#
# numpy array containing the model parameters for every iteration
#
#Usage Notes
# - x,y,sigma must have same number of elements
# - first column of every array is iteration number

def mcmc(x,y,sigma,limits,initial,niterations=1000,func='brokenlinear'):

    import random as rd

    #----Set Up Main Array and Initial Values---

    #-get number of used (including fixed) parameters-
    used = 0
    for par in limits:
#        if par != [0,0]:
        if cmp(list(par),[0.0,0.0]) != 0:
            used = used + 1

    #-make sure limits is numpy array-
    limits = np.array(limits)

    #-set up arrays-
    model_params = np.zeros((niterations,used+1))
    #(first column is iteration number)

    model_params[:,0] = range(niterations)

    for par in xrange(used):
        model_params[0,par+1] = initial[par+1]

    #model_params[0,1] = initial[1]
    #model_params[0,2] = initial[2]
    #model_params[0,3] = initial[3]
    #model_params[0,4] = initial[4]

    #-1D 5-element array to keep the most up-to-date parameters-
    #(first element keeps the current value of i)
    current_params = model_params[0,:]
    current_L = likelihood(x=x,y=y,sigma=sigma,params=current_params,func=func)

    #---------------------------------------
    #              Run MCMC
    #---------------------------------------

    #----Main Iteration Loop----

    for i in model_params[1:,0]: #--skip 0th iteration--

    #-note: i is the iteration that is in the process of being updated
   
        current_params[0] = i

        if i % 100000 == 0: #print progress
            print "Iteration "+str(int(i))+" of "+str(niterations)

        #--choose order of parameters by shuffling column indices--
        #order = [1,2,3,4]
        order = [p+1 for p in xrange(used)]
        rd.shuffle(order)
    
        #--loop through parameters to get new value--
        for par in order:

            #--check if parameter is fixed--
            if limits[par,1] != limits[par,0]:

                #--choose candidate value from uniform distribution--   
                candidate = rd.random()*(limits[par,1]-limits[par,0]) + limits[par,0]

                #--calculate likelihoods--
                candidate_params = np.copy(current_params)            
                candidate_params[par] = candidate
                candidate_L = likelihood(x,y,sigma,candidate_params,func=func)

                #--accept or reject candidate--
                # (Metropolis-Hastings)

                if candidate_L >= current_L: #accept
                    current_params[par] = candidate
                    current_L = candidate_L
                    model_params[i,par] = candidate
                else: 
                    #-check for random acceptance-
                    if rd.random() < candidate_L/current_L: #accept
                        current_params[par] = candidate
                        current_L = candidate_L
                        model_params[i,par] = candidate
                    else: #reject -> repeat previous value in chain
                        model_params[i,par] = model_params[i-1,par]

            #-end if not fixed loop-
            
            else: #-parameter is fixed, repeat intial value-
                model_params[i,par] = model_params[i-1,par]

        #-end for par loop-

    #-end for i loop


    #-----Retrun Final Array----

    return model_params

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#
# plot_posteriors()
#
#Author: Kari A. Frank
#Date: July 17, 2014
#Purpose: Plot the priors output from mcmc()
#Usage: plot_posteriors(model_params[burnin:,:],plotfile)
#
#Input:
#  
# model_params = numpy array as output from mcmc(), optionally 
#                with the burnin rows excluded
#
# plotfile = already open file which will be plotted to
#
# nbins = number of bins in the posterior plots. default=30
#
# names = optional string list of parameter names, to be used as x-axis titles
#
#Output:
#
# plots (to already open file) of the posterior of each fit parameter
#
#Usage Notes

def plot_posteriors(model_params,plotfile,nbins=30,names=None):

    import matplotlib
    import matplotlib.pyplot as plt

    #----get array dimensions (number of parameters)----
    npar = model_params.shape[1]

    #----loop through parameters----
    for p in xrange(npar):
        
        #--skip iteration column or if parameter is fixed--
        if (max(model_params[:,p]) != min(model_params[:,p])) and p != 0:
            
            y = model_params[:,p]        
            y_hist, x_bin = np.histogram(y,bins=nbins)
            fig,ax=plt.subplots()
            plt.bar(x_bin[:-1],y_hist,width=x_bin[1]-x_bin[0])
            if names != None:
                plt.xlabel(names[p])
            ymin,ymax = ax.get_ylim()

            #-plot median and 1sigma range-
            med = np.median(y)
            sig = np.std(y)
            plt.plot([med,med],[ymin,ymax],color='red',linestyle='-')
            plt.fill_betweenx([0.0,ymax],[med-sig,med-sig],[med+sig,med+sig],color='red',alpha=0.2)

            #-save plot-
            plotfile.savefig()
            plt.close()


#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#
# mcmc_align_img()
#
#Author: Kari A. Frank
#Date: November 6, 2014
#Purpose: Fit the given 2-dimensional (e.g. image) data to another 
#Usage: result = mcmc_align_img(img,refimg,parmeans,parsigmas,niterations=1000)
#
#Input:
#  
# img = 2d array of that is to be translated to match refimg
# refimg = 2d array (of same dimensions as img) used as reference image
# parmeans = [dx0,dy0] center of the prior for each parameter 
#            (numpy float array)
# parsigmas = [sigma_dx0,sigma_dy0] standard deviation of each parameter
#             prior. setting an element equal to zero will cause the 
#             parameter to be fixed at the given mean. (numpy float array)
# niterations = number of iterations to perform (integer)
#
#Output:
#
# numpy array containing the model parameters for every iteration
#
#Usage Notes
# - img and refimg must have the same dimensions
# - assumes img and refimg are normalized to 1 
# - assumes gaussian priors (may add capability for uniform priors later)

def mcmc_align_img(img,refimg,parmeans,parsigmas,niterations=1000):

    import random as rd
    import numpy as np

    #----Set Up Main Array and Initial Values---

    #-get number of used (including fixed) parameters-
    used = 0
    for par in parsigmas:
        if par != 0.0:
            used = used + 1

    # convert parameters to integers
    parmeans = [int(p) for p in parmeans]

    #-set up arrays-
    model_params = np.zeros((niterations,used+1))
    #(first column is iteration number)
    model_params[:,0] = range(niterations)
    # set first guess = mean
    for par in xrange(used):
        model_params[0,par+1] = parmeans[par]

    #-1D array to keep the most up-to-date parameters-
    #(first element keeps the current value of i)
    current_params = np.copy(model_params[0,:])
    current_L = imgdiff_likelihood(img,refimg,params=current_params[1:,])

    #---------------------------------------
    #              Run MCMC
    #---------------------------------------

    #----Main Iteration Loop----

    for i in model_params[1:,0]: #--skip 0th iteration--

    # Note: i is the iteration that is in the process of being updated
   
        current_params[0] = i

        # print progress
        if i % 100000 == 0: 
            print "Iteration "+str(int(i))+" of "+str(niterations)
        
        #--choose order of parameters by shuffling column indices--
        #order = [1,2,3,4]
        order = [p for p in xrange(used)]
        rd.shuffle(order)
    
        #--loop through parameters to get new value--
        for par in order:

            #--check if parameter is fixed--
            if parsigmas[par] != 0:

                #--choose candidate value from gaussian distribution--   
                candidate = int(rd.gauss(parmeans[par],parsigmas[par]))

                #--calculate likelihoods--
                candidate_params = np.copy(current_params)            
                candidate_params[par+1] = candidate
                candidate_L = imgdiff_likelihood(img,refimg,params=candidate_params[1:])

                #--accept or reject candidate--
                # (Metropolis-Hastings)

                if candidate_L >= current_L: #accept
                    current_params[par+1] = candidate
                    current_L = candidate_L
                    model_params[i,par+1] = candidate
                else: 
                    #-check for random acceptance-
                    if rd.random() < candidate_L/current_L: #accept
                        current_params[par+1] = candidate
                        current_L = candidate_L
                        model_params[i,par+1] = candidate
                    else: #reject -> repeat previous value in chain
                        model_params[i,par+1] = model_params[i-1,par+1]

            #-end if not fixed loop-1
            
            else: #-parameter is fixed, repeat intial value-
                model_params[i,par+1] = model_params[i-1,par+1]

        #-end for par loop-

    #-end for i loop


    #-----Return Final Array----

    return model_params


#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#Author: Kari A. Frank
#Date: November 7, 2014
#Purpose: Calculate the likelihood comparing two images
#Usage: 
#
#Input:
#  
# img = 2d array of that is to be translated to match refimg
# refimg = 2d array (of same dimensions as img) used as reference image
# 
# params = numpy array or list of the transformation parameters [dx,dy]
# 
#Output:
#
# float value of the likelihood
#
#Usage Notes
# - img and refimg must have the same dimensions

def imgdiff_likelihood(img,refimg,params):

#    import ciao_contrib.runtool as crt 
    import numpy as np
    import pyfits as fits

    # translate img
    newimg = np.roll(img,int(params[0]),axis=0)
    newimg = np.roll(newimg,int(params[1]),axis=1)

    # create difference^2 image
    diffimg = (newimg - refimg)**2.0
    
    # calculate chi2
    chi2 = np.sum(diffimg)
    
    # calculate likelihood
    N = img.shape[0]*img.shape[1] #number of pixels (i.e. measurements)
    Li = (2.0*np.pi)**(-N/2.0)*np.exp(-1.0*chi2/2.0)

    return Li
