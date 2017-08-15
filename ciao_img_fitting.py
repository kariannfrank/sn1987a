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
    current_L = imgdiff_likelihood(img,refimg,params=current_params)

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
                candidate = rd.gauss(parmeans[par],parsigmas[par])

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

    import ciao_contrib.runtool as crt 
    import numpy as np
    import pyfits as fits

    # create temporary copy of img
    tempimg = img+'_trans'    
    crt.dmcopy.punlearn()
    crt.dmcopy.infile = img
    crt.dmcopy.outfile = tempimg
    crt.dmcopy.clobber = 'yes'
    crt.dmcopy()
    
    # shift image
    crt.wcs_update.punlearn()
    crt.wcs_update.infile = tempimg
    crt.wcs_update.wcsfile = tempimg
    crt.wcs_update.deltax = str(params[0])
    crt.wcs_update.deltay = str(params[1])
    crt.wcs_update.rotang = 0.0
    crt.wcs_update.scalefac = 1.0
    crt.wcs_update.clobber = 'yes'
    crt.wcs_update()

    # create difference^2 image

    # difference image
    diffimg = img+'_diff'
    crt.dmimgcalc.punlearn()
    crt.dmimgcalc.infile = tempimg
    crt.dmimgcalc.infile2 = refimg
    crt.dmimgcalc.outfile = diffimg
    crt.dmimgcalc.op = 'sub'
    crt.dmimgcalc.clobber = 'yes'
    crt.dmimgcalc()

    # square difference image
    diffsqimg = img+'_diffsq'
    crt.dmimgcalc.punlearn()
    crt.dmimgcalc.infile = diffimg
    crt.dmimgcalc.infile2 = 'none'
    crt.dmimgcalc.outfile = diffsqimg
    crt.dmimgcalc.op = 'imgout=img1*img1'
    crt.dmimgcalc.clobber = 'yes'
    crt.dmimgcalc()


    # calculate chi2

    # read in diff image as array and get total
    imgfits = fits.open(diffimg)
    img = imgfits[0].data
    imgfits.close()
    chi2 = np.sum(img)
    
    #-calculate likelihood-
    N = img.shape[0]*img.shape[1] #number of pixels (i.e. measurements)
    Li = (2.0*np.pi)**(-N/2.0)*np.exp(-1.0*chi2/2.0)

    return Li
