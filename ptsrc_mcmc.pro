;; Name: ptsrc_mcmc.pro
;;
;; Date: February 5, 2015
;; Author: Kari A. Frank
;;        
;; Purpose:
;;   Run a MCMC to determine upper limit on potential central point
;;   source flux (in counts) given an input SN 1987A image.
;;   
;; 
;; Calling Sequence:
;;   ptsrc_mcmc,infile=infile,niter=niter,burnin=burnin,log=log,verbose=verbose,priortype=priortype,ltype=ltype,convolve=convolve,psffile=psffile,xystep=2.0,xysigma=0.5
;;
;; Parameters:
;;
;; inimage -- input 'raw' (not deconvolved or smoothed) image
;;            file. to properly add convolved point source, this 
;;            needs to have same binning as psf.fits image (e.g.
;;            as output by deconvolve.pro in the usual SN1987A 
;;            pipeline, typically name is <evtfilename>.image)
;; niter -- number of iterations to perform
;; burnin -- number of iterations to skip in final calculated values
;; log -- switch to use log likelihoods (default=1)
;; verbose -- set verbosity (0=print nothing, 1=print progress,
;;            2=print more stuff. Default=1).
;; priortype -- set type of distribution for counts prior. must be
;;              'gaussian', 'poisson', or 'uniform'. Default is 'uniform'.
;; ltype -- set type of likelihood. must be 'chi2' or
;;          'poisson'. Default is 'chi2'
;; psffile -- name of the file to use for psf. not required
;;            if convolve='no' (default='psf.fits')
;; convolve -- optionally specify that the input image is
;;             deconvolved, and therefore the point source
;;             should not be convolved with the psf prior 
;;             to adding (default=='yes')
;;
;;Notes:
;; - Must modify priors manually
;;   (Currently assume the point source position and width are fixed). 
;; - poisson prior will run, but seems to accept almost everything
;;   input_img -- subpixel SN 1987A image, binned to 0.125" pixels
;;
;; Example:
;;   
;;


;-------------------------------------------
;Function to Calculate Likelihood for Two Images
;-------------------------------------------
FUNCTION likelihood,data,model,ltype=ltype,log=log


IF NOT KEYWORD_SET(ltype) THEN ltype='chi2'
IF NOT KEYWORD_SET(log) THEN log=0

nonzero = WHERE(data NE 0.0)

IF ltype EQ 'chi2' THEN BEGIN
;; chi2 likelihood
;; errors assumed to be poisson = sqrt(N)
   diffimg = ((data[nonzero] - model[nonzero])/(data[nonzero])^0.5)^2.0

   chi2 = total(diffimg)
   ;print, 'chi2 = ',chi2

   shape = SIZE(diffimg)
;N = shape[1]*shape[2] ;; number of pixels (data points)
;lnLi = (-1.0*N/2.0)*ALOG(2.0*!PI)-chi2/2.0
;Li = (2.0*!PI)^(-1.0*N/2.0)*exp(-1.0*chi2/2.0)
   Li = exp(-1.0*chi2/2.0)
   logLi = -1.0*chi2/2.0
;PRINT, "Li = ",Li
ENDIF
IF ltype EQ 'poisson' THEN BEGIN
   logLi = TOTAL(model[nonzero]^data[nonzero]*EXP(-1.0*model[nonzero])/data[nonzero])
   Li = EXP(logLi)
   PRINT, 'Li = ',Li
ENDIF

;el = (data[nonzero]^model[nonzero])*EXP(data[nonzero])/FACTORIAL(data[nonzero])
;PRINT, "max,min L_k = ",max(el),min(el)

IF log EQ 0 THEN RETURN, Li ELSE RETURN, logLi

END

;-------------------------------------------
;  Function to add point source to an image
;  (optionally convolving it with psf)
;-------------------------------------------
FUNCTION add_ptsrc,img,cts,x,y,width=width,psfimg=psfimg

;-- default width is true point source (<= 1 pixel) --
IF NOT KEYWORD_SET(width) THEN width = 0.5

;-- set minimum allowed pixel value (for zero for very small numbers)
thresh = 10.0^(-7.0)
  
;-- get image size information --
shape = SIZE(img)

;-- create ptsrc-only image array --
;print, 'shape[1],width,x,y = ',shape[1],width,x,y
ptimg = psf_gaussian(NPIXEL=shape[1],NDIMENSION=2,FWHM=width*2.0,CENTROID=[x,y]);,/NORMALIZE)

;-- convolve ptsrc with psf --
IF KEYWORD_SET(psfimg) NE 0 THEN BEGIN
   ;print, 'convolving'
   ptimg_conv = DOUBLE(convolve(double(ptimg),double(psfimg),FT_PSF=psf_ft))
;PRINT, "max ptimg_conv = ",max(ptimg_conv)
;PRINT, "total ptimg_conv = ",total(ptimg_conv)
   ptimg = ptimg_conv
ENDIF 

ptimg[WHERE(ptimg LE thresh)] = 0.0
ptimg = cts*ptimg/total(ptimg)
ptimg[WHERE(ptimg LE thresh)] = 0.0

;-- add (convolved) ptsrc to image --
RETURN, img+ptimg;_conv
;RETURN, img+psfimg

END

;-------------------------------------------
;Function to find confidence interval
;-------------------------------------------
FUNCTION confidence_interval,param,level=level

  ;; default is 90% confidence interval
  IF KEYWORD_SET(level) EQ 0 THEN level = 0.9

;  truncated_params = FULL_SIGRANGE(param,FRACTION=level,RANGE=interval)
  stats = credible_interval(param,level=level,binsize=1)
  

  RETURN, stats

END

;-------------------------------------------
;Function to randomly shuffle an array
;-------------------------------------------
FUNCTION shuffle,X,seed=seed

;  ind = INDGEN(npar)
  rand = RANDOMU(seed,N_ELEMENTS(X))

  RETURN, X[SORT(rand)]

END

;-------------------------------------------
;Function to Plot Posterior to Open File
;-------------------------------------------
FUNCTION plot_posterior,param,binsize=binsize,parname=parname,prior=prior

;; prior is optional vector specifying the prior parameters
;;    prior = [ptype,par1,par2]
;;    ptype = 0,1 (0=uniform, 1=poisson)
;;    for poisson, par1=mean
;;    for uniform, par1=min,par2=max
;;    for gaussian, par1=mean,par2=sigma

IF KEYWORD_SET(parname) EQ 0 THEN parname = ''

;nbins=max(param)/20.0
nbins = 40.0
binsize=(max(param)-min(param))/nbins
;binsize = 1
;nbins = max(param)-min(param)

;; note that autobin is ignored if binsize is defined
;;PLOTHIST,param,BIN=binsize,/AUTOBIN,
plotted = histogram_plot(param,XTITLE=parname,binsize=binsize)

; set up prior for overplotting
IF KEYWORD_SET(prior) THEN BEGIN
   CASE prior[0] OF
      0: BEGIN 
;         priorx = FINDGEN(prior[2]-prior[1])+prior[1]
         inc = (prior[2]-prior[1])/N_ELEMENTS(param)
         priorx = FINDGEN(N_ELEMENTS(param))/N_ELEMENTS(param)*prior[2]+prior[1]
         priorval = N_ELEMENTS(param)/nbins
         priory = REPLICATE(priorval,N_ELEMENTS(param))
         OPLOT,priorx,priory,LINESTYLE=2,COLOR=128
      END

      1: BEGIN
         priorarr = RANDOMU(seed,N_ELEMENTS(param),POISSON=prior[1])
         PLOTHIST,priorarr,/OVERPLOT,BIN=binsize,COLOR=128,LINESTYLE=2
         ;priorx = (FINDGEN(500)*(binsize*nbins)+min(param))/500.0
         ;priory = prior[1]^priorx*EXP(-1.0*prior[1])/FACTORIAL(priorx)
         ;OPLOT,priorx,priory,LINESTYLE=2,COLOR=128
      END

      ELSE: PRINT, "WARNING: Invalid prior number (must be 0 or 1). Not plotting prior." ;may want to implement gaussian in the future
   ENDCASE
ENDIF

stats = confidence_interval(param) ;assumes only one mode
ymax = 2.0*MAX(HISTOGRAM(param))
OPLOT,[stats[0],stats[0]],[0.0,ymax],THICK=2,LINESTYLE=3 ;plot mode
OPLOT,[stats[1],stats[1]],[0.0,ymax],THICK=2,LINESTYLE=1 ;plot lower 90% bound
OPLOT,[stats[2],stats[2]],[0.0,ymax],THICK=2,LINESTYLE=1 ;plot upper 90% bound

RETURN, 1

END

;-------------------------------------------
;       Function to Make Plots
;-------------------------------------------

FUNCTION mcmc_plots,plotfile,mparams,burnin=burnin,maxiter=maxiter,log=log,prior=prior

IF NOT KEYWORD_SET(burnin) THEN burnin=0.0
IF NOT KEYWORD_SET(log) THEN log = 0
IF NOT KEYWORD_SET(maxiter) THEN maxiter = N_ELEMENTS(mparams[*,0])-1
maxiter = ULONG(maxiter)
burnin = ULONG(burnin)

SET_PLOT,'ps'
DEVICE, /COLOR
DEVICE, FILENAME=plotfile

IF burnin LT maxiter THEN BEGIN
   countmedian=MEDIAN(mparams[burnin:maxiter,2])
   limedian=MEDIAN(mparams[burnin:maxiter,1])
   void = MAX(HISTOGRAM(mparams[burnin:maxiter,1],OMIN=mn),wmode)
   countmode = mn+wmode
ENDIF ELSE BEGIN
   countmedian=MEDIAN(mparams[1:maxiter,2])
   limedian=MEDIAN(mparams[1:maxiter,1])
   void = MAX(HISTOGRAM(mparams[1:maxiter,1],OMIN=mn),wmode)
   countmode = mn+wmode
;   countmode=MODE(mparams[1:maxiter,2])
ENDELSE

;-- counts prior --
IF burnin LT maxiter THEN plotted = plot_posterior(mparams[burnin:maxiter,2],parname='Counts',prior=*prior[2]) ELSE plotted = plot_posterior(mparams[1:maxiter,2],parname='Counts',prior=*prior[2])
;-- x prior --
IF burnin LT maxiter THEN plotted = plot_posterior(mparams[burnin:maxiter,3],parname='X',prior=*prior[3]) ELSE plotted = plot_posterior(mparams[1:maxiter,3],parname='Counts',prior=*prior[3])
;-- y prior --
IF burnin LT maxiter THEN plotted = plot_posterior(mparams[burnin:maxiter,4],parname='Y',prior=*prior[4]) ELSE plotted = plot_posterior(mparams[1:maxiter,4],parname='Counts',prior=*prior[4])

;-- plot counts and mode vs iteration --
PLOT, mparams[1:maxiter,0],mparams[1:maxiter,2],THICK=2,XTHICK=2,YTHICK=2,CHARTHICK=2,XTITLE='Iteration',YTITLE='Counts'
IF burnin LT maxiter THEN OPLOT,[burnin,burnin],[0.0,MAX(mparams[1:maxiter,2])],LINESTYLE=2
modes = FLTARR(maxiter+1)
FOR i=1,maxiter DO BEGIN
   void = MAX(HISTOGRAM(mparams[1:i,2],OMIN=mn),wmode)
   modes[i] = mn+wmode
ENDFOR
PLOT,mparams[1:maxiter,0],modes,THICK=2,XTHICK=2,YTHICK=2,CHARTHICK=2,XTITLE='Iteration',YTITLE='Mode'
OPLOT,mparams[1:maxiter,0],[countmode,countmode],THICK=2,LINESTYLE=0

;-- plot likelihood vs iteration --
IF log EQ 1 THEN BEGIN
   lit='LogLikelihood' 
   nsig = [3,0.1]
ENDIF ELSE BEGIN
   lit='Likelihood'   
   nsig = [3,3]
ENDELSE
ymin = MEDIAN(mparams[1:maxiter,1]) - nsig[0]*STDDEV(mparams[1:maxiter,1])
ymax = MEDIAN(mparams[1:maxiter,1]) + nsig[1]*STDDEV(mparams[1:maxiter,1])
PLOT, mparams[1:maxiter,0],mparams[1:maxiter,1],THICK=2,XTHICK=2,YTHICK=2,CHARTHICK=2,XTITLE='Iteration',YTITLE=lit,YRANGE=[ymin,ymax]
IF burnin LT maxiter THEN  OPLOT,[burnin,burnin],[0.0,MAX(mparams[1:maxiter,1])],LINESTYLE=2
OPLOT,[0,maxiter],[limedian,limedian],LINESTYLE=3

;-- plot counts vs likelihood --
PLOT, mparams[1:maxiter,2],mparams[1:maxiter,1],THICK=2,XTHICK=2,YTHICK=2,CHARTHICK=2,XTITLE='Counts',YTITLE=lit,PSYM=3

DEVICE, /CLOSE
SET_PLOT,'X'

RETURN, 1
END

;-------------------------------------------
;-------------------------------------------
;           MAIN FUNCTION
;-------------------------------------------
;-------------------------------------------



PRO ptsrc_mcmc,infile=infile,niter=niter,burnin=burnin,log=log,verbose=verbose,priortype=priortype,ltype=ltype,outsuffix=outsuffix,psffile=psffile,convolve=convolve,xystep=xystep,xysigma=xysigma,sigmastep=sigmastep;input image file, psf image file, niter

;-------------------------------------------
;     Parse Arguments and Set Defaults
;-------------------------------------------

;IF N_PARAMS() LT 1 THEN MESSAGE, "ERROR: Minimum usage is 'process, input_files'"
IF NOT KEYWORD_SET(niter) THEN niter = 500 ;number of mcmc iterations
IF NOT KEYWORD_SET(burnin) THEN burnin = 100
IF N_ELEMENTS(verbose) EQ 0 THEN verbose=1
IF NOT KEYWORD_SET(priortype) THEN priortype = 'uniform'
IF NOT KEYWORD_SET(ltype) THEN ltype = 'chi2'
IF NOT KEYWORD_SET(log) THEN log=1
IF N_ELEMENTS(outsuffix) EQ 0 THEN outsuffix='' ELSE outsuffix = '_'+outsuffix
IF N_ELEMENTS(psffile) EQ 0 THEN psffile = 'psf.fits'
IF N_ELEMENTS(convolve) EQ 0 THEN convolve = 'yes'
IF NOT KEYWORD_SET(xystep) THEN xystep = 2.0
IF NOT KEYWORD_SET(xysigma) THEN xysigma = 0.5
IF NOT KEYWORD_SET(sigmastep) THEN sigmastep = 0.0
niter = ULONG(niter)
burnin = ULONG(burnin)

;-------------------------------------------
;        Set Up Files and Images
;-------------------------------------------

;; eventually make these arguments (except model_img)
;---- Input Image ----
;in_file ='acisf15810_repro_evt2_subpix_1200-8000.fits.image'
;in_file = 'acisf15810_repro_evt2_subpix_1200-8000_20countptsrc.fits.image'
fit_file = infile+outsuffix+'.ptsrcfit'
in_img = READFITS(infile,in_head)

thresh = 10.0^(-5.0)
in_img[WHERE(in_img LE thresh)] = 0.0

;---- PSF File ----
;psffile = 'psf.fits'
IF convolve EQ 'yes' THEN psf_img = READFITS(psffile,psf_head)

;---- Posterior Plot File ----
post_plot_file = fit_file+'_posteriors'+outsuffix+'.ps'

;---- Chain Parameter File ----
chain_file = fit_file+'_chains'+outsuffix+'.txt'

;---- Create Model Image Array ----
dims = SIZE(in_img)
model_img = FLTARR(dims[1],dims[2])

;-------------------------------------------
;             Set Up MCMC
;-------------------------------------------

;---- parameter priors ---- 

;-- prior for counts --
count_min = 0
count_max = 0.2*total(in_img);4.0*max(in_img)
count_median = MEDIAN(in_img[WHERE(in_img GT 0.0)])
count_avg = MEAN(in_img[WHERE(in_img GT 0.0)])
count_stdev = STDDEV(in_img[WHERE(in_img GT 0.0)])
IF verbose GT 0 THEN PRINT, 'count_max,median,total = ',count_max,count_median,total(in_img)

CASE priortype OF
   'uniform': countsprior = [0,count_min,count_max]
   'poisson': countsprior = [1,4.0*count_median,0.0]
   'gaussian': countsprior = [2,count_avg,count_stdev]
ENDCASE

;-- gaussian prior for position --
x0 = dims[1]/2-1.0 ;; center of image
y0 = dims[2]/2-1.0
IF verbose GE 1.0 THEN PRINT, "x0,y0 initial = ",x0,y0
xstep = xystep
ystep = xystep
;xysigma = 0.5 ;; in pixels
;xysigma = 0.0 ;; in pixels
;xysigmastep = 0.0

CASE priortype OF
   'uniform': xprior = [0,x0-xstep,x0+xstep]
   'poisson': xprior = [1,4.0*x0,x0]
   'gaussian': xprior = [2,x0,xstep]
ENDCASE
CASE priortype OF
   'uniform': yprior = [0,y0-ystep,y0+ystep]
   'poisson': yprior = [1,4.0*y0,y0]
   'gaussian': yprior = [2,y0,ystep]
ENDCASE
CASE priortype OF
   'uniform': xysigmaprior = [0,xysigma-sigmastep,xysigma+sigmastep]
   'poisson': xysigmaprior = [1,4.0*sigmastep,xysigma]
   'gaussian': xysigmaprior = [2,xysigma,sigmastep]
ENDCASE

;; first two pointers are just dummies to make the parameter indices
;; match those in model_params
prior = [PTR_NEW(/Allocate_Heap),PTR_NEW(/Allocate_Heap),PTR_NEW([countsprior],/NO_COPY),PTR_NEW([xprior],/NO_COPY),PTR_NEW([yprior],/NO_COPY),PTR_NEW([xysigmaprior],/NO_COPY)]
IF sigmastep EQ 0.0 THEN free = [2,3,4] ELSE free = [2,3,4,5] ;; edit this manually -- prior indices which correspond to free parameters
;IF verbose GT 0 THEN PRINT, 'prior = ',*prior

;---- parameter array ----
;; model parameter columns are iteration#, likelihood, counts, x, y, sigma
model_params = FLTARR(niter,6)
model_params[*,0] = FINDGEN(niter) ;; fill iteration numbers
model_params[*,3] = x0 ;; fill x and y (for fixed x,y)
model_params[*,4] = y0
model_params[*,5] = xysigma

current_params = model_params[0,*] ;; 1D array to hold the most 
                                   ;;  up-to-date parameters
IF convolve EQ 'yes' THEN BEGIN
   model_img = add_ptsrc(in_img,count_max,x0,y0,width=xysigma,psfimg=psf_img) ;; create 0th model image
ENDIF ELSE BEGIN
   model_img = add_ptsrc(in_img,count_max,x0,y0,width=xysigma) ;; create 0th model image
ENDELSE
;writefits,fit_file,model_img,in_head
current_L = likelihood(in_img,model_img,log=log,ltype=ltype) ;; hold current likelihood
current_params[1] = current_L
current_params[2] = count_max
current_params[3] = x0
current_params[4] = y0
current_params[5] = xysigma

naccepted=0.0
nrejected=0.0

print, 'starting mcmc loop'
;-------------------------------------------
;             MCMC Loop
;-------------------------------------------
;mini=100000UL
;FOR i=mini,niter-1 DO BEGIN
FOR i=1UL,niter-1 DO BEGIN

   ;; note: i is the iteration that is in the process of being updated
   current_params[0] = i

   ;-- choose order of parameters by shuffling column indices --
;; not implemented (only matters if more than one free parameter)
   IF N_ELEMENTS(free) GT 1 THEN BEGIN
;      order = INDGEN(free) ;; array to hold parameter (prior row) indices
      order = shuffle(free) ;; array to hold parameter (prior row) indices
   ENDIF ELSE BEGIN
      order = [2]
   ENDELSE 

   ;-- loop through parameters to get new value --
   FOR pi=0,N_ELEMENTS(free)-1 DO BEGIN
      p = order[pi];+2 ;skip first 2 model_params columns (iterations,likelihood)
      piprior = *prior[p]

      ;-- choose candidate value from distribution --
;;      IF (p EQ 2) OR (p EQ 3) THEN gauss = 1
      CASE priortype OF
         'uniform': candidate = (piprior[2]-piprior[1])*RANDOMU(seed,1) + piprior[1]
         'poisson': candidate = RANDOMU(seed,1,POISSON=piprior[1])
         'gaussian': candidate = RANDOMN(seed,1)*piprior[2]+piprior[1]
      ENDCASE
      if (p eq 2) and (verbose GT 1) then print, 'candidate = ',candidate
;      IF (p EQ 2) AND (priortype EQ 'uniform') THEN candidate = FLOOR(candidate)

      ;-- calculate likelihood --
      candidate_params = current_params
      candidate_params[p] = candidate ;update only parameter p
;      PRINT, "candidate_params = ",candidate_params
      cand_model_img = add_ptsrc(in_img,candidate_params[2], $
                                 candidate_params[3],candidate_params[4], $
                                 width=candidate_params[5],psfimg=psf_img)

;      PRINT, "Max,Min model_img = ",max(cand_model_img),min(cand_model_img)
      candidate_L = likelihood(in_img,cand_model_img,log=log,ltype=ltype)
;      IF verbose GT 2 THEN PRINT, "candidate_L = ",candidate_L

      ;-- accept or reject candidate --
      ;; Metropolis-Hastings
      IF candidate_L GE current_L THEN BEGIN ;; accept
         IF verbose GT 1 THEN PRINT, 'Accepted A: ',i,candidate_L,current_L
         current_params[p] = candidate
         current_L = candidate_L
         model_params[i,p] = candidate
         model_params[i,1] = candidate_L
         naccepted=naccepted+1
      ENDIF ELSE BEGIN
         ;- check for random acceptance - 
         IF log EQ 0 THEN ratio = candidate_L/current_L ELSE ratio = EXP(candidate_L-current_L)
         rnd = RANDOMU(seed,1)
         IF rnd LT ratio THEN BEGIN 
            IF verbose GT 1 THEN PRINT, 'Accepted B: ',i,candidate_L,current_L, ratio, rnd
            ;; randomly accepted
            current_params[p] = candidate
            current_L = candidate_L
            model_params[i,p] = candidate
            model_params[i,1] = candidate_L
            naccepted=naccepted+1
         ENDIF ELSE BEGIN ;; not randomly accepted, repeat previous value
            IF verbose GT 1 THEN PRINT, 'Rejected: ',i,candidate_L,current_L,ratio,rnd
            model_params[i,p] = model_params[i-1,p]
            model_params[i,1] = model_params[i-1,1]
            nrejected=nrejected+1
         ENDELSE
      ENDELSE 
      
   ;-- print progress and update plots --
   IF i MOD (niter/20) EQ 0 THEN BEGIN
      IF verbose GT 0 THEN BEGIN
         PRINT, "Iteration ",i 
         PRINT, "  Current Params:"
         PRINT, "    ",model_params[i,*]
         countmedian = MEDIAN(model_params[1:i,2])
         PRINT, '  Median counts: ',countmedian
         IF i GT burnin THEN imin = burnin ELSE imin = 1
         stats = confidence_interval(model_params[imin:i,2]$
                                     ,level=0.9)
         PRINT, '  Mode Counts = ',stats[*,0]
         PRINT, '  90% Confidence: ',stats[*,1],stats[*,2]
      ENDIF 
      plotted = mcmc_plots(post_plot_file,model_params,burnin=burnin, $
                           maxiter=i,log=log,prior=prior)
      
   ENDIF 

   ENDFOR ;; endfor loop over free parameters
ENDFOR ;; end main mcmc loop

;-------------------------------------------
;         Get Posterior Median(s),
;           Confidence Intervals, and
;           Best Fit Image
;-------------------------------------------

;--Print Final Results to Screen--
PRINT, "Final Results:"
PRINT, ""

;-- print median, stdev, accepted/rejected statistics --
countmedian = MEDIAN(model_params[burnin:-1,2])
stats = confidence_interval(model_params[burnin:-1,2],level=0.9)
xmedian = MEDIAN(model_params[burnin:-1,3])
ymedian = MEDIAN(model_params[burnin:-1,4])
sigmamedian = MEDIAN(model_params[burnin:-1,5])
limedian = MEDIAN(model_params[burnin:-1,1])
PRINT, "median x0,y0,sigma = ",xmedian,ymedian,sigmamedian
PRINT, "total counts in input image = ",total(in_img)
PRINT, "median counts = ",countmedian
PRINT, "stdev counts = ",STDDEV(model_params[burnin:-1,2])
PRINT, "Mode counts = ",stats[*,0]
PRINT, "Accepted/Rejected = ",FLOAT(naccepted)/FLOAT(nrejected)
PRINT, "Accepted/Total    = ",FLOAT(naccepted)/FLOAT(nrejected+naccepted);FLOAT(niter)

;-- get 90% interval -- 
conf_interval = stats[*,1:2]
PRINT, "90% confidence = ",conf_interval

;-- save best fit image --
best_model_img = add_ptsrc(in_img,countmedian,xmedian,ymedian,width=sigmamedian,psfimg=psf_img)
writefits,fit_file,best_model_img,in_head

;-------------------------------------------
;    Plot Posterior(s) and Other Plots
;-------------------------------------------

plotted = mcmc_plots(post_plot_file,model_params,burnin=burnin,maxiter=niter-1,log=log,prior=prior)

;-------------------------------------------
;       Save model_params array
;-------------------------------------------

IF log THEN lhead = 'LogLikelihood' ELSE lhead = 'Likelihood'
header = '#Iteration     '+lhead+'       Counts     X0      Y0      SIGMA'
FORPRINT, model_params[*,0],model_params[*,1],model_params[*,2],model_params[*,3],model_params[*,4],model_params[*,5],textout=chain_file,COMMENT=header

END
