;Use MPFIT to fit the light curve


PRO light_curve_fit,infile=infile,band=band,mincounts=mincounts,func=func

;----Parse arguments and set defaults----
IF NOT KEYWORD_SET(infile) THEN infile = 'spectra465_frank2015a_fits_official.txt_fluxes.txt'
IF NOT KEYWORD_SET(band) THEN band = 'soft'
IF NOT KEYWORD_SET(mincounts) THEN mincounts = 0
IF NOT KEYWORD_SET(func) THEN func = 'brokenlogline'
  
;----Set Useful Constants----
sn1987a_dir = '/data/nasdrive/main/kaf33/sn1987a_monitoring/comparison_dir/spectra_current/'

;----Read Data----
READCOL,/SILENT,infile,obsid,age,grating,counts,soft,softlow,softhigh,soft12,soft12low,soft12high,hard,hardlow,hardhigh,broad,broadlow,broadhigh,SKIPLINE=1,COMMENT='#',FORMAT='A,D,A,I,D,D,D,D,D,D,D,D,D,D,D,D'

CASE band OF
   'hard': BEGIN 
      flux = hard
      fluxlow = hardlow
      fluxhigh = hardhigh
   END
   'soft': BEGIN
      flux = soft
      fluxlow = softlow
      fluxhigh = softhigh
   END
   'soft12': BEGIN
      flux = soft12
      fluxlow = soft12low
      fluxhigh = softhigh12
   END
   'broad': BEGIN
      flux = broad
      fluxlow = broadlow
      fluxhigh = broadhigh12
   END
ENDCASE

;----Convert to Log and Calculate Error Bars----
; force symmetric errors, just use whichever is larger for each point, upper or lower error
fluxerr = FLTARR(N_ELEMENTS(flux))
logflux = FLTARR(N_ELEMENTS(flux))
logfluxerr = FLTARR(N_ELEMENTS(flux))
FOR j=0,N_ELEMENTS(flux)-1 DO BEGIN
   logflux[j] = alog10(flux[j])
   loghigherr = alog10(fluxhigh[j])-alog10(flux[j])
   loglowerr = alog10(flux[j])-alog10(fluxlow[j])
   logfluxerr[j] = max([loglowerr,loghigherr])
;   fluxerr[j] = max([fluxlow[j],fluxhigh[j]])
ENDFOR
;print, 'logfluxerr = ',logfluxerr

;-----------------Broken Linear 2 (two changepoints)-----------------
IF func EQ 'brokenlogline' THEN BEGIN

;----Initial Values and Limits----

;parameter array should be
;[intercept,seg1slope,changepoint1,seg2slope,changepoint2,seg3slope]

;--create parinfo structure--
parinf = REPLICATE({fixed:0,limited:[0,0],limits:[0.D,0.D]},6)

;set intercept limits
parinf[0].limited[0] = 1
parinf[0].limits[0] = -3.
parinf[0].limited[1] = 1
parinf[0].limits[1] = 3.

;initial estimates 
intercept = -2.
slope0 = 0.004 ;log10 e-13 erg/cm^2/s per day
slope1 = 0.0006
slope2 = 0.0002
changepoint1 = 8000.0
changepoint2 = 10000.0

;set slope0 limits
;parinf[1].limited[0] = 1
;parinf[1].limits[0] = slope0;*0.2
;parinf[1].limited[1] = 1
;parinf[1].limits[1] = slope0;*5.0

;set slope1 limits
parinf[3].limited[0] = 1
parinf[3].limits[0] = 0.0
parinf[3].limited[1] = 0.5*slope0
parinf[3].limits[1] = 0;slope1*2.0

;set slope2 limits
parinf[5].limited[0] = 0
;parinf[5].limits[0] = slope0
parinf[5].limited[1] = 0
;parinf[5].limits[1] = 10.0*slope2

;set changepoint1 limits
parinf[2].limited[0] = 1
parinf[2].limits[0] = FLOAT(MIN(age));5000.0
parinf[2].limited[1] = 1
parinf[2].limits[1] = FLOAT(MAX(age));11000.0

;set changepoint2 limits
parinf[4].limited[0] = 1
parinf[4].limits[0] = FLOAT(MIN(age));5000.0
parinf[4].limited[1] = 1
parinf[4].limits[1] = FLOAT(MAX(age));11000.0

;--set starting values--
p0 = [DOUBLE(intercept),DOUBLE(slope0),DOUBLE(changepoint1),DOUBLE(slope1),DOUBLE(changepoint2),DOUBLE(slope2)]
print, 'p0 = ',p0
ENDIF 

;-----------------2 changepoints-----------------
IF func EQ 'brokenline' THEN BEGIN
;----Initial Values and Limits----

;parameter array should be
;;P0 = early slope, P1 = late slope, P2 = intercept, P3 = change point

;--create parinfo structure--
parinf = REPLICATE({fixed:0,limited:[0,0],limits:[0.D,0.D]},4)

;set intercept limits
;parinf[2].limited[0] = alog10(1)
;parinf[2].limits[0] = -0.5
;parinf[2].limited[1] = 1
;parinf[2].limits[1] = 0.5

;initial estimates 
intercept = -2.
slope0 = 0.006 ;log10 e-13 erg/cm^2/s per day
slope1 = 0.0002
changepoint1 = 6000.0

;set slope0 limits
parinf[0].limited[0] = 1
parinf[0].limits[0] = slope0*0.2
parinf[0].limited[1] = 1
parinf[0].limits[1] = slope0*5.0

;set slope1 limits
parinf[1].limited[0] = 1
parinf[1].limits[0] = slope1*0.2
parinf[1].limited[1] = 1
parinf[1].limits[1] = slope1*2.0

;set changepoint limits
parinf[3].limited[0] = 1
parinf[3].limits[0] = 6000.
parinf[3].limited[1] = 1
parinf[3].limits[1] = 10000.

;--set starting values--
p0 = [DOUBLE(slope0),DOUBLE(slope1),DOUBLE(intercept),DOUBLE(changepoint1)]
ENDIF 

;----Fit with MPFIT----
rosatage = DOUBLE(1200)

thefit = MPFITFUN(func,age,logflux,logfluxerr,p0,YFIT=logmodely,PERROR=param_errors,BEST_FJAC=jcb,COVAR=cvr,PFREE_INDEX=pfree,/CALC_FJAC,/QUIET,BESTNORM=chi2,DOF=dof,MAXITER=10000,NITER=niter)
modely = 10.0^logmodely
print, 'niter = ',niter

;----Print the Best Fit Parameters----
;PRINT, thefit[0]*arcsecdays_to_kms,
;thefit[1]*arcsecdays_to_kms,thefit[2],thefit[3]

PRINT, 'model fluxes = ',modely
PRINT, 'best fit params = ',thefit
PRINT, 'best fit param errors = ',param_errors
PRINT, 'chi2, dof, chi2/dof = ',chi2,dof,chi2/dof

;----Save Best Fit to File----
;outfile = 'mpfit_results_3.txt'
outfile = FILE_BASENAME(infile,'.txt')+'_'+band+'_mpfit_'+func+'_params.txt'

OPENW,lun,outfile,/GET_LUN
IF func EQ 'brokenlogline' THEN BEGIN
   PRINTF, lun,'stat'+string(9B)+'intercept'+string(9B)+'slope0'+string(9B)+'changepoint1'+string(9B)+'slope1'+string(9B)+'changepoint2'+string(9B)+'slope3'+string(9B)+'chi2'+string(9B)+'dof'
   PRINTF, lun, FORMAT=('(%"%S\t%F\t%F\t%F\t%F\t%F\t%F\t%F\t%F")'),'bestfit',thefit[0],thefit[1],thefit[2],thefit[3],thefit[4],thefit[5],chi2,dof
   PRINTF, lun,FORMAT=('(%"%S\t%F\t%F\t%F\t%F\t%F\t%F")'),'stdev',param_errors[0],param_errors[1],param_errors[2],param_errors[3],param_errors[4],param_errors[5]
ENDIF ELSE BEGIN
   PRINTF, lun,'stat'+string(9B)+'intercept'+string(9B)+'slope0'+string(9B)+'changepoint1'+string(9B)+'slope1'+string(9B)+'chi2'+string(9B)+'dof'
   PRINTF, lun, FORMAT=('(%"%S\t%F\t%F\t%F\t%F\t%F\t%F\t%F\t%F")'),'bestfit',thefit[2],thefit[0],thefit[3],thefit[1],chi2,dof
   PRINTF, lun,FORMAT=('(%"%S\t%F\t%F\t%F\t%F\t%F\t%F")'),'stdev',param_errors[2],param_errors[0],param_errors[3],param_errors[1]
ENDELSE
FREE_LUN,lun

;----Save fitted y values with errors----
;(similar format as the input file)

;--get fitted y values--
;fittedy = brokenline(r0,thefit)

;--get fitted y errors--
logmodely_err = (MPPROPERR(jcb,cvr,pfree_index=pfree,/DIAG))^0.5
modely_err = 10.0^(logmodely) - 10.0^(logmodely - logmodely_err)

;--print fitted observation fluxes to file--
;fitfile = 'sn1987a_mpfit_radii_fits.txt'
fitfile = FILE_BASENAME(infile,'.txt')+'_'+band+'_mpfit_'+func+'_fluxes.txt'

OPENW,lun,fitfile,/GET_LUN
PRINTF,lun,'obsid'+string(9B)+'age'+string(9B)+'flux'+string(9B)+'fluxlow'+string(9B)+'fluxhigh'+string(9B)+'counts'
FOR d=0,N_ELEMENTS(flux)-1 DO BEGIN
   PRINTF,lun,FORMAT='(%"%S\t%I\t%F\t%F\t%F\t%I")',obsid[d],age[d],modely[d],modely_err[d],modely_err[d],counts[d]
ENDFOR
FREE_LUN,lun

;----Save age and fluxes for 'continuous' best fit line----
;(similar format as the input file)

;--get fitted y values (no errors)--
IF func EQ 'brokenlogline' THEN fittedx = [4000.0,thefit[2],thefit[4],11000.0] ELSE fittedx = [4000.0,thefit[3],11000.0]
IF func EQ 'brokenlogline' THEN fittedy = brokenlogline(fittedx,thefit) ELSE fittedy = brokenline(fittedx,thefit)

;--print to file--
linefile = FILE_BASENAME(infile,'.txt')+'_'+band+'_mpfit_'+func+'_line.txt'

OPENW,lun,linefile,/GET_LUN
PRINTF,lun,'age'+string(9B)+'flux'
FOR d=0,N_ELEMENTS(fittedx)-1 DO BEGIN
   PRINTF,lun,FORMAT='(%"%F\t%F")',fittedx[d],10.0^fittedy[d]
ENDFOR
FREE_LUN,lun

END
