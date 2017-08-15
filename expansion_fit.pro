;Use MPFIT to fit the expansion curve


PRO expansion_fit,infile=infile,band=band,mincounts=mincounts

;----Parse arguments and set defaults----
IF NOT KEYWORD_SET(infile) THEN infile = 'radii_fits_official.txt'
IF NOT KEYWORD_SET(band) THEN band = '300-8000'
IF NOT KEYWORD_SET(mincounts) THEN mincounts = 0
  
;----Set Useful Constants----

;--convert fitted radii and errors to arcsec--
r_to_arcsec = 0.0147928866

;--unit conversions--
kpc_to_km = 1000.0*3.0857*10.0^13.0
distance_kpc = 51.4 ;Panagia2003
distance_km = distance_kpc*kpc_to_km
arcsec_to_radian = !PI/(3600.0*180.0)
arcsec_to_km = arcsec_to_radian*distance_km
days_to_s = 24.0*3600.0
arcsecdays_to_kms = arcsec_to_km/days_to_s


;----Read Data----

;indir = '/export/bulk/rice3/kaf33/Chandra_Observations/SN1987A/comparison_dir/imaging/'
;sn1987a_dir = '/media/backup/kaf33/main/Chandra_Observations/SN1987A/'
;infile = indir+'sn1987a_radii_fits.txt'

;READCOL,/SILENT,infile,obsid,age,r0,r0err_low_raw,r0err_upp_raw,counts,deconcounts,SKIPLINE=1,COMMENT='#';,FORMAT='I,F,F,F,F,I,I'

; file as formatted by get_radius_results.py
; radius units already in arscec
READCOL,/SILENT,infile,obsid,age,rband,r0,r0err_low,r0err_upp,sigr,sigrerr_low,sigrerr_upp,counts,deconcounts,skyns,skynserr_low,skynserr_upp,SKIPLINE=1,COMMENT='#',FORMAT='I,F,A,F,F,F,F,F,F,I,I,F,F,F'


;----Calculate Errors----
;r0err_low = SQRT( r0err_low_raw^2.0 + r0^2.0/counts  )
;r0err_upp = SQRT( r0err_upp_raw^2.0 + r0^2.0/counts  )

; force symmetric errors, just use whichever is larger for each point, upper or lower error
r0err_max = FLTARR(N_ELEMENTS(r0))
FOR j=0,N_ELEMENTS(r0)-1 DO BEGIN
   r0err_max[j] = max([r0err_low[j],r0err_upp[j]])
ENDFOR

;--convert radii to arcsec--
;r0 = r0*r_to_arcsec
;r0err_low = r0err_low*r_to_arcsec
;r0err_upp = r0err_upp*r_to_arcsec
;r0err_max = r0err_max*r_to_arcsec


;----Select Only The Specified Band----
bandi = WHERE((rband EQ band) AND (deconcounts GE mincounts))
;print, 'bandi = ',bandi
PRINT, obsid[bandi]

;----Initial Values and Limits----

;parameter array should be
;[earlyslope,lateslope,intercept,changepoint]

;--create parinfo structure--
parinf = REPLICATE({fixed:0,limited:[0,0],limits:[0.D,0.D]},4)

;set early slope limits
;parinf[0].limited[0] = 1
;parinf[0].limits[0] = 1000.0/arcsecdays_to_kms
;parinf[0].limited[1] = 1
;parinf[0].limits[1] = 15000.0/arcsecdays_to_kms

;set late slope limits
parinf[1].limited[0] = 0
parinf[1].limits[0] = 0.0/arcsecdays_to_kms
;parinf[1].limited[1] = 1
;parinf[1].limits[1] = 15000.0/arcsecdays_to_kms

;set intercept limits
;parinf[2].limited[0] = 1
;parinf[2].limits[0] = -1.0
;parinf[2].limited[1] = 1
;parinf[2].limits[1] = 1.0

;set change point limits
parinf[3].limited[0] = 1
parinf[3].limited[0] = FLOAT(MIN(age[bandi]));5000.0
;parinf[3].limited[1] = 1
;parinf[3].limited[1] = FLOAT(MAX(age[bandi]));11000.0

;--set starting values--
p0 = [10000.0/arcsecdays_to_kms,1800.0/arcsecdays_to_kms,0.2,6000.0];MEAN(FLOAT(age[bandi]))];6000.0]

;----Fit with MPFIT----

;PRINT, brokenline(age,p0)

thefit = MPFITFUN('brokenline',age[bandi],r0[bandi],r0err_max[bandi],p0,YFIT=modely,PERROR=param_errors,BEST_FJAC=jcb,COVAR=cvr,PFREE_INDEX=pfree,/CALC_FJAC,/QUIET,PARINFO=parinf)

;----Print the Best Fit Parameters----
;PRINT, thefit[0]*arcsecdays_to_kms,
;thefit[1]*arcsecdays_to_kms,thefit[2],thefit[3]

;----Save Best Fit to File----
;outfile = 'mpfit_results_3.txt'
outfile = FILE_BASENAME(infile,'.txt')+'_'+band+'_mpfit_params.txt'

PRINT, 'model R0 = ',modely
PRINT, 'best fit params = ',thefit
PRINT, 'best fit param errors = ',param_errors

OPENW,lun,outfile,/GET_LUN
PRINTF, lun,'stat'+string(9B)+'v_early'+string(9B)+'v_late'+string(9B)+'intercept'+string(9B)+'changepoint'
PRINTF, lun,'best',thefit[0]*arcsecdays_to_kms,thefit[1]*arcsecdays_to_kms,thefit[2],thefit[3]
PRINTF, lun,'stdev',param_errors[0]*arcsecdays_to_kms,param_errors[1]*arcsecdays_to_kms,param_errors[2],param_errors[3]
FREE_LUN,lun

;----Save fitted y values with errors----
;(similar format as the input file)

;--get fitted y values--
;fittedy = brokenline(r0,thefit)

;--get fitted y errors--
modely_err = (MPPROPERR(jcb,cvr,pfree_index=pfree,/DIAG))^0.5

;--print to file--
;fitfile = 'sn1987a_mpfit_radii_fits.txt'
fitfile = FILE_BASENAME(infile,'.txt')+'_'+band+'_mpfit_radii.txt'
PRINT, 'fitfile = ',fitfile

OPENW,lun,fitfile,/GET_LUN
;PRINTF,lun,FORMAT='(%"%S\t%S\t%S\t%S\t%S\t%S\t%S\t")','obsid','age','radius','radiuslow','radiushigh','counts','deconcounts'
PRINTF,lun,'obsid'+string(9B)+'age'+string(9B)+'radius'+string(9B)+'radiuslow'+string(9B)+'radiushigh'+string(9B)+'counts'+string(9B)+'deconcounts'
FOR d=0,N_ELEMENTS(bandi)-1 DO BEGIN
   PRINTF,lun,FORMAT='(%"%S\t%I\t%F\t%F\t%F\t%I\t%I")',obsid[bandi[d]],age[bandi[d]],modely[d],modely_err[d],modely_err[d],counts[bandi[d]],deconcounts[bandi[d]]
ENDFOR
FREE_LUN,lun

;----Save age and radii for 'continuous' best fit line----
;(similar format as the input file)

;--get fitted y values (no errors)--
fittedx = [4000.0,thefit[3],11000.0]
fittedy = brokenline(fittedx,thefit)

;--print to file--
linefile = FILE_BASENAME(infile,'.txt')+'_'+band+'_mpfit_line.txt'

OPENW,lun,linefile,/GET_LUN
PRINTF,lun,'age'+string(9B)+'radius'
FOR d=0,N_ELEMENTS(fittedx)-1 DO BEGIN
   PRINTF,lun,FORMAT='(%"%F\t%F")',fittedx[d],fittedy[d]
ENDFOR
FREE_LUN,lun

END
