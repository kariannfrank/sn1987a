;; Name: estimate_pileup_correction.pro
;;
;; Date: March 5, 2014
;; Author: Kari A. Frank
;;
;; Purpose: Use information from Helder2012 Table 1 to
;;           estimate pileup correction factor for a 
;;           given observation.  Currently, makes plots
;;           of pileup factor vs other parameters.
;;
;; Calling Sequence:
;;          estimate_pileup_correction
;;
;; Input:
;;          Requires text file with Helder Table 1 info.
;;   
;; Output:
;;          File containing plots of pileup correction factor
;;            vs other observation parameters.
;;
;; Usage Notes:
;;
;;
;;
;; Example:
;;
;;

PRO estimate_pileup_correction;, arg1, key1=key1
;;FUNCTION estimate_pileup_correction, arg1, arg2, key1=key1

;-------------------------------------------
;     Parse Arguments and Set Defaults
;-------------------------------------------


;-------------------------------------------
;             Set Up File Paths
;-------------------------------------------

sn1987a_path = '~/Dropbox/Research/SN1987A/'
table_file = sn1987a_path+'Helder_table1.txt'
out_plots = sn1987a_path + 'pileup_correction_plots.ps'


;-------------------------------------------
;             Read in Table
;-------------------------------------------


;----Read file----
;fluxes are 0.5-2.0 keV pileup-corrected fluxes in units of 10^-13 erg/s/cm^2
READCOL,table_file,/SILENT,COMMENT='#',FORMAT='A,A,F,F,F,F,F',obsids_str,instr,flux,frame,corrections,counts,exposures,DELIMITER=' '

;----Calculate counts/frametime----
count_rates = counts/(1000.0*exposures) ;count rates in counts/second
frame_counts = count_rates*frame ;counts/frame

;----Separate ACIS-S/NONE from ACIS-S/HETG----
hetg_subs = WHERE(instr EQ 'HETG')
s3_subs = WHERE(instr EQ 'S3')

;----Print out some stats----
hetg_mean = MEAN(corrections[hetg_subs])
hetg_median = MEDIAN(corrections[hetg_subs])
hetg_stdev = STDDEV(corrections[hetg_subs])

PRINT, ''
PRINT, 'mean, median, stdev, 100*stdev/mean of hetg corrections: '
PRINT, hetg_mean,hetg_median, hetg_stdev,100.*hetg_stdev/hetg_mean
PRINT, ''

;-------------------------------------------
;                Make Plots  
;-------------------------------------------

;----Set Up Plot Device----
SET_PLOT,'ps'
DEVICE, /COLOR
DEVICE, FILENAME=out_plots
LOADCT,38

color_levels = INDGEN(8)*256/8

syms = INTARR(N_ELEMENTS(flux))
syms[s3_subs] = 2
syms[hetg_subs] = 4

;--set dummy error array--
errors = FLTARR(N_ELEMENTS(syms))

;----Plot correction vs frametime----
num = N_ELEMENTS(syms)
colors = INDGEN(num)*256/num

plotted = fancy_plot(frame,corrections,errors,errors,PSYMS=syms,TITLE='',XTITLE='Frame Time (s)',YTITLE='Flux Pileup Correction Factor',COLORS=colors,XRANGE=[0.1,3.5])

LEGEND,['ACIS-S3','HETG'],PSYM=[2,4],/TOP,/LEFT,CHARTHICK=3

;----Plot correction vs soft flux----

;--set separate colors for different frametimes--
frame_labels = ['3.2','3.1','1.5','1.1','1.0','0.4','0.2']
frame_colors = color_levels[1:7]
colors = INTARR(N_ELEMENTS(syms))
colors[WHERE(frame EQ 3.2)] = color_levels[1]
colors[WHERE(frame EQ 3.1)] = color_levels[2]
colors[WHERE(frame EQ 1.5)] = color_levels[3]
colors[WHERE(frame EQ 1.1)] = color_levels[4]
colors[WHERE(frame EQ 1.0)] = color_levels[5]
colors[WHERE(frame EQ 0.4)] = color_levels[6]
colors[WHERE(frame EQ 0.2)] = color_levels[7]

;--plot--
plotted = fancy_plot(flux,corrections,errors,errors,PSYMS=syms,TITLE='',XTITLE='0.5-2.0 keV flux',YTITLE='Flux Pileup Correction Factor',COLORS=colors,XRANGE=[1.0,85.])

LEGEND,['ACIS-S3','HETG'],PSYM=[2,4],/TOP,/LEFT,CHARTHICK=3
LEGEND,frame_labels,TEXTCOLORS=frame_colors,/TOP,/RIGHT,CHARTHICK=3

;----Plot correction vs counts/frame----

;--plot--
plotted = fancy_plot(frame_counts,corrections,errors,errors,PSYMS=syms,TITLE='',XTITLE='Counts/frame',YTITLE='Flux Pileup Correction Factor',COLORS=colors,XRANGE=[0.0,0.9])

LEGEND,['ACIS-S3','HETG'],PSYM=[2,4],/TOP,/LEFT,CHARTHICK=3
LEGEND,frame_labels,TEXTCOLORS=frame_colors,/TOP,/RIGHT,CHARTHICK=3

;--fit line--
;fitted = LINEFIT_OPLOT(frame_counts,corrections,A=[50.0,50.0])
f_err = frame_counts*0.1
c_err = corrections*0.1
fitted = LINMIX_FIT(frame_counts,corrections,f_err,c_err)
PRINT, fitted

fitx=FINDGEN(10)/10.0
fity=fitted[0]*fitx+fitted[1]
fity_mhigh=(fitted[0]+fitted[3])*fitx+fitted[1]-fitted[4]
fity_mlow=(fitted[0]-fitted[3])*fitx+fitted[1]+fitted[4]

OPLOT,fitx,fity,LINESTYLE=0
OPLOT,fitx,fity_mhigh,LINESTYLE=1
OPLOT,fitx,fity_mlow,LINESTYLE=2

;----Close Plot File----
DEVICE, /CLOSE
SET_PLOT,'X'

END
