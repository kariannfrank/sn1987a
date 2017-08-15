;; Name: point_soure_limit.pro
;;
;; Date: 2015-03-23
;; Author: Kari A. Frank
;;
;; Purpose: Iteratively add point source until the lobes-only fit
;;          is worse at the 90% (or other) level to determine
;;          upper limit on point source flux
;;
;;
;; Calling Sequence:
;;          point_source_limit,smoothed_image,exptime=exptime,fparsfile=fparsfile
;;
;; Input:
;;
;;    smoothed_image -- standard smoothed image created by process.pro
;;
;;    exptime -- optionally include observation (filtered) exposure
;;               used for the initial count range.
;;    
;;    fparsfile -- if not given, expects to find this file in the
;;                 relevant lobes/ directory with the default name
;; 
;;
;; Output:
;;
;;    outimage -- copy of smoothed_image with a point source added at
;;                the upper flux limit.  named <smoothed_image>.ptsrclimit
;;
;;    outtext -- text file with the quantitative results
;;
;; Usage Notes:
;;
;;
;;
;; Example:
;;
;;

PRO point_source_limit,smoothed_image,exptime=exptime,fparsfile=fparsfile
;;FUNCTION template, arg1, arg2, key1=key1

;-------------------------------------------
;     Parse Arguments and Set Defaults
;-------------------------------------------

modeldir='lobes/'
IF NOT KEYWORD_SET(exptime) THEN exptime=48000
IF NOT KEYWORD_SET(fparsfile) THEN fparsfile=modeldir+smoothed_image+'_lobes_fpars.fits'

;-------------------------------------------
;             Set Up File Paths
;-------------------------------------------

outimage = modeldir+smoothed_image+'.ptsrclimit'
outtext = modeldir+smoothed_image+'.ptsrclimit.txt'

;-------------------------------------------
;           
;-------------------------------------------


;----Subsection 1----


;-Subsubsection 1-




END
