;; Name: deconvolve_process.pro
;;
;; Date: May 21, 2014
;; Author: Kari A. Frank
;;
;; Purpose: run deconvolve.pro and process.pro on all subband event files
;;
;; Calling Sequence:
;;         
;;         deconvolve_process,source_x,source_y
;;
;; Input:
;;       
;;         source_x,source_y -- source position (in physical coords)
;;
;;         must be at least one subpixeled event file with standard
;;         name in the current folder
;;   
;; Output:
;;
;;        all usual output from deconvolve.pro and process.pro for
;;        subpixeled event file
;;
;; Usage Notes:
;; 
;;
;;
;; Example:
;;
;;

PRO deconvolve_process,source_x,source_y

;-------------------------------------------
;     Parse Arguments and Set Defaults
;-------------------------------------------






;-------------------------------------------
;        Set Up Files and Paths
;-------------------------------------------


;-------------------------------------------
;           Run deconvolution
;-------------------------------------------

;----Create file list----

SPAWN,'ls acisf*filtered*_subpix_*.fits >& subband_evtfiles.lis'
batch_deconvolve,'subband_evtfiles.lis',source_x,source_y

;-------------------------------------------
;          Run process (smoothing)
;-------------------------------------------

SPAWN,'ls acisf*filtered*_subpix*deconvolved > & deconvolved_images.lis'


END
