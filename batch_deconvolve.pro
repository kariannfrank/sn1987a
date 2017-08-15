;; Name: batch_deconvolve.pro
;;
;; Date: May 21, 2014
;; Author: Kari A. Frank
;;
;; Purpose:
;;   Read in SN1987A fits event file or list of event files
;;   and run deconvolve.pro on each of them
;; 
;; Calling Sequence:
;;   batch_deconvolve,input_files,source_x,source_y[,image_size=image_size]
;;
;; Input:
;;
;;   input_files -- [REQUIRED STRING] path (including file name) 
;;                   to the input event fits files. 
;;                   - if the file name ends with '.lis', then will run
;;                     in batch mode, reading file names from the 
;;                     .lis file (similar to ciao stacks).
;;                   - if input_files = 'auto', will use (or create
;;                     and use) a .lis file with all 
;;                     filtered*_subpix*.fits files
;;                     and then proceed in batch mode
;;
;;   image_size --[OPTIONAL FLOAT] side length of the input image
;;                 in pixels. default=40.0
;;
;;   source_x,source_y -- [REQUIRED FLOATS] position of the source
;;                         in physical coords. 
;;
;; Output:
;;   - fits file(s) containing the deconvolved image(s). 
;;   - if ran in batchmode, decon_images.lis file containing 
;;     list of deconvolved images 
;;
;; Usage Notes:
;;    - If using batch mode, all images must have the same
;;      source coordinates
;;    - If writing to deconvolved_images.lis, will append to the file if
;;      it already exists 
;;
;; Example:
;;   
;;

PRO batch_deconvolve,input_files,source_x,source_y,image_size=image_size

;-------------------------------------------
;     Parse Arguments and Set Defaults
;-------------------------------------------

IF N_PARAMS() LT 3 THEN MESSAGE, "ERROR: Minimum usage is 'batch_deconvolve, input_files,source_x,source_y'"

IF N_ELEMENTS(image_size) EQ 0 THEN image_size = 40

;----Check for batch mode (.lis file)----
batchmode = 0

IF (input_files EQ 'auto' OR input_files EQ 'AUTO') THEN BEGIN
   input_files = 'subband_evtfiles.lis'
; create list file if it doesn't exist
   IF FILE_TEST(input_files) EQ 0 THEN BEGIN
      SPAWN,'ls *evt2*_subpix_*-*.fits >& subband_evtfiles.lis'
   ENDIF; ELSE BEGIN
;      input_files = 'subband_evtfiles.lis'
;   ENDELSE
ENDIF

parts = STRSPLIT(input_files,'.',/EXTRACT)
IF N_ELEMENTS(parts) GT 1 THEN BEGIN
   IF parts[N_ELEMENTS(parts)-1] EQ 'lis' THEN batchmode = 1
ENDIF 

;-------------------------------------------
;        Set Up File Paths and Names
;-------------------------------------------

;----Store input file names image parameters----
IF batchmode EQ 0 THEN BEGIN

   IF N_ELEMENTS(image_size) EQ 0 THEN image_size = 40

   input_evtfiles = [input_files] 
   image_sizes = [image_size]
   sources_x = [source_x]
   sources_y = [source_y]
   n_files = 1
ENDIF ELSE BEGIN
;   IF (N_ELEMENTS(image_size) EQ 0) THEN BEGIN
;      ;if no image parameter arguments provided, read them from file
;      READCOL, input_files,FORMAT='A,F,F,F,F',input_evtfiles,input_image_sizes,output_image_sizes,sources_x,sources_y,/SILENT 
;      n_images = N_ELEMENTS(input_evtfiles)     
;   ENDIF ELSE BEGIN
      ;if image parameter arguments provided, create arrays
      READCOL, input_files,FORMAT='A',input_evtfiles,/SILENT      
      n_files = N_ELEMENTS(input_evtfiles)
      image_sizes = FLTARR(n_files)+image_size
      sources_x = FLTARR(n_files)+source_x
      sources_y = FLTARR(n_files)+source_y

;   ENDELSE
ENDELSE ;endelse batchmode

; open list file for appending new file names
OPENW,unit,'deconvolved_images.lis',/APPEND,/GET_LUN

;set size of new image (number of pixels in each dimension after rebinning)
;new_image_size = 338 

;-------------------------------------------
;            Loop Over Images
;-------------------------------------------

FOR i=0,n_files - 1 DO BEGIN

   ;----Read in FITS Image----
   evt_file = input_evtfiles[i]

   ;----Run deconvolve.pro----
   PRINT, evt_file,sources_x[i],sources_y[i],image_sizes[i]
   deconvolve,infile=evt_file,source_x=sources_x[i],source_y=sources_y[i],img_size=image_sizes[i]

   ;----Remove files to prevent it asking for removal confirmation----
   SPAWN,'rm -f idl*.ps marx.par'

   ;----Add new deconvolved image file name to list----
   PRINTF,unit,evt_file+'.deconvolved'

ENDFOR ;endfor i evt files

;----close list file----
CLOSE,unit

END
