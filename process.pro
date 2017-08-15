;; Name: process.pro
;;
;; Date: July 12, 2013
;; Author: Kari A. Frank
;;         Rewritten version of the script by S. Park/E. Helder
;;
;;
;; Purpose:
;;   Read in SN1987A fits image file or list of fits image files
;;   produced with deconvolve.pro and produce a smoothed image
;;   file suitable for overlays with HST and radio images.
;; 
;; Calling Sequence:
;;   process,
;;input_files[,input_image_size=image_size][,output_image_size=output_image_size][,source_x=source_x,source_y=source_y][,suffix=suffix][gauss=gauss]
;;
;; Input:
;;
;;   input_files -- [REQUIRED STRING] path (including file name) 
;;                   to the input image fits files. if the file
;;                   name ends with '.lis', then will run in 
;;                   batch mode, reading files names from the 
;;                   .lis file (similar to ciao stacks).  see 
;;                   usage notes for formatting.
;;                
;;   input_image_size --[OPTIONAL FLOAT] desired side length of the 
;;                      input image subarray in pixels. default=40.0
;;
;;   output_image_size -- [OPTIONAL FLOAT] side length of the output
;;                        in pixels. must be smaller than 
;;                        input_image_size*rebin_factor (default
;;                        is 200).
;;
;;   source_x,source_y -- [OPTIONAL FLOATS] position of the source
;;                         in pixels.  must provide either neither
;;                         or both.
;;   rebin_factor -- [OPTIONAL FLOAT] rebinning factor for the image.
;;                   should be greater than 1 to create an image with 
;;                   more pixels than the original. new_image_size is 
;;                   calculated from this parameter as
;;                   input_image_size*rebin_factor.  setting rebin_factor=1
;;                   will skip the rebin command. default = 8.45.
;;
;;   suffix -- [OPTIONAL STRING] if provided, will be appended to end
;;             of the output file base name (after '_smoothed').
;;
;;   gauss -- [OPTIONAL FLOAT] provide the size of the gaussian kernel 
;;            for smoothing, in pixels. default=new_image_size/48.
;;
;; Output:
;;   - fits file containing the smoothed image(s).  naming convention
;;     root of the input file plus '_smoothed.img'
;;   - smoothed_image.lis text file with list of smoothed images (if
;;     already exists, will be overwritten)
;;
;; Usage Notes:
;;    - If using batch mode, all images must either have the same
;;      input and output sizes with the same source position, all 
;;      of which are provided as arguments (required for 'auto'), 
;;      or have the image sizes (side length) and source position 
;;      (x,y) provided as four columns in the .lis file. format for 
;;      each line should be:
;;      "file_name  input_side_length  output_side_length source_x  source_y"
;;    - Assumes images are square.
;; 
;; Example:
;;   process,'auto',source_x=80.32,source_y=79.22,input_image_size=40,output_image_size=200
;;   
;;

PRO process, input_files,input_image_size=input_image_size,source_x=source_x,source_y=source_y,rebin_factor = rebin_factor,suffix=suffix,gauss=gauss,output_image_size=output_image_size

;-------------------------------------------
;     Parse Arguments and Set Defaults
;-------------------------------------------

IF N_PARAMS() LT 1 THEN MESSAGE, "ERROR: Minimum usage is 'process, input_files'"

;----Check for batch mode (.lis file)----
batchmode = 0

IF (input_files EQ 'auto' OR input_files EQ 'AUTO') THEN BEGIN
   input_files = 'deconvolved_images.lis'
; create list file if it doesn't exist
   IF FILE_TEST(input_files) EQ 0 THEN BEGIN
      SPAWN,'ls *evt2*_subpix_*-*.deconvolved >& deconvolved_images.lis'
   ENDIF; ELSE BEGIN
ENDIF

parts = STRSPLIT(input_files,'.',/EXTRACT)
IF N_ELEMENTS(parts) GT 1 THEN BEGIN
   IF parts[-1] EQ 'lis' THEN batchmode = 1
ENDIF 
IF N_ELEMENTS(rebin_factor) EQ 0 THEN rebin_factor = 8.45
IF N_ELEMENTS(suffix) EQ 0 THEN suffix=''

;-------------------------------------------
;        Set Up File Paths and Names
;-------------------------------------------

;----Store input file names image parameters----
IF batchmode EQ 0 THEN BEGIN

   IF N_ELEMENTS(input_image_size) EQ 0 THEN input_image_size = 40
   IF N_ELEMENTS(output_image_size) EQ 0 THEN output_image_size = 200

   input_images = [input_files] 
   input_image_sizes = [input_image_size]
   output_image_sizes = [output_image_size]
   sources_x = [source_x]
   sources_y = [source_y]
   n_images = 1
ENDIF ELSE BEGIN
   IF (N_ELEMENTS(input_image_size) EQ 0) THEN BEGIN
      ;if no image parameter arguments provided, read them from file
      READCOL, input_files,FORMAT='A,F,F,F,F',input_images,input_image_sizes,output_image_sizes,sources_x,sources_y,/SILENT 
      n_images = N_ELEMENTS(input_images)     
   ENDIF ELSE BEGIN
      ;if image parameter arguments provided, create arrays
      READCOL, input_files,FORMAT='A',input_images,/SILENT      
      n_images = N_ELEMENTS(input_images)
      input_image_sizes = FLTARR(n_images)+input_image_size
      output_image_sizes = FLTARR(n_images)+output_image_size
      sources_x = FLTARR(n_images)+source_x
      sources_y = FLTARR(n_images)+source_y
   ENDELSE
ENDELSE ;endelse batchmode

;-open list file-
OPENW,unit,'smoothed_images.lis',/GET_LUN

;set size of new image (number of pixels in each dimension after rebinning)
;new_image_size = 338 

;-------------------------------------------
;            Loop Over Images
;-------------------------------------------

FOR i=0,n_images - 1 DO BEGIN

   ;----Read in FITS Image----
   image_file = input_images[i]
   image = READFITS(image_file,hdr)

   ;----Set image subarray and smoothing parameters----
   new_image_size = input_image_sizes[i]*rebin_factor   
   IF N_ELEMENTS(gauss) EQ 0 THEN gauss=new_image_size/48 ;for default values, this is 7 pixels ((40*8.45)/48)
   PRINT,'Beginning '+image_file
   PRINT, 'new_image_size = ',new_image_size

   ;----Construct output file name----
   parts = STRSPLIT(image_file,'.',/EXTRACT)
   n_parts = N_ELEMENTS(parts)
   out_image_file_str = ''
   FOR p = 0,n_parts-2 DO BEGIN
      out_image_file_str = out_image_file_str + parts[p]
   ENDFOR 
   IF suffix EQ '' THEN out_image_file = out_image_file_str + '_smoothed.img' ELSE out_image_file = out_image_file_str + '_smoothed_'+suffix+'.img'
   PRINT, out_image_file
   PRINTF,unit,out_image_file

   IF suffix EQ '' THEN unsm_out_image_file = out_image_file_str + '_unsmoothed.img' ELSE unsm_out_image_file = out_image_file_str + '_unsmoothed_'+suffix+'.img'
   
   ;----Set Image Parameters----
   image_radius = input_image_sizes[i]/2.0
   x0 = sources_x[i] - image_radius
   x1 = sources_x[i] + image_radius - 1.0
   y0 = sources_y[i] - image_radius
   y1 = sources_y[i] + image_radius - 1.0

   ;----Extract input image subarray----
   HEXTRACT,image,hdr,image1,image1hdr,x0,x1,y0,y1 ;40x40pixels, pixel=0.125"
   PRINT, 'x0,x1,y0,y1 = ',x0,x1,y0,y1

   ;----Rebin Image----

   IF rebin_factor NE 1.0 THEN BEGIN
     ;How to determine what the output
     ;sizes should be (338,160-100,etc)?
      HREBIN,image1,image1hdr,image2,image2hdr,new_image_size,new_image_size,/SAMPLE ;338x338pixels (40*8.45), pixel=40*0.125"/338=0.01479"
     ; ratio between pixel size and resolution?

   ;----Extract subarray of rebinned image----
      new_source_x = new_image_size/2
      new_source_y = new_source_x
;      x0 = new_source_x - new_image_size/2
;      x1 = new_source_x + new_image_size/2 - 1.0
;      y0 = new_source_y - new_image_size/2
;      y1 = new_source_y + new_image_size/2 - 1.0
      out_radius = output_image_size/2
      x0 = new_source_x - out_radius
      x1 = new_source_x + out_radius - 1.0
      y0 = new_source_y - out_radius
      y1 = new_source_y + out_radius - 1.0

      HEXTRACT,image2,image2hdr,image3,image3hdr,x0,x1,y0,y1 ;same pixel size as image2 (pix=0.01479",but image size is trimmed to 200x200pixels)
;      HEXTRACT,image2,image2hdr,image3,image3hdr,(160-100),(160+99), (165-100),(165+99)

   ENDIF ELSE BEGIN ;skip rebinning if rebin_factor = 1
      image3=image1
      image3hdr=image1hdr
   ENDELSE 
         
   ;----Smooth Rebinned Image with 7 pixel FWHM Gaussian----
;   image4 = FILTER_IMAGE(image3,FWHM=7)
   image4 = FILTER_IMAGE(image3,FWHM=gauss)
   image4hdr = image3hdr

   ;----Scale by binning factor to recover units of counts----
   ;unsmoothed image
   image3_scl = image3*(float(input_image_size)/float(new_image_size))^2.0

   ;smoothed image
   image5 = image4*(float(input_image_size)/float(new_image_size))^2.0
   image5hdr = image4hdr

   ;----Write New Smoothed Image to File----
   SXADDHIST,'Smoothed with '+STRING(gauss,FORMAT='(I0)')+' pixel FWHM Gaussian',image5hdr
;   WRITEFITS,out_image_file,image4,image4hdr
   WRITEFITS,out_image_file,image5,image5hdr

   ;----Write Unsmoothed Image to File----
   WRITEFITS,unsm_out_image_file,image3_scl,image3hdr

ENDFOR ;endfor i images

;-close list file-
CLOSE,unit

END
