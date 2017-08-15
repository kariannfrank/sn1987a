@conf_error2d

pro new_sn1987a_image_analysis_parallel,input_files,incounts=incounts,lobes=lobes,ring=ring,ptsrc=ptsrc,bilat=bilat, lobsrc=lobsrc ,all=all,band=band,serial=serial,nomin=nomin,fixedsigr=fixedsigr

;-------------------------------------------
;     Parse Arguments and Set Defaults
;-------------------------------------------

IF N_PARAMS() LT 1 THEN MESSAGE, "ERROR: Minimum usage is 'new_sn1987a_image_analysis_parallel,input_files"

;--set up 'auto' parameters--
IF (input_files EQ 'auto' OR input_files EQ 'AUTO') THEN BEGIN
   input_files = 'smoothed_images.lis'
; create list file if it doesn't exist
   IF FILE_TEST(input_files) EQ 0 THEN BEGIN
      SPAWN,'ls *evt2*_subpix_*-*_smoothed.img >& smoothed_images.lis'
   ENDIF
   incounts = 'deconvolved_images_counts.txt'
ENDIF
IF N_ELEMENTS(serial) EQ 0 THEN serial = 0

;----Check for batch mode (.lis file)----
batchmode = 0
parts = STRSPLIT(input_files,'.',/EXTRACT)
IF N_ELEMENTS(parts) GT 1 THEN BEGIN
   IF parts[-1] EQ 'lis' THEN batchmode = 1
ENDIF 

IF N_ELEMENTS(incounts) EQ 0 THEN MESSAGE,'ERROR: no incounts parameter provided. Quitting.'

PRINT,'input_files = ',input_files
PRINT,'incounts = ',incounts

;--create array of file names and counts--
IF batchmode EQ 0 THEN BEGIN
   files = [input_files]
   counts = [incounts]
ENDIF ELSE BEGIN
   READCOL,input_files,/SILENT,files,FORMAT='A'
   READCOL,incounts,/SILENT,countfiles,counts,FORMAT='A,I',DELIMITER=' '
ENDELSE

IF NOT keyword_set(lobes) THEN lobes=0 ;& outdir = 'lobes/'
IF NOT keyword_set(ring) THEN ring=0 ;& outdir = 'ring/'
IF NOT keyword_set(ptsrc) THEN ptsrc=0 ;& outdir = 'ptsrc/'
IF NOT keyword_set(bilat) THEN bilat=0 ;& outdir = 'bilat/'
IF NOT keyword_set(lobsrc) THEN lobsrc=0 ;& outdir = 'lobsrc/'
IF NOT keyword_set(all) THEN all=0

IF keyword_set(fixedsigr) THEN suffix='_fixedwidth' ELSE suffix=''
IF keyword_set(lobes) THEN outdir = 'lobes'+suffix+'/'
IF keyword_set(ring) THEN outdir = 'ring'+suffix+'/'
IF keyword_set(ptsrc) THEN outdir = 'ptsrc/'
IF keyword_set(bilat) THEN outdir = 'bilat'+suffix+'/'
IF keyword_set(lobsrc) THEN outdir = 'lobsrc'+suffix+'/'

IF NOT keyword_set(nomin) THEN nomin = 0
IF NOT keyword_set(fixedsigr) THEN fixedsigr = 0

nf=n_elements(files)

;-------------------------------------------
;  Run (parallel) Loop of Image Fitting
;-------------------------------------------

  ;regular for loop
;  for k=0,nf-1 do begin 
;     file=files[k]
;     count=counts[k]

     ;call image fitting
;     new_sn1987a_image_analysis_single,file,count,lobes=lobes,ring=ring,ptsrc=ptsrc,bilat=bilat,lobsrc=lobsrc,all=all

;  endfor 
  ;return

IF batchmode EQ 1 THEN BEGIN

   IF serial EQ 0 THEN BEGIN
  ; parallel loop
      split_for, 0,nf-1, commands=[$
                'file=files[i]',$
                'count=counts[i]',$

     ;call image fitting
                'new_sn1987a_image_analysis_single,file,count,lobes=lobes,ring=ring,ptsrc=ptsrc,bilat=bilat,lobsrc=lobsrc,all=all,nomin=nomin,fixedsigr=fixedsigr'],$
             
                varnames=['files','counts','lobes','ring','ptsrc','bilat','lobsrc','all','nomin','fixedsigr'],$
                nsplit=7
   ENDIF ELSE BEGIN
      ; serial loop
      FOR f=0,N_ELEMENTS(files)-1 DO BEGIN
         new_sn1987a_image_analysis_single,files[f],counts[f],lobes=lobes,ring=ring,ptsrc=ptsrc,bilat=bilat,lobsrc=lobsrc,all=all,nomin=nomin,fixedsigr=fixedsigr
      ENDFOR
   ENDELSE
ENDIF ELSE BEGIN
     new_sn1987a_image_analysis_single,files[0],counts[0],lobes=lobes,ring=ring,ptsrc=ptsrc,bilat=bilat,lobsrc=lobsrc,all=all,nomin=nomin,fixedsigr=fixedsigr
ENDELSE
;-------------------------------------------
; Gather all results in a single text file
;-------------------------------------------

;SPAWN,'anaconda'
SPAWN,'anaconda & $pydir/sn1987a/get_radii_results.py --clobber yes --results_dir '+outdir

end
