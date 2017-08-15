;-------------------------------------------
; Function to calculate chi2 for 2 images
;-------------------------------------------
FUNCTION imgchi2,img1,img2

; assumes the image was fit with the lobes model
; typically img1 should be the observed image, img2 the model fit

err = sqrt(img1) ; poisson error
ddof = 19 ; number of degrees of freedom lobes+ptsrc model
w=WHERE(img1 ne 0 and img2 ne 0,nw); and err ne 0,nw)
dof=nw-ddof
;PRINT, 'dof = ',dof
chisq=total((img1[w]-img2[w])^2/err[w]^2);err[w]^2)

RETURN, chisq

END

;-------------------------------------------
;  Function to add point source to an image
;  (optionally convolving it with psf)
;-------------------------------------------
FUNCTION add_ptsrc,img,cts,x,y,width=width,psfimg=psfimg

IF NOT KEYWORD_SET(width) THEN width = 2.0*13.52

;-- set minimum allowed pixel value (for zero for very small numbers)
thresh = 10.0^(-7.0)
  
;-- get image size information --
shape = SIZE(img)

;-- create ptsrc-only image array --
ptimg = psf_gaussian(NPIXEL=[shape[1],shape[2]],FWHM=width,CENTROID=[x,y]);,/NORMALIZE)

;-- convolve ptsrc with psf --
IF N_ELEMENTS(psfimg) NE 0 THEN BEGIN
   ptimg_conv = DOUBLE(convolve(double(ptimg),double(psfimg),FT_PSF=psf_ft))
   ptimg_conv = cts*ptimg_conv/total(ptimg_conv)
;PRINT, "max ptimg_conv = ",max(ptimg_conv)
;PRINT, "total ptimg_conv = ",total(ptimg_conv)
   ptimg = ptimg_conv
ENDIF 

ptimg[WHERE(ptimg LE thresh)] = 0.0
ptimg = ptimg*cts/total(ptimg)
ptimg[WHERE(ptimg LE thresh)] = 0.0

;-- add (convolved) ptsrc to image --
RETURN, img+ptimg;_conv
;RETURN, img+psfimg

END

;-------------------------------------------
;              MAIN FUNCTION
;-------------------------------------------
PRO ptsrc_chi2,fitimgfile,level=level,x0=x0,y0=y0,sigma=sigma,verbose=verbose,countstep=countstep,maxcounts=maxcounts,detect=detect

; KAF 2016-04-22
; reads in the observed image and the best fit image (as used
;  in the usual image fitting with new_sn1987a_image_analysis)
;  and iteratively adds a point source of increasing counts to 
;  the best fit model until the chi2 increases the designated 
;  amount, to obtain an upper limit on point source counts

IF NOT KEYWORD_SET(level) THEN level=2.706 ; 90% confidence
IF NOT KEYWORD_SET(sigma) THEN sigma = 2.0*13.52
;IF NOT KEYWORD_SET(sigma) THEN sigma = 0.5; gaussian width of point source
                                           ;  in pixels (~1 deconvolved pixel)
IF NOT KEYWORD_SET(verbose) THEN verbose=1
IF NOT KEYWORD_SET(detect) THEN detect='no'

thresh = 10.0^(-7.0)

; read in deprojected observed image
obsimg = mrdfits(fitimgfile,0)
obsimg[WHERE(obsimg LT thresh)] = 0.0

; read in best fit image
modimg = mrdfits(fitimgfile,1)
;modimg = obsimg
modimg[WHERE(modimg LT thresh)] = 0.0

; calculate count increment
IF NOT KEYWORD_SET(maxcounts) THEN maxcounts = 10.0*max(obsimg)
;countstep = 
IF NOT KEYWORD_SET(countstep) THEN countstep = MAX([0.0000001,maxcounts/10000.0])
PRINT, 'maxcounts,countstep = ',maxcounts,countstep
;counts = MIN(obsimg[WHERE(obsimg NE 0.0)]);countstep
counts = 0.0

; determine x0,y0 in image coordinates
parsfile = REPSTR(fitimgfile,'fitim','fpars')
print, 'fpars file = ',parsfile
fpars = mrdfits(parsfile,1,range=1,columns=['X0','Y0'])
shape = SIZE(modimg)
IF NOT KEYWORD_SET(x0) THEN x0 = shape[1]/2.0 + fpars.X0
IF NOT KEYWORD_SET(y0) THEN y0 = shape[1]/2.0 + fpars.Y0
PRINT, 'x0, y0 = ',x0,y0

; initial chi2 value
;modimg = add_ptsrc(modimg,counts,x0,y0,sigma)
chisq0 = imgchi2(obsimg,modimg)
delchisq = 0.0
delchisq_check = delchisq
IF verbose GT 0 THEN PRINT, 'chisq initial = ',chisq0

; loop
WHILE delchisq_check LT level DO BEGIN

; choose new counts
counts = counts + countstep
IF counts GE maxcounts THEN BEGIN
   MESSAGE, 'Level not reached up to counts = '+STRING(counts)
ENDIF

; add point source to best fit image
ptsrcimg = add_ptsrc(modimg,counts,x0,y0,width=sigma)

; calculate new chi2
newchisq = imgchi2(obsimg,ptsrcimg)
delchisq = newchisq - chisq0
IF verbose EQ 1 THEN BEGIN
   IF counts MOD 0.1 EQ 0 THEN PRINT,'counts = ',counts
ENDIF
IF verbose GT 1 THEN BEGIN
   PRINT, 'counts = ',counts
   PRINT, 'newchisq = ',newchisq
   PRINT, 'delchisq = ',delchisq
   WRITEFITS,fitimgfile+'_ptsrcfit',ptsrcimg
ENDIF

IF detect EQ 'no' THEN delchisq_check = delchisq ELSE delchisq_check = ABS(delchisq)

ENDWHILE
WRITEFITS,fitimgfile+'_ptsrcfit',ptsrcimg
PRINT, 'newchisq = ',newchisq
PRINT, 'delchisq = ',delchisq

; print results
PRINT, 'count upper limit = ', counts

END
