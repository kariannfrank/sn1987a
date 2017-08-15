;	this idl program extracts images
;	from one .evt file 
;	and deconvolves the image using a PSF from Marx.
;	Note that this version smooths the image (3x3) before 
;	deconvolving it.
;
;	12/31/99: based on George Chartas' image.pro program, which
;	combined several .evt files and deconvolved the images, then
;	combined the deconvolved images.
;
;       7/20/2011: This version uses Marx 4.4 (EAH)
;
;
;       05/21/2014: Updated to be handle providing source position
;                   as argument (non-interactive).  To run
;                   non-interactively, must provide all arguments.
;                   It is also possible to provide some arguments and
;                   not others.  Unspecified arguments will be
;                   prompted for. (KAF)
;                     
;

@eb_parameter
@eb_property
@eb_derived_property
@domain_dataset
	

PRO my_deconvolve,infile=infile,source_x=source_x,source_y=source_y,img_size=img_size

;	input parameters


;dn2evfile = '/bulk/pkg/xray/idl_lib/pro/acis.dn2ev'
split_threshold = 14

IF N_ELEMENTS(infile) EQ 0 THEN BEGIN
   evtfile = pickfile(/READ, TITLE='Choose events file...', /NOCONFIRM)
   ;print,'Input file: ', evtfile

   print,'Starting ds9 so you can select desired subimage ...'
   spawn,'ciao; ds9 ' + evtfile + '&'
ENDIF ELSE BEGIN
   evtfile = infile
ENDELSE

;	define rectangular event filter region 

;;original code
;read,'Enter the physical X, Y coordinates for the lower left corner: ', xllc, yllc
;read,'Enter the physical X, Y coordinates for the upper right corner: ', xurc, yurc

; get coordinates if not specified
IF N_ELEMENTS(source_x) EQ 0 THEN BEGIN
  read,'Enter the physical X coordinate of source: ', source_x
ENDIF
IF N_ELEMENTS(source_y) EQ 0 THEN BEGIN
  read,'Enter the physical Y coordinate of source: ', source_y
ENDIF
IF N_ELEMENTS(img_size) EQ 0 THEN BEGIN
  read,'Enter side length of output image: ', img_size
ENDIF
; define box
xllc = source_x - img_size/2.0
yllc = source_y - img_size/2.0
xurc = source_x + img_size/2.0
yurc = source_y + img_size/2.0



pixel_size = 0.125 	; output pixel size in arcseconds

dx = pixel_size/0.492	; pixel width for deconvolution in units of original pix
dy = dx	


; Get grating from event file header
evdata = readfits(infile,head,exten_no=1)
grat = fxpar(head,'GRATING')


outfile = evtfile + '.image'
;read,'Enter name for output raw image file: ', outfile

; Run dmcopy to make an image of the desired portion of the events file
; convert box size and steps into strings
sxllc = strtrim(string(xllc),2) 
sxurc = strtrim(string(xurc),2)
sdx = strtrim(string(dx),2)
syllc = strtrim(string(yllc),2)
syurc = strtrim(string(yurc),2)
sdy = strtrim(string(dy),2)

print,'Spawning dmcopy to create the raw image FITS file ...'

spawn, 'ciao; dmcopy "' + evtfile + '[EVENTS][grade=0,2,3,4,6]' $
 	+ '[bin x=' + sxllc + ':' + sxurc + ':' + sdx $
 	+ ', y=' + syllc + ':' + syurc + ':' + sdy + ']" ' + $
 	outfile + ' option=image clobber=yes'

print, 'ciao; dmcopy "' + evtfile + '[EVENTS][grade=0,2,3,4,6]' $
 	+ '[bin x=' + sxllc + ':' + sxurc + ':' + sdx $
 	+ ', y=' + syllc + ':' + syurc + ':' + sdy + ']" ' + $
 	outfile + ' option=image clobber=yes'

; werkt vanaf hier met .image. Dus.. als ik zelf een puntbron toevoeg
; aan image, vanaf hier inlezen.. (google translation = "works
; from here with .image. So .. if I do add a point source; to image, from reading here")

print,'Reading and displaying raw image ...'
im0 = readfits(outfile, header)
im = smooth(im0,3)
s = size(im)
syz1 = s(1)
syz2 = s(2)
repl = fix(512/syz2) + 1
tvscl, rebin(im, syz1*repl, syz2*repl, /sample)



;	Determine psf's using MARX 

;       first remove current marx parameter file and cp original one
;       spawn, '/usr/bin/rm marx.par'
;       spawn, 'cp /bulk/pkg/cxc/marx/marx.par .'
;       spawn, 'cp /usr/astro/marx/marx.par .'
;        spawn, 'cp /usr/astro/marx-dist-4.3.0/share/marx/pfiles/marx.par .'
        spawn, 'cp /usr/astro/marx/share/marx/pfiles/marx.par .'

;       Format the marx command line we need, and execute it:

        f='("marx ExposureTime=",F0, " OutputDir=", A, " GratingType=", A, " DetectorType=", A, " DetOffsetX=",F0, " DetOffsetY=",F0, " DetOffsetZ=",F0, " DetIdeal=", A, " SourceFlux=",F0, " SpectrumType=", A, " SpectrumFile=", A, " SourceType=", A, " SourceOffsetZ=",F0, " SourceOffsetY=",F0, "  S-GaussSigma=",F0)'


;ExposureTime = 100. ; exposure time  sec - 0 ==> run MARX in ray
;                                           generation mode ->
;                                           apparently no longer used...(EAH)
;OutputDir ="/bulk/draco2/burrows/marxout1"
OutputDir ="./"
;GratingType = grat
GratingType = "NONE"
DetectorType ="ACIS-S"
;	dox = 0.0
;	doy = 0.0
;	doz = 0.0
        dox = 0.001444942
        doy = 1.0826158009
        doz = 0.996576645
	sy = 0.0
	sz = 0.0
DetOffsetX = dox
DetOffsetY = doy
DetOffsetZ = doz
DetIdeal="no"
SourceFlux = -1
SpectrumType = "FILE"
;SpectrumFile = "/bulk/draco2/burrows/knotmod.dat"   ; wabs pow = 0.09 1.7 1
SpectrumFile = "~/IDL_Programs/sn1987a/knotmod.dat"
;SpecctrumFile = "../obs15810_spectrum_chart.dat"
SourceType="GAUSS"
;SourceType="SAOSAC"
SourceOffsetZ =sz
SourceOffsetY =sy
SGaussSigma = 0.17

	nl=1
	source_circle_xc = dblarr(nl)
	source_circle_yc = dblarr(nl)

	source_circle_xc(*) = 0
	source_circle_yc(*) = 0
xc0 = source_circle_xc(0)
yc0 = source_circle_yc(0)

	
;	 Get psf for current obsid
;ypsf = read_marx_file('/bulk/draco2/burrows/marxout1/xpixel.dat')
;zpsf = read_marx_file('/bulk/draco2/burrows/marxout1/ypixel.dat')
;ypsf = read_marx_file('./xpixel.dat')
;zpsf = read_marx_file('./ypixel.dat')
;ypsf = read_marx_file('/astro/research/kaf33/IDL_Programs/sn1987a/xpixel.dat')
;zpsf = read_marx_file('/astro/research/kaf33/IDL_Programs/sn1987a/ypixel.dat')
ypsf = read_marx_file('xpixel.dat')
zpsf = read_marx_file('ypixel.dat')

;	find centroid of psf and create a 499x499 psf array with the peak of
;	the psf center on pixel 249,249 

;	center for zero offset is at 511, 250
ypsf_m = mean(ypsf)
zpsf_m = mean(zpsf)

;
;	PSF image is 80 x 80;  each pixel is 0.492*dx arcsec
;
;psf = make_image(ypsf,zpsf,xbinsize=dx,ybinsize=dy,xrange=[ypsf_m-10,ypsf_m+10],yrange=[zpsf_m-10,zpsf_m+10])
psf = make_image(ypsf,zpsf,xbinsize=dx,ybinsize=dy,xrange=[ypsf_m-40,ypsf_m+40],yrange=[zpsf_m-40,zpsf_m+40])


;
;	Now, select the 21 x 21 subarray containing the PSF peak
;
psf_max = max(psf,ind)
psf_index = ind
ymax= ind mod n_elements(psf(*,0))
zmax= ind/n_elements(psf(0,*))

psf = psf(ymax-10.:ymax+10.,zmax-10.:zmax+10.)
;psf = psf(ymax-20.:ymax+20.,zmax-20.:zmax+20.)

print, psf[where(psf NE 0)]

title='MARX-generated PSF used for image deconvolution!C!Cfor '$
		+ strtrim(string(total(psf)),2) + ' photons'
xlabel = strtrim(string(pixel_size),2) + ' arcsecond pixels'
ylabel = xlabel
contour,psf,nlevels=10,title=title, xtitle=xlabel, ytitle=ylabel
surface,psf,title=title, xtitle=xlabel, ytitle=ylabel
;stop

; print out plots
prt = ''
;read,'Do you want printed output (Y/N)?', prt
prt = 'Y'
if (prt eq 'y' or prt eq 'Y') then begin
	set_plot,'ps'
	device,/landscape
	contour,psf,nlevels=10,title=title, xtitle=xlabel, ytitle=ylabel
	surface,psf,title=title, xtitle=xlabel, ytitle=ylabel
	device,/close
;	spawn,'lpr -r idl.ps'
	spawn,'mv idl.ps idl1.ps'
	set_plot,'X'
endif
	

;	normalize psf
print,'Normalizing PSF ...'
norm = total(psf)
psf = double(psf/norm)

;psf = smooth(psf,3)

;   deconvolve image with maximum likelihood
print,'Deconvolving image with max likelihood ...'
deconv = im
s = size(deconv)
xdim = s(1)
ydim = s(2)

tvscl,rebin(deconv, xdim*repl, ydim*repl, /sample)
xyouts,1,1,'Iteration 0'

Niter = 20
cube = fltarr(xdim,ydim,Niter+1)
cube[*,*,0] = im
for i=1,Niter do begin
	print,'Starting iteration #',i
	Max_Likelihood, im, psf, deconv, FT_PSF=psf_ft
	tvscl,rebin(deconv, xdim*repl, ydim*repl, /sample)
	xyouts,1,1,'Iteration '+strtrim(string(i),2)
	cube[*,*,i] = deconv
endfor
writefits,'deconvolved_images_cube.fits',cube

; Write out the deconvolved image as a FITS file
outfile2 = evtfile + '.deconvolved'
;read,'Enter output file name for deconvolved image: ', outfile2

writefits,outfile2,deconv, header

; print out plots
prt = ''
;read,'Do you want printed output (Y/N)?', prt
prt = 'Y'
if (prt eq 'y' or prt eq 'Y') then begin
	set_plot,'ps'
	device,/landscape
	contour,im,nlevels=10,$
		title='SN1987A: raw 0th order image',xtitle=xlabel,ytitle=ylabel
	contour,deconv,nlevels=10,$
		title='SN1987A: deconvolved 0th order image',$
		xtitle=xlabel,ytitle=ylabel
	device,/close
;	spawn,'lpr -r idl.ps'
        spawn,'mv idl.ps idl2.ps'
	set_plot,'X'
endif

d1 = n_elements(im(*,0))
d2 = n_elements(im(0,*))
projx = dblarr(d1)
projy = dblarr(d2)
imx = dblarr(d1)
imy = dblarr(d2)

for i =0,d1-1 do projx(i) = total(deconv(i,*))
for i =0,d2-1 do projy(i) = total(deconv(*,i))

for i =0,d1-1 do imx(i) = total(im(i,*))
for i =0,d2-1 do imy(i) = total(im(*,i))


print,'Plotting X projection of deconvolved image and original image...'
xvals = findgen(d1)*dx + xllc
;function_1d, id, xvals, projx, /block, dataset_name='Deconvolved',$
;	title='Projections onto X axis', xlabel='0.123 arcsecond pixels'
;function_1d, id, xvals, imx, /block, dataset_name='Data',linestyle=2
;xmanager
	
plot, xvals, projx,$
	title='Projections onto X axis!CDeconvolved (solid) and Data (dashed)',$
	xtitle=xtitle
oplot,xvals,imx,linestyle=2
prt = ''
;read,'Do you want printed output (Y/N)?', prt
prt = 'Y'
if (prt eq 'y' or prt eq 'Y') then begin
	set_plot,'ps'
	device,/landscape
;; Show numbers HERE
        FOR pp=0,(d1-1) DO BEGIN
        if (projx[pp] NE 0.) then print, xvals[pp], projx[pp]
        ENDFOR
	plot, xvals, projx,$
		title='Projections onto X axis!CDeconvolved (solid) and Data (dashed)',$
		xtitle=xtitle
	oplot,xvals,imx,linestyle=2
	device,/close
;	spawn,'lpr -r idl.ps'
        spawn,'mv idl.ps idl3.ps'
	set_plot,'X'
endif

;plot,projy
;oplot,imy

ss_im = im
ss_deconv = deconv


	end






