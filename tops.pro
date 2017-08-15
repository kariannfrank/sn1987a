pro tops,aspect=aspect,landscape=landscape, $
    print=print,filename=filename,margin=margin, $
    full=full, columns=columns, font=font, $
    charthick=charthick,charscale=charscale,tickscale=tickscale,$
    standardsize=standardsize,thickps=thickps, $
    screen=screen,yscale=yscale;,xsize=xsize,ysize=ysize

; think something up
; if keyword_set(screen) then defsysv,'!screen',1
; 
; if keyword_set(!screen) then begin
;  print,'Out kept to screen'
;  goto, SKIPWHOLE
; endif


; future ideas
; x and ymargin

; Default Keyword settings
cm=1      ; standard use of cm's
if not keyword_set(filename) then psfile='plot.ps' else psfile=filename
if not keyword_set(thickps) then thickps=3 ; ADAPT PostScript thickness (2-5)
if not keyword_set(charthick) then charthick=1.
if not keyword_set(charscale) then charscale=1.
if not keyword_set(standardsize) then standardsize=10

if !d.name eq 'X' then begin

;***************************** PSWINDOW INSERT

 ; save old settings
 defsysv,'!oldp',!p
 defsysv,'!oldx',!x
 defsysv,'!oldy',!y
 defsysv,'!oldz',!z
 
  IF KEYWORD_SET(landscape) THEN BEGIN
      pageXsize = 11.0 * 2.54
      pageYsize = 8.5 * 2.54
      inches = 0
  ENDIF ELSE BEGIN
      pageXsize = 8.5 * 2.54
      pageYsize=11.0 * 2.54
      inches = 0
  ENDELSE

;    !p.charsize=charscale*0.378*sqrt(pageYsize/2.54) ; RR size trick
     ; But I prefer font_size
    if keyword_set(columns) then charscale=charscale*columns
    font_size=standardsize*charscale
;    !p.charthick=charthick
    !p.thick=thickps
    !x.thick=thickps
    !y.thick=thickps
;    if keyword_set(font) then !p.font = font  ;doesnt work ;(


  ; Determine the margin of the window.
  ; Thijs Addin!!!

  if not keyword_set(margin) then margin=1.27 ; .5 inch
  IF margin LT 1.27 OR margin GE (pageXsize > pageYsize) THEN margin = 1.27 

   ; Get the aspect ratio of the current display window.

  aspectRatio = FLOAT(!D.Y_VSIZE) / !D.X_VSIZE
  ;JMK
  if keyword_set(aspect) then AspectRatio=aspect
  if keyword_set(full) then AspectRatio=(pageYsize-margin)/(pageXsize-margin)
        ; Fit to the largest plot possible.

  ; JMK
  xsize = pageXsize - margin
  ysize = xsize * aspectRatio
  IF ysize GT pageYsize or keyword_set(yscale) THEN BEGIN
   print,'RESCALING IMAGE'
   ysize = pageYsize - margin
   xsize = ysize / aspectRatio
  ENDIF
  
;   if not keyword_set(tickscale) then  begin
;   !x.ticklen=0.003*(xsize/2.54)*charscale
;   !y.ticklen=0.003*(ysize/2.54)*charscale
;   endif

  if keyword_set(tickscale) then  begin
  !x.ticklen=0.002*(xsize/2.54)*charscale*tickscale
  !y.ticklen=0.002*(xsize/2.54)*charscale*tickscale*aspectratio
  endif


   ; Calculate the offsets.

  xoffset = (pageXsize - xsize) / 2.0
  yoffset = (pageYsize - ysize) / 2.0
  IF KEYWORD_SET(landscape) THEN yoffset = pageYSize - yoffset

;***************************** END PSWINDOW INSERT
 set_plot,'ps', /copy, /interpolate
 device, font_size=font_size, XSIZE=xsize, YSIZE=ysize, $ 
        XOFFSET=xoffset, YOFFSET=yoffset, INCHES=inches, $
	  landscape=keyword_set(landscape),filename=psfile,/color
 print,'%% Output now to   : '+psfile
 print,'%% with AspectRatio: '+strtrim(string(aspectratio),2)
endif else begin
 device,/close 
 set_plot,'x'
 print,'%% Output back to X-window'
 if keyword_set(print) then begin
   print,'%% Sending to Printer'
   spawn,'lpr -h '+psfile
   spawn,'lpq'
 endif  
 ;restore old settings
 !p= !oldp
 !x= !oldx
 !y= !oldy
 !z= !oldz

endelse

SKIPWHOLE:
; simple way of avoiding tops when needed

end
