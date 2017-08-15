function modelpileup, xbin, ybin, a,model=model;

if not keyword_set(model) then model='lobes'
;create model with ring, lobes and pointsource. Then convolve this
;model with the PSF to get what it would like as raw data. 
  x0 = a[2]
  y0 = a[3]
  counts = a[22]
  xsrcf = a[23] 
  ysrcf = a[24]
  wsrcf = a[25]
  fsrcf = a[26]  ;add pointsource AFTER projecting 
   
  common grid_com, xpos, ypos

  if model eq 'lobes' then begin
     bla = sn1987a_model_lobes(xpos, ypos, a[0:23])
  endif else begin
     if model eq 'lobsrc' then bla = sn1987a_model_lobes_ptsrc(xpos, ypos, a[0:26])
  endelse 
     
  bigmod = dblarr(398,398);zodat het uiteindelijke beter matcht ("so that the final better matches")
  bigmod[99:298,99:298]  = bla
  proj = project2(bigmod, 81.2d, 43d, [x0,y0])   ;project model PA = 81.2 and inc is 43 graden
  rsrc = sqrt((xpos - xsrcf)^2d + (ypos - ysrcf)^2d)
  imsrc = fsrcf[0]*exp(-rsrc^2d/(2d*wsrcf^2d))/(total(exp(-rsrc^2d/(2d*wsrcf^2d))))
 

  ;dum     = (float(indgen(398)) - 199.0) - 0.5
  ;xposhal = rebin(dum,398,398,/sample)
  ;yposhal = transpose(xposhal) 
  ;rihal = 40 ;inner ring to halo
  ;rc    = sqrt((xposhal - xhal)^2d + (yposhal - yhal)^2d)
  ;stop
  ;emiss = dblarr(max(rc))   ;emissivity function
  ;emiss[rihal:rohal] = 1    ;constant emissivity
  ;r     = dindgen(max(rc))  ;radius
  ;surf  = (halo2(emiss, 2d*!dpi)/(2d*r+ 1d)) ;for cartesian grid
  ;halo  = interpol(surf, r, rc) ;interpol from 1d to 2d
  ;halo  = halo*fhal/total(halo) ;normalising
;  stop
  
 ; proj = proj*(counts - fhal - fsrcf)/total(proj);to keep total # of counts
 proj = proj;*(counts - fhal - fsrcf)/total(proj)
 if model eq 'lobes' then begin
    proj[99:298,99:298] = proj[99:298,99:298] ; + imsrc ;use this line for lobes-only image
 endif else begin
    proj[99:298,99:298] = proj[99:298,99:298] + imsrc 
 endelse

  return, proj
 
 end
