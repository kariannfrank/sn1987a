function project2,image,ang,inc, cen
;function to project the model of 1987A to the desired
;inclination. Assumes the center of 1987A is at the center of the
;image. Largely stolen from gal_flat
;help, cen
;stop
 
  if ( N_params() lt 3 ) then begin  
      print,'Syntax - result = gal_flat( image, ang, inc, [ cen, /INTERP ])'
      print,'ANG - Position Angle of major axis (degrees)'                
      print,'INC - Inclination of galaxy (degrees)'
      return, -1
  endif 

  if not keyword_set( INTERP ) then interp = 0

  angr = (ang+90)/!RADEG
 
  tanang = tan(angr)
  cosang = cos(angr)
  cosinc = cos(inc/!RADEG)
;                                    Parameters of image
  dims = SIZE(image)

  if N_elements(cen) NE 2 then begin 

      xcen = dims[1]/2.0                  ;Center
      ycen = dims[2]/2.0
  ;    if not !QUIET then message, 'center 1987A assumed in image center',/CONT

  endif else begin

      xcen = dims[1]/2.0 + cen[0]
      ycen = dims[2]/2.0 + cen[1]

  endelse
;                                    Equation of rotation axis
  b = ycen - xcen*tanang
;                                    Fiducial grid (as in ROT_INT)   
  gridx = xcen + [ [-1,1], [-1,1] ] * dims[1]/6.0
  gridy = ycen + [ [-1,-1], [1,1] ] * dims[2]/6.0      
;                                    Distorted version of grid
  yprime = gridx*tanang + b            ;Equation of major axis
  r0 = (gridy-yprime)*cos(angr)        ;Dist of control pts to major axis
  delr = -r0*(1d -cosinc)               ;;Correct distance for inclination ->
                                       ;EAH: * -> / to get inverse
 
  dely = -delr*cos(angr)               
  delx =  delr*sin(angr)
  distx = gridx + delx
  disty = gridy + dely
;                                    Parameters of undistorted grid
  x0 = dims[1]/3.0
  y0 = dims[2]/3.0
  dx = x0                              ;In this case only
  dy = y0

                                 
  polywarp, distx, disty, gridx, gridy, 1, kx, ky
  RETURN,poly_2d( image, kx, ky, interp, MISSING = 0) 
 

end
