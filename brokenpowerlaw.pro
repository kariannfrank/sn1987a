; function for a broken line, to be used by MPFIT

FUNCTION brokenline, X, P

  ;P0 = first power law coefficient
  ;P1 = first power law index
  ;P2 = first break point
  ;P3 = second power law coefficient
  ;P4 = second power law index
  ;P5 = second breakpoint
  ;P6 = third power law coefficient
  ;P7 = third power law index

  P1 = [P[0],P[1]] ;params for first powerlaw
  P2 = [P[3],P[4]]
  P3 = [P[6],P[7]]

  b1 = powerlaw(P[2],P1) ;y-value at first break point
  b2 = powerlaw(P[5],P2)+b1 ;y-value at second break point

;  print, n_elements(X)

  Y = FLTARR(N_ELEMENTS(X))
  ;print, n_elements(X)

  FOR i=0,N_ELEMENTS(X)-1 DO BEGIN

     IF X[i] LT P[2] THEN BEGIN
        Y[i] = powerlaw(X[i],P1)
     ENDIF ELSE BEGIN
        IF X[i] LT P[5] THEN BEGIN
           Y[i] = powerlaw(X[i],P2)+b1
        ENDIF ELSE BEGIN
           Y[i] = powerlaw(X[i],P3)+b2
        ENDELSE 
     ENDELSE
  ENDFOR

  RETURN, Y

END
