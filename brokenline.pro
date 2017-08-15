; function for a broken line, to be used by MPFIT

FUNCTION brokenline, X, P

  ;P0 = early slope, P1 = late slope, P2 = intercept, P3 = change point

  ;late intercept
  b2 = P[3]*(P[0]-P[1])+P[2]

;  print, n_elements(X)

  Y = FLTARR(N_ELEMENTS(X))
  ;print, n_elements(X)

  FOR i=0,N_ELEMENTS(X)-1 DO BEGIN

     IF X[i] LT P[3] THEN BEGIN
        Y[i] = P[0]*X[i]+P[2]
     ENDIF ELSE BEGIN
        Y[i] = P[1]*X[i]+b2
     ENDELSE
  ENDFOR

  RETURN, Y

END
