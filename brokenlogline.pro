; function for a broken line, to be used by MPFIT

FUNCTION brokenlogline, X, P

  ;P0 = intercept
  ;P1 = segment1 slope
  ;P2 = changepoint1
  ;P3 = segment2 slope
  ;P4 = changepoint2
  ;P5 = segment3 slope 

  ; intercept2 and 3
  b2 = P[2]*(P[1]-P[3])+P[0]
  b3 = P[4]*(P[3]-P[5])+b2

  Y = FLTARR(N_ELEMENTS(X))
  ;print, n_elements(X)

  FOR i=0,N_ELEMENTS(X)-1 DO BEGIN
;     IF X[i] LE -1.0*P[0]/P[1] THEN BEGIN ;because log can't be negative
;        Y[i] = 0.0
;     ENDIF ELSE BEGIN
        IF X[i] LT P[2] THEN BEGIN
           Y[i] = P[1]*X[i]+P[0]
        ENDIF ELSE BEGIN
           IF X[i] LT P[4] THEN BEGIN
              Y[i] = P[3]*X[i]+b2
           ENDIF ELSE BEGIN
              Y[i] = P[5]*X[i]+b3
           ENDELSE 
        ENDELSE
 ;    ENDELSE
  ENDFOR

  RETURN, Y

END
