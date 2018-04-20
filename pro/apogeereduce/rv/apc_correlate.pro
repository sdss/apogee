; $Id: //depot/idl/IDL_70/idldir/lib/c_correlate.pro#2 $
;
; Copyright (c) 1995-2008, ITT Visual Information Solutions. All
;       rights reserved. Unauthorized reproduction is prohibited.
;+
; NAME:
;       C_CORRELATE
;
; PURPOSE:
;       This function computes the cross correlation Pxy(L) or cross
;       covariance Rxy(L) of two sample populations X and Y as a function
;       of the lag (L).
;
; CATEGORY:
;       Statistics.
;
; CALLING SEQUENCE:
;       Result = C_correlate(X, Y, Lag)
;
; INPUTS:
;       X:    An n-element vector of type integer, float or double.
;
;       Y:    An n-element vector of type integer, float or double.
;
;     LAG:    A scalar or n-element vector, in the interval [-(n-2), (n-2)],
;             of type integer that specifies the absolute distance(s) between
;             indexed elements of X.
;
; KEYWORD PARAMETERS:
;       COVARIANCE:    If set to a non-zero value, the sample cross
;                      covariance is computed.
;
;       DOUBLE:        If set to a non-zero value, computations are done in
;                      double precision arithmetic.
;
; EXAMPLE
;       Define two n-element sample populations.
;         x = [3.73, 3.67, 3.77, 3.83, 4.67, 5.87, 6.70, 6.97, 6.40, 5.57]
;         y = [2.31, 2.76, 3.02, 3.13, 3.72, 3.88, 3.97, 4.39, 4.34, 3.95]
;
;       Compute the cross correlation of X and Y for LAG = -5, 0, 1, 5, 6, 7
;         lag = [-5, 0, 1, 5, 6, 7]
;         result = c_correlate(x, y, lag)
;
;       The result should be:
;         [-0.428246, 0.914755, 0.674547, -0.405140, -0.403100, -0.339685]
;
; PROCEDURE:
;       See computational formula published in IDL manual.
;
; REFERENCE:
;       INTRODUCTION TO STATISTICAL TIME SERIES
;       Wayne A. Fuller
;       ISBN 0-471-28715-6
;
; MODIFICATION HISTORY:
;       Written by:  GGS, RSI, October 1994
;       Modified:    GGS, RSI, August 1995
;                    Corrected a condition which excluded the last term of the
;                    time-series.
;       	     - GGS, RSI, April 1996
;                    Simplified CROSS_COV function. Added DOUBLE keyword.
;                    Modified keyword checking and use of double precision.
;       	     - W. Biagiotti,  Advanced Testing Technologies
;       	     Inc., Hauppauge, NY, July 1997, Moved all
;       	     constant calculations out of main loop for
;       	     greatly reduced processing time.
;   CT, RSI, September 2002. Further speed improvements, per W. Biagiotti.
;                Now handles large vectors and complex inputs.
;-
function APC_Correlate, X, Y, Lag, Covariance = Covariance, Double = doubleIn, nomean=nomean

    compile_opt idl2

    ; Compute the sample cross correlation or cross covariance of
    ; (Xt, Xt+l) and (Yt, Yt+l) as a function of the lag (l).

    ON_ERROR, 2

    typeX = SIZE(X, /TYPE)
    typeY = SIZE(Y, /TYPE)
    nX = N_ELEMENTS(x)

    if (nX ne N_ELEMENTS(y)) then $
        MESSAGE, "X and Y arrays must have the same number of elements."

    ;Check length.
    if (nX lt 2) then $
        MESSAGE, "X and Y arrays must contain 2 or more elements."

    isComplex = (typeX eq 6) or (typeX eq 9) or $
        (typeY eq 6) or (typeY eq 9)

    ;If the DOUBLE keyword is not set then the internal precision and
    ;result are identical to the type of input.
    useDouble = (N_ELEMENTS(doubleIn) eq 1) ? $
        KEYWORD_SET(doubleIn) : $
        (typeX eq 5 or typeY eq 5) or (typeX eq 9 or typeY eq 9)

    fx = finite(X)
    fy = finite(Y)

    ; This will now be in double precision if Double is set.
    ;Xd = x - TOTAL(X, Double = useDouble) / nX ;Deviations
    ;Yd = y - TOTAL(Y, Double = useDouble) / nX
    if not keyword_set(nomean) then begin
      xmn = TOTAL(X, Double=useDouble,/nan) / TOTAL(fx)
      ymn = TOTAL(Y, Double=useDouble,/nan) / TOTAL(fy)
      Xd = x - xmn
      Yd = y - ymn
    endif else begin
      Xd = x
      Yd = y
    endelse

    gdx = where(fx,ngdx,comp=bdx,ncomp=nbdx)
    if nbdx gt 0 then Xd[bdx]=0.0
    gdy = where(fy,ngdy,comp=bdy,ncomp=nbdy)
    if nbdy gt 0 then Yd[bdy]=0.0

    nLag = N_ELEMENTS(Lag)

    Cross = useDouble ? $
        (isComplex ? DCOMPLEXARR(nLag) : DBLARR(nLag)) : $
        (isComplex ? COMPLEXARR(nLag) : FLTARR(nLag))

    num = lonarr(nlag)  ; number of "good" points at this lag
    for k = 0L, nLag-1 do begin
        ; Note the reversal of the variables for negative lags.
        Cross[k] = (Lag[k] ge 0) ? $
            TOTAL(Xd[0:nX - Lag[k] - 1L] * Yd[Lag[k]:*]) : $
            TOTAL(Yd[0:nX + Lag[k] - 1L] * Xd[-Lag[k]:*])
        num[k] = (Lag[k] ge 0) ? $
            TOTAL(fx[0:nX - Lag[k] - 1L] * fy[Lag[k]:*]) : $
            TOTAL(fy[0:nX + Lag[k] - 1L] * fx[-Lag[k]:*])
    endfor
    Cross *= max(num)/(float(num)>1)  ; normalize by number of "good" points

    ; Divide by N for covariance, or divide by variance for correlation.
    ;Cross = TEMPORARY(Cross) / $
    ;    (KEYWORD_SET(Covariance) ? nX : SQRT(TOTAL(Xd^2)*TOTAL(Yd^2)))
    if ngdx gt 2 then rmsx = SQRT(TOTAL(Xd[gdx]^2)) else rmsx=1.0
    if rmsx eq 0.0 then rmsx=1.0
    if ngdy gt 2 then rmsy = SQRT(TOTAL(Yd[gdy]^2)) else rmsy=1.0
    if rmsy eq 0.0 then rmsy=1.0
    Cross = TEMPORARY(Cross) / $
        (KEYWORD_SET(Covariance) ? nX : rmsx*rmsy)
    ;    (KEYWORD_SET(Covariance) ? nX : SQRT(TOTAL(Xd[gdx]^2)*TOTAL(Yd[gdy]^2)))

    return, useDouble ? Cross : $
        (isComplex ? COMPLEX(Cross) : FLOAT(Cross))

end

