;+
; NAME:
;    sincint
; PURPOSE: (one line)
;    Sinc interpolation of a 1-D vector of data.
; DESCRIPTION:
;
; CATEGORY:
;    Numerical
; CALLING SEQUENCE:
;    result = sincint( x, nres, f )
; INPUTS:
;    x  : Independent variable values for which f is to be interpolated.
;    nres : number of pixels per resolution element (2 = critical sampling)
;    f  : Vector of function values (dependent variable).
;
; OPTIONAL INPUT PARAMETERS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;    Interpolated value(s).
; COMMON BLOCKS:
;
; SIDE EFFECTS:
;
;  WARNING!!!  The output type of the function is derived from the rank,
;      size, and type of the "x" input vector.  It is NOT controlled in any
;      way by the type of "f".
;
; RESTRICTIONS:
;
; PROCEDURE:
;
; MODIFICATION HISTORY:
;    Written by Doug Loucks, Lowell Observatory, September, 1993.
;    Adapted from the IDL function sshift.pro written by Marc Buie.
;    01/14/94, DWL, Documentation update.
;   2/12: Jon Holtzman, modified to allow independent variable input
;      (e.g. for significantly oversampled data)
;-
FUNCTION sincint, x, nres, f, ferr

dampfac = 3.25 * nres/2.
ksize   = fix(21 * nres/2.)
if (ksize mod 2) eq 0 then ksize+=1
nhalf = ksize/2

nx = N_ELEMENTS( x )
nf = N_ELEMENTS( f )

ix = FIX( x )
fx = x - ix
i = WHERE( finite(fx) eq 1, counti )
;z = WHERE( fx EQ 0, countz )
;i = WHERE( fx NE 0, counti )

r = x * 0

if n_elements(ferr) gt 0 then begin
  fvar=ferr*ferr
  err = x * 0
endif
;
; with windowing for oversampled data, want to reconstruct even at integral pixel locations
;IF countz NE 0 THEN BEGIN
;   ;There are integral values of the independent variable. Select the function
;   ;values directly for these points.
;   r[ z ] = f[ ix[z] ]
;ENDIF

;IF counti NE 0 THEN BEGIN
   ;Use sinc interpolation for the points having fractional values.
   FOR point=0, counti-1 DO BEGIN
      xkernel = ( FINDGEN( ksize ) - nhalf ) - fx[ i[ point ] ]
      xkernel /= (nres/2.)
      u1 = xkernel / dampfac
      u2 = !pi * xkernel
      sinc = EXP( -( u1*u1 ) ) * SIN( u2 ) / u2
      sinc /= (nres/2.)
;if point eq 0 then begin
;plot,sinc
;stop
;endif
      lobe = ( INDGEN( ksize ) - nhalf ) + ix[ i[point] ]
      vals = FLTARR( ksize )
      vars = FLTARR( ksize )
      w = WHERE( lobe LT 0, count )
      ;IF count NE 0 THEN vals[w] = f[0]
      IF count NE 0 THEN begin
        vals[w] = !values.f_nan
        if n_elements(ferr) gt 0 then vars[w] = !values.f_nan
      endif
      w = WHERE( lobe GE 0 AND lobe LT nf, count )
      IF count NE 0 THEN begin
        vals[w] = f[ lobe[w] ]
        if n_elements(ferr) gt 0 then vars[w] = fvar[ lobe[w] ]
      endif
      w = WHERE( lobe GE nf, count )
      ;IF count NE 0 THEN vals[w] = f[ nf-1 ]
      IF count NE 0 THEN begin
        vals[w] = !values.f_nan
        if n_elements(ferr) gt 0 then vars[w] = !values.f_nan
      endif
      r[ i[ point ] ] = TOTAL( sinc * vals )
      if n_elements(ferr) gt 0 then err[ i[ point ] ] = TOTAL( sinc * sinc * vars )
   ENDFOR
;ENDIF
   if n_elements(ferr) gt 0 then ferr=sqrt(err)

RETURN, r

END
