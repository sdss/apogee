function sincinterlaced,arr1,arr2,shift,outshift,err1=err1,err2=err2,errout=err,tp=stp

;+
;
; SINCINTERLACED
;
; This does SINC interpolation of two undersampled arrays, but with
; the sampling shifted between the two.  This is considered
; "interlaced sampling" and the equations from Bracewell pgs.201-202
; are used to do the sinc interpolation.
;
; INPUTS:
;  arr1      First undersampled array of data
;  arr2      Second undersampled array of data (with the sampling
;              shifted by "shift").
;  shift     The shift between where the two arrays were sampled.
;              A positive shift means that arr1 is to the right of
;              arr2.
;  outshift  At which points to do the interpolation.  This is assumed
;              to be the array which is on the LEFT.
;
; OUTPUTS:
;  outarr  The SINC interpolated array shifted by "outshift"
;            from arr1. 
;             
; USAGE:
;  IDL>out = sincinterlaced(arr1,arr2,0.55,0.5)
;
; By D.Nidever  March 2010
;-

narr1 = n_elements(arr1)
narr2 = n_elements(arr2)
nshift = n_elements(shift)
noutshift = n_elements(outshift)

; Not enough inputs
if narr1 eq 0 or narr2 eq 0 or nshift eq 0 or noutshift eq 0 then begin
  print,'Syntax - out = sincinterlaced(arr1,arr2,shift,outshift)'
  return,-1
endif

; No interpolation needed
if outshift eq shift then begin
  if shift gt 0 then begin
    if keyword_set(err1) then err=err1*err1
    return,arr1 
  endif else begin
    if keyword_set(err2) then err=err2*err2
    return,arr2
  endelse
endif
if outshift eq 0.0 then begin
  if shift gt 0 then begin
    if keyword_set(err2) then err=err2*err2
    return,arr2 
  endif else begin
    if keyword_set(err1) then err=err1*err1
    return,arr1 
  endelse
endif

; Interlaced sampling
; See Bracewell "The Fourier Transform and Its Applications" pg.201-202
; f(x) = a(x) * (IIIf) + b(x) * (IIIaf)
;  where * is convolution, and IIIf is the "normally" sampled data
;  and IIIaf is the data sampled at an offset of "a"
;  a(x) = sinc(2x) - (pi*cot(a*pi))*x*(sinc(x))^2
;  b(x) = a(-x)
;  this is for critical sampling of 0.5 pixels
;  NEED to use 1-shift in COT!!!

; ARR2 is the reference frame (on the left)
if shift gt 0 then begin
  leftarr = arr2
  rightarr = arr1
  if keyword_set(err1) and keyword_set(err2) then begin
    lefterr = err2*err2
    righterr = err1*err1
  endif

; ARR1 is the reference frame (on the left)
endif else begin
  leftarr = arr1
  rightarr = arr2
  if keyword_set(err1) and keyword_set(err2) then begin
    lefterr = err1*err1
    righterr = err2*err2
  endif
endelse

; Interpolation of LEFTARR
;   outshift pixels to the RIGHT in the LEFTARR frame
dampfac = 3.25
ksize = 21
xkernel1 = ( DINDGEN( ksize ) - ksize/2 ) - outshift
sincx1 = EXP( -( xkernel1 / dampfac )^2 ) * SIN( !dpi*xkernel1 ) / ( !dpi*xkernel1 )
sinc2x1 = EXP( -( xkernel1 / dampfac )^2 ) * SIN( 2*!dpi*xkernel1 ) / ( 2*!dpi*xkernel1 )
;sincx1 = SIN( !dpi*xkernel1 ) / ( !dpi*xkernel1 )
;sinc2x1 = SIN( 2*!dpi*xkernel1 ) / ( 2*!dpi*xkernel1 )
afunc = sinc2x1 - (!dpi/tan((1-shift)*!dpi))*xkernel1*(sincx1)^2
;afunc = EXP( -( xkernel1 / dampfac )^2 ) * afunc
apart = CONVOL(leftarr,afunc,/center,/edge_truncate,/nan)  ;/normalize
if keyword_set(err1) then aerr = CONVOL(lefterr,afunc^2,/center,/edge_truncate,/nan)  ;/normalize

; Interpolation of RIGHTARR
;   outshift-shift to the RIGHT in the RIGHTARR frame
xkernel2 = ( DINDGEN( ksize ) - ksize/2 ) - (outshift-shift)
xkernel2 = -xkernel2    ; b(x) = a(-x)
sincx2 = EXP( -( xkernel2 / dampfac )^2 ) * SIN( !dpi*xkernel2 ) / ( !dpi*xkernel2 )
sinc2x2 = EXP( -( xkernel2 / dampfac )^2 ) * SIN( 2*!dpi*xkernel2 ) / ( 2*!dpi*xkernel2 )
;sincx2 = SIN( !dpi*xkernel2 ) / ( !dpi*xkernel2 )
;sinc2x2 = SIN( 2*!dpi*xkernel2 ) / ( 2*!dpi*xkernel2 )
bfunc = sinc2x2 - (!dpi/tan((1-shift)*!dpi))*xkernel2*(sincx2)^2
;bfunc =  EXP( -( xkernel2 / dampfac )^2 ) * bfunc
bpart = CONVOL(rightarr,bfunc,/center,/edge_truncate,/nan)  ;/normalize
if keyword_set(err2) then berr = CONVOL(righterr,bfunc^2,/center,/edge_truncate,/nan)  ;/normalize

; Now add the two parts
out = apart + bpart
if keyword_set(err1) and keyword_set(err2) then err = aerr + berr

; Are the edges okay?

if keyword_set(stp) then stop

return,out

end
