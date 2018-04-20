function gfunc,x,par,dp,noderiv=noderiv

;+
;
; GFUNC
;
;  This function makes gaussians
;
; each gaussian has 3 parameters height, center, width
; one extra for a constant offset
;
;  INPUT
;   x         Array of X values
;   par       Gaussian parameters
;   /noderiv  Don't return the derivative
;
;  OUTPUT
;   th        Theoretical Y values
;   dp        Derivative of gaussian with the parameters (optional)
;
; When there are problems in this program it returns:
;  th = x*0.+999999.
;  dp = x*0.+999999.
;
; Created by David Nidever April 2005
;-

npar = n_elements(par)
npts = n_elements(x)
ngauss = npar/3
th = 0.
der = 0.

; Bad Input Values
if (n_params() eq 0) or (npts eq 0) or (npar eq 0) or (ngauss eq 0) then begin
  print,'Syntax - dev=gfunc(x,par,dp,/noderiv)'
  if npts gt 0 then th = x*0.+999999. else th = 999999.
  if not keyword_set(noderive) then dp = th
  return,th
endif

; Looping through gaussians
for i=0.,(ngauss-1) do begin
    ipar = par(i*3:i*3+2.)
    ;if ipar(0) lt 0 then ipar(0)=0.
    g = ipar(0)*exp(-0.5*((x-ipar(1))/ipar(2))^2.)
    th = th + g
  ;  der = der + g*((ipar(1)-x)/ipar(2)^2.)
endfor

; Adding power terms
npow = npar-3.*ngauss
if npow gt 0 then for i=0,npow-1 do th = th + par(i+3*ngauss)*x^float(i)
;if npow gt 0 then for i=0,npow-1 do th = th + par(i+3)*x^float(i)
;if npar-ngauss*3 gt 0 then th = th + par(ngauss*3)

; Computing partial derivatives with respect to the parameters
if not keyword_set(noderiv) and arg_present(dp) then begin
  dp = fltarr(npts,npar)
  for i=0.,ngauss-1 do begin
    ipar = par(i*3:i*3+2.)
    ; height
    dp[*,i*3] = exp(-0.5*((x-ipar(1))/ipar(2))^2.)
    ; center
    dp[*,i*3+1.] = dp(*,i*3.)*ipar(0)*(x-ipar(1))/ipar(2)^2.
    ;  width.
    dp[*,i*3+2.] = dp(*,i*3)*ipar(0)*((x-ipar(1))^2.)/ipar(2)^3.
  end ; for i
  ; derivatives for polynomial terms, at most constant+linear terms
  if npow gt 0 then begin
    dp[*,3*ngauss] = 1.0  ; constant
    if npow gt 1 then dp[*,1+3*ngauss]=x
  endif
endif

;stop

;dum = where(finite([th,(dp)(*)]) ne 1,ndum)
;if ndum ne 0 then stop

return,th

end
