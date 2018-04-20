function gaussbin,x,par,dx=dx

;+
;
; GAUSSBIN
;
; This function returns a binned Gaussian
;  par = [height, center, sigma]
;
; INPUTS:
;  x       The array of X-values.
;  par     The Gaussian parameters. par=[height,center,sigma,yoffset]
;  =dx     The width of each "pixel" (scalar).  The default is 1.0.
;
; OUTPUTS:
;  geval   The binned Gaussian in the pixel
;
; USAGE:
;  IDL>f=gaussbin(x,par)
;
; By D.Nidever  March 2010
;-

ht = par[0]
cen = par[1]
sig = par[2]
if n_elements(par) gt 3 then yoffset=par[3] else yoffset=0.0
if n_elements(dx) eq 0 then dx=1.0d0  ; the width of the bins/pixels

xcen = x-cen           ; relative to the center
x1cen = xcen - 0.5d0*dx  ; left side of bin
x2cen = xcen + 0.5d0*dx  ; right side of bin

t1cen = x1cen/(sqrt(2.0d0)*sig)  ; scale to a unitless Gaussian
t2cen = x2cen/(sqrt(2.0d0)*sig)

; For each value we need to calculate two integrals
;  one on the left side and one on the right side

; Evaluate each point
;   ERF = 2/sqrt(pi) * Integral(t=0-z) exp(-t^2) dt
;   negative for negative z
geval_lower = ERF(t1cen)
geval_upper = ERF(t2cen)

geval = ht*sqrt(2.0d0)*sig * sqrt(!dpi)/2.0 * ( geval_upper - geval_lower )

; Add YOFFSET
geval = geval + yoffset

return,geval

end
