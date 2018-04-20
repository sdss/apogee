function pix2wave,x,par

;+
;
; PIX2WAVE
;
; This function returns the wavelength solution given
; an array of pixel values and the wavelength coefficients
;
; INPUTS:
;  x      The array of pixel values along the spectrum.
;           This should be relative to the beginning of the CHIP.
;           So the first pixel of the chip should have X=0.
;  par    The wavelength coefficients:
;           Xoffset
;           4 sine parameters
;           6 poly parameters
;           OLD SCHEME -> 7 poly parameters (first one is a zero-point offset)
;
; OUTPUTS:
;  wave   The wavelength array
;
; USAGE:
;  IDL>wave = pix2wave(x,wcoef)
;
; By D.Nidever  March 2010
;-

radeg = 180.0d0/!dpi

; Not enough inputs
if n_elements(x) eq 0 or n_elements(par) eq 0 then begin
  print,'Syntax - wave=pix2wave(x,par)'
  return,-1
endif

xoffset = par[0]
sinpars = par[1:4]
polyoffset = par[5]
polypars = par[6:*]
;polypars = par[5:*]

; Apply the X offset
XB = x+xoffset

; Sine term
sinout = sinpars[0]*( SIN( (XB+sinpars[1])/sinpars[2]/radeg ) + sinpars[3])

; Polynomial term
;  add polynomial X-offset and divide by 3000 to get the
;  values close to unity
out = sinout + POLY( (XB+polyoffset)/3000.,polypars)
;out = sinout + POLY(XB,polypars)

return,out
end
