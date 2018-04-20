function wave2pix,wave,wave0,norev=norev
 
;+
;
; WAVE2PIX
;
; This function returns the pixels corresponding to 
;   an array of wavelength values
; Currently does this using the wavelength array, but
;    would like to do this with inverse coefficients
;
; INPUTS:
;  wave    The array of wavelenth values along the spectrum.
;  wave0   The array of wavelenths corresponding to pixel numbers
;
; OUTPUTS:
;  pix   The pixel array
;
; USAGE:
;  IDL>pix = wave2pix(wave,wave0)
;
; By J. Holtzman  Feb 2012
;-

; The pixel array
pix0=indgen(n_elements(wave0))*1.d0

; Put in ascending order
if not keyword_set(norev) then begin
  wave0=reverse(wave0)
  pix0=reverse(pix0)
endif
nwave0=n_elements(wave0)

; Only interpolate to get pixel values within the array, NaNs elsewhere
gd=where(wave ge wave0[0] and wave le wave0[nwave0-1])
out=SPLINE(wave0,pix0,wave[gd],/double)
pix=dblarr(n_elements(wave))
pix[*]=!values.f_nan
pix[gd]=out

return,pix

end
