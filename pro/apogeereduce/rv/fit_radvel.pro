function fit_radvel,x,par,synwave=synwave,synspec=synspec,lsf=lsf,minwave=minwave,maxwave=maxwave,$
                    wext=wext, next=next,gdpts=gdpts

; This program calculates the model spectrum (using the synthetic
; spectrum) to compare to the observed spectrum.
;
; x         The observed wavelength array
; par       The parameters: relative radial velocity
; =synwave  The synthetic wavelength array
; =synspec  The synthetic spectrum array
; =minwave  Minimum wavelength of the observed spectrum
; =maxwave  Maximum wavelength of the observed spectrum
; =wext     The "extended" observed wavlength array
; =next     Number of points "xext" is extended on both edges

cspeed = 2.99792458d5  ; speed of light in km/s

npix = n_elements(x)

; 1.) Correct for the radial velocity
w = synwave*(1.0d0 + par[0]/cspeed)

; 2.) Interpolate onto the observed wavelength grid
interp_spec = SPLINE(w,synspec,wext)

; 3.) Convolve with the LSF
;model = CONVOL(interp_spec,lsf,/center,/edge_truncate,/normalize,/nan) 
model = interp_spec     ; DON'T CONVOLVE, need unsmoothed synspec

; 4.) Now trim the edges
fmodel = model[next:next+npix-1]

; Only use good points
fmodel = fmodel[gdpts]

;stop

return,fmodel

end
