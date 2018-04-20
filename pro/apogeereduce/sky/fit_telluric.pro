function fit_telluric,x,par,mspec=mspec,ind=ind

; mspec = [NPix,3] with the telluric model spectra
; par = [scale1, scale2, scale3, cont_poly_coef]
; where cont_poly_coef are polynomial coefficients for
; the continuum
; =ind the indices used for X

mspec1 = reform(mspec[*,0])
mspec2 = reform(mspec[*,1])
mspec3 = reform(mspec[*,2])

;relspec = ( (par[0]*(mspec1-1)+1) + (par[1]*(mspec2-1)+1) + (par[2]*(mspec3-1)+1) )/3.0
relspec = (par[0]*(mspec1-1)+1.0) * (par[1]*(mspec2-1)+1) * (par[2]*(mspec3-1)+1)
if n_elements(ind) gt 0 then relspec=relspec[ind]
cont_coef = par[3:*]
cont = POLY(x,cont_coef)
spec = cont*relspec
;spec = relspec

;stop

return,spec

end
