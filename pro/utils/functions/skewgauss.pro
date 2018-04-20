function skewgauss,x,par

; This is a skewed Gaussian
;  See http://en.wikipedia.org/wiki/Skew_normal_distribution

; par = [height, center, sigma, alpha]
; alpha is the skewness of the Gaussian

ht = par[0]
cen = par[1]
sig = par[2]
alpha = par[3]

x2 = (x-cen)/sig  ; rescale 

;y = 2/sqrt(2*!dpi) * exp(-0.5d0 * x2^2 ) * 0.5*(1 + erf(alpha*x2/sqrt(2)))
y = ht * exp(-0.5d0 * x2^2 ) * 0.5*(1 + erf(alpha*x2/sqrt(2)))

;print,par
;dum = where(finite(y) eq 0,nbad)
;if nbad gt 0 then stop

return,y
end
