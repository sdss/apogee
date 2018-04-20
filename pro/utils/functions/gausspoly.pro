function gausspoly,x,par

; par = [heigh, center, sigma, poly_coef0, poly_coef1, ...]

gpar = par[0:2]
ppar = par[3:*]

y = gpar[0]*exp(-0.5d0 * ((x-gpar[1])/gpar[2])^2 ) + poly(x,ppar)

return,y
end
