function skyfit_lsf_gh_single,x,par,coef=coef

; This helps fit the LSF to a line
; COEF  the Gauss-Hermite LSF parameters

area = par[0]
center = par[1]

y = area * LSF_GH(x,center,coef)

;stop

return,y

end
