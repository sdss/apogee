function func_multi_poly,x,par,group=group,chipnum=chipnum,xb=xb

;; Use chipgaps to get XB
chipgap1 = par[0]
chipgap2 = par[1]
chipxoff = [-1023.5-2048-chipgap1, -1023.5, -1023.5+2048+chipgap2]
xb = x + chipxoff[chipnum-1]

;; Make small corrections for each chip/group
;;  always additive corrections to xb
ngroups = max(group)-min(group)+1
grpchipxoff = par[2:ngroups*3+1]

; The green chip offset is for all three chips of that group
; add this value to the blue and red chips
grpchipxoff[indgen(ngroups)*3] += grpchipxoff[indgen(ngroups)*3+1]
grpchipxoff[indgen(ngroups)*3+2] += grpchipxoff[indgen(ngroups)*3+1]

; Now add the corrections
; [grp1_chip1, grp1_chip2, grp1_chip3, grp2_chip1, ...]
xb += grpchipxoff[(group-1)*3+chipnum-1]


; Polynomial coefficients
polypars = par[ngroups*3+2:*]

return,poly(xb,polypars)
end
