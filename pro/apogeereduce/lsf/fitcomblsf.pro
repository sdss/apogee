function fitcomblsf,x,par,xlsf=xlsf,mask=mask,gapbeg=gapbeg,gapend=gapend

; This function creates a model of the 2D LSF
; array.  This helps with the 2D LSF fitting
; done by apvisitcomb.pro

xcenter = reform(x[*,0])
lsf2d = LSF_GH(xlsf,xcenter,par,dp,/globalderiv)

; Fix the gaps
;lsf2d[gapbeg[0]:gapend[0],*] = 0.0
;lsf2d[gapbeg[1]:gapend[1],*] = 0.0
for i=0,n_elements(gapend)-1 do lsf2d[gapbeg[i]:gapend[i],*] = 0.0

if n_elements(mask) gt 0 then lsf2d*=mask

return,lsf2d

end
