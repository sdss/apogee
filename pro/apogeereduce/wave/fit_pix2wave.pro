function fit_pix2wave,x,par,chipnum=chipnum,xb=xb

radeg = 180.0d0/!dpi

; The parameters
; 1 xoffset
; 4 sine parameters
; 2 chip gaps
; 6 poly parameters (first one is a zero-point offset)
; The chip number must also be input (1, 2 or 3)

; Break up the input parameters
;sinpars = par[0:3]
;chipgap1 = par[4]
;chipgap2 = par[5]
;polypars = par[6:*]
xoffset = par[0]
sinpars = par[1:4]
chipgap1 = par[5]
chipgap2 = par[6]
polypars = par[7:*]

; Construct XB, chipnum is 1, 2, or 3
; Fix the zero-point at the center of the middle chip
XB = X + (chipnum eq 1)*(-1023.5-2048-chipgap1) + $
         (chipnum eq 2)*(-1023.5) + $
         (chipnum eq 3)*(-1023.5+2048+chipgap2)

;; Construct XB, chipnum is 1, 2, or 3
;XB = X + (chipnum-1)*2048 + $
;         (chipnum eq 2)*chipgap1 + $
;         (chipnum eq 3)*(chipgap1+chipgap2)

; Construct pix2wave parameters
;   set xoffset to zero, since we are inputting XB
;pars2 = [0.0, sinpars, polypars]
pars2 = [xoffset, sinpars, polypars]

; Run pix2wave
out = PIX2WAVE(XB,pars2)

return,out
end
