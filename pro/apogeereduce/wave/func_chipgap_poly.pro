function func_chipgap_poly,x,par,chipnum=chipnum,xb=xb

;XB = X + (chipnum eq 1)*(-1023.5-2048-chipgap1) + $
;         (chipnum eq 2)*(-1023.5) + $
;         (chipnum eq 3)*(-1023.5+2048+chipgap2)
chipgap1 = par[0]
chipgap2 = par[1]
XB = X + (chipnum eq 1)*(-1023.5-2048-chipgap1) + $
         (chipnum eq 2)*(-1023.5) + $
         (chipnum eq 3)*(-1023.5+2048+chipgap2)

return,poly(xb,par[2:*])
end
