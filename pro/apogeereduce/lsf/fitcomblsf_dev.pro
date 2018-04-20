function fitcomblsf_dev,par,dp,x=x,y=y,err=err,lsfx=lsfx,mask=mask,porder=porder,$
                               wporder=wporder,wproftype=wproftype,binsize=binsize,$
                               offsetx=offsetx,inpar=inpar,model=lsf2d

; This function creates a model of the 2D LSF
; array.  This helps with the 2D LSF fitting
; done by apvisitcomb.pro

; Construct the proper parameter array
; [binsize, xoffset, Horder, Porder, GHcoefs, wproftype, nWpar, WPorder, Wcoefs]
Horder = n_elements(porder)-1
nGHcoefs = total(Porder+1)
GHcoefs = par[0:nGHcoefs-1]
nWcoefs = total(WPorder+1)
Wcoefs = par[nGHcoefs:*]
nWpar = n_elements(WPorder)
inpar = [binsize, offsetx, Horder, Porder, GHcoefs, wproftype, nWpar, WPorder, Wcoefs]

xcenter = reform(x[*,0])
lsf2d = LSF_GH(lsfx,xcenter,inpar,dlsf,/globalderiv)

; Fix the gaps
;lsf2d[gapbeg[0]:gapend[0],*] = 0.0
;lsf2d[gapbeg[1]:gapend[1],*] = 0.0
;for i=0,n_elements(gapend)-1 do lsf2d[gapbeg[i]:gapend[i],*] = 0.0

if n_elements(mask) gt 0 then lsf2d*=mask

; Put DP in the right form
;  dlsf is [height, center, GHcoefs, Wcoefs]
sz = size(y)
npar = n_elements(par)
dp = dblarr(sz[1]*sz[2],npar)
for i=0,npar-1 do dp[*,i]=dlsf[*,*,2+i]
;dp = dlsf[*,*,2:*]

;apgundef,dp
;stop

dev = (y-lsf2d)/err
dev = (dev)(*)    ; make it 1D

;bd = where(finite(dev) eq 0,nbd)
;if nbd gt 0 then stop
;bd2 = where(finite(dp) eq 0,nbd2)
;if nbd2 gt 0 then stop

return,dev

end
