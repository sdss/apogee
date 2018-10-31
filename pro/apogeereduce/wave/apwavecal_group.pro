pro apwavecal_group,mlinestr,flinestr,parstr,coefstr,npoly=npoly,nfibers=nfibers,stp=stp

;; This program fits wavelength solutions for a single group of
;; adjacent exposures for apwavecal_multi.pro.

if n_elements(npoly) eq 0 then npoly=8
if n_elements(nfibers) eq 0 then nfibers=300

  ; Step 1 - Get the chipgaps by fitting a polynomial to each fiber
  ;------------------------------------------------------------------
  resid = mlinestr.model_wave*0
  npoly = 8 ;6 ; 4
  parstr = replicate({group:0,fiber:0,pars:fltarr(npoly+2)+999999.0,perror:fltarr(npoly+2)+999999.0,ypos:999999.0,sig:999999.0,rms:999999.0},nfibers)
  parstr.group = mlinestr[0].group
  for i=0,nfibers-1 do begin
    ind = where(mlinestr.fiber eq i,nind)
    if nind gt 0 then begin
      xx = mlinestr[ind].x
      yy = mlinestr[ind].model_wave
      ;err = yy*0.0+1.0
      err = mlinestr[ind].gfit_perror[1] > 0.1 
      chipnum = mlinestr[ind].chipnum
      initpars = [140.0d, 150.,fltarr(npoly)]
      fa = {chipnum:chipnum}
      pars1 = MPFITFUN('func_chipgap_poly',xx,yy,err,initpars,status=status,dof=dof,$
                         functargs=fa,bestnorm=chisq,perror=perror,yfit=yfit,/quiet)  
      yfit1 = func_chipgap_poly(xx,pars1,chipnum=chipnum,xb=xb)
  
      ; Remove outliers and refit
      diff = yy-yfit1
      sig = mad(diff)
      gd = where(abs(diff) lt 2.5*sig,nbd)
  
      fa2 = {chipnum:chipnum[gd]}
      pars2 = MPFITFUN('func_chipgap_poly',xx[gd],yy[gd],err[gd],pars1,status=status2,dof=dof2,$
                       functargs=fa2,bestnorm=chisq2,perror=perror2,yfit=yfit2,/quiet)      
      yfit2 = func_chipgap_poly(xx,pars2,chipnum=chipnum)
      rchisq2 = chisq2/dof2
      diff2 = mlinestr[ind].model_wave-yfit2
      sig2 = mad(diff2)
  
      resid[ind] = yy-yfit2
      parstr[i].fiber = i
      parstr[i].pars = pars2
      parstr[i].perror = perror2
      parstr[i].ypos = median(mlinestr[ind].ypos)
      parstr[i].sig = mad(yy-yfit2)
      parstr[i].rms = stddev(yy[gd]-yfit2[gd])
    endif else begin
    ;; no lines for this fiber
      parstr[i].fiber = i
      parstr[i].pars = 999999.
      parstr[i].perror = 999999.
      ;parstr[i].ypos = median(mlinestr[ind].ypos)
      parstr[i].sig = 999999.
    endelse
  endfor
  
  ; fit chipgaps with YPOS
  gdfib = where(parstr.sig lt 1000,ngdfib)
  chipgap1_coef = ap_robust_poly_fit(parstr[gdfib].ypos,parstr[gdfib].pars[0],1)
  chipgap1 = poly(parstr.ypos,chipgap1_coef)
  chipgap2_coef = ap_robust_poly_fit(parstr[gdfib].ypos,parstr[gdfib].pars[1],1)
  chipgap2 = poly(parstr.ypos,chipgap2_coef)
  
  ; Get XB using the linear fits to the chipgaps
  for i=0,nfibers-1 do begin
    ind = where(mlinestr.fiber eq i,nind)
    xx = mlinestr[ind].x
    chipnum = mlinestr[ind].chipnum
    pars = [chipgap1[i],chipgap2[i],fltarr(5)]
    yfit = func_chipgap_poly(xx,pars,chipnum=chipnum,xb=xb)
    mlinestr[ind].xb = xb
  endfor
  
  ; These residuals are already at 0.033A
  ;  even if we use npoly=4
  print,'Initial poly fits to each fiber. Sig = ',stringize(mad(resid),ndec=4),' A'
    
  ; Step 2 - Fit polynomial for each fiber the constraining the chipgaps
  ;---------------------------------------------------------------------
  resid = mlinestr.model_wave*0
  coefstr = replicate({group:0,fiber:0,chipgap1:0.0,chipgap2:0.0,coef:fltarr(npoly),perror:fltarr(npoly),ypos:0.0,sig:0.0,rms:0.0},nfibers)
  coefstr.group = mlinestr[0].group
  flinestr = mlinestr
  for i=0,nfibers-1 do begin
    ind = where(mlinestr.fiber eq i,nind)
    if nind gt 0 then begin
      xx = mlinestr[ind].x
      yy = mlinestr[ind].model_wave
      ;err = yy*0.0+1.0
      err = mlinestr[ind].gfit_perror[1] > 0.1 
      chipnum = mlinestr[ind].chipnum
      initpars = [chipgap1[i],chipgap2[i],parstr[i].pars[2:*]]
      fa = {chipnum:chipnum}
      parinfo = replicate({limited:[0,0],limits:[0.0,0.0],fixed:0},n_elements(initpars))
      parinfo[0:1].fixed = 1
      ;parinfo[0].limited = 1 & parinfo[0].limits = [-1.0,1.0]+chipgap1[i]
      ;parinfo[1].limited = 1 & parinfo[1].limits = [-1.0,1.0]+chipgap2[i]
      pars1 = MPFITFUN('func_chipgap_poly',xx,yy,err,initpars,status=status,dof=dof,$
                       functargs=fa,bestnorm=chisq,parinfo=parinfo,perror=perror,yfit=yfit,/quiet)  
      yfit1 = func_chipgap_poly(xx,pars1,chipnum=chipnum,xb=xb)
  
      ; Remove outliers and refit
      diff = yy-yfit1
      sig = mad(diff)
      gd = where(abs(diff) lt 2.5*sig,nbd)
  
      fa2 = {chipnum:chipnum[gd]}
      pars2 = MPFITFUN('func_chipgap_poly',xx[gd],yy[gd],err[gd],pars1,status=status2,dof=dof2,$
                       functargs=fa2,bestnorm=chisq2,parinfo=parinfo,perror=perror2,yfit=yfit2,/quiet)      
      yfit2 = func_chipgap_poly(xx,pars2,chipnum=chipnum)
      rchisq2 = chisq2/dof2
      diff2 = mlinestr[ind].model_wave-yfit2
      sig2 = mad(diff2)
  
      ;plot,xb,diff2,ps=8,xs=1,ys=1,tit='Fiber = '+strtrim(i,2)
  
      flinestr[ind].wave_fit = yfit2
  
      resid[ind] = yy-yfit2
      coefstr[i].fiber = i
      coefstr[i].chipgap1 = chipgap1[i]
      coefstr[i].chipgap2 = chipgap2[i]
      coefstr[i].coef = pars2[2:*]
      coefstr[i].perror = perror2[2:*]
      coefstr[i].ypos = median(mlinestr[ind].ypos)
      coefstr[i].sig = sig2
      coefstr[i].rms = stddev(diff2)
      ;stop
   endif else begin
   ;; No lines for this fiber
      coefstr[i].fiber = i
      coefstr[i].chipgap1 = chipgap1[i]
      coefstr[i].chipgap2 = chipgap2[i]
      coefstr[i].coef = 999999.
      coefstr[i].perror = 999999.
      ;coefstr[i].ypos = median(mlinestr[ind].ypos)
      coefstr[i].sig = 999999.
      coefstr[i].rms = 999999.
    endelse  
  endfor
  
  print,'Poly fits to each fiber. Sig = ',stringize(mad(resid),ndec=4),' A'

  if keyword_set(stp) then stop

;stop

end
