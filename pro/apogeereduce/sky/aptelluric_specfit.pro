pro aptelluric_specfit,frame,starind,telmodelstr,outstr,specfitopt=specfitopt,pl=pl,stp=stp,$
                       convolved_telluric=convolved_telluric,template=template,bestmod=bestmod,test=test
;+
;
; APTELLURIC_SPECFIT
;
; This fits the telluric spectra to a single stellar spectrum.
;
; INPUTS:
;  frame       A structure with the header/data information for an
;                  undersampled frame that has been wavelength
;                  calibrated and airglow line subtracted (not essential)
;  starind     The IDL index of the star to use.
;  telmodelstr A structure with the telluric model spectrum
;  =specfitopt The spectral fitting option. 1) Chi-square minimization,
;               2) RMS minimization.  Default telspecopt=2.
;              Note that (1) fits all species simultaneously, (2) fits one species at a time
;              The latter is needed to separately select the best of multiple models for a given species
;  convolved_telluric  Gives the pre-convolved model telluric spectra
;  template    If given, divide the stellar spectra by the template before doing telluric fitting
;  bestmod     If given, gives the array of the best model index for each species
;  /pl         Make some plots
;  /stp        Stop at the end of the program.
;
; OUTPUTS:
;  outstr      The output structure with various parameters.
;  =error      The error message if one occurred.
;
; USAGE:
;  IDL>aptelluric_specfit,frame,starind,plugmap,telmodelstr,telstr,specfitopt=specfitopt
;
; By D. Nidever  December 2011
; Modified by Holtz. April 2015
;-

; Not enough inputs
if n_elements(frame) eq 0 or n_elements(starind) eq 0 then begin
  print,'Syntax - aptelluric_specfit,frame,starind,telmodelstr,telstr,specfitopt=specfitopt,pl=pl,stp=stp'
  return
endif

if n_elements(telmodelstr) eq 0 and not keyword_set(convolved_telluric) then begin
  print,' Must give a valid telluric model structure if you dont give a pre-convolved model'
  return
endif

if n_elements(specfitopt) eq 0 then specfitopt=2
; specfitopt=1  Chi-squared minimzation of telluric-dominated pixels
; specfitopt=2  RMS minimzation of large regions

sz = size(frame.chipa.flux)
npix = sz[1]
pix = findgen(npix)
nfibers = sz[2]

; Getting the spectrum, concatenate the chips
spec = [frame.(0).flux[*,starind], frame.(1).flux[*,starind], $
        frame.(2).flux[*,starind] ]
errspec = [frame.(0).err[*,starind], frame.(1).err[*,starind], $
           frame.(2).err[*,starind] ]
errspec = errspec > 1
mask = [frame.(0).mask[*,starind], frame.(1).mask[*,starind], $
        frame.(2).mask[*,starind] ]
wave = [frame.(0).wavelength[*,starind], frame.(1).wavelength[*,starind], $
        frame.(2).wavelength[*,starind] ]
skyspec = [frame.(0).sky[*,starind], frame.(1).sky[*,starind], $
        frame.(2).sky[*,starind] ]

; Initialize the output structure
outstr = {specfitopt:0,fiber:-1,snr:0.0,par:fltarr(3),cont_coef:fltarr(7),bestmod:[-1,-1,-1],$
          status:0,rchisq:99.99,rchisq2:99.99,x:fltarr(3*npix),wave:dblarr(3*npix),$
          spec:fltarr(3*npix),nspec:fltarr(3*npix),cont:fltarr(3*npix),telluric:fltarr(3*npix)}
outstr.wave = wave
outstr.spec = spec

; Remove continuum from skyspec
skycont = skyspec*0
for j=0,2 do skycont[j*npix:(j+1)*npix-1]=MEDFILT1D(skyspec[j*npix:(j+1)*npix-1], 51,/edge)
skyspec_lines = skyspec-skycont
; we will avoid fitting on pixels where sky makes a significant contribution
badskythresh=0.5

; Signal/Noise
gdpix = where((mask and badmask()) eq 0,ngdpix)
if ngdpix gt npix/2 then snr = median(spec[gdpix]/errspec[gdpix]) else snr=0
outstr.snr = snr

; Bad spectrum
if snr lt 10 or median(spec) lt 0.0 then begin
  print,'Bad spectrum.  Skipping'
  goto,BOMB1
endif

; Convolved model spectrrum
;--------------------------

; Pre-convolved
if keyword_set(convolved_telluric) then begin
  x = dblarr(3*npix)
  for j=0,2 do begin
    wcoef1 = reform( frame.(0).wcoef[starind,*] )
    wcoef = frame.(j).wcoef[starind,*]
    xoff = wcoef[0]-wcoef1[0]
    x[j*npix:j*npix+npix-1] = dindgen(npix)+xoff
  endfor
  mspec = reform(convolved_telluric[*,starind,*,*]) 
endif else begin
  ; On-the-fly convolution
  stop,'On-the-fly convolution not implemented for multi-scale telluric models!'
  apgundef,lsfpars,wcoef
  for j=0,2 do PUSH,lsfpars,frame.(j).lsfcoef[starind,*]
  for j=0,2 do PUSH,wcoef,frame.(j).wcoef[starind,*]
  APTELLURIC_SINGLECONVOLVE,lsfpars,wcoef,telmodelstr,x,mspec
endelse

; do we have multiple models for each species? If so we will try all of them, unless
;   we have been given the one to use
sz=size(mspec)
if sz[0] eq 2 then nscale=1 else $
  if n_elements(bestmod) gt 0 then nscale=1 else nscale=sz[2]

; remove stellar spectrum?
if keyword_set(template) then spec/=template

; Different method options
CASE specfitopt of


;-------------------------------------------------------------------------
; OPTION 1: CHI-SQUARED Fitting to telluric-dominated pixels
;-------------------------------------------------------------------------
1: begin

  ; choose the appropriate telluric sub-model which must have been previously determined (e.g., with specfitopt=2)
  if n_elements(bestmod) gt 0 then $
    bestmspec=[[mspec[*,0,bestmod[0]]],[mspec[*,1,bestmod[1]]],[mspec[*,2,bestmod[2]]]] $
  else if sz[0] eq  2 then $
    bestmspec=mspec $
  else $
    stop,'need to specify bestmod array for multi-scale telluric models'

  ; First cut telluric estimate
  telnorm = fit_telluric(x,[1.0,1.0,1.0,1.0,0.0],mspec=bestmspec)

  ; Initialize the fitting arrays
  initpar = [1.0,1.0,1.0, 1.0, 0.0, 0.0, 0.0]
  npar = n_elements(initpar)
  parinfo = REPLICATE({limited:[0,0],limits:[0.0,0.0],fixed:0},npar)
  parinfo[0:2].limited = 1
  parinfo[0:2].limits = [0.0,10.0]

  ; Iterate until we converge
  flag = 0
  count = 0
  lastpar = initpar
  smbin = 200  
  WHILE (flag ne 1) do begin

    ; Get the continuum and make the normalized spectrum
    smspec = MEDFILT1D(spec/telnorm,smbin,/edge) > 1
    specnorm = spec / smspec
    errspecnorm = errspec / smspec  > 1.0/median(smspec)
    ; set a floor of 0.02 for error to recognize systematic errors in telluric correction
    j=where(errspecnorm lt 0.02,nj)
    if nj gt 0 then errspecnorm[j]=0.02

    ; Get the pixels to fit
    gd = where(telnorm lt 0.99 and specnorm gt 0.0 and $
           skyspec_lines lt badskythresh*smspec and (mask and badmask()) eq 0,ngd)
    ; mask first 4 H Brackett lines in APOGEE region
    if ~keyword_set(template) then begin
      for ir=11,14 do begin
       wmask=rydberg(4,ir)
       bd=where(abs(wave[gd]-wmask) lt 10,nbd,comp=gd2,ncomp=ngd)
       gd=gd[gd2]
      endfor
    endif
    fa = {mspec:bestmspec,ind:gd}
    if ngd eq 0 then begin
      print,'No good pixels for fitting'
      goto,BOMB1
    endif

    ; Fit the relative scalings of the three models to a normalized
    ; spectrum with MPFITFUN.PRO, even though it is a simple linear problem
    ;---------------------------------------------------------------
    par = MPFITFUN('fit_telluric',x[gd],specnorm[gd],errspecnorm[gd],lastpar,yfit=yfit1,functargs=fa,$
                   parinfo=parinfo,status=status,perror=perror1,dof=dof1,bestnorm=chisq,/quiet)
    rchisq = chisq/dof1

    ; Current fit
    telnorm = fit_telluric(x,[par[0:2],1.0,0.0],mspec=bestmspec)

    ; Are we done?
    if count gt 0 then begin
      diffpar = 100*(par-lastpar)/abs(par)
      maxdiffpar = max(diffpar)
      if maxdiffpar lt 0.1 then flag=1
      if ngd ne nlastgd then flag=0  ; keep iterating if Ngd has changed (unless count>10)
      if count gt 10 then flag=1
    endif
    lastpar = par
    nlastgd = ngd

    ; Plotting
    if keyword_set(pl) or keyword_set(test) then begin
      !p.multi=[0,1,2]
      ;plot,x,spec,xtit='Pixels',ytit='Counts',xs=1,tit='Count='+strtrim(count,2)+' ChiSq='+strtrim(rchisq,2)
      ;oplot,x,telnorm*smspecnorm,co=250,linestyle=2
      ;oplot,x,spec/telnorm,co=150
      ;oplot,x[gd],spec[gd],ps=1
      ;legend,['Original','Telluric','Corrected'],textcolor=[255,250,150],/bottom,/left
      plot,x,specnorm/telnorm,xtit='Pixels',ytit='corrected',xs=1,tit='Count='+strtrim(count,2)+' ChiSq='+strtrim(rchisq,2),yr=[0.7,1.3],ystyle=1

      plot,x,specnorm,xtit='Pixels',ytit='Counts',xs=1,tit='Count='+strtrim(count,2)+' ChiSq='+strtrim(rchisq,2)+' starind='+strtrim(starind,2),yr=[0,1.2]
      oplot,x,telnorm,co=250;,linestyle=2
      oplot,x,specnorm/telnorm,co=150
      oplot,x[gd],specnorm[gd],ps=1,symsize=0.25
      print,par
      print,ngd
      if flag eq 1 then stop
    endif

    count++

  ENDWHILE

  ; Final parameters
  smspec = MEDFILT1D(spec/telnorm,100,/edge)

  ; Plug the data into the structure
  outstr.specfitopt = 1
  outstr.fiber = starind
  outstr.par = par[0:2]
  outstr.cont_coef = par[3:*]
  outstr.status = status
  outstr.rchisq = rchisq

End  ; chi-squared fitting method



;-------------------------------------------------------------------------
; OPTION 2: RMS Fitting to telluric-dominated pixels
;-------------------------------------------------------------------------
2: begin

  ; Continuum fit
  ;------------------
  smspec1 = MEDFILT1D(spec,100,/edge)

  spec2 = spec / smspec1
  errspec2 = errspec / smspec1  > 1.0/median(smspec1)

  ;cont_coef_pix = ROBUST_POLY_FIT(x,smspec,6)  ; with pixels

  ; Apply the correction method
  ;-----------------------------

  telnorm0 = fit_telluric(x,[1.0,1.0,1.0,1.0,0.0],mspec=mspec)
  ;gd = where(telnorm0 lt 0.99 and spec2 lt 1.5 and spec2 gt 0.0 and $
  ;gd = where(spec2 lt 1.5 and spec2 gt 0.0 and $
  gd = where(spec2 lt 1.5 and spec2 gt 0.5 and $
;             abs(spec/smspec1/telnorm0-1) lt 5*errspec2 and $
             skyspec_lines lt badskythresh*smspec1 and (mask and badmask())  eq 0,ngd)
  telmask0 = mask*0
  if ngd gt 0 then telmask0[gd] = 1

  ; Normalize the telluric spectra the same way
  ;mspec_norm = mspec*0
  ;for j=0,2 do begin
  ;  tspec = mspec[*,j]
  ;  tspec += randomn(seed,3*npix)*errspec2
  ;  smtspec = MEDFILT1D(tspec,100,/edge)
  ;  mspec_norm[*,j] = mspec[*,j] / smtspec
  ;endfor
  ;  THIS IS WRONG.  NEED TO DO IT ITERATIVELY!!!
  ;  NOT BEFORE THE TELLURIC SPECTRA ARE SCALED, AFTERWARDS!!

  ; Iterate until we converge
  initpar = [1.0,1.0,1.0, 1.0, 0.0, 0.0, 0.0]
  flag = 0
  count = 0
  telnorm1 = telnorm0
  lastpar = initpar
  smbin = 200  ; 100
  WHILE (flag ne 1) do begin

    ; Get the continuum and make the normalized spectrum
    smspec2 = MEDFILT1D(spec/telnorm1,smbin,/edge) > 1
    ; remove bad pixels
    ;smspec2_sig = mad(spec/telnorm1 - smspec2)
    ;bdpix = where(abs(spec/telnorm1-smspec2) gt 3*smspec2_sig,nbdpix)
    ;temp = spec/telnorm1
    ;if nbdpix gt 0 then temp[bdpix]=smspec2[bdpix]
    ;smspec2 = MEDFILT1D(temp,smbin,/edge) > 1

    ;cont_coef = ap_robust_poly_fit(x,spec/telnorm1,5)
    ;smspec2 = poly(x,cont_coef) > 1
    spec2 = spec / smspec2
    errspec2 = errspec / smspec2  > 1.0/median(smspec2)

    ; Loop through the three telluric spectra
    telnormstr = replicate({scale:0.0,std:0.0,best:-1},3)
    best_scale=fltarr(nscale)
    best_std=fltarr(nscale)
    For j=0,2 do begin

      ; find the pixels where tellurics dominate
      norm1 = [0.0, 0.0, 0.0, 1.0, 0.0]
      norm1[j] = 1.0
      telspec0 = fit_telluric(x,norm1,mspec=mspec)
      telmask = long(telspec0 lt 0.99)
      ngrow = 5 ;12 ; 3 ;5  ; 5-12 same for 3500
      telmask = convol(telmask,indgen(2*ngrow+1)+1,/center)
      telmask = telmask/(telmask>1)
      telpix = where(telmask eq 1 and telmask0 eq 1,ntelpix)
      ;telpix = where(telspec0 lt 0.92,ntelpix)
      if ntelpix lt 10 then begin
        print,'No good pixels for fitting'
        goto,BOMB1
      endif

      ; Measure the scale!
      ;fa = {x:x,mspec:mspec_norm,telpix:telpix}
      for iscale=0,nscale-1 do begin
        if n_elements(bestmod) eq 0 then jscale=iscale else jscale=bestmod[j] 
        tmspec = mspec[*,*,jscale]*0
        tmspec[*,0] = mspec[*,j,jscale]
        fa = {x:x,mspec:tmspec,telpix:telpix}
        best_scale[iscale] = CHECK_TELLURIC_CORRECTION_FUNC(spec2,std,_extra=fa)
        best_std[iscale]=std
      endfor
      best=min([best_std],imin)
      telnormstr[j].scale = best_scale[imin]
      telnormstr[j].std = best_std[imin]
      if n_elements(bestmod) eq 0 then telnormstr[j].best=imin else $
        telnormstr[j].best = jscale

    Endfor ; species loop

    par1 = [telnormstr.scale,1.0,0.0]
    tmspec=reform(mspec[*,*,0])*0
    for j=0,2 do tmspec[*,j]=mspec[*,j,telnormstr[j].best]
    telnorm1 = fit_telluric(x,par1,mspec=tmspec)

    ;print,telnormstr.scale
    ;
    ;plot,spec2,yr=[0,2]
    ;oplot,telnorm1,co=250

    ; IS IT POSSIBLE TO GET AN ERROR ESTIMATE???

    ; Are we done
    if count gt 0 then begin
      diffpar = 100*(par1-lastpar)/abs(par1)
      maxdiffpar = max(diffpar)
      if maxdiffpar lt 0.1 then flag=1
      ;if ngd ne nlastgd then flag=0  ; keep iterating if Ngd has changed (unless count>10)
      if count gt 10 then flag=1
    endif
    lastpar = par1
    ;nlastgd = ngd

    count++

    ;stop

  ENDWHILE


  ; Final parameters
  smspec = MEDFILT1D(spec/telnorm1,100,/edge)
  telnorm = telnorm1
  par = par1


  ; Measure shifts in X-dimension
  ;xmask = long(telnorm lt 0.99)
  ;ngrow = 5
  ;xmask = convol(xmask,indgen(2*ngrow+1)+1,/center)
  ;xmask = xmask/(xmask>1)
  ;xmask = long(xmask eq 1 and skyspec_lines lt 2.0*smspec)
  ;XCORLB,spec/smspec,telnorm,5,xshft,mask=xmask
  ;print,'Xshift = ',stringize(xshft,ndec=3)
  ;add_tag,outstr,'XSHIFT',xshft,outstr

  ; Plug the data into the structure
  outstr.specfitopt = 2
  outstr.fiber = starind
  outstr.par = par[0:2]
  outstr.cont_coef = par[3:*]
  outstr.status = 1
  outstr.bestmod = telnormstr.best
  ; Getting final chisq value
  ;  spec2/errspec2, normalized spectrum and errors
  ;  telnorm, final telluric spectrum
  telnorm0 = fit_telluric(x,[1.0,1.0,1.0,1.0,0.0],mspec=mspec)
  gdpix = where(spec2 lt 1.5 and spec2 gt 0.0 and $
                telnorm0 lt 0.99 and $
                skyspec_lines lt badskythresh*smspec1 and (mask and badmask()) eq 0,ngdpix)
  rchisq = total( (spec2[gdpix]-telnorm[gdpix])^2 / errspec2[gdpix]^2) / ngdpix
  outstr.rchisq2 = rchisq

  ;stop

End ; RMS-minimization method

else: begin
   print,specfitopt+' Not supported'
   return
end

ENDCASE

; Put stuff in the output structure
outstr.x = x
outstr.nspec = spec/smspec
outstr.cont = smspec
outstr.telluric = telnorm

print,'Normalization=',strtrim(par[0:2],2)

; Plotting
if keyword_set(pl) then begin
  set_plot,'X'
  ifiber = 300-starind

  co = 255
  ;if keyword_set(save) then begin
  ;  ;psfile1 = plots_dir+'aptelluric_'+expname+'_hotstarfit_fiber'+strtrim(ifiber,2)
  ;  psfile1 = plots_dir+'aptelluric_'+expname+'_starfit_fiber'+strtrim(ifiber,2)
  ;  PUSH,psfiles,psfile1
  ;  ps_open,psfile1,thick=4,/color,/encap
  ;  device,/inches,xsize=14,ysize=7
  ;  loadct,39,/silent
  ;  co = 0
  ;endif

  medspec = median(spec)
  yr = [(medspec-0.5*abs(medspec))<0, medspec*1.5]
  xr = minmax(x)

  plot,x,spec,xtit='Pixels',ytit='Counts',xs=1,ys=1,xr=xr,yr=yr,tit='Fiber '+strtrim(ifiber,2)
  ;oplot,x[gd],yfit1,co=250,linestyle=2
  ;oplot,x,yfit1*smspec,co=250,linestyle=2
  oplot,x,telnorm*smspec,co=250,linestyle=2
  oplot,x,spec/telnorm,co=150
  legend,['Original','Telluric','Corrected'],textcolor=[co,250,150],/bottom,/left
  xyouts,mean(xr),yr[1]-0.05*range(yr),'N ormalization='+strjoin(stringize(par[0:2],ndec=4),' '),align=0.5,charsize=1.5,charthick=4

  ;if keyword_set(save) then ps_close

endif

BOMB1:

if keyword_set(stp) then stop

end
