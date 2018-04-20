pro apxcorr,wave,tempspec,obsspec,obserr,outstr,auto=auto,pl=pl,sum=sum,maxlag=maxlag,nofit=nofit,$
            pixlim=pixlim,errccf=errccf, fail=fail

;+
;
; APXCORR
;
; This program measures the cross-correlation shift between
; a template spectrum (can be synthetic or observed) and
; an observed spectrum (or multiple spectra) on the same
; logarithmic wavelength scale.
;
; INPUTS:
;  wave     wavelength array
;  tempspec template spectrum: normalized and on log-lambda scale
;  obsspec  observed spectra: normalized and sampled on tempspec scale
;  obserr   observed error; normalized and sampled on tempspec scale
;  /sum     set to sum the multiple cross-correlation functions; used
;           if we have three chips of data that are stacked
;  /pl      Plot the cross-correlation function and best-fit
;            values to the peak.
;  /fail   Return default structure
; 
; OUTPUTS:
;  outstr  The output structure of the final derived RVs and
;            errors.
;  =auto   The auto-correlation function of the template
;  
;
; USAGE:
;  IDL>apxcorr,wave,tempspec,spec,err,outstr
;
; By D.Nidever  June 2012
; By J.Holtzman  October 2012
;-

pl=0
cspeed = 2.99792458d5  ; speed of light in km/s

apgundef,outstr

nwave = n_elements(wave)
ntempspec = n_elements(tempspec)
nspec = n_elements(obsspec)
nerr = n_elements(obserr)

; Not enough inputs
if nwave eq 0 or ntempspec eq 0 or nspec eq 0 or nerr eq 0 then begin
  print,'Syntax - apspecshift,wave,tempspec,spec,err,outstr,auto=auto,pl=pl'
  return
endif

sz = size(obsspec)
if sz[0] eq 2 then nspec = sz[2] else nspec = 1
nsum=1
if keyword_set(sum) then begin
  nsum=nspec
  nspec=1
endif
if keyword_set(pixlim) then nsum = n_elements(pixlim[0,*])

; Set up the cross-correlation parameters
;  this only gives +/-450 km/s with 2048 pixels, maybe use larger range
nlag = 401 ;321 ;161
if n_elements(maxlag) gt 0 then nlag=2*round(abs(maxlag[0]))+1
if nlag mod 2 eq 0 then nlag++  ; make sure nlag is even
dlag = 1
minlag = -nlag/2
lag = findgen(nlag)*dlag+minlag

outstr = replicate({xshift0:0.0,ccp0:0.0,xshift:0.0,xshifterr:0.0,$
                    xshift_interp:0.0,ccf:fltarr(nlag),$
                    ccferr:fltarr(nlag),cclag:lag,ccpeak:0.0,$
                    ccpfwhm:0.0,ccp_pars:fltarr(5),$
                    ccp_perror:fltarr(5),ccp_polycoef:fltarr(4),vrel:0.0,$
                    vrelerr:0.0,w0:0.0d0,dw:0.0d0,chisq:0.0},nspec)
outstr.xshift = !values.f_nan
outstr.xshifterr = !values.f_nan
outstr.vrel = !values.f_nan
outstr.vrelerr = !values.f_nan
outstr.chisq = !values.f_nan

IF keyword_set(fail) THEN return

; Spectrum cross-correlation loop
;---------------------------------
for i=0,nspec-1 do begin
  pl = 0
  if keyword_set(pl) then begin
    !p.multi=[0,2,4]
    ;smcolor
  endif

  nsumgood = 0
  for ichip=0,nsum-1 do begin
    ; Interpolate the visit spectrum onto the synth wavelength grid
    ;---------------------------------------------------------------
    wobs = reform(wave)
    nw=n_elements(wobs)
    if keyword_set(pixlim) then begin
      IF pixlim[0,ichip] LT pixlim[1,ichip] THEN BEGIN
       spec = reform(obsspec[pixlim[0,ichip]:pixlim[1,ichip],i])
       err = reform(obserr[pixlim[0,ichip]:pixlim[1,ichip],i])
       template = tempspec[pixlim[0,ichip]:pixlim[1,ichip]]
      ENDIF ELSE BEGIN
       spec = 0.
       err = 0.
       template = 0.
      ENDELSE
      nw = n_elements(spec)
    endif else begin
      spec = reform(obsspec[*,i+ichip])
      err = reform(obserr[*,i+ichip])
      template=tempspec
    endelse

    ; find first and last non-zero pixels in observed spec
    gd=where(spec ne 0.0,ngd)
    if ngd eq 0 then begin
      ccf = float(lag)*0.0
      if keyword_set(errccf) or not keyword_set(nofit) then ccferr=ccf
      lo = 0
      hi = nw-1
      goto,BOMB1
    endif
    goodlo=gd[0] > 0
    goodhi=gd[ngd-1] < (nw-1)

    ; fix any pixels within good range
    ;fix=where(spec[goodlo:goodhi] le 0.01,nfix)
    ;if nfix gt 0 then spec[fix+goodlo]=0.0
    ;fix=where(template[goodlo:goodhi] le 0.01,nfix)
    ;if nfix gt 0 then template[fix+goodlo]=0.0

    ; mask bad pixels, set to NAN
    fix=where(spec le 0.01,nfix)
    if nfix gt 0 then spec[fix]=!values.f_nan
    fix=where(template le 0.01,nfix)
    if nfix gt 0 then template[fix]=!values.f_nan

    ; How masking pixels affects the cross-correlation
    ; -mean
    ; -total(x*y)
    ; -1/rmsx*rmsy
    ; apc_correlate deals with all this appropriately if bad pixels are NAN

    ; set cross-corrlation window to be good range + nlag
    lo=(gd[0]-nlag) > 0
    hi=(gd[ngd-1]+nlag) < (nw-1)

    indobs = where(finite(spec) eq 1,nindobs)  ; only finite values, in case any NAN
    indtemp = where(finite(template) eq 1,nindtemp)
    if (hi-lo+1) gt nlag and nindobs gt 0 and nindtemp gt 0 then begin

      ; Cross-Correlation
      ;------------------
      ccf = APC_CORRELATE(template[lo:hi],spec[lo:hi],lag)
      ;ccf = C_CORRELATE(template[lo:hi],spec[lo:hi],lag)

      if keyword_set(pl) then begin
        plot,spec[lo:hi],yrange=[-0.5,0.2]+1
        oplot,template[lo:hi],color=250 ;2
        plot,ccf
      endif

      ; Calculate the CCF uncertainties using propagation of errors
      ;if keyword_set(errccf) then begin
      if keyword_set(errccf) or not keyword_set(nofit) then begin
        ; Make median filtered error array
        ;   high error values give crazy values in ccferr
        obserr1 = err[lo:hi]
        bderr = where(obserr1 gt 1 or obserr1 le 0.0,nbderr,comp=gderr,ncomp=ngderr)
        if nbderr gt 0 and ngderr gt 1 then obserr1[bderr]=median([obserr1[gderr]])
        obserr1 = medfilt1d(obserr1,51,/edge)
        ccferr = C_CORRELATE_ERROR(template[lo:hi],spec[lo:hi],obserr1,lag)
      endif

      nsumgood++
    endif else begin ; no good pixels
      ccf = float(lag)*0.0
      if keyword_set(errccf) or not keyword_set(nofit) then ccferr=ccf
    endelse

    BOMB1:

    if ichip eq 0 then begin
      xcorr=ccf 
      if keyword_set(ccferr) then xcorrerr=ccferr^2
      tout=template[lo:hi] ;+1
      sout=spec[lo:hi]     ;+1
      errout=err[lo:hi]
    endif else begin
      xcorr+=ccf
      if keyword_set(ccferr) then xcorrerr+=ccferr^2  ; add in quadrature
      tout=[tout,template[lo:hi]] ;+1
      sout=[sout,spec[lo:hi]]     ; +1
      errout=[errout,err[lo:hi]]
    endelse
  endfor ; chip loop

  ; No good cross-correlation
  if nsumgood eq 0 then begin
    error = 'No good cross-correlation'
    ;print,error
    goto,BOMB
  endif
  xcorr /= (nsumgood>1) ; want to average XCORR arrays from multiple chips
  if keyword_set(xcorrerr) then xcorrerr=sqrt(xcorrerr)/(nsumgood>1)

  ; Remove the median
  xcorr -=median(xcorr)
  if keyword_set(pl) then plot,xcorr

  ; Best shift
  best_shiftind0 = first_el(maxloc(xcorr))
  best_xshift0 = lag[best_shiftind0]
  temp = shift( tout, best_xshift0)

  ; Find Chisq for each synthetic spectrum
  gdpix = where( finite(sout) eq 1 and finite(temp) eq 1 and sout gt 0.0 and errout gt 0.0 and errout lt 1e5,ngdpix)
  if ngdpix eq 0 then begin
    error = 'Bad spectrum'
    print,error
    goto,BOMB
  endif

  ;chisq = sqrt( total( (spec[gdpix]-temp[gdpix])^2/err[gdpix]^2 )/ngdpix )
  chisq = sqrt( total( (sout[gdpix]-temp[gdpix])^2/errout[gdpix]^2 )/ngdpix )
  if keyword_set(pl) then begin
    plot,sout[gdpix],yrange=[-0.5,2]
    oplot,temp[gdpix],color=2
  endif

  outstr[i].chisq = chisq
  outstr[i].ccf = xcorr
  if keyword_set(errccf) then outstr[i].ccferr = xcorrerr

  if keyword_set(nofit) then return

  ; Remove smooth background at large scales
  cont = GSMOOTH(xcorr,100)
  xcorr_diff = xcorr-cont

  ;func = 'gfunc'
  func = 'gaussbin'

  ; Fit Xcorr peak with a Gaussian plus a line
  ;---------------------------------------------
  ; Some CCF peaks are SOOO wide that they span the whole width
  ; do the first one without background subtraction
  parinfo0 = REPLICATE({limited:[0,0],limits:[0.0,0.0],fixed:0},4)
  parinfo0[0].limited[0]=1 & parinfo0[0].limits[0]=1e-3
  parinfo0[1].limited=1 & parinfo0[1].limits=[min(lag),max(lag)]
  parinfo0[2].limited[0]=1 & parinfo0[2].limits[0]=0.1
  estimates0 = [xcorr_diff[best_shiftind0],best_xshift0,4.0,0.0]
  pars0 = MPFITFUN(func,lag,xcorr_diff,xcorrerr,estimates0,yfit=yfit0,$
                   bestnorm=chisq0,dof=dof0,parinfo=parinfo0,perror=perror0,yerror=yerror0,status=status0,/quiet)
  if status0 eq 0 then goto,BOMB

  ; Find peak in the background subtracted XCOR
  best_shiftind = first_el(maxloc(xcorr_diff))
  best_xshift = lag[best_shiftind]
  
  ; fit the width
  ;  keep height, center and constant constrained
  estimates = pars0
  estimates[1] = best_xshift
  parinfo = REPLICATE({limited:[0,0],limits:[0.0,0.0],fixed:0},n_elements(estimates))
  parinfo[0].limited=1 & parinfo[0].limits=[0.5,1.5]*estimates[0]
  parinfo[1].limited = 1
  parinfo[1].limits = [-4,4]+best_xshift
  parinfo[2].limited=1 & parinfo[2].limits=[0.3,1.5]*pars0[2]
  parinfo[3].limited = 1
  parinfo[3].limits = [(min(xcorr_diff)<0) < (pars0[3]-0.1),max(xcorr_diff)*0.5 > (pars0[3]+0.1)]
  lo1 = floor( best_shiftind-(pars0[2]*2>5) ) > 0
  hi1 = ceil( best_shiftind+(pars0[2]*2>5) ) < (n_elements(lag)-1)
  pars1 = MPFITFUN(func,lag[lo1:hi1],xcorr_diff[lo1:hi1],xcorrerr[lo1:hi1],estimates,yfit=yfit1,$
                   bestnorm=chisq1,dof=dof1,parinfo=parinfo,perror=perror1,yerror=yerror1,status=status1,/quiet)
  if status1 eq 0 then goto,BOMB


  ; refit and let constant vary more, keep width constrained
  lo2 = floor( best_shiftind-(pars1[2]*2>5) ) > 0
  hi2 = ceil( best_shiftind+(pars1[2]*2>5) ) < (n_elements(lag)-1)
  est2 = pars1
  est2[3] = median(xcorr_diff[lo1:hi1]-yfit1)   + pars1[3]
  parinfo2 = parinfo
  parinfo2[3].fixed = 0  ; let it float now
  parinfo2[0].limited=1 & parinfo2[0].limits=[0.5,1.5]*pars1[0]
  parinfo2[1].limits = [min(lag) > (best_xshift-(pars1[2]>1)) < (pars1[1]-1),$
                        max(lag) < (best_xshift+(pars1[2]>1)) > (pars1[1]+1)]
  parinfo2[2].limited=1 & parinfo2[2].limits=[0.3,1.5]*pars1[2]
  parinfo2[3].limits = [(min(xcorr_diff)<0) < (est2[3]-0.1),max(xcorr_diff)*0.5 > (est2[3]+0.1)]
  pars2 = MPFITFUN(func,lag[lo2:hi2],xcorr_diff[lo2:hi2],xcorrerr[lo2:hi2],est2,yfit=yfit2,$
                   bestnorm=chisq2,dof=dof2,parinfo=parinfo2,perror=perror2,yerror=yerror2,status=status2,/quiet)
  if status2 lt 1 then goto,BOMB


  ; refit with even narrower range
  lo3 = floor( best_shiftind-(pars2[2]>5) ) > 0
  hi3 = ceil( best_shiftind+(pars2[2]>5) ) < (n_elements(lag)-1)
  parinfo2[0].limited=1 & parinfo2[0].limits=[0.5,1.5]*pars2[0]
  parinfo2[1].limits = [min(lag) > (best_xshift-(pars2[2]>1)) < (pars2[1]-1),$
                        max(lag) < (best_xshift+(pars2[2]>1)) > (pars2[1]+1)]
  parinfo2[2].limited=1 & parinfo2[2].limits=[0.3,1.5]*pars2[2]
  parinfo2[3].limits = [(min(xcorr_diff)<0) < (pars2[3]-0.1),max(xcorr_diff)*0.5 > (pars2[3]+0.1)]
  pars3 = MPFITFUN(func,lag[lo3:hi3],xcorr_diff[lo3:hi3],xcorrerr[lo3:hi3],pars2,yfit=yfit3,$
                   bestnorm=chisq3,dof=dof3,parinfo=parinfo2,perror=perror3,yerror=yerror3,status=status3,/quiet)
  if status3 lt 1 then goto,BOMB
  ; this seems to fix high shift/sigma errors
  if perror3[0] gt 10 or perror3[1] gt 10 then begin
    dparinfo = parinfo2
    dparinfo[1].limits = [-10,10]+pars2[1]
    dparinfo[2].limits = [0.01,2*pars2[2]]
    dparinfo[[0,3]].fixed = 1
    dpars2 = MPFITFUN(func,lag[lo3:hi3],xcorr_diff[lo3:hi3],xcorrerr[lo3:hi3],pars2,parinfo=dparinfo,$
                     perror=perror3,/quiet)
  endif

  ; Final parameters
  pars = pars3
  perror = perror3
  xshift = pars[1]
  xshifterr = perror[1]
  ccpfwhm_pix = pars[2]*2.35482  ; ccp fwhm in pixels
  ; v = (10^(delta log(wave))-1)*c
  dwlog = median(slope(alog10(wave)))
  ccpfwhm = ( 10.d0^(ccpfwhm_pix*dwlog)-1.0d0 )*cspeed  ; in km/s

  ; Convert pixel shift to velocity
  ;---------------------------------
  ; delta log(wave) = log(v/c+1)
  ; v = (10^(delta log(wave))-1)*c
  dwlog = median(slope(alog10(wave)))
  vrel = ( 10.d0^(xshift*dwlog)-1.0d0 )*cspeed
  ; Vrel uncertainty
  dvreldshift = alog(10.0)*(10.d0^(xshift*dwlog))*dwlog*cspeed  ; derivative wrt shift
  vrelerr = dvreldshift * xshifterr

  ; Make XCORR structure and add to STR
  ;------------------------------------
  outstr[i].xshift0 = best_xshift
  outstr[i].ccp0 = max(xcorr)
  outstr[i].xshift = xshift
  outstr[i].xshifterr = xshifterr
  ;outstr[i].xshift_interp = xshift_interp
  outstr[i].ccpeak = pars[0] 
  outstr[i].ccpfwhm = ccpfwhm  ; in km/s
  outstr[i].ccp_pars = pars
  outstr[i].ccp_perror = perror
  ;outstr[i].ccp_polycoef = polycoef
  outstr[i].vrel = vrel
  outstr[i].vrelerr = vrelerr
  outstr[i].w0 = min(wave)
  outstr[i].dw = dwlog
  ;outstr[i].chisq = chisq

  ; Printing
  ;----------
  ;print,' Best Xshift = ',stringize(best_xshift,ndec=2),' pixels'
  ;print,' Best Gaussian fit Xshift = ',stringize(xshift,ndec=2),' +/- ',stringize(xshifterr,ndec=2),' pixels'
  ;print,' Vrel = ',stringize(vrel,ndec=2),' +/- ',stringize(vrelerr,ndec=2),' km/s'

  ; Plotting
  if keyword_set(pl) then begin
    plot,lag,xcorr_diff,ps=-1,xs=1 ;xr=xr
    oplot,lag[lo1:hi1],xcorr_diff[lo1:hi1],color=3,ps=-1
    oplot,lag[lo1:hi1],yfit1,color=4
    stop
    plot,lag,xcorr_diff,ps=-1,xs=1 ;xr=xr
    oplot,lag[lo2:hi2],xcorr_diff[lo2:hi2],color=3,ps=-1
    oplot,lag[lo2:hi2],yfit2,color=4
    stop
    plot,lag,xcorr_diff,ps=-1,xs=1 ;xr=xr
    oplot,lag[lo3:hi3],xcorr_diff[lo3:hi3],color=3,ps=-1
    oplot,lag[lo3:hi3],yfit3,color=4
    stop

    xr = [-50,50]
    plot,lag,xcorr_diff,ps=-1,xs=1,xr=xr
    oplot,[0,0]+xshift,[-10,10],co=250
    gd=where(finite(yfit3) eq 1,ngd)
    if ngd gt 1 then oplot,lag[lo3:hi3],yfit3,co=150
    !p.multi=0
    print,pars
    print,perror
    wait,0.5
    stop
  endif

  BOMB:

  ;wait,0.7
  ;stop

Endfor  ; spectrum loop

; Auto-correlation
;auto = C_CORRELATE(template,template,lag)

;stop

end
