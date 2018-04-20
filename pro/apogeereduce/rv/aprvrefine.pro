pro aprvrefine,allstr,allvisits,starstr,sinc=sinc,log=log,holtz=holtz,plotdir=plotdir,bccomb=bccomb,synth=synth, trimgrid = trimgrid, allsmooth = allsmooth, nosmooth = nosmooth, chisq=chisq, fixmasked = fixmasked

; KEYWORDS:
;  bccomb - Set to: 1 to do initial combination using baricentric correction, 2 to combine using visit-level RVs. Uses this combined spectrum in the first iteration rather than highest S/N visit.
;  synth - Allows intermediate and final visit combination using synth-template velocities IF they produce less scatter than the relative observed-template RVs. 
;  trimgrid - Set to restrict synthetic grids that can be used based on the J-K color and Washington photometry of the star
;  allsmooth - Set to have smoothing occur at every iteration (with lower S/N threshold (default 10) after 2 iterations         
;  nosmoth - Set to have smoothing never occur despite S/N of the spectrum
; fixmasked - Interpolate over bad pixels rather than just setting them to nan or 0. 
; chisq - Set to use aprv_template_nofile rather than apxcorr to derive RVs (aprv_template_nofile includes option to use direct-pixel chisq method as well as xcorr). Possibly useful for low S/N spectra. - NOT FULLY IMPLEMENTED YET(partially for observed spectra not for synthetic) - Needed auxillary functions not committed yet. - DO NOT SET                                                                             

; Original written by D.Nidever. Still available in v3_31 (DR12 tagged version)
; Heavy modifications by N.Troup  6/2016 - Added functionality to iterate using synthetic spectrum 

pl = 0           ;Debug mode off/on
;aprvprep parameters:  
IF n_elements(fixmasked) EQ 0 then fixmasked = 1 ;If not set, default to DR13 value (See below) 
          ; DR13 default = 1(interpolating over the bad pixels in aprvprep, 
          ; 0= fixmasked keyword not set in aprvprep. - This is the best with synth keyword on - reccomended for DR14.
          ; 2 = use setnan override keyword in apfixmasked,- NOT COMMITTED!

nvisits=n_elements(allvisits)

IF nvisits GT 1 THEN BEGIN
; VREL started with barycentric correction
  ;   Subtracting the mean here is VITAL otherwise the final apStar
  ;   RV will be WRONG
 IF keyword_set(bccomb) && (abs(bccomb) EQ 2) THEN allvisits.synthvrel = allvisits.vrel ELSE BEGIN
   ; start with barycentric correction
  allvisits.vrel = -allvisits.bc+mean(allvisits.bc)
  allvisits.synthvrel = -allvisits.bc+mean(allvisits.bc)
 ENDELSE
ENDIF ELSE BEGIN
  ref = 0 ;Need to set ref to 0 for single-visit spectra
  allvisits.vrel = 0.0
  allvisits.vrelerr = 0.0
  IF keyword_set(synth) THEN BEGIN
    allvisits.synthvrel = 0.0
    allvisits.synthvrelerr = 0.0
  ENDIF
ENDELSE
 

 apvisitcomb,allstr,allvisits,starstr,/nolsffit,sinc=sinc,log=log,/quick
 ;GET SYNTHETIC GRID TO MAKE TEMPLATES 
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 synthfile = 'apg_rvsynthgrid_v2.fits'
 apgetgrid,synthfile,grid=grid,wave=starstr.wave
 IF keyword_set(trimgrid) THEN grid = aptrimgrid(grid, allvisits[0])
 ;DOES THIS NEED TO BE INSIDE THE WHILE? Don't think the pixel limits update. 
 if tag_exist(grid,'PIXLIM') eq 0 then add_tag,grid,'PIXLIM',starstr.pixlim_overlap,grid else begin
    gridpixlim = reverse(grid.pixlim,2)
    ind = where(starstr.pixlim LT 0)
    pixlim = gridpixlim*0
    pixlim[0,*] = gridpixlim[0,*] > starstr.pixlim[0,*]  ; want pixel range where both are non-zero
    pixlim[1,*] = gridpixlim[1,*] < starstr.pixlim[1,*]
    ; Check for negative pixlim 
    IF total(ind) GE 0 THEN pixlim[ind] = gridpixlim[ind]
    grid.pixlim = pixlim	
 endelse


; IF nvisits > 1 then perform iterations. Otherwise skip to the end and produce final velocities
IF nvisits GT 1 THEN BEGIN

 ; SET UP INITIAL PARAMETERS FOR ITERATIONS
 allvisits.vhelio=0
 allvisits.synthvhelio = 0
 niter=0 & maxiter=10 ;1 = TEMPORARY FOR DEBUGGING
 last_vrel=allvisits.vrel
 last_synthvrel = allvisits.synthvrel
 last_vscatter  = !Values.F_nan
 last_synthvscatter  = !Values.F_nan

 ; use best S/N spectrum as reference for first iteration unless bccomb is set
 bestsnr=max(allvisits.snr,ibestsnr)
 endflag = 0
 if keyword_set(bccomb) then begin 
    ref=0 
    IF bccomb LT 0 THEN endflag = 1
 endif else ref=2+ibestsnr
 print,'S/N:',allvisits.snr,format='(a4,10f8.2)'
 print,'ITER      Relative RVs      RV errors'
 ;while niter lt maxiter do begin

  while (endflag ne 1) do begin

    if niter gt 0 then reject=1 else reject=0  ; reject outliers once we have comb spec

    ; Smooth low-S/N spectra for early iterations
    IF niter LT 2 THEN snrlim = 20 ELSE snrlim = 10
    IF keyword_set(nosmooth) THEN snrlim = 0
    if niter lt 2 and (min(allvisits.snr) lt snrlim) then begin
      lowsnr = where(allvisits.snr lt snrlim,nlowsnr)
      smlen = 3.0 ;1.5
      for k=0,nlowsnr-1 do starstr.spec[*,2+lowsnr[k]]=gsmooth(starstr.spec[*,2+lowsnr[k]],smlen)
    endif else begin
     IF keyword_set(allsmooth) AND (min(allvisits.snr) lt snrlim) then begin
      lowsnr = where(allvisits.snr lt snrlim,nlowsnr)
      smlen = 3.0 ;1.5
      for k=0,nlowsnr-1 do starstr.spec[*,2+lowsnr[k]]=gsmooth(starstr.spec[*,2+lowsnr[k]],smlen)
     END
    endelse
    ;Normalize and mask observed spectra
    aprvprep,starstr.wave,starstr,specout,errout,/gaps, fixmasked=fixmasked,reject=reject,pl=pl
    add_tag, starstr, 'nspec', specout, starstr
    add_tag, starstr, 'nerr', errout, starstr    

    ;DERIVE RVs using observed spectrum as template
    apgundef,maxlag
    if niter gt 2 then if max(vout.ccp_pars[2]*2.35) lt 75 then maxlag=75
    IF keyword_set(chisq) THEN BEGIN
      ;print, 'Using chisq mode - NOT IMPLEMENTED YET'
      stop
      ;vout = aprv_template_nofile(allvisits,starstr, ref = ref, maxlag=maxlag, dir_plots = plotdir, pl=pl, /apstar)   
    ENDIF ELSE BEGIN
     apxcorr,starstr.wave,specout[*,ref],specout[*,2:nvisits+1],$
            errout[*,2:nvisits+1],vout,maxlag=maxlag,pixlim=starstr.pixlim_overlap
    ENDELSE
    print,niter,'  RESIDRV',vout.vrel,format='(i4,A-10,20f9.3)'
    print,'','  ERR',vout.vrelerr,format='(A4,A-10,20f9.3)'

    ;Derive RVs using synthetic template if synth keyword is set
    IF keyword_set(synth) THEN BEGIN
     aprvprep,grid.wave,starstr,synth_specout,synth_errout,/gaps,fixmasked=fixmasked, reject=reject,pl=pl
     synth_specout = specout
     synth_errout	= errout

     ; Choose best grid from current ref spectrum
     apxcorr_newgrid,grid,synth_specout[*,ref],synth_errout[*,ref],synthvout,bestgrid

     IF keyword_set(chisq) THEN BEGIN
       ; TO BE IMPLEMENTED!
      print, 'This feature not available yet for synthetic templates. Using apxcorr' 
      apxcorr,grid.wave,grid.ndata[bestgrid,*],synth_specout,synth_errout,synthvout,sum=sum,/errccf,pixlim=grid.pixlim ; TEMPORARY until this is implemented.
      ;synthvout = aprv_template_nofile(allvisits,starstr, maxlag=maxlag, dir_plots = plotdir, pl=pl, /apstar)
     ENDIF ELSE BEGIN
       apxcorr,grid.wave,grid.ndata[bestgrid,*],synth_specout,synth_errout,synthvout,sum=sum,/errccf,pixlim=grid.pixlim
     ENDELSE
     refvout = synthvout[ref]
     synthvout = synthvout[2:nvisits+1]
     synthvout.vrel = synthvout.vrel - refvout.vrel ;Put on same scale as vrel with ref spectrum at 0 velocity
    ENDIF ELSE usesynth = 0 ; DEFAULT: Don't combine with synthvrel
   
    ;plot CCFs
    ;pl = 1
    if keyword_set(pl) then begin
      ccf = vout.ccf
      mx = max(ccf,dim=1)
      nlag = n_elements(ccf[*,0])
      ccf /= replicate(1,nlag)#mx
      plot,vout.cclag,ccf[*,0],yr=[-0.1,1.1],xs=1
      for j=0,nvisits-1 do oplot,vout.cclag,ccf[*,j],co=50+j*15
      ;stop
    endif

    if keyword_set(holtz) then begin
      for ivisit=0,n_elements(allvisits)-1 do begin
        allvisits[ivisit].vrel+=vout[ivisit].vrel
        allvisits[ivisit].vhelio+=vout[ivisit].vrel
      endfor
      IF ~keyword_set(synth) THEN BEGIN
       for ivisit=0,n_elements(allvisits)-1 do begin
        allvisits[ivisit].synthvrel+=synthvout[ivisit].vrel
        allvisits[ivisit].synthvhelio+=synthvout[ivisit].vrel
       endfor	
      ENDIF
      apvisitcomb,allstr,allvisits,starstr,/nolsffit,sinc=sinc,log=log, synth = usesynth
      if max(abs(vout.vrel)) lt 0.05 then niter=maxiter else niter+=1
    endif else begin 

     ; Remove mean from resid_vrel so we don't drift over many iterations 
     resid_vrel = vout.vrel
     gdrv = where(allvisits.snr gt 10 and finite(resid_vrel) eq 1,ngdrv)      ; exclude bad RVs
     if ngdrv eq 0 then gdrv=where(allvisits.snr gt 7 and finite(resid_vrel) eq 1,ngdrv)
     if ngdrv eq 0 then gdrv=where(finite(resid_vrel) eq 1,ngdrv)
     if ngdrv gt 0 then mnresid=median(resid_vrel[gdrv]) else mnresid=0.0
     vrel = allvisits.vrel + resid_vrel - mnresid
     vhelio = vrel + allvisits.bc
     vrelerr = vout.vrelerr


     ; Do the same thing with synthvrel if synth keyword is set.
     IF keyword_set(synth) THEN BEGIN
      resid_vrel = synthvout.vrel
      gdrv = where(allvisits.snr gt 10 and finite(synthvout.vrel) eq 1,ngdrv)      ; exclude bad RVs
      if ngdrv eq 0 then gdrv=where(allvisits.snr gt 7 and finite(synthvout.vrel) eq 1,ngdrv)
      if ngdrv eq 0 then gdrv=where(finite(synthvout.vrel) eq 1,ngdrv)
      if ngdrv gt 0 then mnresid=median(resid_vrel[gdrv],/even) else mnresid=0.0
      synthvrel = allvisits.synthvrel + resid_vrel - mnresid
      synthvhelio = synthvrel + allvisits.bc
      synthvrelerr  = synthvout.vrelerr
     ENDIF

	 ; Set maximum allowable |RV| value. Visits with values above this will be set to nan and excluded in future iterations.
      vmax = 1000. ;Maximum |vrel| offset allowed in a single iteration
      sigmax = 10. ;Maximum number of median absolute devations allowed before rejecting RV. -- BE CAREFUL WITH THIS; Don't want to erase real variability.
      vheliosig = mad(vhelio)
	 IF vheliosig LT 0 THEN vheliosig = 999999. ;If mad returns an error value, then make it so that RVs cannot be rejected via sigma clipping
      print,'','  VHELIO',vhelio,format='(A4,A-10,20f9.3)'
      bdrv = where(finite(vrelerr) eq 0,nbdrv,comp=gdrv,ncomp=ngdrv)
      if nbdrv gt 0 then begin 
       vrel[bdrv]=!values.f_nan
       vrelerr[bdrv]=!values.f_nan
       vhelio[bdrv]=!values.f_nan
      endif
      bdrv = where((abs(vrel) gt vmax) OR (abs(vhelio-median(vhelio,/even)) GT max([4.,sigmax*vheliosig])),nbdrv,comp=gdrv,ncomp=ngdrv)
       if nbdrv gt 0 then begin 
       vrel[bdrv]=!values.f_nan
       vrelerr[bdrv]=!values.f_nan
       vhelio[bdrv]=!values.f_nan
      endif
      print,'','  VHELIO',vhelio,format='(A4,A-10,20f9.3)'

     IF keyword_set(synth) THEN BEGIN
       synthvheliosig = mad(synthvhelio)
        IF synthvheliosig LT 0 THEN synthvheliosig = 999999. ;If mad returns an error value, then make it so that RVs cannot be rejected via sigma clipping
       print,'','  SYNTHVHELIO',synthvhelio,format='(A4,A-10,20f9.3)'
       bdrv = where(finite(synthvrelerr) eq 0,nbdrv,comp=gdrv,ncomp=ngdrv)
       if nbdrv gt 0 then begin 
        synthvrel[bdrv]=!values.f_nan
        synthvrelerr[bdrv]=!values.f_nan
        synthvhelio[bdrv]=!values.f_nan
       endif
       bdrv = where((abs(synthvrel) gt vmax) OR (abs(synthvhelio-median(synthvhelio,/even)) GT max([4.,sigmax*synthvheliosig])),nbdrv,comp=gdrv,ncomp=ngdrv)
       if nbdrv gt 0 then begin 
        synthvrel[bdrv]=!values.f_nan
        synthvrelerr[bdrv]=!values.f_nan
        synthvhelio[bdrv]=!values.f_nan
       endif
       print,'','  SYNTHVHELIO',synthvhelio,format='(A4,A-10,20f9.3)'

      ;Choose between synthvrel and vrel based on which produces the smallest scatter
      synthvscatter = stddev(synthvhelio-median(synthvhelio,/even),/nan)
      vscatter = stddev(vhelio-median(vhelio,/even),/nan)
      IF (synthvscatter LT vscatter) THEN usesynth=1 ELSE usesynth=0

      ; Need to add chosen velocity to BOTH for them to remain comparable. - Issue: This muddies the definition of vrel and synthvrel....
      IF (usesynth EQ 1) THEN BEGIN
       allvisits.vrel = synthvrel
       allvisits.vhelio = synthvhelio
       allvisits.vrelerr = synthvrelerr
       allvisits.synthvrel = synthvrel
       allvisits.synthvhelio = synthvhelio
       allvisits.synthvrelerr = synthvrelerr
      ENDIF ELSE BEGIN
       allvisits.vrel = vrel
       allvisits.vhelio = vhelio
       allvisits.vrelerr = vrelerr
       allvisits.synthvrel = vrel
       allvisits.synthvhelio = vhelio
       allvisits.synthvrelerr = vrelerr
      ENDELSE
      print, 'vscatter = ', vscatter
      print, 'synthvscatter = ', synthvscatter
      IF usesynth THEN print, 'Using Synthetic Template' ELSE print, 'Using observed template'
     ENDIF ELSE BEGIN 
      usesynth = 0
      allvisits.vrel =  vrel
      allvisits.vhelio = vhelio
      allvisits.vrelerr = vout.vrelerr
     ENDELSE ; keyword_set(synth)
     gdrv = where(finite(allvisits.vrel),ngdrv)
     
     ; Combine spectrum using velocity method with the lowest scatter
     apvisitcomb,allstr,allvisits,starstr,/nolsffit,sinc=sinc,log=log,/quick, synth = usesynth
     npixcomb0 = n_elements(starstr.mask[*,0])
     ngdcomb0 = n_elements(where(starstr.mask[*,0] EQ 0))
    	IF n_elements(starstr.mask[0,*]) GT 1 THEN BEGIN 
         npixcomb1 = n_elements(starstr.mask[*,1])
         ngdcomb1 = n_elements(where(starstr.mask[*,1] EQ 0))
     ENDIF ELSE ngdcomb1 = 0
     print, "Unmasked pixels in combined spectra", ngdcomb0, ngdcomb1

     ; check for bad RV errors, sometimes MPFIT spits out 0.0
     bderror = where(vout.vrelerr lt 1e-5,nbderror)
     if nbderror gt 0 and niter gt 0 then vout[bderror].vrelerr=last_vout[bderror].vrelerr

     IF keyword_set(synth) THEN BEGIN
      bderror = where(synthvout.vrelerr lt 1e-5,nbderror)
      if nbderror gt 0 and niter gt 0 then synthvout[bderror].vrelerr=last_synthvout[bderror].vrelerr
     ENDIF

     ; Have we converged?
     diff_vrel = allvisits.vrel-last_vrel
     diff_synthvrel = allvisits.synthvrel - last_synthvrel
     maxdiff = max(abs(diff_vrel))
     meddiff = median(abs(diff_vrel))
     if (maxdiff lt 0.1 and meddiff lt 0.1 and niter gt 3) or niter gt maxiter then endflag=1
     if ngdrv eq 0 then endflag=1  ; all RVs are bad
     last_vrel = allvisits.vrel
     last_vout = vout
     IF keyword_set(synth) THEN last_synthvout = synthvout
     niter++
    endelse ;keyword_set(holtz)
    
    IF keyword_set(pl) THEN BEGIN
     !p.multi = 0
     window,0
     plotc, allvisits.JD-2450000, vhelio, allvisits.bc, psym = 5, yrange = [-25,25]
     window,1
     plotc, allvisits.JD-2450000, synthvhelio, allvisits.bc, psym = 5, yrange = [-25,25]
     window,2
     cgplot, starstr.wave, specout[*,ref], xr = [1.51,1.58]*(10^4)
     cgplot, starstr.wave, specout[*,2], /overplot,color='red'
    ENDIF
    
    IF (ngdcomb0 LT 0.5*npixcomb0) THEN ref=1 ELSE ref = 0 ;use combined spectrum for reference in subsequent cross correlations
                                                           ;If the OR masked combined spectrum (index 0) has more the half of the spectrum masked out, then use the AND maksed spectrum (index 1)
                                                           ;NOTE the combined spectrum index 1 uses a different combination scheme (global vs. pixel-by-pixel snr weighting) than the index 0
  endwhile

  IF ~keyword_set(bccomb)||(bccomb GE 0) THEN vout0 = vout ELSE synth = 0
  print,'Final relative Velocities:'
  print,allvisits.vrel,format='(20f9.3)'
  print,allvisits.vrelerr,format='(20f9.3)'
endif; Nvisits=1

;FINALIZE RVs and put them on absolute scale.
; TO DO: Convert vhelio, etc. to final chosen velocity, and move old contents of vhelio,etc. to obsvhelio -- MAYBE DO THIS IS APSTAR?
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if not tag_exist(allvisits,'OBSVREL') then begin
 ADD_TAG, allvisits, 'obsvrel', 999999.0, allvisits
 ADD_TAG, allvisits, 'obsvrelerr', 999999.0, allvisits
 ADD_TAG, allvisits, 'obsvhelio', 999999.0, allvisits
endif
; All RVs are bad
gdrv = where(finite(allvisits.vrel) eq 1,ngdrv,comp=bdrv,ncomp=nbdrv)
if ngdrv eq 0 then begin
  print,'All RVs bad'
  allvisits.vtype = -1
  allvisits.vrel = 999999.0
  allvisits.vrelerr = 999999.0
  allvisits.vhelio = 999999.0
  allvisits.vlsr = 999999.0
  allvisits.vgsr = 999999.0
  allvisits.chisq = 999999.0
  allvisits.obsvrel = 999999.0
  allvisits.obsvrelerr = 999999.0
  allvisits.obsvhelio = 999999.0
  allvisits.synthvrel = 999999.0
  allvisits.synthvrelerr = 999999.0
  allvisits.synthvhelio = 999999.0
  allvisits.synthfile = ''
  allvisits.rv_teff = 999999.0
  allvisits.rv_logg = 999999.0
  allvisits.rv_feh = 999999.0
  allvisits.rv_alpha = 999999.0
  allvisits.rv_carb = 999999.0
  ; Setup starstr with needed output values, initially set to their error/default values.
  ADD_TAG,starstr,'RV_TEFF',999999.0,starstr
  ADD_TAG,starstr,'RV_LOGG',999999.0,starstr
  ADD_TAG,starstr,'RV_FEH',999999.0,starstr
  ADD_TAG,starstr,'RV_ALPHA',999999.0,starstr
  ADD_TAG,starstr,'RV_CARB',999999.0,starstr
  ADD_TAG,starstr,'CCPFWHM',999999.0,starstr
  ADD_TAG,starstr,'AUTOFWHM',999999.0,starstr
  ADD_TAG,starstr,'XSHIFT',fltarr(nvisits)+999999.0,starstr
  ADD_TAG,starstr,'VREL',fltarr(nvisits)+999999.0,starstr
  ADD_TAG,starstr,'VRELERR',fltarr(nvisits)+999999.0,starstr
  ADD_TAG,starstr,'CCPEAK',fltarr(nvisits)+999999.0,starstr
  ADD_TAG,starstr,'CHISQ',fltarr(nvisits)+999999.0,starstr
  ADD_TAG,starstr,'CCF',fltarr(401,nvisits)+999999.0,starstr
  ADD_TAG,starstr,'CCFERR',fltarr(401,nvisits)+999999.0,starstr
  ADD_TAG,starstr,'CCFLAG',lindgen(401)-200,starstr
  ADD_TAG,starstr,'CCFW0',999999.0,starstr
  ADD_TAG,starstr,'CCFDW',999999.0,starstr  ; log dw
  ADD_TAG,starstr,'AUTOCF',fltarr(401)+999999.0,starstr
  return
endif

; now cross-correlate combined spectrum with grid to get absolute RV
print,'Cross-correlating with minigrid: '

aprvprep,grid.wave,starstr,specout,errout,/gaps,fixmasked=fixmasked, reject = reject, pl=pl
 
; new routine
apxcorr_newgrid,grid,specout[*,ref],errout[*,ref],vout,bestgrid
print,'Teff: ',grid.teff[bestgrid],' log g: ',grid.logg[bestgrid],' [Fe/H]',grid.metals[bestgrid]
; old routine
;print,'old routine'
;apxcorr_template,grid,specou[t*,0],errout[*,0],vout,bestgrid

; correct individual visit RVs
  allvisits.synthfile=synthfile
  allvisits.rv_teff=float(grid.teff[bestgrid])
  allvisits.rv_logg=float(grid.logg[bestgrid])
  allvisits.rv_feh=float(grid.metals[bestgrid])
  if tag_exist(grid,'ALPHA') then allvisits.rv_alpha=float(grid.alpha[bestgrid]) else allvisits.rv_alpha=99.99
  if tag_exist(grid,'CARBON') then allvisits.rv_carb=float(grid.carbon[bestgrid]) else allvisits.rv_carb=99.99
  allvisits.vrel+=vout.vrel
  allvisits.vrelerr = sqrt(allvisits.vrelerr^2+vout.vrelerr^2)  ; add errors
  ; VHELIO = VRAD + BC
  allvisits.vhelio = allvisits.vrel + allvisits.bc
  IF keyword_set(synth) THEN BEGIN
    allvisits.synthvrel+=vout.vrel
    allvisits.synthvrelerr = sqrt(allvisits.synthvrelerr^2+vout.vrelerr^2)
    allvisits.synthvhelio = allvisits.synthvrel + allvisits.bc
  ENDIF
  gdvis = where(finite(allvisits.vrel),ngdvis, comp=bdvis)
  allvisits.vtype = -1
  IF ngdvis GT 0 THEN allvisits[gdvis].vtype = 3

; if we got a bad RV from the grid, give up
if finite(vout.vrel) eq 0 then goto,outrv
  allvisits.obsvrel = allvisits.vrel
  allvisits.obsvrelerr = allvisits.vrelerr
  allvisits.obsvhelio = allvisits.vhelio
  IF keyword_set(synth) THEN BEGIN
    IF nvisits GT 1 THEN BEGIN 
     vscatter = stddev(allvisits.vhelio - median(allvisits.vhelio, /even),/nan)    
     synthvscatter = stddev(allvisits.synthvhelio - median(allvisits.synthvhelio, /even),/nan)
    ENDIF ELSE BEGIN
     vscatter = allvisits.vrelerr
     synthvscatter = allvisits.synthvrelerr
     IF synthvscatter LT vscatter THEN usesynth = 1 ELSE usesynth = 0
    ENDELSE
    IF usesynth THEN BEGIN 
     bestvscatter = synthvscatter
     allvisits.vrel = allvisits.synthvrel
     allvisits.vrelerr = allvisits.synthvrelerr
     allvisits.vhelio = allvisits.synthvhelio
     print, 'Final Combination Using synth'
    ENDIF ELSE BEGIN
     bestvscatter=vscatter
     ;vhelio already contains obsvhelio.
     print, 'Final Combination Using obs'
    ENDELSE
  ENDIF ELSE usesynth = 0 ; DEFAULT: Use vrel for combination

; one final combination with the final velocities to put on best zero-velocity scale
; this one includes the LSF!
apvisitcomb,allstr,allvisits,starstr,sinc=sinc,log=log,synth = usesynth
  
; Cross-correlate each visit spectrum with the best-fit template one
; last time.  This might be better than the relative RVs for hot stars
  aprvprep,grid.wave,starstr,specout,errout,/gaps,fixmasked=fixmasked, reject=reject,pl=pl
  if nvisits gt 1 then nspec=nvisits+2 else nspec=1
  ccf = fltarr(n_elements(vout.ccf),nspec)
  ccferr = ccf
  rvchisq = fltarr(nspec)
  synth_vrel = fltarr(nspec)
  synth_vrelerr = fltarr(nspec)
  for i=0,nspec-1 do begin
    apxcorr,grid.wave,grid.ndata[bestgrid,*],specout[*,i],errout[*,i],rvout,sum=sum,/errccf,pixlim=grid.pixlim
    ccf[*,i] = rvout.ccf  
    ccferr[*,i] = rvout.ccferr
    rvchisq[i] = rvout.chisq
    synth_vrel[i] = rvout.vrel
    synth_vrelerr[i] = rvout.vrelerr
  endfor
  ; the visit spectra were doppler shifted by allvisits.vrel
  ; so to get the "real" doppler shift we need to do:
  ; synth_vrel + allvisits.vrel
  if nvisits gt 1 then begin
    allvisits.synthvrel = synth_vrel[2:*]+allvisits.vrel
    allvisits.synthvrelerr = synth_vrelerr[2:*]
    allvisits.synthvhelio = allvisits.synthvrel+allvisits.bc
    allvisits.chisq = rvchisq[2:*]
  endif else begin
    allvisits.synthvrel = synth_vrel+allvisits.vrel
    allvisits.synthvrelerr = synth_vrelerr
    allvisits.synthvhelio = allvisits.synthvrel+allvisits.bc
    allvisits.chisq = rvchisq
  endelse
  IF keyword_set(synth) THEN BEGIN
    ; To be fair in synth mode we need to do one more run with observed template (if there is at least two visits)
   IF nvisits GT 1 THEN BEGIN
    ;Normalize and mask observed spectra
    aprvprep,starstr.wave,starstr,specout,errout,/gaps, fixmasked=fixmasked,reject=reject,pl=pl
    add_tag, starstr, 'nspec', specout, starstr
    add_tag, starstr, 'nerr', errout, starstr

    ;DERIVE RVs using observed spectrum as template
    ;apgundef,maxlag	
    ;if niter gt 2 then if max(vout.ccp_pars[2]*2.35) lt 75 n maxlag=75
    IF keyword_set(chisq) THEN BEGIN
     print, 'Using chisq mode-NOT IMPLEMENTED YET'
     stop
     ;obsvout = aprv_template_nofile(allvisits,starstr, ref = 0, maxlag=maxlag, dir_plots = plotdir, pl=pl, /apstar)   
    ENDIF ELSE BEGIN
     apxcorr,starstr.wave,specout[*,ref],specout[*,2:nvisits+1],$
     errout[*,2:nvisits+1],obsvout,maxlag=maxlag,pixlim=starstr.pixlim_overlap
    ENDELSE
    allvisits.obsvrel = obsvout.vrel + allvisits.vrel
    allvisits.obsvrelerr = obsvout.vrelerr
    allvisits.obsvhelio = allvisits.obsvrel+allvisits.bc
    newobsvscatter = stddev(allvisits.obsvhelio,/nan)
    newsynthvscatter = stddev(allvisits.synthvhelio,/nan)
   ENDIF ELSE BEGIN 
    newobsvscatter = 999999.
    newsynthvscatter = allvisits.synthvrelerr
   ENDELSE ;nvisits GT 1

    ; Check to see if the re-derivation of synthvhelio or vhelio is better than the current best velocity.
     ; If final values produce the smallest scatter, combine the spectra one last time using these velocities, and use that spectrum as the final combined spectrum
   IF (newobsvscatter LT bestvscatter) OR (newsynthvscatter LT bestvscatter) THEN BEGIN
     ; Pick the best of the two
     IF newobsvscatter LT newsynthvscatter THEN BEGIN 
      newvscatter = newobsvscatter 
      allvisits.vrel = allvisits.obsvrel
      allvisits.vrelerr = allvisits.obsvrelerr
      allvisits.vhelio = allvisits.obsvhelio
      gdvis = where(finite(allvisits.vrel),ngdvis, comp=bdvis)
      allvisits.vtype = -1
      IF ngdvis GT 0 THEN allvisits[gdvis].vtype = 4
      usesynth = 0
      print, 'Using final obs velocities'
     ENDIF ELSE BEGIN 
      newvscatter = newsynthvscatter
      allvisits.vrel = allvisits.synthvrel
      allvisits.vrelerr = allvisits.synthvrelerr
      allvisits.vhelio = allvisits.synthvhelio
      gdvis = where(finite(allvisits.vrel),ngdvis, comp=bdvis)
      allvisits.vtype = -1
      IF ngdvis GT 0 THEN allvisits[gdvis].vtype = 5
      usesynth = 1
      print, 'Using final synth velocities'
     ENDELSE

     apvisitcomb,allstr,allvisits,starstr,sinc=sinc,log=log,synth=usesynth
    if usesynth then aprvprep,grid.wave,starstr,specout,errout,/gaps,fixmask=fixmask, reject=reject else aprvprep,starstr.wave,starstr,specout,errout,/gaps,fixmask=fixmask, reject=reject
   ENDIF; (newobsvscatter LT bestvscatter) OR (newsynthvscatter LT bestvscatter)
  ENDIF ;keyword_set(synth)

; Get autocorrelation of best-fit template, copied from apxcorr_newgrid.pro
err = fltarr(n_elements(grid.ndata[bestgrid,*]))+0.001
apxcorr,grid.wave,grid.ndata[bestgrid,*],reform(grid.ndata[bestgrid,*]),err,autoout,sum=sum,pixlim=grid.pixlim

; add grid parameters to starstr after it is created
ADD_TAG,starstr,'RV_TEFF',float(grid.teff[bestgrid]),starstr
ADD_TAG,starstr,'RV_LOGG',float(grid.logg[bestgrid]),starstr
ADD_TAG,starstr,'RV_FEH',float(grid.metals[bestgrid]),starstr
ADD_TAG,starstr,'RV_ALPHA',999999.0,starstr
ADD_TAG,starstr,'RV_CARB',999999.0,starstr
if tag_exist(grid,'ALPHA') then starstr.rv_alpha=float(grid.alpha[bestgrid])
if tag_exist(grid,'CARBON') then starstr.rv_carb=float(grid.carbon[bestgrid])

; Add CCF peak, fwhm, chisq, etc.
ADD_TAG,starstr,'XSHIFT', vout.xshift,starstr
ADD_TAG,starstr,'VREL', vout.vrel,starstr
ADD_TAG,starstr,'VRELERR', vout.vrelerr,starstr
ADD_TAG,starstr,'CCPEAK', vout.ccpeak,starstr
ADD_TAG,starstr,'CCPFWHM', vout.ccpfwhm,starstr
ADD_TAG,starstr,'CHISQ', vout.chisq,starstr
ADD_TAG,starstr,'AUTOFWHM', vout.autofwhm,starstr

; Add CCF for all spectra
ADD_TAG,starstr,'CCF',ccf,starstr
ADD_TAG,starstr,'CCFERR',ccferr,starstr
ADD_TAG,starstr,'CCFLAG',rvout.cclag,starstr
ADD_TAG,starstr,'CCFW0',vout.w0,starstr
ADD_TAG,starstr,'CCFDW',vout.dw,starstr  ; log dw
ADD_TAG,starstr,'AUTOCF',autoout.ccf,starstr

; Fix output values for visits with BAD RVs
outrv:
gdrv = where(finite(allvisits.vrel) eq 1,ngdrv,comp=bdrv,ncomp=nbdrv)
if nbdrv gt 0 then begin
  allvisits[bdrv].vrel = 999999.0
  allvisits[bdrv].vrelerr = 999999.0
  allvisits[bdrv].vhelio = 999999.0
  allvisits[bdrv].vlsr = 999999.0
  allvisits[bdrv].vgsr = 999999.0
  allvisits[bdrv].chisq = 999999.0
  allvisits[bdrv].rv_teff = 999999.0
  allvisits[bdrv].rv_logg = 999999.0
  allvisits[bdrv].rv_feh = 999999.0
  allvisits[bdrv].rv_alpha = 999999.0
  allvisits[bdrv].rv_carb = 999999.0
  allvisits[bdrv].synthvrel = 999999.0
  allvisits[bdrv].synthvrelerr = 999999.0
  allvisits[bdrv].synthvhelio = 999999.0
  allvisits[bdrv].obsvrel = 999999.0
  allvisits[bdrv].obsvrelerr = 999999.0
  allvisits[bdrv].obsvhelio = 999999.0
  ; CCF stuff in starstr is already ZERO
endif

; Make the RV and CCF plot
IF keyword_set(bccomb) && (bccomb LT 0) THEN ngdrv = 0
if ngdrv gt 0  and keyword_set(plotdir) then aprvrefine_plot,allvisits,starstr,grid,bestgrid,vout0,plotdir

end