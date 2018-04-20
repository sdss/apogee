pro aprvprep,wref,allstr,specout,errout,ncorder=ncorder,gaps=gaps,masked=masked,$
             reject=reject,fixmasked=fixmasked,pl=pl

;+
;
; APRVPREP
;
; Prepares one or more spectra for cross correlation by normalizing and
; fixing bad pixels, then resmpling to reference wavelength grid
; Bad pixels are set to ZERO/Nan
;
; INPUTS:
;  wave     wavelength array
;  allstr   structure with input data
;  specout  output normalized and fixed spectrum
;  error    output normalized and fixed error
;  ncorder= order to use for continuum removal
;  /gaps    if specified, use pixlim data to fit continuum in chunks
;  /reject  Reject outliers using the combined spectrum (first one).
;             This is to reject residual sky lines.
; fixmasked= Run apfixmasked to interpolate over masked pixels. If set to 2 then, masked pixels are set to NaN.
; OUTPUTS:
; 
;-

;APNORMSPEC params
noerrcorr = 1; Default = 1 (Don't correct for error array)
IF (keyword_set(fixmasked) && (fixmasked EQ 2)) THEN setnan = 1
; determine number of spectra to fix
sz=size(allstr.spec)
if sz[0] eq 2 then nspec=sz[2] else nspec=1

; determine value for pixels marked as bad;
badpixvalue = !Values.F_nan ;0.0

; determine input tags
tags=tag_names(allstr)
specind = where(tags eq 'SPEC',nerrind)
waveind = where(tags eq 'WAVE',nerrind)
errind = where(tags eq 'ERR',nerrind)
maskind = where(tags eq 'MASK',nmaskind)
skyind = where(tags eq 'SKY',nskyind)
; add continuum tag if needed
contind= where(tags eq 'CONTINUUM',ncontind)
if ncontind eq 0 then ADD_TAG,allstr,'CONTINUUM',allstr.(specind[0])*0,allstr
contind= where(tags eq 'CONTINUUM',ncontind)

; number of output wavelengths
nwave=n_elements(wref)
; declare output arrays
specout=fltarr(nwave,nspec)
errout=fltarr(nwave,nspec)
masked=fltarr(nwave,nspec)


; loop over each spectrum
for i=0,nspec-1 do begin
  good=intarr(nwave)

  ; do we need to do spectrum in chunks?
  if keyword_set(gaps) then begin
     pixlim=allstr.pixlim_overlap
     sz=size(pixlim,/dim)
     nchip=sz[1]
     if not keyword_set(ncorder) then ncorder=3
  endif else nchip=1

  if keyword_set(pl) then !p.multi=[0,1,4]
  for ichip=0,nchip-1 do begin
  
    if keyword_set(gaps) then begin
       lo=pixlim[0,ichip]
       hi=pixlim[1,ichip]
    endif else begin
      lo=0
      hi=n_elements(allstr.spec[*,i])-1
    endelse

    IF lo GE hi THEN continue ;if any pixlim has lower limit higher than upper limit, then the chip is bad. Skip to the next chip, and use default values for specout and errout.

    ; extract this spectrum, with chunk if needed, and normalize
    sz=size(allstr.wave)
    ; if we are only given one wavelength array, then assume all spectra are on same scale
    if sz[0] eq 1 then w=allstr.wave else w=allstr.wave[*,i]
    str = {spec:allstr.spec[lo:hi,i],err:allstr.err[lo:hi,i],wave:w[lo:hi]}
    if nmaskind gt 0 then str=CREATE_STRUCT(str,'MASK',allstr.mask[lo:hi,i])
    if nskyind gt 0 then str=CREATE_STRUCT(str,'SKY',allstr.sky[lo:hi,i])

    if keyword_set(pl) then plot,str.spec, background = cgcolor('white'), color = cgcolor('black')

    if keyword_set(fixmasked) then begin
      fixmask = maskval('SIG_SKYLINE') or maskval('SIG_TELLURIC') or maskval('LITTROW_GHOST') or badmask()
      new=apfixmasked(str.spec,(str.mask and fixmask) gt 0, setnan=setnan)
      str.spec=new
    endif
    if keyword_set(pl) then oplot,str.spec,color=2
    ;apnormspec,str,ncorder=ncorder,/fixbadpix,nsky=4,growsky=5,/noerrcorr
    apnormspec,str,ncorder=ncorder,noerrcorr=noerrcorr,normtype=1,nsky=4,growsky=5;, /silent,pl=pl; DEFAULT: nsky = 4, growsky =5 ; Used nsky = 8 in one run
    ;apnormspec,str,ncorder=ncorder,/noerrcorr,normtype=1,/silent
    allstr.continuum[lo:hi,i]=str.continuum 
;    ; Mask Littrow ghost
;    littrowpix = where(str.mask AND maskval('LITTROW_GHOST'),nlittrowpix)
;    if nlittrowpix gt 0 then begin
;      str.nspec[littrowpix] = 0
;      str.masked[littrowpix] = 1
;    endif

    ; resample to reference wavelength scale
    resamp = 1
    if n_elements(wref) eq n_elements(w) then if max(abs(wref-w)) lt 1e-10 then resamp=0
    if keyword_set(resamp) then begin
      wsi=sort(str.wave)
      nw=n_elements(str.wave)
      gd=where(wref gt str.wave[wsi[20]] and wref lt str.wave[wsi[nw-20]],ngd)
      specout[gd,i] = SPLINE(str.wave[wsi],str.nspec[wsi],wref[gd],/double)
      errout[gd,i] = SPLINE(str.wave[wsi],str.err[wsi]/str.continuum[wsi],wref[gd],/double)
;      masked1 = SPLINE(str.wave[wsi],float(str.masked[wsi]),wref[gd],/double)
;      ; set "bad" pixels to zero
;      bd = where(masked1 gt 0.1,nbd)
;      if nbd gt 0 then begin
;        specout[gd[bd],i] = badpixvalue;0.0
;        masked[gd[bd],i] = 1  ; bad
;      endif
      good[gd]=1
     endif else begin
      ; set "bad" pixels to zero
;      bd = where(str.masked eq 1,nbd)
;      if nbd gt 0 then str.nspec[bd]= badpixvalue;0.0
      specout[lo:hi,i] = str.nspec
      errout[lo:hi,i] = str.err/str.continuum
      masked[lo:hi,i] = str.masked
      good[lo:hi]=1
    endelse
  endfor
  ; Reject outliers using combined spectrum
  ;   this is to reject residual sky lines
  if keyword_set(reject) and i ge 2 then begin
    gg = where(specout[*,i] gt 0.0 and specout[*,0] gt 0.0,ngg)
    if ngg gt 5 then begin
      frac = median(specout[gg,i]/specout[gg,0])
      diff = specout[*,i]-specout[*,0]*frac
      diff *= (specout[*,i] gt 0.0 and specout[*,0] gt 0.0)  ; set masked pixels to zero
      diff /= (errout[*,i] > 1e-5)
      diff -= medfilt1d(diff,201,/edge)
      bd = where(abs(diff) gt 5 and abs(specout[*,i]-median(specout[gg,i])) gt 0.1,nbd)
      if nbd gt 0 then begin  ; only allow up to 50 pixels to be masked this way
       si = reverse(sort(abs(diff[bd])))  ; sort in Nsig
        bd = bd[si[0:49<(nbd-1)]]   
        specout[bd,i] = badpixvalue; 0.0
        masked[bd,i] = 1
      endif
    endif
  endif

  if keyword_set(pl) then plot,specout[*,i], color = cgcolor('black')

  if not keyword_set(fixmasked) then begin
    ; Set high error or negative error pixels to zero
    bd = where(errout[*,i] gt 10 or errout[*,i] lt -0.01,nbd)
    if nbd gt 0 then begin
      specout[bd,i] = badpixvalue;0.0
      masked[bd,i] = 1
    endif
    if keyword_set(pl) then oplot,specout[*,i],color= cgcolor('red')
    ; Set very high and very low pixels to bad pixel value (0 or Fnan)
    gd = where(specout[*,i] gt 0.0,ngd)
    if ngd gt 50 then begin
      sig = mad(specout[gd,i]) > 0.01
      med = median(specout[gd,i])
      bd = where((specout[*,i]-med gt 5*sig) OR (specout[*,i] LT 0),nbd)
      if nbd gt 0 then begin
        specout[bd,i] = badpixvalue ;0.0
        masked[bd,i] = 1
      endif
    endif else begin
      bd = where(specout[*,i] gt 3,nbd)
      if nbd gt 0 then begin
        specout[bd,i] = badpixvalue;0.0
        masked[bd,i] = 1
      endif
    endelse
    if keyword_set(pl) then oplot,specout[*,i],color=cgcolor('green')
    
    ; set all unfilled pixels to zero/badpixvalue
    bd=where(good eq 0, nbd)
    if nbd gt 0 then begin
      specout[bd,i]= badpixvalue;0.
      errout[bd,i]=badpixvalue ;0.
    endif
    if keyword_set(pl) then oplot,specout[*,i],color=cgcolor('blue') 
  endif


  ;if keyword_set(pl) then stop
endfor

end
