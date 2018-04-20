pro apradvel_pieces,str,error=error,pl=pl,silent=silent,verbose=verbose,stp=stp

;+
;
; APRADVEL_PIECES
;
; This computes the radial velocity with ChiSq minimizations
; of an APOGEE Visit spectrum with a synthetic spectrum.
; The comparison is done in small pieces of the spectrum.
;
; INPUTS:
;  str      An input visit spectrum structure.  The structure
;              should have already been run through APXCORR_GRID.PRO
;              and have the best-fit radial velocity and synthetic
;              spectrum in it (syn_wave and syn_spec).
;              STR gets updated with the output paramaters.
;  /pl       Plot the spectra.
;  /silent   Don't print anything to the screen.
;  /verbose  Print lots of information to the screen.
;  /stp      Stop at the end of the program.
;
; OUTPUTS:
;  The "str" structure is updated with the output parameters
;  =error  The error message if there was one.
;
; USAGE:
;  IDL>apradvel_pieces,str
;
; By D.Nidever  July 2010
;-

cspeed = 2.99792458d5  ; speed of light in km/s

; Not enough inputs
if n_elements(str) eq 0 then begin
  print,'Syntax - apradvel_pieces,str,error=error,pl=pl,silent=silent,verbose=verbose,stp=stp'
  error = 'Not enough inputs'
  return
endif

; Best-fit synthetic spectrum not in STR
if tag_exist(str,'syn_wave') eq 0 or tag_exist(str,'syn_spec') eq 0 then begin
  error = 'SYN_WAVE/SYN_SPEC NOT FOUND IN INPUT STRUCTURE'
  print,error
  return
endif


; Fitting the synthetic spectra to the observed chunks
;------------------------------------------------------

npix = n_elements(str.spec[*,0])
nchunks = 15
nchpix = npix/nchunks

chunkstr = REPLICATE({chip:0L,chunk:0L,midwave:0.0,snr:0.0,ew:0.0,vrel:0.0,vrelerr:0.0,status:-1L,$
                      chisq:0.0,dof:0L,niter:0L},nchunks*3)
chunknum = 0L

if keyword_set(verbose) then begin
  print,' Chip Chunk SNR    EW    RVEL  RVERR  CHISQ  NITER'
endif

; Loop through the chips
For i=0,2 do begin

  ; Loop through the chunks
  For j=2,nchunks-1 do begin

   ; SKIP the first TWO chunks, the synth wavelengths don't cover it!!!!

    ; Get the spectrum
    x = findgen(nchpix)+j*nchpix
    wobs = str.wave[j*nchpix:(j+1)*nchpix-1,i]       ; wavelength
    ;sobs = str.nspec[j*nchpix:(j+1)*nchpix-1,i]      ; normalized spectrum
    spec = str.spec[j*nchpix:(j+1)*nchpix-1,i]      ; normalized spectrum
    cont = str.continuum[j*nchpix:(j+1)*nchpix-1,i]  ; continuum
    sobs = spec/cont         ; renomalize to make sure no pixels were set to continuum
    ;invobs = str.var[j*nchpix:(j+1)*nchpix-1,i]      ; inverse variance
    ;eobs = sqrt(1.0/invobs)/cont                     ; normalized error spectrum
    eobs = str.err[j*nchpix:(j+1)*nchpix-1,i]        ; error
    eobs = eobs/cont                                 ; normalized error spectrum
    mask = str.mask[j*nchpix:(j+1)*nchpix-1,i]       ; mask
    sky = str.sky[j*nchpix:(j+1)*nchpix-1,i]         ; sky
    snr = median(sobs/eobs)

    ; Sort by wavelength
    wsi = sort(wobs)
    wobs = wobs[wsi]
    sobs = sobs[wsi]
    cont = cont[wsi]
    eobs = eobs[wsi]
    mask = mask[wsi]
    sky = sky[wsi]

    ; Extend the observed wavelength range on the edges for the model
    wcoef = str.wcoef[*,i]
    next = 20
    xext = findgen(nchpix+2*next)+j*nchpix-next
    wobs_ext = PIX2WAVE(xext,wcoef)
    wsi2 = sort(wobs_ext)
    wobs_ext = wobs_ext[wsi2]

    ; Get the LSF coefficients at the center of this chunk
    lcoef = str.lcoef[*,i]
    xcen = j*nchpix+nchpix/2
    xpix = findgen(21)-10 + xcen
    lsf = LSF_GH(xpix,xcen,lcoef)

  ; I SHOULD CONVOLVE WITH THE LSF BEFORE THE FITTING, SO IT ONLY
  ; HAS TO BE DONE ONCE!!

    ; Fit the synthetic spectrum to the observed spectrum
    ; Get synthetic spectral indices with the wavelength range
    ind = where(str.syn_wave ge min(wobs)-40.0 and str.syn_wave le max(wobs)+40.0,nind)
    if nind eq 0 then begin
      if keyword_set(verbose) then print,'Synthetic spectrum not in the observed wavelength range'
      goto,BOMB1
    endif
    ; Only want goood observed pixels
    gdpts = where(sobs gt 0.0 and (mask and badmask()) eq 0 and sky lt 10*cont,ngdpts)   ; only want good points to fit
    if ngdpts lt 10 then begin
      if keyword_set(verbose) then print,'Not enough good pixels'
      goto,BOMB1
    endif
    ; Check that the flux isn't negative
    if median(sobs) lt 0.0 then begin
      if keyword_set(verbose) then print,'Continuum is negative'
      goto,BOMB1
    endif

    ; Fit the Radial Velocity
    fa = {synwave:str.syn_wave[ind], synspec:str.syn_spec[ind], lsf:lsf, minwave:min(wobs), maxwave:max(wobs),$
          wext:wobs_ext, next:next,gdpts:gdpts}
    initpar = [str.vrel]
    par = MPFITFUN('fit_radvel',wobs[gdpts],sobs[gdpts],eobs[gdpts],initpar,functargs=fa,yfit=yfit1,parinfo=parinfo,$
                   status=status,perror=perror,dof=dof,bestnorm=bestnorm,niter=niter,/quiet)
    if status lt 0 then begin
      print,'ERROR: mpfitfun failed in apradvel_pieces...'
      goto,BOMB1
    endif
    pcerror = perror * sqrt(bestnorm/dof)
    chisq = sqrt(bestnorm/dof)

    ; Get shifted synthetic spectrum for the whole thing
    yfit = fit_radvel(wobs,par,synwave=str.syn_wave[ind],synspec=str.syn_spec[ind],lsf=lsf,minwave=min(wobs),$
                      maxwave=max(wobs),wext=wobs_ext,next=next,gdpts=indgen(n_elements(wobs)))

    ; Calculate the equivalent width of the spectrum
    ;   This gives a measure of the number of lines there are
    ew = TOTAL(1.0-sobs)*median(slope(wobs))

    ; Print out the values
    if keyword_set(verbose) then begin
      fmt = '(I4,I5,I6,2F7.2,F7.2,F7.2,I5)'
      print,format=fmt,i+1,j+1,snr,ew,par[0],pcerror[0],chisq,niter
    end

    ; Plug values into the structure
    chunkstr[chunknum].chip = i
    chunkstr[chunknum].chunk = j
    chunkstr[chunknum].midwave = mean(minmax(wobs))
    chunkstr[chunknum].snr = snr
    chunkstr[chunknum].ew = ew
    chunkstr[chunknum].vrel = par[0]
    chunkstr[chunknum].vrelerr = pcerror[0]
    chunkstr[chunknum].status = status
    chunkstr[chunknum].chisq = chisq
    chunkstr[chunknum].dof = dof
    chunkstr[chunknum].niter = niter

    ; Plotting
    ;pl=1 ;0
    if keyword_set(pl) then begin
      xr = minmax(wobs)
      yr = minmax([sobs,yfit])
      yr = [yr[0]-0.1*range(yr),yr[1]+0.1*range(yr)]
      plot,wobs,sobs,xr=xr,yr=yr,xs=1,ys=1,tit='Chip '+strtrim(i+1,2)+' Chunk '+strtrim(j+1,2)+'/'+strtrim(nchunks,2)
      xyouts,mean(xr),yr[1]-0.05*range(yr),'VREL='+stringize(par[0],ndec=3)+'+/-'+stringize(pcerror[0],ndec=3)+$
             ' km/s  ChiSq='+stringize(chisq,ndec=2),align=0.5
      oplot,wobs,yfit,co=250
      wait,1
      ;stop
    endif

    BOMB1:

    chunknum++  ; increment the chunk counter

    ;stop

  End ; chunk loop

End ; chip loop

;stop

;=========================
; Analyze the results
;=========================
gdchunk = where(chunkstr.status gt 0,ngdchunk)
if ngdchunk lt 3 then begin
  error = 'Not enough good chunks'
  return
endif
chunkstr2 = chunkstr[gdchunk]

; Initial Estimates
;------------------
; median
medvrel0 = median([chunkstr2.vrel])

; Sigma clipping
;sig0 = MAD(chunkstr2.vrel-medvrel0,/zero) > 0.05
;gd = where(abs(chunkstr2.vrel-medvrel0) lt sig0*5.0,ngd)
ROBUST_MEAN,chunkstr2.vrel,robmean,robsig,ind=gd,numused=ngd

; Final Values
;--------------
medvrel = median(chunkstr2[gd].vrel)

; weighted mean, use SNR, EW, VRELERR, and CHISQ
;wt = (chunkstr2[gd].snr*chunkstr2[gd].ew/(chunkstr2[gd].vrelerr*chunkstr2[gd].chisq))^2
wt = (chunkstr2[gd].snr*chunkstr2[gd].ew/chunkstr2[gd].chisq)^2
wtvrel = TOTAL(wt*chunkstr2[gd].vrel)/TOTAL(wt)
; Sigma scatter in the observations
;sigvrel = MAD(chunkstr2[gd].vrel-wtvrel,/zero)
sigvrel = sqrt( total( ((chunkstr2[gd].vrel-wtvrel)^2.)*wt ) * ngd / (( ngd-1.) * total(wt)) )
; Final error in the mean
vrelerr = sigvrel/sqrt(ngd)

; Final chisq
;--------------
; rchisq = sqrt(total(diff^2)/dof)
;chisq1 = sqrt( total( (chunkstr2.chisq^2)*chunkstr2.dof ) / total(chunkstr2.dof) )

synwave = str.syn_wave*(1.0d0 + wtvrel/cspeed)
synspec = str.syn_spec

; Loop through the chips
diffsq = fltarr(npix,3)
numpts = 0
for i=0,2 do begin
  ; Get the spectrum
  wobs = str.wave[*,i]       ; wavelength
  si = sort(wobs)             ; sort by wavelength
  wobs = wobs[si]
  spec = str.spec[si,i]
  ;nspec = str.nspec[si,i]      ; normalized spectrum
  ;spec = str.spec[si,*]       ; original spectrum
  cont = str.continuum[si,i]  ; continuum
  nspec = spec/cont
  mask = str.mask[si,i]
  err = str.err[si,i]        ; error
  err = err/cont                                 ; normalized error spectrum
  sky = str.sky[si,i]         ; sky

  ; Interpolate the synthetic spectrum
  ind = where(wobs ge min(synwave) and wobs le max(synwave) and $
              nspec gt 0.0 and (mask and badmask()) eq 0 and sky lt 10*cont,nind)
  medcont = median(spec)         ; make sure the continuum is positive
  if nind gt 10 and medcont gt 0.0 then begin
    synspec1 = SPLINE(synwave,synspec,wobs[ind])

    ; Measure the difference
    diffsq[ind,i] = (synspec1-nspec[ind])^2/err[ind]^2
    numpts += nind    ; number of points used
  endif

endfor

chisq = sqrt( total(diffsq) / numpts )


; Printing
;----------
if not keyword_set(silent) then begin
  print,' Vrel = ',stringize(wtvrel,ndec=3),' km/s'
  print,' Vrelerr = ',stringize(vrelerr,ndec=3),' km/s'
  print,' Scatter = ',stringize(sigvrel,ndec=3),' km/s'
  snr = median(str.spec/str.err)
  print,' S/N = ',stringize(snr,ndec=1)
  ;print,' Chisq = ',stringize(chisq,ndec=2)
endif

; Add information to the STR structure
;--------------------------------------
radvel = {vrel:wtvrel, vrelerr:vrelerr, sigvrel:sigvrel, medvrel:medvrel, wtvrel:wtvrel, nused:ngd, $
          chisq:chisq, chunkstr:chunkstr}
str = CREATE_STRUCT(str,'radvel',radvel)

;; Update Vrel, Vrelerr, Vhelio in STR
;if (str.radvel.vrelerr lt str.xcorr.vrelerr) then begin
;  print,'Using Chisq Velocity values'
;  str.vrel = wtvrel
;  str.vrelerr = vrelerr
;; Using the XCORR value
;endif else begin
;  print,'Using the XCORR Velocity values'
;endelse
;str.vhelio = str.vrel+str.bc

;if str.xcorr.teff gt 8000. then stop

;stop

if keyword_set(stp) then stop

end
