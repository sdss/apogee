pro apvisitcomb,allvisits,visitstr,starstr,globalwt=globalwt,stp=stp,nolsffit=nolsffit,log=log,sinc=sinc,quick=quick,synth=synth

;+
;
; APVISITCOMB
;
; Combine all the visits into one spectrum and output
; the apStar file
;
; INPUTS:
;  allvisits   An array of pointers to the individual visit data
;  visitstr    Tables of the visit parameters
;  starstr     Output combined spectra
;  /synth      Use the synthvrel to combine spectra
;  /globalwt    When combining the visit spectra use a single
;                 weight for each spectrum (NOT pixel-by-pixel).
;  /stp        Stop at the end of the program.
;
; OUTPUTS:
;  starstr     The structure of the final spectrum
;
; USAGE:
;  IDL>apvisitcomb,allvisits,starstr
;
; By D.Nidever  July 2010
; Major mods J. Holtzman 3/12
; New functionality to combine using synthvrel vs. vrel added by J. Holtzman/N.Troup 5/2016
;-

cspeed = 2.99792458d5  ; speed of light in km/s

apgundef,starstr

; Not enough inputs
if n_elements(allvisits) eq 0 then begin
  print,'Syntax - apvisitcomb,allvisits,starstr,globalwt=globalwt,stp=stp'
  return
endif

nvisits = n_elements(allvisits)

; Choose whether to combine using vrel or synthvrel. Save combtype in starstr
; Idea for later: make combtype a bitmask to save all the settings one can toggle in combination (e.g. weighting, LSFfitting, quick mode on/off)
if keyword_set(synth) then begin 
  vrel=visitstr.synthvrel 
  combtype = 1
endif else begin 
  vrel=visitstr.vrel
  combtype = 0
endelse

if keyword_set(log) then begin

  ;npix = 7117
  ;w0=4.1808101d0
  ;dw=6.0831101e-6
  ;npix = 8850
  ;dw=4./cspeed*alog10(exp(1))
  w0=4.179d0
  npix=8575
  dw=6.d-6
  wave_final = w0+indgen(npix)*dw
  wave_final= 10^wave_final
  p2w_coef = [w0,dw]

endif else begin

  ; We need to shift the spectra to REST wavelengths
  nord = 4 ;3 ;4  ; the order for the wavelength polynomial fit
  ; 3-good to within ~10-20%
  ; 4-good to within ~1%
  
  npix = n_elements((*allvisits[0]).spec[*,0])
  
  ; Loop through the visits
  For i=0,nvisits-1 do begin
  
    ; Get the rest wavelength range to cover
    str = (*allvisits[i])
    ;wave = str.wave
    wave = [str.wave[*,0], str.wave[*,1], str.wave[*,2]]
    wave_rest = wave/(1.0d0+vrel[i]/cspeed)
  
    wcoef = str.wcoef
  
    ; sort by wavelength
    wsi = sort(wave_rest)
    wave = wave[wsi]
    wave_rest = wave_rest[wsi]
  
    ; Create the pixel array
    xoffset = wcoef[0,1]+1023.5
    x = [ dindgen(npix)+wcoef[0,0], dindgen(npix)+wcoef[0,1], dindgen(npix)+wcoef[0,2] ]
    x -= xoffset
    ; Want the X-values to start at ONE
    x = x-min(x)+1
  
    ; Fit wavelength wrt pixels (starting with 1)
    ;   This fits to within about 1% of a pixel (if nord=4)
    coef = POLY_FIT((x)(*),(wave_rest)(*),nord)
    f = POLY(x,coef)
  
    ; We need to use the pixel to wavelength polynomial
    ; fit for the visit that has the LOWEST rest wavelength
  
    ; First visit
    if i eq 0 then begin
      minwave = min(wave_rest)
      maxwave = max(wave_rest)
      minx = min(x)
      maxx = max(x)
  
      p2w_coef = coef
  
    ; Other visits
    endif else begin
      minwave = minwave < min(wave_rest)
      maxwave = maxwave > max(wave_rest)
      minx = minx < min(x)
      maxx = maxx >max(x)
  
      if min(wave_rest) eq minwave then p2w_coef = coef
  
    endelse
  
  endfor

  ; How many pixels do we need to cover the entire wavelength range
  xbig = dindgen(ceil(maxx)+200L)+1  ; start with 1
  wbig = POLY(xbig,p2w_coef)
  lastpix = first_el(where(wbig le maxwave),/last)  ; last pixel to use
  npix = lastpix+1    ; number of pixels in the final array
  
  ; Final wavelength array
  wave_final = POLY(dindgen(npix)+1,p2w_coef)

endelse

; VISIT stack structure
;  all of the interpolated visit arrays go in here
one = fltarr(nvisits,npix)
oneint=intarr(nvisits,npix)
vstackstr = {spec:one,err:one,specnorm:one,errnorm:one,skynorm:one,cont0:one,$
             cont:one,mask:oneint,sky:one,skyerr:one,telluric:one,telerr:one}
vstackstr.spec = !values.f_nan  ; all bad for now
vstackstr.err = !values.f_nan  ; all bad for now
vstackstr.specnorm = !values.f_nan  ; all bad for now
vstackstr.errnorm = !values.f_nan  ; all bad for now
vstackstr.skynorm = !values.f_nan  ; all bad for now
vstackstr.mask = 0

; Initialize the LSF arrays
nLSFpix = 15
totlsfwt = fltarr(npix,nLSFpix)
;totwt = fltarr(npix)
totwt = fltarr(npix,nLSFpix)

; Set up the final pixel array
; Get the minimum pixel value

pixlim=intarr(2,3)
pixlim[0,*]=npix-1
pixlim[1,*]=0
pixlim_overlap=intarr(2,3)
pixlim_overlap[0,*]=0
pixlim_overlap[1,*]=npix-1
y = findgen(npix)

if keyword_set(sinc) and n_elements(sinc) lt 3 then begin
  print,'Must specify 3 elements of sinc for filter widths for each chip'
  stop
endif

; setup for masking
getmaskvals,flag,badflag,maskcontrib
maskarr=lonarr(n_elements(flag))
for ibit=0,n_elements(flag)-1 do maskarr[ibit]=2L^ibit
if keyword_set(quick) then maskarr = maskarr and (badmask() or maskval('SIG_SKYLINE') or maskval('SIG_TELLURIC'))

missingChipCount = [0,0,0]; Counts number of visits with a bad chip. Used to set pixlim to erorr value if needed.
; Loop through the visits
For i=0,nvisits-1 do begin

  str = (*allvisits[i])

;  fiber = long(first_el(strsplit(file_basename(str.file,'.fits'),'-',/extract),/last))

  ; Check that we have a doppler shift
  if finite(vrel[i]) eq 0 then begin
    vstackstr.mask[i,*] = 1
    vstackstr.cont[i,*] = !values.f_nan
    goto,BOMB
  endif

  wcoef = str.wcoef
  ; Create the pixel array
  xoffset = wcoef[0,1]+1023.5
  xpix = dindgen(npix)-xoffset

  lsf2d = fltarr(npix,nLSFpix)

  ; Loop through the chips
  ; npix per resolution element for each chip (~*0.9)
  For j=0,2 do begin

    ; SINC interpolation
    if keyword_set(sinc) then begin
      wave=wave_final*(1.0d0+vrel[i]/cspeed)
      pix=wave2pix(wave,str.wave[*,j])
      ind=where(finite(pix))
      err_interp=str.err[*,j]
 
      ; set flux and errors in bad pixels to zero so that they won't unduly affect nearby regions
      ; later we will mask any pixels in which the bad pixel had a significant contribution
      ; but if bad values are very large, they can affect more distant pixels even if they
      ; are not significant contributors
      tmp_interp=str.spec[*,j]
      bd=where(str.mask[*,j] and  badmask(),nbd)
      if nbd gt 0 then begin
        ;tmp_interp[bd]=0.
        ;err_interp[bd]=0.
        tmp_interp[bd]=!values.f_nan
        err_interp[bd]=!values.f_nan
        tmpflux=smooth(medfilt1d(tmp_interp,501,edge=2),100,/nan)
        tmperr=smooth(medfilt1d(err_interp,501,edge=2),100,/nan)
        tmp_interp[bd]=tmpflux[bd]
        err_interp[bd]=tmperr[bd]
      endif
      spec_interp = sincint( pix[ind], sinc[j], tmp_interp ,err_interp)

      ; get output pixel limits for input chip limits
      ilim=[min(where(finite(spec_interp))),max(where(finite(spec_interp)))]
      xlim=ind[ilim]
      ; pixlim gives maximum range from ANY input spectra
      ; pixlim_overlap gives overl range that occurs in ALL input spectra  
      pixlim[0,j]=min([pixlim[0,j],floor(xlim[0])])
      pixlim[1,j]=max([pixlim[1,j],ceil(xlim[1])])
      pixlim_overlap[0,j]=max([pixlim_overlap[0,j],ceil(xlim[0])])
      pixlim_overlap[1,j]=min([pixlim_overlap[1,j],floor(xlim[1])])

      ;olderr = sincint( pix[ind], sinc[j], str.err[*,j] )
      sky_interp =  sincint( pix[ind], sinc[j], str.sky[*,j] )
      if keyword_set(quick) then begin
        skyerr_interp=fltarr(n_elements(spec_interp))
        telluric_interp=fltarr(n_elements(spec_interp))
        telerr_interp=fltarr(n_elements(spec_interp))
      endif else begin
        skyerr_interp =  sincint( pix[ind], sinc[j], str.skyerr[*,j] )
        telluric_interp =  sincint( pix[ind], sinc[j], str.telluric[*,j] )
        telerr_interp =  sincint( pix[ind], sinc[j], str.telerr[*,j] )
      endelse
      ; special handling for propagation of mask bits
      mask_interp=intarr(n_elements(spec_interp))
      for ibit=0,n_elements(maskarr)-1 do begin
       if maskarr[ibit] gt 0 then begin
        mask1=str.mask[*,j] and maskarr[ibit]
        junk=where(mask1 gt 0,nmask1)
        ; only continue if we have any of these bits set!
        if nmask1 gt 0 then begin
          tmp = sincint( pix[ind], sinc[j], float(mask1)<1.)
          bd = where(abs(tmp) gt maskcontrib[ibit],nbd)
          if nbd gt 0 then mask_interp[bd] = mask_interp[bd] or maskarr[ibit]
        endif
       endif
      endfor
      ; set bad pixels to NaN
      ;bd=where(mask_interp and badmask(),nbd)
      ;if nbd gt 0 then begin
      ;  ;spec_interp[bd]=0.
      ;  err_interp[bd]=baderr()
      ;endif

    ; SPLINE interpolation
    endif else begin
 
      ; Rest wavelengths
      wave = str.wave[*,j]
      wave_rest = wave/(1.0d0+vrel[i]/cspeed)

      ; sort by wavelength
      wsi = sort(wave_rest)
      xpix = xpix[wsi]
      wave = wave[wsi]
      wave_rest = wave_rest[wsi]
      spec = str.spec[wsi,j]
      err = str.err[wsi,j]
      mask = str.mask[wsi,j]
      sky = str.sky[wsi,j]
      skyerr = str.skyerr[wsi,j]
      telluric = str.telluric[wsi,j]
      telerr = str.telerr[wsi,j]

      ; Trim the edge pixels
      lo = 1
      hi = n_elements(wave)-2 ;4094

      ; What "final" wavelengths are within our range
      ind = where(wave_final ge min(wave_rest[lo:hi]) and wave_final le max(wave_rest[lo:hi]),nind)

      ; Now interpolate onto the final wavelength grid
      spec_interp = SPLINE(wave_rest[lo:hi],spec[lo:hi],wave_final[ind],/double)
      err_interp = SPLINE(wave_rest[lo:hi],err[lo:hi],wave_final[ind],/double)
      mask_interp = SPLINE(wave_rest[lo:hi],mask[lo:hi],wave_final[ind],/double)
      mask_interp = ( round(mask_interp) > 0) < 1
      sky_interp = SPLINE(wave_rest[lo:hi],sky[lo:hi],wave_final[ind],/double)
      skyerr_interp = SPLINE(wave_rest[lo:hi],skyerr[lo:hi],wave_final[ind],/double)
      telluric_interp = SPLINE(wave_rest[lo:hi],telluric[lo:hi],wave_final[ind],/double)
      telerr_interp = SPLINE(wave_rest[lo:hi],telerr[lo:hi],wave_final[ind],/double)

      ; Create the LSF array for this spectrum
      ; use lsfgh2d
      ; get the x-values of the final array (in terms of this one)
      ;  so we can get the correct LSF for these pixels
      ;xcenter = spline(wave_rest,xpix[*,j],wave_final[ind],/double)
      xcenter = spline(wave_rest,xpix,wave_final[ind],/double)
      nnew = n_elements(xcenter)
      ; Make 2D LSF array
      lsfpars = str.lcoef[*,j]
      dx = 1.0
      xlsf = REPLICATE(1.0d0,nnew)#(dindgen(nLSFpix)-nLSFpix/2)*dx
      xlsf += xcenter#REPLICATE(1.0d0,nLSFpix)
      lsf2d = LSF_GH(xlsf,xcenter,lsfpars)

      ; Add to TOTLSFWT which is the running sum of the
      ;  lsf times the weight for each pixel which is just 1/err^2

      ; Global weight per spectrum
      ;---------------------------
      if keyword_set(globalwt) then begin
        mederr = median(err_interp)   ; median error for entire specrum
        totlsfwt[ind,*] += lsf2d * (1.0/mederr^2)
        ; Add to TOTWT which is just the running sum of the weights (1/err^2)
        totwt[ind] += (1.0/mederr^2)
      ; Pixel-by-pixel weighting
      ;-------------------------
      endif else begin
        ; Weight every spectrum separately
        totlsfwt[ind,*] += lsf2d * ((1.0/err_interp^2)#replicate(1.0,nLSFpix))
        ; Add to TOTWT which is just the running sum of the weights (1/err^2)
        totwt[ind] += (1.0/err_interp^2)
      endelse

    endelse ; sinc/spline interpolation


    if not keyword_set(nolsffit) then begin
      ; Create the LSF array for this spectrum
      ; get the x-values of the final array (in terms of this one)
      ;  so we can get the correct LSF for these pixels
      ;xcenter = spline(wave_rest,xpix,wave_final[ind],/double)
      ;nnew = n_elements(xcenter)
      ; Make 2D LSF array
      lsfpars = str.lcoef[*,j]
      dx = slope(pix[ind]) & dx=[dx[0],dx]
      ;lsfpars[0] = abs(median(dx))  ; change binsize
      xlsf = REPLICATE(1.0d0,n_elements(ind))#(dindgen(nLSFpix)-nLSFpix/2)
      xlsf *= dx#replicate(1,nLSFpix)
      xlsf += pix[ind]#REPLICATE(1.0d0,nLSFpix)
      if lsfpars[0] eq 1 then lsf2d[ind,*]=LSF_GH(xlsf/2.,pix[ind]/2.,lsfpars) else $
        lsf2d[ind,*]=LSF_GH(xlsf,pix[ind],lsfpars)

      ; Make sure there are NO negative LSF values
      lsf2d >= 0

      ;lsf2d[ind,*] = reverse(lsf2d[ind,*],1)

      ; BUT WHAT PIXEL SCALE ARE THE INDIVIDUAL LSFs HERE ON?
      ;  the dither-combined pixel scale, the final pixel scale????
      ; the pixel steps in the 2nd dimension need to be on the *ORIGINAL*
      ; scale!!!!
      ;
      ; and I think I need to flip the LSFs in the 2nd dimension
      ; since now wavlength is INCREASING with pixels and before
      ; it was DECREASING.
      ;
      ; How about BINSIZE in LSFCOEF?? does that need to change too??


      ;; Interpolate the LSF from file
      ;;  this is a double-check that I'm doing it correctly
      ;if j eq 0 then lsf3 = lsf2d*0
      ;flsf = reform(lsfstr.(j).lsf[*,300-fiber,*])
      ;fwave = wavestr.(j).wave[*,300-fiber]
      ;dfwave = slope(fwave) & dfwave=[dfwave[0],dfwave]
      ;dw = slope(wave) & dw=[dw[0],dw]
      ;for k=0,n_elements(ind)-1 do begin
      ;  flsf1 = reform(flsf[pix[ind[k]]/2.,*])
      ;  fwave0 = reform(fwave[pix[ind[k]]/2.,*])
      ;  fwave1 = (dindgen(21)-21/2)*dfwave[pix[ind[k]]/2.] ;+fwave0[0]
      ;  si = sort(fwave1)
      ;  fwave1=fwave1[si] & flsf1=flsf1[si]
      ;  ;wave2 = (dindgen(21)-21/2)*dw[k]+wave[ind[k]]
      ;  wave2 = wave[ind[k]-21/2:ind[k]+21/2]-wave[ind[k]]
      ;  new = spline(fwave1,flsf1,wave2)
      ;  lsf3[ind[k],*] = new
      ;endfor

    endif

    ; Stuff it in the VSTACKSTR
    vstackstr.spec[i,ind] = spec_interp
    vstackstr.err[i,ind] = err_interp
    vstackstr.mask[i,ind] = mask_interp
    vstackstr.sky[i,ind] = sky_interp
    vstackstr.skyerr[i,ind] = skyerr_interp
    vstackstr.telluric[i,ind] = telluric_interp
    vstackstr.telerr[i,ind] = telerr_interp

  Endfor  ; chip loop

  ; Remove the median so that missing/bad values in one or more spectra match
  bad = where(vstackstr.mask[i,*] and badmask(),nbad,complement=gddata)    ; only want good values
  tmpspec=vstackstr.spec[i,*]
  if nbad gt 0 then tmpspec[bad] = !values.f_nan
  cont = smooth(medfilt1d(tmpspec,501,edge=2),100,/nan)
  ;if nbad gt 0 then vstackstr.spec[i,bad] = !values.f_nan
  ;cont = medfilt1d(reform(vstackstr.spec[i,*]),501,edge=2)
  vstackstr.cont[i,*] = cont
  vstackstr.specnorm[i,gddata] = vstackstr.spec[i,gddata]/cont[gddata]       ; normalize spectrum
  vstackstr.errnorm[i,gddata] = vstackstr.err[i,gddata]/cont[gddata]         ; normalize error
  vstackstr.skynorm[i,gddata] = vstackstr.sky[i,gddata]/cont[gddata]         ; normalize error

;  ; Fit a low-order polynomial to the data
;  ;;  bad/missing values at the ends have flag=NAN
;  gddata = where(finite(vstackstr.mask[i,*]) eq 1,ngddata)    ; only want good values
;  ;yy = cgscalevector(findgen(npix),-1,1)
;  scalecoef = ROBUST_POLY_FITQ(y[gddata],vstackstr.spec[i,gddata],5)
;  ;scalecoef = ROBUST_POLY_FIT(y,spec,5)
;  cont = POLY(y,scalecoef)

;  vstackstr.cont[i,*] = cont
;  vstackstr.specnorm[i,gddata] = vstackstr.spec[i,gddata]/cont[gddata]       ; normalize spectrum
;  vstackstr.errnorm[i,gddata] = vstackstr.err[i,gddata]/cont[gddata]         ; normalize error

  if not keyword_set(nolsffit) then begin
    ; Make sure the LSF is normalized
    ;  this 2D LSF is not normalized at *ALL*, not sure why
    ;  maybe because of the different pixel scale
    lsftot = total(lsf2d,2)
    bd = where(lsftot lt 0.01,nbd)
    if nbd gt 0 then lsftot[bd]=1
    lsf2d /= lsftot#replicate(1,nLSFpix)  ; make sure they are normalized

    ; Add to total weight LSF array
    mederr = MEDIAN(vstackstr.errnorm[i,*])              ; median err per spectrum, gaps are NAN
    totlsfwt += lsf2d * (1.0/mederr^2)
    ; Add to TOTWT which is just the running sum of the weights (1/err^2)
    ;  to this at the pixel level since each visit covers different pixels
    totwt += (lsf2d gt 1e-8)*(1.0/mederr^2)
  endif

  BOMB:

Endfor  ; visit loop
; IF all of the visits had a bad chip, then set that chip's pixlim to the error value.
allvisitbad = where(missingChipCount EQ nvisits, nallvisitbad)
IF nallvisitbad GT 0 THEN pixlim[*,allvisitbad] = make_array(2,value = -9999)

; All spectra bad, skip combination and proceed to end
visitmask = fix(TOTAL((vstackstr.mask AND badmask()) eq 0,2) gt 0)
visitmask2d = visitmask#replicate(1,npix)
ngoodvisits = total(visitmask) 
if ngoodvisits eq 0 then begin
  combspec = fltarr(npix)
  comberr = fltarr(npix)+baderr()
  combmask = intarr(npix)+1
  combsky = fltarr(npix)
  combskyerr = fltarr(npix)
  combtel = fltarr(npix)
  combtelerr = fltarr(npix)
  comblsf = fltarr(nlsfpix,npix)
  flsfpar = fltarr(4)
  goto,BOMBEND
endif


;; Get a better normalization by using the median spectrum
;;---------------------------------------------------------
;;  This doesn't do that much.  Mainly helps on the ends
;medspec0 = MEDIAN(vstackstr.specnorm,dim=1)   ; median spectrum
;; Try to make the median spectrum more normalized by dividing
;;  by the mode in the histogram
;f = histogram(medspec0,bin=0.01,min=0,max=2,locations=loc)
;maxind = first_el(maxloc(f))
;maxval = loc[maxind]
;medspec0 = medspec0/maxval
;; Set gaps to 1.0
;bd = where(finite(medspec0) eq 0,nbd)
;if nbd gt 0 then medspec0[bd]=1.0  ; set gaps to 1.0

;; Visit loop
;for i=0,nvisits-1 do begin
;  specfrac = vstackstr.spec[i,*]/medspec0  ; this should basically be the continuum
;  gddata = where(finite(specfrac) eq 1,ngddata)    ; only want good values
;  scalecoef = ROBUST_POLY_FITQ(y[gddata],specfrac[gddata],5)
;  fcont = POLY(y,scalecoef)
;  vstackstr.cont[i,*] = fcont
;  vstackstr.specnorm[i,gddata] = vstackstr.spec[i,gddata]/fcont[gddata]       ; normalize spectrum
;  vstackstr.errnorm[i,gddata] = vstackstr.err[i,gddata]/fcont[gddata]         ; normalize error
;endfor

; The average continuum
;  should we weight this by average snr?

;avgcont = TOTAL(vstackstr.cont,1)/nvisits
avgcont = TOTAL(vstackstr.cont,1,/nan)/(TOTAL(finite(vstackstr.cont),1)>1)

; don't include PERSIST visits in average, if there are non-PERSIST visits
npersist=0
igd=0
for i=0,nvisits-1 do begin
  if ( (visitstr[i].starflag and starflagval('PERSIST_HIGH')) gt 0 or $
       (visitstr[i].starflag and starflagval('PERSIST_MED')) gt 0 or $
       (visitstr[i].starflag and starflagval('PERSIST_LOW')) gt 0 ) then begin
     npersist+=1 
  endif else begin
     if igd eq 0 then gd=i else gd=[gd,i]
     igd+=1
  endelse
endfor
if npersist gt 0 and npersist lt nvisits then begin
  avgcont0=avgcont
  avgcont = TOTAL(vstackstr.cont[gd,*],1,/nan)/(TOTAL(finite(vstackstr.cont[gd,*]),1)>1)
endif

; Now combine the spectra and errors
;--------------------------------------
dataspec = vstackstr.specnorm[*,*]
dataerr = vstackstr.errnorm[*,*]
datamask = vstackstr.mask[*,*]
datasky = vstackstr.skynorm[*,*]

; Do outlier rejection
nbdpix=0
;medspec = MEDIAN(dataspec,dim=1,/even)   ; median spectrum
;if nvisits eq 1 then medspec=reform(dataspec)
;mnspec = TOTAL(dataspec,1,/nan)/(TOTAL(finite(dataspec),1) > 1)  ; divide by number of "good" pts
;diffspec_avg = dataspec - replicate(1.0,nvisits)#mnspec
;diffspec_med = dataspec - replicate(1.0,nvisits)#medspec
;medsig = MAD(diffspec_avg,dim=1)  ; the real sigma between them, need to use avg
;                            ; or there could be lots of zeros
;                            ; where the med is equal to the actual
;                            ; data point
;if nvisits eq 1 then medsig2d=mad(medspec) else $
;                     medsig2d=replicate(1,nvisits)#medsig

;  Use the median diff for outlier rejection because then there will
;  always be one data point that can't be rejected
;bdpix = where(abs(diffspec_med/medsig2d) gt 5 OR $
;              finite(diffspec_med) eq 0 OR $
;              datamask eq 1 ,nbdpix)  ; 5sigma outliers
bdpix= where(datamask and badmask(),nbdpix)
maskspec = dataspec
maskerr = dataerr
if nbdpix gt 0 then maskspec[bdpix] = !values.f_nan  ; set them to NAN
if nbdpix gt 0 then maskerr[bdpix] = !values.f_nan   ; set them to NAN
; Number of "good" points for each pixel
masknvisits = replicate(1L,nvisits,npix)
if nbdpix gt 0 then masknvisits[bdpix] = 0
nvisits2d = TOTAL(masknvisits,1) 

; Take a weighted mean of the spectrum
;  wt = 1/error^2
;  weighted mean = Sum(wt*X)/Sum(wt)

; Global weight per spectrum
;---------------------------
;if keyword_set(globalwt) then begin

  ; Create the error array using same values for all
  ;   pixels of a spectrum
  ;mederr = MEDIAN(dataerr,dim=2)              ; median err per spectrum, gaps are NAN
  ;meddataerr = mederr#replicate(1.0,npix)     ; make 2d array
  ; mask bad pixels
  ;if nbdpix gt 0 then meddataerr[bdpix] = !values.f_nan
  ; Now do a weighted mean (global) while ignoring the outliers (NANs)
  ;gcombspec = TOTAL(maskspec/meddataerr^2,1,/NAN) / TOTAL(1.0/meddataerr^2,1,/NAN)
  ;gcomberr = sqrt( TOTAL(maskerr^2,1,/NAN)/nvisits2d^2 )   ; ignore NANs

  ; use median filtered error unless given pixel has significantly deviant error
  newerr=zap(dataerr,[1,100])

  ; for "global" combination, deweight areas around sky lines
  bd=where(datamask and maskval('SIG_SKYLINE'),nbd)
  if nbd gt 0 then newerr[bd]*=100

  ; get enhanced error around sky lines
;  newflux=zap(dataspec,[1,100])
;  eerr=dataerr*0.
;  for ivisit=0,nvisits-1 do $
;   eerr[ivisit,*]=enhancederr(reform(newflux[ivisit,*]),reform(dataerr[ivisit,*]),reform(datasky[ivisit,*]),$
;        sig=3,skyfact=1000)
;  ; use median error + enhanced error
;  newerr+=(eerr-dataerr)
;  bigerr=where(dataerr/newerr gt 1.5,nbigerr)
;  if nbigerr gt 0 then newerr[bigerr]=dataerr[bigerr]

  meddataerr=newerr

  if nbdpix gt 0 then meddataerr[bdpix] = !values.f_nan

  gcombspec = TOTAL(maskspec/meddataerr^2,1,/NAN) / TOTAL(1.0/meddataerr^2,1,/NAN)
  gcomberr = sqrt( 1./TOTAL(1./meddataerr^2,1,/NAN) )   ; ignore NANs

; Pixel-by-pixel weighting
;-------------------------
;endif else begin

  persistpix = where(datamask and (maskval('PERSIST_HIGH') or maskval('PERSIST_MED') or maskval('PERSIST_LOW')) gt 0,np)
  phi = where((datamask and maskval('PERSIST_HIGH')) gt 0, np)
  if np gt 0 then maskerr[phi] *= sqrt(5.)
  pmed = where((datamask and maskval('PERSIST_HIGH')) eq 0 and $
               (datamask and maskval('PERSIST_MED')) gt 0, np)
  if np gt 0 then maskerr[pmed] *= sqrt(4.)
  plo = where((datamask and maskval('PERSIST_HIGH')) eq 0 and $
               (datamask and maskval('PERSIST_MED')) eq 0 and $
               (datamask and maskval('PERSIST_LOW')) gt 0, np)
  if np gt 0 then maskerr[plo] *= sqrt(3.)

  ; downweight sky lines
  psky = where((datamask and maskval('SIG_SKYLINE')) gt 0, np)
  if np gt 0 then maskerr[psky] *= sqrt(100.)

  ; Now do a weighted mean (pixel-by-pixel) while ignoring the outliers (NANs)
  pcombspec = TOTAL(maskspec/maskerr^2,1,/NAN) / TOTAl(1.0/maskerr^2,1,/NAN)
  pcombsky = TOTAL(datasky/maskerr^2,1,/NAN) / TOTAl(1.0/maskerr^2,1,/NAN)
  ; Combine the errors for mean and ignore outlier points
  ;pcomberr = sqrt( TOTAL(maskerr^2,1,/NAN)/nvisits2d^2 )   ; ignore NANs
  pcomberr = sqrt( 1./TOTAL(1./maskerr^2,1,/NAN) )   ; ignore NANs

;endelse

; Rescale spec/error using the average continuum
gcombspec *= avgcont
pcombspec *= avgcont
pcomberr *= avgcont
gcomberr *= avgcont
pcombsky *= avgcont

; for the mask, bitwise combine them, but only for frames that contributed data
combmask = fix(pcombspec*0)
combandmask = reform(vstackstr.mask[0,*])
for i=0,nvisits-1 do begin
  gd=where(finite(vstackstr.spec[i,*]) eq 1,ngd)
  if ngd gt 0 then begin
    combmask[gd]=(combmask[gd] OR vstackstr.mask[i,gd])
    combandmask[gd]=(combandmask[gd] AND vstackstr.mask[i,gd])
  endif
endfor

; For multiple visits, save 2 combined spectra (pixel-weighted and global 
;   weighted), plus individual visits, in separate rows
; For single visit, just save the one spectrum
if nvisits gt 1 then begin
  combspec=[[pcombspec],[gcombspec]]
  combspec=[[combspec],[transpose(vstackstr.spec)]]
  comberr=[[pcomberr],[gcomberr]]
  comberr=[[comberr],[transpose(vstackstr.err)]]
  ; combmask=TOTAL(vstackstr.mask,1)
  combmask = [[combmask],[combandmask],[transpose(vstackstr.mask)]]
  c = TOTAL(vstackstr.sky*visitmask2d,1)/ngoodvisits               ; sum sky
  ;combsky = [[c],[c],[transpose(vstackstr.sky)]]
  combsky = [[pcombsky],[c],[transpose(vstackstr.sky)]]
  if keyword_set(quick) then begin
    combskyerr = combsky*0.
    combtel = combskyerr
    combtelerr = combskyerr
  endif else begin
    c = sqrt(TOTAL((vstackstr.skyerr*visitmask2d)^2,1))/ngoodvisits         ; sum sky error
    ;combskyerr = [[c],[c],[transpose(vstackstr.skyerr)]]
    combskyerr = [[pcomberr],[c],[transpose(vstackstr.skyerr)]]
    c = TOTAL(vstackstr.telluric*visitmask2d,1)/ngoodvisits  ; average telluric
    combtel = [[c],[c],[transpose(vstackstr.telluric)]]
    c = sqrt(TOTAL((vstackstr.telerr*visitmask2d)^2,1))/ngoodvisits         ; sum telluric error
    combtelerr = [[c],[c],[transpose(vstackstr.telerr)]]
  endelse
endif else begin
  combspec = transpose(vstackstr.spec)
  comberr = transpose(vstackstr.err)
  combmask = transpose(vstackstr.mask)   
  combsky = transpose(vstackstr.sky)
  combskyerr = transpose(vstackstr.skyerr)
  combtel = transpose(vstackstr.telluric)
  combtelerr = transpose(vstackstr.telerr)
endelse

; Set pixels with no flux to ZERO
bdpix = where(finite(combspec) eq 0,nbdpix)
if nbdpix gt 0 then begin
  combspec[bdpix]=0.
  comberr[bdpix]=baderr()
  combmask[bdpix]=1
endif

pl=0
if keyword_set(pl) then begin
  smcolor
  !p.multi=[0,1,2]
  xs=15000 & xe=17000 & ys=0 & ye=2*median(combspec[*,0])

  ; Plot of spectra with telluric and sky
  plot,wave_final,combtel,color=3,xrange=[xs,xe]
  !p.multi=[0,1,2]
  plot,wave_final,combspec[*,0],yrange=[0,2*median(combspec[*,0])],xrange=[xs,xe],/noerase
  oplot,wave_final,combsky,color=3
  for i=2,n_elements(combspec[0,*])-1 do oplot,wave_final,combspec[*,i],color=(2+(i mod 6))
  oplot,wave_final,combspec[*,1],color=3
  oplot,wave_final,combspec[*,0]
 
  ; S/N plot of summed and individual spectra 
  !p.multi=[1,1,2]
  plot,wave_final,combspec[*,0]/comberr[*,0],yrange=[0,2*median(combspec[*,0]/comberr[*,0])]
  ;oplot,wave_final,comberr2,color=2
  for i=2,n_elements(comberr[0,*])-1 do oplot,wave_final,combspec[*,i]/comberr[*,i],color=(2+(i mod 6))
  oplot,wave_final,combspec[*,1]/comberr[*,0],color=3
  oplot,wave_final,combspec[*,0]/comberr[*,0]
  print,'S/N:'
  for i=0,n_elements(combspec[0,*])-1 do print,median(combspec[*,i]/comberr[*,i])
  stop
endif

;stop


; Make the final combined LSF array
;-----------------------------------
; divide LSF*WT array by the total weights array
;comblsf = totlsfwt / (totwt#replicate(1.0,nLSFpix) > 1.0)
;comblsf = totlsfwt / (totwt#replicate(1.0,nLSFpix) > 1e-8)
comblsf = totlsfwt / (totwt > 1e-8)

; Fit the combined LSF array
;----------------------------
; use MFITFUN and LSFGH2D
if not keyword_set(nolsffit) then begin

  ; Make the 2D X array
  dx = 1.0
  xlsf = REPLICATE(1.0d0,npix)#(dindgen(nLSFpix)-nLSFpix/2)*dx
  xlsf += dindgen(npix)#REPLICATE(1.0d0,nLSFpix)
  ;lsf2d = LSF_GH(xlsf,xcenter,lsfpars)

  ; Figure out the gap locations
  ;tempmask = total(comblsf,2) ge 0.1
  ;gapbeg = where(tempmask eq 0 and shift(tempmask,1) eq 1)
  ;gapend = where(tempmask eq 0 and shift(tempmask,-1) eq 1)
  lsfmask = total(comblsf,2) ge 0.1

  ;;fa = {xlsf:xlsf,gapbeg:gapbeg,gapend:gapend}
  ;fa = {xlsf:xlsf,mask:lsfmask#replicate(1,nLSFpix)}
  ;initpar = lsfpars
  ;initpar = str.lcoef[*,1]
  ;initpar[1] = -npix/2   ; center of the array
  ;horder = initpar[2]
  ;npar = n_elements(initpar)
  ;parinfo = REPLICATE({limited:[0,0],limits:[0.0,0.0],fixed:0},npar)
  ;;parinfo[0:horder+3].fixed = 1   ; # fixed parameters
  ;parinfo.fixed = 1  ; everything fixed to start with

  ; Breaking up the parameters
  ;binsize = lsfpars[0]
  binsize = 1.0     ; force this
  ;Xoffset = lsfpars[1]   ; Additive Xoffset
  Horder = lsfpars[2]
  Porder = lsfpars[3:Horder+3]   ; Horder+1 array
  nGHcoefs = total(Porder+1)
  ;GHcoef_ind = indgen(nGHcoefs)+Horder+4
  ;GHcoefs = lsfpars[GHcoef_ind]
  GHcoefs = lsfpars[Horder+4:Horder+4+nGHcoefs-1]
  GHcoef_ind = indgen(nGHcoefs)
  initpar = GHcoefs
  ;parinfo[GHcoef_ind].fixed = 0   ; allow the parameters to float
  ; Wing parameters
  if n_elements(lsfpars) gt (3+Horder+1+nGHcoefs) then begin
    wpar = lsfpars[3+Horder+1+nGHcoefs:*]
    ; Nwpar     number of W parameters
    ; WPorder   the polynomial order for each
    ; Wing coefficients
    wproftype = wpar[0]
    nWpar = wpar[1]
    wPorder = wpar[2:2+nWpar-1]
    nWcoefs = total(wPorder+1)
    Wcoefs = lsfpars[3+Horder+1+nGHcoefs+2+nWpar:*]
    Wcoef_ind = indgen(nWcoefs)+nGHcoefs
    ;Wcoef_ind = indgen(nWcoefs)+3+Horder+1+nGHcoefs+2+nWpar
    ;Wcoefs = lsfpars[Wcoef_ind]
    ;parinfo[Wcoef_ind].fixed = 0   ; allow the parameters to float

    ; Make both linear
    ;wPorder+=1
    ;Wcoefs = [Wcoefs[0],0.0,Wcoefs[1],0.0]
    ;nWcoefs = total(wPorder+1)
    ;Wcoef_ind = indgen(nWcoefs)+nGHcoefs

    initpar = [initpar,Wcoefs]
  endif
  npar = n_elements(initpar)
  parinfo = REPLICATE({limited:[0,0],limits:[0.0,0.0],fixed:0},npar)
  parinfo.fixed = 0
  ; set limits on Wing parameters
  ;parinfo[wcoef_ind[0]].limited = 1
  ;parinfo[wcoef_ind[0]].limits = [0,1]
  ;parinfo[wcoef_ind[1]].limited = 1
  ;parinfo[wcoef_ind[1]].limits = [0,10]
  ;parinfo[wcoef_ind].limited = [1,0]
  ;parinfo[wcoef_ind].limits = 0.01
  parinfo[wcoef_ind[[0,2]]].limited = [1,0]
  parinfo[wcoef_ind[[0,2]]].limits = 0.01

  ; Set Xoffset to middle of spectrum
  ;xoffset = -npix/2
  xoffset = npix/2

  ; Maybe expand the Wing order to linear??
  ;  right now they are constants and can't fit
  ;  the variation across all three chips

  ; The LSF parameter polynomial fits are going to be
  ; reversed since we flipped the X-axis.
  ; MPFIT will just have to fix it.

  ; Fit the parameters
  print,'Fitting the combined LSF parameters'
  x = dindgen(npix)#REPLICATE(1.0d0,nLSFpix)
  ;lcoef = MPFITFUN('fitcomblsf',x,comblsf,x*0+1,initpar,$
  ;                 parinfo=parinfo,yfit=yfit,dof=dof,status=status,bestnorm=chisq,$
  ;                 perror=perror,functargs=fa,niter=niter,maxiter=15,/quiet)
  fa = {lsfx:xlsf,mask:lsfmask#replicate(1,nLSFpix),x:x,y:comblsf,err:x*0+1,$
        binsize:binsize,offsetx:xoffset,Porder:Porder,wproftype:wproftype,WPorder:WPorder}
  ;t0 = systime(1)
  maxiter = 10 ;50
  lcoef = MPFIT('fitcomblsf_dev',initpar,$
                parinfo=parinfo,dof=dof,status=status,bestnorm=chisq,$
                perror=perror,functargs=fa,niter=niter,maxiter=maxiter,/quiet,autoderivative=0)
  ;fpar = [binsize, xoffset, Horder, Porder, lcoef[0:nGHcoefs-1], wproftype, nWpar, WPorder, lcoef[nGHcoefs:*]]
  dum = fitcomblsf_dev(lcoef,_extra=fa,inpar=flsfpar,model=model)
  yfit = comblsf*0.0
  yfit[*,*] = model   ; put in 2D form
  ;print,systime(1)-t0,chisq

  ;stop

endif else flsfpar=0

BOMBEND:

; Make sure to set ERROR high for bad pixels
bdpix = where(combmask eq 1,nbdpix)
if nbdpix gt 0 then comberr[bdpix] = baderr()

;---------------------------------------------------
; Now put everything in the final STARSTR structure
; Put the wavelenth coefficients in the header
; Put the LSF coefficients in the header
;---------------------------------------------------
one = fltarr(npix)
starstr = {nvisits:nvisits,header:strarr(nvisits,500),$
           spec:combspec,$
           wave:wave_final,$
           err:comberr,$
           mask:combmask,$
           sky:combsky,$
           skyerr:combskyerr,$
           telluric:combtel,$
           telerr:combtelerr,$
           comblsf:comblsf,$
           wcoef:reform(p2w_coef),$
           pixlim:pixlim,$
           pixlim_overlap:pixlim_overlap,$
           lcoef:flsfpar, combtype:combtype}

; Put the visit headers into STARSTR
for i=0,nvisits-1 do begin
  head = (*allvisits[i]).head

  ; Add weighting method (currently including both, so comment these for now!)
;  if keyword_set(globalwt) then wtmethod = 'Using global spectrum weighting' else $ ; weight method
;    wtmethod = 'Using pixel-by-pixel weighting'
;  sxaddhist,'APVISITCOMB: '+wtmethod,head

  nhead = n_elements(head)
  starstr.header[i,0:nhead-1] = head
end

if keyword_set(stp) then stop

end
