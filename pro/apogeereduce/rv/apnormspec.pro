pro apnormspec,str,ncorder=ncorder,fixbadpix=fixbadpix,noerrcorr=noerrcorr,sepchip=sepchip,$
               normtype=normtype,binsize=binsize0,perclevel=perclevel0,growsky=growsky,$
               nsky=nsky,error=error,stp=stp,pl=pl,silent=silent

;+
;
; APNORMSPEC
;
; This program normalizes an APOGEE spectrum
;
; INPUTS:
;  str         An APOGEE spectrum structure.  This at least needs
;                to have a SPEC or FLUX tag and a WAVE tag.
;  /sepchip    Do the normalization chip-by-chip.
;  =normtype   The normalization type to use:
;                1-Polynomial fitting in Nth percentile binned
;                    spectrum (default).
;                2-Gaussian Smooothing
;                3-Fourier Filtering
;  =ncorder    The continuum polynomial order.  The default is 6.
;  /noerrcorr  Do not use a correction for the effects of the errors
;                on the continuum measurement.  The default is to make
;                this correction if errors are included.
;  /fixbadpix  Set bad pixels to the continuum
;  =binsize    The binsize to use (in units of 900A) for determining
;                the Nth percentile spectrum to fit with a polynomial.
;                The default is 0.05.
;  =perclevel  The Nth percentile to use to determine the continuum.
;                The default is 0.95 or 95%.
;  =growsky    Number of pixels by which to grow the sky pixel masking.
;  =nsky       The sky threshold to use for sky pixels to mask in
;                units of the median sky flux.
;  /pl         Plot the spectrum and the continuum fit.
;  /silent     No output to the screen.
;  /stp        Stop at the end of the program.
;
; OUTPUTS:
;  The STR is updated with the normalized spectrum (NSPEC) and the
;  continuum (CONTINUUM).  If a chip normalization failed then NSPEC=NAN
;  and CONTINUUM=1
;  =error     The error message if one occurred.  If /sepchip is set
;               then an error will be returned ONLY if there were
;               errors in all chip normalizations.
;
; USAGE:
;  IDL>apnormspec,str,outstr
;
; By D.Nidever  July 2011
;-

apgundef,error

; Not enough inputs
if n_elements(str) eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax - apnormspec,str,fixbadpix=fixbadpix,ncorder=ncorder,noerrcorc=noerrcorr,normtype=normtype,'
  print,'                    binsize=binsize,perclevel=perclevel,growsky=growksy,error=error,sepchip=sepchip,stp=stp,silent=silent'
  return
endif

; Defaults
if n_elements(ncorder) gt 0 then ncontorder=ncorder[0] else ncontorder = 6          ; polynomial order
if n_elements(binsize0) gt 0 then binsize=binsize0 > 0.01 else binsize=0.05         ; binsize
if n_elements(perclevel0) gt 0 then perclevel=perclevel > 0.01 else perclevel=0.95  ; percent level
if n_elements(normtype) eq 0 then normtype=1                                        ; normalization method
if n_elements(nsky) eq 0 then nsky=5                                                ; sky threshold

tags = tag_names(str)
; Check for SPEC or FLUX tag
specind = where(tags eq 'SPEC',nspecind)
if nspecind eq 0 then specind = where(tags eq 'FLUX',nspecind)
if nspecind eq 0 then begin
  error = 'No SPEC or FLUX tag in structure'
  print,error
  return
endif

; Check for WAVE tag
waveind = where(tags eq 'WAVE',nwaveind)
if nwaveind eq 0 then begin
  error = 'No WAVE tag in structure'
  print,error
  return
endif

; Check for ERR tag
errind = where(tags eq 'ERR',nerrind)
maskind = where(tags eq 'MASK',nmaskind)
skyind = where(tags eq 'SKY',nskyind)
sz = size(str.(specind))



; Do the normalization chip-by-chip
;----------------------------------
if keyword_set(sepchip) and sz[0] eq 2 then begin

  ; Initialize the output tags
  dum = where(tags eq 'NSPEC',num_nspec)
  if num_nspec eq 0 then ADD_TAG,str,'NSPEC',str.(specind[0])*0,str
  dum = where(tags eq 'CONTINUUM',num_continuum)
  if num_continuum eq 0 then ADD_TAG,str,'CONTINUUM',str.(specind[0])*0,str

  ; Use 3rd order for chip normalization
  if n_elements(ncorder) eq 0 then ncorder1=3 else ncorder1=ncorder

  ; Loop through chips
  errorarr = strarr(3)
  for i=0,2 do begin
    ; Make temporary chip structure
    chipstr = {spec:str.(specind)[*,i],err:str.(errind)[*,i],wave:str.(waveind)[*,i]}
    if nmaskind gt 0 then chipstr=CREATE_STRUCT(chipstr,'MASK',str.(maskind)[*,i])
    if nskyind gt 0 then chipstr=CREATE_STRUCt(chipstr,'SKY',str.(skyind)[*,i])

    ; Normalize
    APNORMSPEC,chipstr,ncorder=ncorder1,fixbadpix=fixbadpix,noerrcorr=noerrcorr,$
               binsize=binsize0,perclevel=perclevel0,error=error,stp=stp,pl=pl,normtype=normtype

    ; Now put the normalized spectra in
    if n_elements(error) eq 0 then begin
      str.nspec[*,i] = chipstr.nspec
      str.continuum[*,i] = chipstr.continuum
    endif else begin
      errorarr[i] = error
      str.nspec[*,i] = !values.f_nan
      str.continuum[*,i] = 1.0
    endelse

  endfor ; chip loop

  ; Only return an error if all three chips had problems
  apgundef,error
  chiptag = ['a','b','c']
  bdchip = where(errorarr ne '',nbdchip)
  if nbdchip eq 3 then error = strjoin('CHIP '+chiptag+' '+errorarr,' ')

  return

endif


; Initialize the output structure
dum = where(tags eq 'NSPEC',num_nspec)
if num_nspec eq 0 then ADD_TAG,str,'NSPEC',str.(specind[0])*0,str
dum = where(tags eq 'CONTINUUM',num_continuum)
if num_continuum eq 0 then ADD_TAG,str,'CONTINUUM',str.(specind[0])*0,str
;if keyword_set(fixbadpix) then if tag_exist(str,'MASKED') eq 0 then add_tag,str,'MASKED',str.mask*0,str
if tag_exist(str,'MASKED') eq 0 then add_tag,str,'MASKED',str.spec*0,str
str.nspec = !values.f_nan
str.continuum = 1


; All pixels are bad
if nmaskind gt 0 then begin
  gd = where( (str.mask and badmask()) eq 0 and finite(str.spec) eq 1,ngd,comp=bd,ncomp=nbd)
endif else begin
  gd = where(finite(str.spec) eq 1,ngd,comp=bd,ncomp=nbd)
endelse 
if ngd eq 0 then begin
  str.masked = 1
  return
endif


; Normalization Methods
;----------------------
Case normtype of

  ;----------------------------------------------------------------
  ; Polynomial fitting to Nth binned spectrum with error correction
  ;----------------------------------------------------------------
  1: begin

    ; Continuum Normalize
    ;----------------------
    w = (str.(waveind[0]))(*)
    x = (w-median(w))/range(w*0.5)  ; -1 to +1
    ;x = (w-16000.0)/900.0
    y = (str.(specind[0]))(*)
    if nerrind gt 0 then yerr = (str.(errind[0]))(*)


    ; Get good pixels, and set bad pixels to NAN
    ;--------------------------------------------
    gdmask = (y gt 0.0)  ; need positive fluxes
    ytemp = y

    ; Exclude pixels with mask=bad
    if nmaskind gt 0 then begin
      mask = (str.(maskind[0]))(*)
      gdmask = gdmask AND ((mask and badmask()) eq 0)

      gdmask = gdmask AND ((mask and maskval('SIG_SKYLINE')) eq 0)
      nskyind = 0
    endif
    ; Exclude pixels with bright airglow lines
    if nskyind gt 0 then begin
      sky = (str.(skyind[0]))(*)
      ; Calculate the median-filtered sky spectrum
      sz = size(str.sky)
      if sz[0] eq 2 then nsum=sz[2] else nsum=1
      skymask = long(str.sky*0)
      for i=0,nsum-1 do begin  ; loop ofver 2nd dim if necessary
        sky = str.sky[*,i]
        medsky = MEDFILT1D(str.sky[*,i],201,/edge)
        medcoef = ap_robust_poly_fit(x,medsky/median(medsky),2)
        medsky2 = poly(x,medcoef)*median(medsky)
        skymask1 = (sky gt nsky*medsky2)  ; pixels Nsig above median sky

        ; THIS WORKS PRETTY WELL, BUT I'M WORRIED
        ; THAT IT'S MASKING OUT TOO MANY PIXELS

        ; Find peaks in sky spectrum
        ;skyslp = slope(sky)
        ;skysig1 = mad(sky-medsky2)
        ;cskyind = where(abs(sky-medsky2) lt 5*skysig1,ncskyind)
        ;skysig = mad(sky[cskyind]-medsky2[cskyind])
        ;skymedlocal = medfilt1d(sky,31,/edge)
        ;skypeak = (([-1,skyslp] gt 0.0 and [1,skyslp[1:*],1] le 0.0) and (sky-medsky2) gt 3*skysig)
        ;grwskypeak = convol(skypeak,fltarr(31)+1)  ; grow
        ;grwskypeak /= (grwskypeak>1)
        ;grwskypeak = (grwskypeak and (sky-medsky2) gt 2*skysig)
        ;;pkind = where( ([-1,skyslp] gt 0.0 and [1,skyslp[1:*],1] le 0.0) and sky-skymedlocal gt 3*skysig,nmaxind)
        ;skymask1 = (skymask1 OR grwskypeak)
        ; grow the sky mask
        if keyword_set(growsky) then begin
          skymask1 = convol(skymask1,lonarr(ceil(growsky)>2)+1,/center)
          skymask1 = skymask1/(skymask1>1)
        endif

        skymask[*,i] = skymask1
      endfor ; nsum

      ;gdmask = gdmask AND (sky lt nsky*medsky)
      gdmask = gdmask AND (skymask eq 0)
    endif

    gdpix = where(gdmask eq 1,ngdpix,comp=bdpix,ncomp=nbdpix)
    if nbdpix gt 0 then ytemp[bdpix]=!values.f_nan  ; set bad pixels to NAN for now

    ; First attempt at continuum
    ;----------------------------
    ; Bin the data points
    BINDATA,x,ytemp,xbin,ybin,binsize=binsize,perc=perclevel
    gdbin = where(finite(ybin) eq 1,ngdbin)
    if ngdbin lt ncontorder+1 then begin
      error = 'not halted: Not enough good flux points to fit the continuum'
      if not keyword_set(silent) then print,error
      return
    endif
    ; Fit with robust polynomial
    coef1 = ROBUST_POLY_FITQ(xbin[gdbin],ybin[gdbin],ncontorder,/silent)
    cont1 = POLY(x,coef1)

    ; Subtract smoothed error from it to remove the effects
    ;  of noise on the continuum measurement
    if nerrind gt 0 and not keyword_set(noerrcorr) then begin
      smyerr = MEDFILT1D(yerr,151,/edge)                 ; first median filter
      smyerr = GSMOOTH(smyerr,100)                       ; Gaussian smoothing
      coef_err = ROBUST_POLY_FITQ(x,smyerr,ncontorder,/silent)   ; fit with robust poly
      ;poly_err = poly(x,coef_err)
      ;cont1 -= 2*poly_err   ; is this right????
      med_yerr = median(smyerr)                          ; median error
      cont1 -= 2*med_yerr
    endif

    ; Second iteration
    ;-----------------
    ;  This helps remove some residual structure
    ytemp2 = ytemp/cont1
    BINDATA,x,ytemp2,xbin,ybin2,binsize=binsize,perc=perclevel,gdind=gdind2
    gdbin2 = where(finite(ybin2) eq 1,ngdbin2)
    if ngdbin2 lt ncontorder+1 then begin
      error = 'not halted: Not enough good flux points to fit the continuum'
      if not keyword_set(silent) then print,error
      return
    endif
    ; Fit with robust polynomial
    coef2 = ROBUST_POLY_FITQ(xbin[gdind2],ybin2[gdind2],ncontorder,/silent)
    cont2 = POLY(x,coef2)

    ; Subtract smoothed error again
    if nerrind gt 0 and not keyword_set(noerrcorr) then begin
      ;cont2 -= 2*poly_err/cont1  ; Is this right??
      cont2 -= med_yerr/cont1
    endif

    ; Final continuum
    cont = cont1*cont2  ; final continuum

  end ; normtype=1


  ;-------------------
  ; Gaussian Smoothing
  ;-------------------
  2: begin

    print,'Gaussian Smoothing NOT supported yet'
    return

    sz = size(str.(specind))
    ; Loop over chips
    nloop = (sz[0] eq 2) ? sz[2] : 1
    for i=0,nloop-1 do begin

      gdmask = (y gt 0.0)  ; need positive fluxes
      ytemp = y

      ; Exclude pixels with mask=bad
      if nmaskind gt 0 then begin
        mask = (str.(maskind[0]))(*)
        gdmask = gdmask AND ((mask and badmask()) eq 0)
      endif
      ; Exclude pixels with bright airglow lines
      if nskyind gt 0 then begin
        sky = (str.(skyind[0]))(*)
        ; Calculate the median-filtered sky spectrum
        sz = size(str.sky)
        if sz[0] eq 2 then begin
          medsky = str.sky*0
          for i=0,sz[2]-1 do medsky[*,i]=MEDFILT1D(str.sky[*,i],201,/edge)
          medsky = (medsky)(*)
        endif else medsky=MEDFILT1D(str.sky,201,/edge)
        skymask = (sky gt nsky*medsky)
        ; grow the sky mask
        if keyword_set(growsky) then begin
          skymask = convol(skymask,lonarr(ceil(growsky)>2)+1,/center)
          skymask = skymask/(skymask>1)
        endif
        ;gdmask = gdmask AND (sky lt nsky*medsky)
        gdmask = gdmask AND (skymask eq 0)
      endif

      gdpix = where(gdmask eq 1,ngdpix,comp=bdpix,ncomp=nbdpix)
      if nbdpix gt 0 then ytemp[bdpix]=!values.f_nan  ; set bad pixels to NAN for now
      
    endfor

  end ; normtype=2, Gaussian smoothing


  ;------------------
  ; Fourier Filtering
  ;------------------
  3: begin

    print,'Fourier Filtering NOT supported yet'
    return

  end ; normtype=3, Fourier Filtering


  ; Not supported
  else: begin
    error = 'NORMTYPE='+strtrim(normtype,2)+' NOT SUPPORTED'
    if not keyword_set(silent) then print,error
    return
  end

ENDCASE


; "Fix" bad pixels
if nbdpix gt 0 and keyword_set(fixbadpix) then y[bdpix]=cont[bdpix]


; Stuff in the output
str.nspec = y/cont
str.continuum = cont

; Add "masked" array
if keyword_set(fixbadpix) then begin
  if tag_exist(str,'MASKED') eq 0 then add_tag,str,'MASKED',str.mask*0,str
  if nbdpix gt 0 then str.masked[bdpix]=1
endif

; Plot the spectrum
if keyword_set(pl) then begin
  plot,str.(waveind)/1e4,str.(specind),/xsty,/ysty,xtit='Wavelength (um)',ytit='Spectrum',charsize=1.2
  oplot,str.(waveind)/1e4,str.continuum,co=250
  if nbdpix gt 0 then oplot,w[bdpix],y[bdpix],ps=1,co=200
  legend,['Spectrum','Continuum','Rejected Pixels'],textcolor=[255,250,200],/top,/left,charsize=1.2
endif

;stop

if keyword_set(stp) then stop

end
