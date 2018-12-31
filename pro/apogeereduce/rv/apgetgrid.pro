pro apgetgrid,synthfile,grid=grid,wave=wave,error=error,stp=stp,pl=pl,normalize=normalize

;+
;
; APGETGRID
;
; This loads a synthetic spectral grid. Grid is continuum normalized,
;  and resampled in wavelength if requested
;
; Input
;    wave=wave : if specified, return grid at input (vacuum) wavelengths
;                Otherwise return at native grid wavelengths (log lambda
;                in air). In either case, load the grid wavelengths
;                in vacuum into grid.wave
; Output:
;    grid  : grid structure with normalized spectra and wavelength arrays
;            giving reference wavelengths in vacuum
;
; USAGE:
;  IDL>apgetgrid,grid=grid,wave=wave
;
; By D.Nidever  July 2010
;-

; Not enough inputs
if n_elements(synthfile) eq 0 then begin
  print,'Syntax - apgetgrid,synthfile,grid=grid'
  return
endif

; Get APOGEE directories
dirs=getdir(apogee_dir,cal_dir,spectro_dir,apogee_vers,lib_dir)
synth_dir = lib_dir+'synthgrid/'
if FILE_TEST(synth_dir,/directory) eq 0 then begin
  print,'SYNTHGRID Directory ',synth_dir,' NOT FOUND'
  return
endif

; read synthetic grid, which should be in vacuum waveengths and continuum normalized
; HDU1 has parameter table, HDU2 spectra on apStar grid, HDU3 spectra on higher sampling for visit spectra
; data arrays have wavelength in axis 0 
if n_elements(grid) eq 0 then begin
  head=mrdfits(synth_dir+synthfile,0)
  stepar=mrdfits(synth_dir+synthfile,1)
  if keyword_set(apstar) then data=mrdfits(synth_dir+synthfile,2,headspec) else $
                              data=mrdfits(synth_dir+synthfile,3,headspec) 
  logwave=sxpar(headspec,'CRVAL1')+indgen(sxpar(headspec,'NAXIS1'))*sxpar(headspec,'CDELT1')
  logdw=sxpar(headspec,'CDELT1')
  outwave=10.**logwave
  grid = {file:synth_dir+synthfile,
              data:transpose(data),ndata:transpose(ndata),$
              head:head,metals:stepar.mh,teff:stepar.teff,logg:stepar.logg,$
              wave:outwave,logwave:logwave,logdw:dw,res:res}
endif
return

; =======================================
;  EVERYTHING BELOW HERE FOR OLD RV GRIDS
; =======================================

; Load the grid file if not input
if n_elements(grid) eq 0 then begin

  ; Read in the new RV synthetic spectral grid
  FITS_READ,synth_dir+synthfile,specdata,head1,/no_abort,message=message1,exten=1
  FITS_READ,synth_dir+synthfile,metals,head2,/no_abort,message=message2,exten=2
  FITS_READ,synth_dir+synthfile,teff,head3,/no_abort,message=message3,exten=3
  FITS_READ,synth_dir+synthfile,logg,head4,/no_abort,message=message4,exten=4

  ; New grid has stellar parameters structre in HDU2
  if message1 eq '' and message2 eq '' and message3 ne '' and message4 ne '' then begin
    stepar = mrdfits(synth_dir+synthfile,2,/silent)
    message3 = ''
    message4 = ''
  endif

  if message1+message2+message3+message4 ne '' then begin
    error = 'Problem opening the synthetic spectral grid '+synth_dir+synthfile
    error += message1+message2+message3+message4
    print,error
    return
  endif

  ; load wavelengths which are on log-lambda scale IN AIR WAVELENGTHS
  wtype = sxpar(head1,'CTYPE2',comment=wcom)
  w0 = sxpar(head1,'CRVAL2')   
  dw = sxpar(head1,'CDELT2')
  npix = sxpar(head1,'NAXIS2')
  res = sxpar(head1,'RESOL')
  loggridwave = w0+dindgen(npix)*dw   ; log wave
  gridwave = 10^loggridwave

  ; Convert from AIR to VACUUM wavelengths if necessary
  if stregex(wcom,' air ',/fold_case,/boolean) then AIRTOVAC,gridwave,gridwavevac $
    else gridwavevac=gridwave

  sz = size(specdata)
  nspec=sz[1]

  ; set output values
  if keyword_set(wave) then begin
    data=dblarr(nspec,n_elements(wave))
    ndata=dblarr(nspec,n_elements(wave))
    outwave=wave
    logwave=alog10(outwave)
  endif else begin
    data=dblarr(nspec,npix)
    ndata=dblarr(nspec,npix)
    outwave=gridwavevac
    logwave=alog10(gridwavevac)
  endelse
  dw = median(slope(logwave))

  ; initialize grid structure
  if n_elements(stepar) gt 0 then begin
    grid = {file:synth_dir+synthfile,origdata:specdata,$
            data:data,ndata:ndata,head:head1}
    tags = tag_names(stepar)
    tagstoadd = tags
    bd = where(tags eq 'PIXLIM' or tags eq 'LO' or tags eq 'HI',nbd)
    if nbd gt 0 then REMOVE,bd,tagstoadd
    for i=0,n_elements(tagstoadd)-1 do begin
      ind = where(tags eq tagstoadd[i],nind)
      grid = create_struct(grid,tagstoadd[i],stepar.(ind[0]))
    endfor
    grid = create_struct(grid,'WAVE',outwave,'LOGWAVE',logwave,'LOGDW',dw,'RES',res)
    grid = create_struct(grid,'PIXLIM',lonarr(2,3))
    if tag_exist(stepar,'PIXLIM') then grid.pixlim=stepar.pixlim
    if tag_exist(stepar,'LO') and tag_exist(stepar,'HI') then grid.pixlim=transpose([[stepar.lo],[stepar.hi]])
  endif else begin
    grid = {file:synth_dir+synthfile,origdata:specdata,$
            data:data,ndata:ndata,$
            head:head1,metals:metals,teff:teff,logg:logg,$
            wave:outwave,logwave:logwave,logdw:dw,res:res}
  endelse


  ; Do we need to resample onto input wavelength scale
  ;  resamp=0  no resampling
  ;  resamp=1  resampling
  ;  resamp=2  subset
  resamp = 0
  if keyword_set(wave) then begin
    resamp = 1
    if n_elements(gridwavevac) eq n_elements(wave) then if max(abs(gridwavevac-wave)) lt 1e-10 then resamp=0
    ; Check if input wavelength array is a SUBSET of the grid wavelengths
    if resamp eq 1 then if min(abs(gridwavevac-wave[0])) lt 1e-10 then begin
      s1 = strtrim(string(gridwavevac*1e10,format='(I30)'),2)
      s2 = strtrim(string(wave*1e10,format='(I30)'),2)
      match,s1,s2,subind,ind2,count=nmatch
      if nmatch eq n_elements(wave) then resamp=2   ; subset
    endif
  endif


  ; Continuum Normalize each spectrum and resample if requested
  for i=0,nspec-1 do begin
    temp = {spec:reform(grid.origdata[i,*]),wave:gridwavevac}  ; temporary structure
    if keyword_set(normalize) then begin
      APNORMSPEC,temp
      nspec = temp.nspec
    endif else nspec=reform(grid.origdata[i,*])
    if resamp eq 2 then nspec=nspec[subind]  ; subset
    if resamp eq 1 then begin
      grid.data[i,*]=spline(gridwavevac,double(nspec),wave,/double)
      grid.ndata[i,*]=grid.data[i,*]
    endif else begin
      grid.data[i,*] = nspec
      grid.ndata[i,*] = nspec
    endelse
  endfor

endif

if keyword_set(stp) then stop

end
