function aptelluric_convolve,frame,fiber=fiber,clobber=clobber,convonly=convonly,chip=chip

;+
;
; APTELLURIC_CONVOLVE
;
; This function returns a set of convolved telluric spectra for a given
;    input frame, where the information in the input specifies the LSF and
;    the wavelength scale for the frame
; If a convolved telluric has already been made for this LSF, then it is 
;    read and sampled on the desired wavelength scale
; If a convolved telluric has not been made, then it makes it first
;
; INPUTS:
;  frame       A structure with the header/data information for an
;                  undersampled frame that has been wavelength
;                  calibrated and airglow line subtracted (not essential)

; OUTPUTS:
;     Function returns a chip-concatenated array of dimension 
;           convolved_telluric[3*npix,nfiber,3], where the
;       third dimension is for the three species
;
; USAGE:
;  IDL>conv=aptelluric_convolve(frame)
;
; By J. Holtzman  March 2012
;   based off of some routines by D. Nidever, but with modified convolution
;-
 
; Read in the Telluric model spectra into telluric structure
dirs=getdir(apogee_dir,cal_dir,spectro_dir,apogee_vers,lib_dir)
telluric_dir = lib_dir+'telluric/'
if FILE_TEST(telluric_dir,/directory) eq 0 then begin
  print,'TELLURIC Directory ',telluric_dir,' NOT FOUND'
  stop
  return,0
endif

; Convolved telluric file name
tmp=strsplit(file_basename(frame.chipa.lsffile),'-',/ext)
tmp1=strsplit(tmp[2],'.',/ext)
lsfid=tmp1[0]

tmp=strsplit(file_basename(frame.chipa.wavefile),'-',/ext)
tmp1=strsplit(tmp[2],'.',/ext)
waveid=tmp1[0]

if keyword_set(fiber) then fibers=fiber else fibers=indgen(300)

chips=['a','b','c']
outfile=cal_dir+'telluric/'+dirs.prefix+'Telluric-'+chips+'-'+waveid+'-'+lsfid+'.fits'

; wait if another process is already working on this frame
lockfile=cal_dir+'telluric/'+dirs.prefix+'Telluric-'+waveid+'-'+lsfid+'.lock'
while file_test(lockfile) do apwait,lockfile,10

; does convolved telluric file already exist? If not, make it!
if keyword_set(clobber) or (not file_test(outfile[2])) then begin

  openw,lock,/get_lun,lockfile
  free_lun,lock

  ; we'll make the convolved spectrum on a subsampling of the wavecal grid
  for k=0,2 do begin
    FITS_READ,frame.(k).wavefile,wtmp,whead,exten=1
    if k eq 0 then wcoef=wtmp else wcoef=[[[wcoef]],[[wtmp]]]
  endfor

 ; get the input telluric spectra
  nspecies=3
  FITS_READ,telluric_dir+'CH4.fits',telim1,telhead1,message=message1,/no_abort
  FITS_READ,telluric_dir+'CO2.fits',telim2,telhead2,message=message2,/no_abort
  FITS_READ,telluric_dir+'H2O.fits',telim3,telhead3,message=message3,/no_abort
  if message1+message2+message3 ne '' then begin
    print,'ERROR loading the Telluric spectra'
    stop
    return,0
  endif
  sz=size(telim1)
  if sz[0] eq 3 then begin 
    twave1=sxpar(telhead1,'CRVAL1')+indgen(sxpar(telhead1,'NAXIS1'))*sxpar(telhead1,'CDELT1')
    twave2=sxpar(telhead2,'CRVAL1')+indgen(sxpar(telhead2,'NAXIS1'))*sxpar(telhead2,'CDELT1')
    twave3=sxpar(telhead3,'CRVAL1')+indgen(sxpar(telhead3,'NAXIS1'))*sxpar(telhead3,'CDELT1')
    si1 = sort(twave1)  ; sort with wavelength
    si2 = sort(twave2)  ; sort with wavelength
    si3 = sort(twave3)  ; sort with wavelength
    nair=sz[2]
    nscale=sz[3]
    tspec1_orig = telim1
    tspec2_orig = telim2
    tspec3_orig = telim3
  endif else begin
    twave1 = reform(telim1[2,*])
    tspec1_orig = reform(telim1[1,*],n_elements(telim1[1,*]),1,1)
    si1 = sort(twave1)  ; sort with wavelength
    twave2 = reform(telim2[2,*])
    tspec2_orig = reform(telim2[1,*],n_elements(telim2[1,*]),1,1)
    si2 = sort(twave2)  ; sort with wavelength
    twave3 = reform(telim3[2,*])
    tspec3_orig = reform(telim3[1,*],n_elements(telim3[1,*]),1,1)
    si3 = sort(twave3)  ; sort with wavelength
    nair=1
    nscale=1
  endelse

  ; loop over chips
  if n_elements(chip) gt 0 then begin
    j1=chip
    j2=chip
  endif else begin
    j1=0
    j2=2
  endelse
  for j=j1,j2 do begin

    ; wait if another process is already working on this frame
    while file_test(outfile[j]+'.lock') do apwait,outfile[j]+'.lock',10

    ; Does the output file already exist?
    if not keyword_set(nowrite) and file_test(outfile[j]) eq 1 and not keyword_set(clobber) then begin
      error = outfile[j]+' already exists and CLOBBER=0'
      if not keyword_set(silent) then print,error
      goto,nextchip
    endif

    ; open lock file to prevent multiple processes from working on same frame
    openw,lock,/get_lun,outfile[j]+'.lock'
    free_lun,lock
 
    ; get frame size 
    sz = size(frame.chipa.flux)
    npix = sz[1]
    pix = findgen(npix)
    nfibers = sz[2]
 
    ; Calculate at the same wavelengths for all fibers, allowing for some extra on both ends
    dx = 0.2d   ; 0.1 is close to the native model wavelength sampling
    dx = 0.5d   ; 0.1 is close to the native model wavelength sampling
    osamp = fix(1/dx)
    extend = 20  ; extra on ends
    nx = osamp*(npix+2*extend)
    out = dblarr(nx,300,3,nscale)
    twave = dblarr(nx,300)
    xfine = dindgen(nx)*dx-extend
    wave1 = PIX2WAVE(xfine,reform(wcoef[150,*,j]))
    wsi = sort(wave1)

    nLSFpix = 2*7*osamp + 1   ; +/-7 underampled pixels
    xlsf = REPLICATE(1.0d0,nx)#(dindgen(nLSFpix)-nLSFpix/2)*dx
    xlsf += (dindgen(nx)*dx-extend)#REPLICATE(1.0d0,nLSFpix)

    ; Loop over airmass and scale factors if we have them
    for iair=0,nair-1 do begin
    out*=0.
    for iscale=0,nscale-1 do begin

    ; interpolate telluric spectra to output wavelengths
    tspec1 = dblarr(nx)
    tspec1[wsi] = SPLINE(twave1[si1],tspec1_orig[si1,iair,iscale],wave1[wsi],/double)
    tspec2 = dblarr(nx)
    tspec2[wsi] = SPLINE(twave2[si2],tspec2_orig[si2,iair,iscale],wave1[wsi],/double)
    tspec3 = dblarr(nx)
    tspec3[wsi] = SPLINE(twave3[si3],tspec3_orig[si3,iair,iscale],wave1[wsi],/double)
    ; Loop through the fibers
    for i=0,nfibers-1 do begin
      ;i=fibers(ii)
 
      print,'Convolving fiber index: ', i , iscale, iair, j
      ; Convolve the model spectra with this fiber's LSF
      ;-------------------------------------------------
      x = dblarr(3*npix)
      twave[*,i]=wave1
  
      ; get the LSF for this fiber
      lsfpars = reform( frame.(j).lsfcoef[i,*] )
      lsf2d = LSF_GH(xlsf,xfine,lsfpars)
      lsftot = total(lsf2d,2)
      lsf2d /= lsftot#replicate(1,nLSFpix)  ; make sure the LSF is normalizedS

      ; do the convolution
      nhalf=nLSFpix/2
      for k=nhalf,nx-1-nhalf do begin
        out[k-nhalf:k+nhalf,i,0,iscale] += lsf2d[k,*]*tspec1[k]
        out[k-nhalf:k+nhalf,i,1,iscale] += lsf2d[k,*]*tspec2[k]
        out[k-nhalf:k+nhalf,i,2,iscale] += lsf2d[k,*]*tspec3[k]
      endfor
;
;      ; instead of using a very high sampling we could rebin
;      ; to lower sampling (i.e. summing/binning) and then use
;      ; less LSF pixels.
;      ; but we need the integer pixels in there.
;  
;      ; Step 2-Make 2D arrays of the model spectra
;      nLSFpix = 2*7*osamp + 1   ; +/-7 underampled pixels
;      ;ind0 = REPLICATE(1,npix)#(lindgen(nLSFpix)-nLSFpix/2)
;      ;indcen = (lindgen(npix)*osamp + osamp*extend )#replicate(1,nLSFpix)
;      ind0 = REPLICATE(1,nx)#(lindgen(nLSFpix)-nLSFpix/2)
;      ;indcen = (lindgen(nx) + extend )#replicate(1,nLSFpix)
;      ;ind = ind0+indcen
;      ind = ind0
;      tspec1_2d = tspec1[ind]
;      tspec2_2d = tspec2[ind]
;      tspec3_2d = tspec3[ind]
;    
;      ; Step 3-Make 2D LSF array
;      ;xlsf = REPLICATE(1.0d0,npix)#(dindgen(nLSFpix)-nLSFpix/2)*dx
;      ;xlsf += dindgen(npix)#REPLICATE(1.0d0,nLSFpix)
;      ;xcenter = dindgen(npix)
;      xlsf = REPLICATE(1.0d0,nx)#(dindgen(nLSFpix)-nLSFpix/2)*dx
;      xlsf += dindgen(nx)*dx#REPLICATE(1.0d0,nLSFpix)
;      xcenter = dindgen(nx)*dx
;  
;      lsf2d = LSF_GH(xlsf,xcenter,lsfpars)
;      lsftot = total(lsf2d,2)
;      lsf2d /= lsftot#replicate(1,nLSFpix)/osamp  ; make sure the LSF is normalized
;  
;      ; Step 4-Multiply LSF by spectrum 2D arrays and sum/collapse
;      out[*,i,0] = TOTAL(lsf2d*tspec1_2d,2)/osamp
;      out[*,i,1] = TOTAL(lsf2d*tspec2_2d,2)/osamp
;      out[*,i,2] = TOTAL(lsf2d*tspec3_2d,2)/osamp
    endfor  ;end fibers


    endfor  ; end scale
    ; write this out for subsequent use
    if iair eq 0 then begin
      mkhdr,hdr,twave
      if sxpar(telhead1,'CRVAL2') gt 0 then air0=sxpar(telhead1,'CRVAL2') else air0=1.
      if sxpar(telhead1,'CDELT2') gt 0 then dair=sxpar(telhead1,'CDELT2') else dair=1.
      sxaddpar,hdr,'AIR0',air0
      sxaddpar,hdr,'DAIR',dair
      sxaddpar,hdr,'NSPECIES',nspecies
      sxaddpar,hdr,'NSCALE',nscale
      mwrfits,twave,outfile[j],hdr,/create
    endif
    mwrfits,out,outfile[j]

    endfor  ; end airmass
    ; remove lock file
    file_delete,outfile[j]+'.lock'

    nextchip:
  endfor

  file_delete,lockfile
endif


; now get the convolved telluric at the desired wavelengths for this frame
; these will be returned in a chip-concatenated array

if keyword_set(convonly) then return,0

alt=sxpar(frame.(0).header,'ALT')
; if not ALT card, adopt 60 degrees
if alt lt 5 then alt=60
airmass=1./cos((90-alt)*!pi/180.)
for j=0,2 do begin

  ; read in convolved spectra
  twave=mrdfits(outfile[j],0,hdr)
  if sxpar(hdr,'AIR0') eq 0 then iair=0 else $
    iair=fix((airmass-sxpar(hdr,'AIR0'))/sxpar(hdr,'DAIR'))
print,'airmass: ', airmass,' iair: ', iair, sxpar(hdr,'air0'), sxpar(hdr,'dair')
  ; do go past last tabulated airmass!
  status=-1
  while status ne 0 do begin
    out1=mrdfits(outfile[j],1+iair,status=status)
    if status ne 0 then iair-=1
  endwhile
  ; interpolate in airmass if we have airmass grid and bracket the observed airmass
  out2=mrdfits(outfile[j],2+iair,status=status)
  if status eq 0 then begin
    air1=sxpar(hdr,'AIR0')+iair*sxpar(hdr,'DAIR')
    air2=sxpar(hdr,'AIR0')+(iair+1)*sxpar(hdr,'DAIR')
    out=out1+(airmass-air1)/(air2-air1)*(out2-out1)
  endif else out=out1
  apgundef,out1
  apgundef,out2
  sz=size(out)
  if sz[0] eq 2 then out=reform(out,sz[1],sz[2],1,1)

  sz=size(frame.(j).flux)
  npix=sz[1]

  nspecies=sxpar(hdr,'NSPECIES')
  if nspecies eq 0 then nspecies=3
  nscale=sxpar(hdr,'NSCALE')
  if nscale eq 0 then nscale=1
  pix=indgen(npix)
  if j eq 0 then conv=dblarr(3*2048,300,nspecies,nscale)

  ; Loop over fibers
  for ii=0,n_elements(fibers)-1 do begin
    i=fibers(ii)
    ; Interpolate to desired pixel sampling
    wcoef = reform(frame.(j).wcoef[i,*])
    wave = PIX2WAVE(pix,wcoef)
    wsi=sort(twave[*,i])
    si=sort(wave)
    for k=0,nspecies-1 do begin
     for iscale=0,nscale-1 do begin
      interp=dblarr(npix)
      interp[si] = SPLINE(twave[wsi,i],out[wsi,i,k,iscale],wave[si],/double)
      conv[j*npix:j*npix+npix-1,i,k,iscale] = interp
      ;conv[j*npix:j*npix+npix-1,i,k,iscale] = SPLINE(twave[wsi,i],out[wsi,i,k,iscale],wave[si],/double)
     endfor
    endfor
  endfor
endfor

return, conv

end
