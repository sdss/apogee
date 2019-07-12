pro aspcap_rdobsnew,file,obs,gaps,obspixrang=obspixrang,apvisit=apvisit,endgaps=endgaps,skyerr=skyerr,skyfact=skyfact,visits=visits,persist=persist

; NAME:   
;  aspcap_rdobsnew    
;                    
; PURPOSE:     
;       This procedure reads the APOGEE visit and star files and returns
;       a structure with the fluxes, vacuum wavelengths in Angstrom,
;       the invariance, the sky fluxes and their errors, 
;       the telluric fluxes and their errors and a flag defining
;       the bad pixels. Fluxes, errors and invariance are given 
;       in (10^-17 ergs/s/cm^2/Ang)
;                                                                                                                                                        
; CALLING SEQUENCE:      
;       ASPCAP_RDOBNEW(file)
;
; INPUT:
;        file              -string      name of the apVisit file to be
;                                       read
; OUTPUT:
;        obs               -structure   
;            flux          -fltarr      fluxes TO ASK DAVID? (10^-17 ergs/s/cm^2/Ang)
;            wave          -dblarr      wavelengths
;            invar         -dblarr      square invariance
;            flag          -intarr       1   bad pixels
;                                        2   cosmic ray
;                                        4   saturated 
;                                        8   unfixable
;            sky           -fltarr       sky fluxes           (10^-17 ergs/s/cm^2/Ang)
;            sky_error     -fltarr       sky flux errors      (10^-17 ergs/s/cm^2/Ang)
;            tell          -fltarr       telluric fluxes      (10^-17 ergs/s/cm^2/Ang)
;            tell_erro     -fltarr       telluric flux errors (10^-17 ergs/s/cm^2/Ang)
;            wav_coeff     -fltarr       coefficients for the
;                                        wavelength solution?
;            lsf           -fltarr       coefficient for the a
;                                        Gauss-Hermite LSF
;                                        decomposition
;           gaps           -lonarr      is an array specifying the
;                                       first and last pixel of each
;                                       detector
; KEYWORDS:
;
;       apvisit            -             if set, the procedure assumes
;                                        an apVisit format (apStar by default) 
;       endgaps            -fltarr       the pixel number (with first
;                                        pixel eq 1) for the last
;                                        pixel of each detector gap

; BY ANA ELIA GARCIA PEREZ - MARCH 2010
;    Modified by J. Holtzman, 8/16/12
; 

; get main header 
tmp = readfits(file,head0,ext=0,/silent)

; read apStar or apVisit file to get wavelengths
if n_elements(apvisit) eq 0 then begin 

  ; for apStar file, wavelength scale is in header, use modified error array from sky error
  tmp=readfits(file,head,ext=0,/silent)
  nvisits=sxpar(head,'NVISITS')
  ; with visits keyword, take all visits for visits=1 or min(visits,nvisits) for visits>1
  if keyword_set(visits) and nvisits gt 1 then begin
    if visits gt 1 then nvisits=min([visits,nvisits])
    nspec=nvisits+1 
  endif else nspec=1

  ; flux
  tmp=readfits(file,head,ext=1,/silent)
  sz=size(tmp,/dim)
  npix=sz[0]
  if keyword_set(visits) and nvisits gt 1 then begin
    flux=fltarr(npix,nvisits+1)
    flux[*,0]=tmp[*,0]
    flux[*,1:nvisits]=tmp[*,2:nvisits+1]
  endif else flux=tmp[*,0]

  ; wavelength
  tmp=sxpar(head,'CRVAL1')+indgen(sxpar(head,'NAXIS1'))*sxpar(head,'CDELT1')
  wave=10^tmp

  ; error
  tmp=readfits(file,head,ext=2,/silent)
  if keyword_set(visits) and nvisits gt 1 then begin
    err=fltarr(npix,nvisits+1)
    err[*,0]=tmp[*,0]
    err[*,1:nvisits]=tmp[*,2:nvisits+1]
  endif else err=tmp[*,0]

  ; sky
  tmp=readfits(file,head,ext=4,/silent)
  ntmp=n_elements(tmp)
  if ntmp eq 1 then tmp=intarr(npix)
  if keyword_set(visits) and nvisits gt 1 then begin
    sky=fltarr(npix,nvisits+1)
    sky[*,0]=tmp[*,0]
    sky[*,1:nvisits]=tmp[*,2:nvisits+1]
  endif else sky=tmp[*,0]
  if keyword_set(skyerr) then begin
    eerr=enhancederr(flux,err,sky,sig=skyerr,skyfact=skyfact)
    err=eerr
  endif

  ; mask
  tmp=readfits(file,head,ext=3,/silent)
  ntmp=n_elements(tmp)
  if ntmp eq 1 then tmp=intarr(npix)
  if keyword_set(visits) and nvisits gt 1 then begin
    flag=intarr(npix,nvisits+1)
    flag[*,0]=tmp[*,0]
    flag[*,1:nvisits]=tmp[*,2:nvisits+1]
  endif else flag=tmp[*,0]
  ;gaps=intarr(2)
  ;gaps[0]=0
  ;gaps[1]=n_elements(flux)-1
  gaps=intarr(2,3)
  gaps[0,0]=sxpar(head0,'ROVERMIN')
  gaps[1,0]=sxpar(head0,'ROVERMAX')
  gaps[0,1]=sxpar(head0,'GOVERMIN')
  gaps[1,1]=sxpar(head0,'GOVERMAX')
  gaps[0,2]=sxpar(head0,'BOVERMIN')
  gaps[1,2]=sxpar(head0,'BOVERMAX')

  ; with persist keyword, do a simple adjustment of blue chip level to match red chip level
  if keyword_set(persist) then begin
    bmed=median(flux[gaps[0,2]:gaps[1,2],0])
    gmed=median(flux[gaps[0,1]:gaps[1,1],0])
    jk=sxpar(head0,'J')-sxpar(head0,'K')
    btarg=(1.09151-0.0736803*jk)*gmed
    if file_test(persist) then openw,pfile,persist,/get_lun,/append else openw,pfile,persist,/get_lun
    printf,pfile,format='(5f10.2,6i6)',sxpar(head0,'H'),jk,bmed,gmed,btarg,$
     sxpar(head0,'ANDFLAG') and starflagval('PERSIST_HIGH'),$
     sxpar(head0,'ANDFLAG') and starflagval('PERSIST_MED'),$
     sxpar(head0,'ANDFLAG') and starflagval('PERSIST_LOW'),$
     sxpar(head0,'STARFLAG') and starflagval('PERSIST_HIGH'),$
     sxpar(head0,'STARFLAG') and starflagval('PERSIST_MED'),$
     sxpar(head0,'STARFLAG') and starflagval('PERSIST_LOW')
    free_lun,pfile

    if (sxpar(head0,'ANDFLAG') and starflagval('PERSIST_HIGH')) gt 0 then begin
      oldflux=flux[*,0]
      flux[gaps[0,2]:gaps[1,2],0]+=(btarg-bmed)
      ;smcolor
      ;plot,oldflux,yr=[0,gmed*1.5]
      ;oplot,flux,color=2
      ;stop
    endif
  endif

endif else begin 

  nvisits=1

  ; apVisit files have data in 3 lines (one per chip), in decreasing wavelength
  ; transform these into 1D arrays sorted in increasing wavelength
  wave = readfits(file,head,ext=4,/silent)
  tmp=reform(wave)
  index=sort(tmp)
  wave=tmp[index]
  flux = readfits(file,head,ext=1,/silent)
  tmp=reform(flux)
  flux=tmp[index]
  err = readfits(file,head,ext=2,/silent)
  tmp=reform(err)
  err=tmp[index]
  flag = readfits(file,head,ext=3,/silent)
  tmp=reform(flag)
  flag=tmp[index]

  sky = readfits(file,head,ext=5,/silent)
  tmp=reform(sky)
  sky=tmp[index]
  if keyword_set(skyerr) then begin
    eerr=enhancederr(flux,err,sky,sig=skyerr,skyfact=skyfact)
    err=eerr
  endif

  gaps=intarr(2,3)
  npix=0
  for ichip=0,2 do begin
   gaps[0,ichip]=npix
   npix+=n_elements(flux)/3-1
   gaps[1,ichip]=npix
   npix+=1
  endfor
  wav_coef = readfits(file,head,ext=9,/silent)
  lsf = readfits(file,head,ext=10,/silent)

endelse

;; renormalize in case there is a bad normalization (e.g., screwed-up delta Oph)
gd=where(flux[0:npix-1] gt 0,ngd)
if ngd gt 0 then amed=median(flux[gd]) else amed=1.
flux/=amed
err/=amed


invar=1.d0/err^2
;bad=where(err lt 0 or finite(err) eq 0 or finite(flux) eq 0 or flux lt 0.001)
;invar[bad]=0.
;flux[bad]=0.
sz=size(flux,/dim)
npix=sz[0]
if n_elements(apvisit) eq 0 then $
  obs = {head: head0, wave: wave, flux: reform(flux,npix,nspec), $
         invar: reform(invar,npix,nspec), flag: reform(flag,npix,nspec)} else $
  obs = {head: head0, wave: wave, flux: reform(flux,npix,nspec), $
         invar: reform(invar,npix,nspec), flag: reform(flag,npix,nspec), $
         wav_coef: wav_coef, lsf: lsf}
  
end
