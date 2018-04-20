;+
;
; APPERSISTMODEL
;
; This creates a persistence correction image using a persistence
; model file for a specific "target" image.
;
; INPUTS:
;  file           The "target" file for which to create the
;                   persistence correction.
;  histfile       The history file of 2D image for the night.
;  persistmodelfile  The persistence model file.
;  =bpmfile       Filename of bad pixel mask.
;  /correction    Apply a pre-dark correction to the persistence model.  The default is no correction.
;
; OUTPUTS:
;  pmodelim       The 2D persistence model image.
;  par            The correction factors.
;  =error         The error message if one occurred.
;
; USAGE:
;  IDL>appersistmodel,'apR-c-12690044.apz','apHist-c-56831.fits','apPersist-c-57184.fits',pmodelim,par
;
; By D.Nidever and D.Nguyen  June 2015
;-

pro appersistmodel,file,histfile,persistmodelfile,pmodelim,par,bpmfile=bpmfile,correction=correction,error=error

;apgundef,pmodelim
pmodelim = dblarr(2048,2048)
par = [1.0, 1.0, 0.0]

if n_elements(file) eq 0 or n_elements(histfile) eq 0 or n_elements(persistmodelfile) eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax -  appersistmodel,file,histfile,persistmodelfile,pmodelim,par,bpmfile=bpmfile,correction=correction,error=error'
  return
endif

target_expnum = first_el(strsplit(file,'-',/extract),/last)
target_expnum = first_el(strsplit(target_expnum,'.',/extract))

print,'Creating persistence correction model for ',target_expnum

; Steps:
; 1) Load the apHist file
; 2) Find target file, previous stimulus and pre-dark images in apHist file
; 3) Create persistence models for the pre-darks
; 4) Derive corrections to the persistence model using the pre-darks
; 5) Create persistence model for the "target" file including the
;     correction factors


;--------------------------------------
; 1) Load the apHist file information
;--------------------------------------
; Get number of extensions
next = 0
fits_open,histfile,fcb
next=fcb.nextend
fits_close,fcb
expstr = replicate({exten:0L,expnum:'',exptype:'',exptime:0.0,nreads:0L,dateobs:'',plateid:'',jd:0.0d0,expstart:0.0d0,expend:0.0d0},next)
for i=0,next-1 do begin
  head = headfits(histfile,exten=i+1)
  rawfilename = sxpar(head,'FILENAME')  ; '/data-ql/data/56831/apRaw-12690001-003.fits'
  expnum = strmid(file_basename(rawfilename,'.fits'),6,8)
  expstr[i].exten = i+1
  expstr[i].expnum = expnum
  expstr[i].exptype = strupcase(strtrim(sxpar(head,'EXPTYPE'),2))
  expstr[i].exptime = sxpar(head,'EXPTIME')
  expstr[i].nreads = sxpar(head,'NREAD')
  expstr[i].dateobs = sxpar(head,'DATE-OBS')
  expstr[i].plateid = strtrim(sxpar(head,'PLATEID'),2)
  expstr[i].jd = date2jd(expstr[i].dateobs)
endfor
; Calculating exposure start/end times relative to start of first
;  exposure of the night
minjd = min(expstr.jd)
expstr.expstart = (expstr.jd-minjd)*24*3600 ; convert to sec from start of 1st exp
expstr.expend = expstr.expstart+expstr.exptime


;---------------------------------------------------------------------------
; 2) Find target file, previous stimulus and pre-dark images in apHist file
;---------------------------------------------------------------------------
ind = first_el(where(expstr.expnum eq target_expnum,nind))
if nind eq 0 then begin
  print,target_expnum,' NOT FOUND in ',histfile
  return
endif
texpstr = expstr[ind]
min_stimulus_expstart = texpstr.expstart
; Check for prior object exposure, from different plate and within the last hour
ind_object = first_el(where(expstr.exptype eq 'OBJECT' and (expstr.jd-texpstr.jd lt 0) and (expstr.jd-texpstr.jd gt -3.0/24.0) and $
                            expstr.plateid ne texpstr.plateid,nind_object),/last)
; Check for prior domeflat exposure
ind_domeflat = first_el(where(expstr.exptype eq 'DOMEFLAT' and (expstr.jd-texpstr.jd lt 0) and $
                              (expstr.jd-texpstr.jd gt -3.0/24.0),nind_domeflat),/last)
; Check for prior pre-darks, after the "stimulus" exposures and within
;   the last hours
if nind_object gt 0 then min_stimulus_expstart=min_stimulus_expstart < expstr[ind_object].expstart
if nind_domeflat gt 0 then min_stimulus_expstart=min_stimulus_expstart < expstr[ind_domeflat].expstart
ind_predark = first_el(where(expstr.exptype eq 'DARK' and (expstr.jd-texpstr.jd lt 0) and (expstr.jd-texpstr.jd gt -3.0/24.0) and $
                             expstr.expstart gt min_stimulus_expstart,nind_predark),/last)
; No prior object or domeflat exposures
if nind_object eq 0 and nind_domeflat eq 0 then begin
  print,'No prior OBJECT or DOMEFLAT exposure.  Not producing a persistence correction model.'
  return
endif
print,'Using  object: ', ind_object
print,'Using  dome: ', ind_domeflat
print,'Using  dark: ', ind_predark
   

;------------------------------------------------
; 3) Create persistence models for the pre-darks
;------------------------------------------------
cc_data = APCALCPERSIST_LOADTABLE(persistmodelfile)

; Object exposure
if nind_object gt 0 then begin
  print,'Using  object: ', ind_object, expstr[ind_object].exten
  FITS_READ, histfile, object_flux, exten_no=expstr[ind_object].exten
  object_flux = double(object_flux[0:2047,*])
  object_coeffs = APCALCPERSIST_COEFFS(cc_data, object_flux,/domask)
endif

; Domeflat exposure
if nind_domeflat gt 0 then begin
  print,'Using  dome: ', ind_domeflat, expstr[ind_domeflat].exten
  FITS_READ, histfile, domeflat_flux, exten_no=expstr[ind_domeflat].exten
  domeflat_flux = double(domeflat_flux[0:2047,*])
  domeflat_coeffs = APCALCPERSIST_COEFFS(cc_data, domeflat_flux,/domask)
endif

; Predark
if nind_predark gt 0 then begin
  print,'Using  dark: ', ind_predark, expstr[ind_predark].exten
  FITS_READ, histfile, predark_flux, exten_no=expstr[ind_predark].exten
  predark_flux = double(predark_flux[0:2047,*])

  if nind_object gt 0 then begin
    object_persistence_start = APCALCPERSIST(object_coeffs, expstr[ind_predark].expstart-expstr[ind_object].expend)
    object_persistence_end = APCALCPERSIST(object_coeffs, expstr[ind_predark].expend-expstr[ind_object].expend)
    object_persistence = object_persistence_end - object_persistence_start
  endif else object_persistence=dblarr(2048,2048)
  if nind_domeflat gt 0 then begin
    domeflat_persistence_start = APCALCPERSIST(domeflat_coeffs, expstr[ind_predark].expstart-expstr[ind_domeflat].expend)
    domeflat_persistence_end = APCALCPERSIST(domeflat_coeffs, expstr[ind_predark].expend-expstr[ind_domeflat].expend)
    domeflat_persistence = domeflat_persistence_end - domeflat_persistence_start
  endif else domeflat_persistence=dblarr(2048,2048)
endif
  

;--------------------------------------------------------------------
; 4) Derive corrections to the persistence model using the pre-darks
;--------------------------------------------------------------------

if nind_predark gt 0 and keyword_set(correction) then begin

  ; Fix the third quadrant of the dark
  ;med1 = median(predark_flux[4:511,100:200])
  ;med2 = median(predark_flux[512:2*512-1,100:200])
  ;med3 = median(predark_flux[2*512:3*512-1,100:200])
  ;med4 = median(predark_flux[3*512:2044,100:200])
  ;q3off = med3 - median([med1,med2,med4],/even)
  ;predark_flux[2*512:3*512-1,*] -= q3off

  sz = size(object_flux)
  x = findgen(sz[1])
  y = findgen(sz[2])
  z = predark_flux > 0
  ;err = sqrt(predark_flux>1)
  err = predark_flux*0+1  ; unweighted gives better results

  im1 = object_persistence
  bd1 = where(finite(im1) eq 0,nbd1)
  if nbd1 gt 0 then im1[bd1]=0
  im2 = domeflat_persistence
  bd2 = where(finite(im2) eq 0,nbd2)
  if nbd2 gt 0 then im2[bd2]=0

  ; mask the pixels
  mask = cc_data[*,*,0]
  z*=mask & z[0:3,*]=0 & z[2044:2047,*]=0 & z[*,2020:2047]=0
  im1*=mask & im1[0:3,*]=0 & im1[2044:2047,*]=0 & im1[*,2020:2047]=0
  im2*=mask & im2[0:3,*]=0 & im2[2044:2047,*]=0 & im2[*,2020:2047]=0
  err=err*mask + (1-mask)*1e20 & err[0:3,*]=1e20 & err[2044:2047,*]=1e20 & err[*,2020:2047]=1e20

  ; Mask bad pixels using the bad pixel mask file
  if keyword_set(bpmfile) then begin
     ;fits_read,'~/apogee/spectro/r5/cal/bpm/apBPM-c-05560001.fits',bpmim
     fits_read,bpmfile,bpmim
     bpm = (bpmim gt 0)
     err = err*(1-bpm) + bpm*1e20
     z *= (1-bpm)
  endif

  ;z = z < 2000
  ;im1 = im1 < 2000
  ;im2 = im2 < 2000
  ;inp1 = im42
  ;inp2 = im43
  ;inp1[*,0:1600]=0 & inp1[0:3,*]=0 & inp1[2044:2047,*]=0 & inp1[*,2020:2047]=0
  ;inp2[*,0:1600]=0 & inp2[0:3,*]=0 & inp2[2044:2047,*]=0 & inp2[*,2020:2047]=0

; TRY TO FIT BOTH PRE-DARKS SIMULTANEOUSLY!
  
  fa = {im1:im1, im2:im2}
  initpar = [1.0, 1.0, 0.0]
  parinfo = replicate({limited:[0,0],limits:[0.0,0.0],fixed:0},n_elements(initpar))
  parinfo[0:1].limited[0]=1
  par = MPFIT2DFUN('func_lincomb', X, Y, Z, ERR, initpar, parinfo=parinfo, functargs=fa, status=status, bestnorm=chisq, /quiet)
  model = func_lincomb(x,y,par,_extra=fa)

  print,'Correction factors = ',par
  
Endif else par = [1.0, 1.0, 0.0]  ; no correction


;------------------------------------------------------------------
; 5) Create persistence model for the "target" file including the
;     correction factors
;------------------------------------------------------------------

target_exposure_start = (expstr[ind].jd-min(expstr.jd))*24*3600   ; convert to sec from start of 1st exp
target_exposure_duration = expstr[ind].exptime
target_exposure_end = target_exposure_start + target_exposure_duration

if nind_object gt 0 then begin
  target_object_persistence_start = APCALCPERSIST(object_coeffs, texpstr.expstart-expstr[ind_object].expend)
  target_object_persistence_end = APCALCPERSIST(object_coeffs, texpstr.expend-expstr[ind_object].expend)
  target_object_persistence = target_object_persistence_end - target_object_persistence_start
endif else target_object_persistence=dblarr(2048,2048)
if nind_domeflat gt 0 then begin
  target_domeflat_persistence_start = APCALCPERSIST(domeflat_coeffs, texpstr.expstart-expstr[ind_domeflat].expend)
  target_domeflat_persistence_end = APCALCPERSIST(domeflat_coeffs, texpstr.expend-expstr[ind_domeflat].expend)
  target_domeflat_persistence = target_domeflat_persistence_end - target_domeflat_persistence_start
endif else target_domeflat_persistence=dblarr(2048,2048)
  
; Apply the correction
pmodelim = func_lincomb(x,y,par,im1=target_object_persistence,im2=target_domeflat_persistence)
pmodelim = float(pmodelim)
pmodelim = pmodelim > 0  ; make sure it's positive

;stop

end
