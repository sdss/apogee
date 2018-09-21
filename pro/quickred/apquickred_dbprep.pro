;+
;
; APQUICKRED_DBPREP
;
; This bins up the 2D image to be output to the database
; Also puts the extracted spectra and other info in STR
;
; INPUTS:
;  output2d    The 2D image output from ap3dquickred.pro
;  dbstr       The database structure created in apquickred.pro
;
; OUTPUTS:
;  The binned 2D images and 1D spectra (when available) are
;  updated in the DBSTR structure.
;  =error      The error message if one occurred.
;
; USAGE:
;  IDL>apquickred_dbprep,output2d,dbstr
;
; By D.Nidever  May 2011
;-

pro apquickred_dbprep,output2d,dbstr,error=error,stp=stp

; Error Handling
;------------------
; Establish error handler. When errors occur, the index of the
; error is returned in the variable Error_status:  
;CATCH, Error_status 
;
;This statement begins the error handler:  
;if (Error_status ne 0) then begin 
;   error = !ERROR_STATE.MSG  
;   if not keyword_set(silent) then print,error
;   CATCH, /CANCEL 
;   return
;endif

; Not enough inputs
if n_elements(output2d) eq 0 or n_elements(dbstr) eq 0 then begin
  error = 'Not enough inputs'
  if not keyword_set(silent) then print,'Syntax - apquickred_dbprep,output2d,dbstr,error=error,stp=stp'
  return
endif

; Parameters
npix = 2048L
nchips = 3L
nbin = 8
npixbin = npix/nbin
chipgap = 300
tags = tag_names(dbstr)
dum = where(tags eq 'FRAME',ntags_frame)

;if ntags_frame gt 0 then nfibers = n_elements(dbstr.frame[0].flux[0,*]) else nfibers=0
if ntags_frame gt 0 then nfibers = n_elements(dbstr.frame.(0).flux[0,*]) else nfibers=0

; sub regions
;yloarr = [500,1000,1500]
;yhiarr = yloarr+99
;nsub = n_elements(yloarr)
;
;header = output1d[0].header
;
;; Create the structure
;str = {frameid:info.fid8,exptime:info.exptime,header:output1d[0].header,plate:info.plate,$
;       date:info.dateobs,mjd:info.mjd5,arraydisplay_nbin:nbin,$
;       arraydisplay:{data:bytarr(npixbin,npixbin*nchips),bscale:0.0,bzero:0.0,zscale:[0.0,0.0]},$
;       arraydisplay_sub:REPLICATE({data:PTR_NEW(),bscale:0.0,bzero:0.0,zscale:[0.0,0.0],yrange:[0L,0L]},nsub),$
;       spec:fltarr(nfibers,npix,nchips)}

; Bin the entire array
nbin = dbstr.arraydisplay_nbin
im = output2d[*,*,0]
npixbin = npix/nbin
;binim = REBIN(im,npixbin,npixbin*nchips)  ; averages values
binim = REBIN(im,npixbin*nchips,npixbin)  ; averages values
;binim = transpose(binim)          ; flip it
zscale,binim,z1,z2,contrast=0.10
 ;0.25; scale it
minim = min(binim)
maxim = max(binim)
; real_values = byte_values * bscale + bzero
bzero = minim
bscale = (maxim-minim)/255. > 1  ; don't scale up
byte_binim = ( binim - bzero )/bscale
byte_binim = byte( round( byte_binim ) ) ; round and make byte type
; stuff into the structure
dbstr.arraydisplay.data = byte_binim
dbstr.arraydisplay.bscale = bscale
dbstr.arraydisplay.bzero = bzero
dbstr.arraydisplay.zscale = [z1,z2]

nsub = n_elements(dbstr.arraydisplay_sub)

; Loop through the subranges
FOR i=0,nsub-1 DO BEGIN

  ;ylo = yloarr[i]
  ;yhi = yhiarr[i]
  ylo = dbstr.arraydisplay_sub[i].yrange[0]
  yhi = dbstr.arraydisplay_sub[i].yrange[1]
  nypix = yhi-ylo+1

  im1 = im[*,ylo:yhi]
     
  ; Only bin in Y-dimension
  npixbin = npix/nbin
  nypixbin = nypix/nbin
  ;binim = REBIN(im,nxpix,npixbin*nchips)  ; averages values
  binim = REBIN(im1,npixbin*nchips,nypix)  ; averages values
  zscale,binim,z1,z2,contrast=0.10 ;0.25
  ; scale it
  minim = min(binim)
  maxim = max(binim)
  ; real_values = byte_values * bscale + bzero
  bzero = minim
  bscale = (maxim-minim)/255. > 1  ; don't scale up
  byte_binim = ( binim - bzero )/bscale
  byte_binim = byte( round( byte_binim ) ) ; round and make byte type
  ; stuff into the structure
  ;if PTR_VALID(dbstr.arraydisplay_sub[i].data) then PTR_FREE,dbstr.arraydisplay_sub[i].data
  ;dbstr.arraydisplay_sub[i].data = PTR_NEW(byte_binim,/no_copy)
  dbstr.arraydisplay_sub[i].data = byte_binim
  dbstr.arraydisplay_sub[i].bscale = bscale
  dbstr.arraydisplay_sub[i].bzero = bzero
  ;dbstr.arraydisplay_sub[i].yrange = [ylo,yhi]
  dbstr.arraydisplay_sub[i].zscale = [z1,z2]

  ;stop

ENDFOR


; Extracted spectra
if nfibers gt 0 then begin

  ; Loop through the fibers
  for i=0L,nfibers-1 do begin

    flux = fltarr(npix*nchips)
    err = fltarr(npix*nchips)
    ; for j=0,2 do flux[j*npix:j*npix+npix-1]=dbstr.frame[j].flux[*,i]
    ; for j=0,2 do err[j*npix:j*npix+npix-1]=dbstr.frame[j].err[*,i]
    for j=0,2 do flux[j*npix:j*npix+npix-1]=dbstr.frame.(j).flux[*,i]
    for j=0,2 do err[j*npix:j*npix+npix-1]=dbstr.frame.(j).err[*,i]

    ; Scale the data
    minflux = min(flux)
    maxflux = max(flux)
    ; real_values = byte_values * bscale + bzero
    bzero = minflux
    ; UINT goes from 0-65535
    ; integer 
    bscale = (maxflux-minflux)/65535. > 1  ; don't scale up
    uint_flux = ( flux - bzero )/bscale
    uint_flux = uint( round( uint_flux ) ) ; round and make uint type

    dbstr.qr_spectrum[i].spectrum = uint_flux
    dbstr.qr_spectrum[i].bscale = bscale
    dbstr.qr_spectrum[i].bzero = bzero
    dbstr.qr_spectrum[i].fiberid = 300-i   ; fiberid = 300-index
    dbstr.qr_spectrum[i].medsnr = MEDIAN(flux/err) 

  end ; fiber loop

endif

if keyword_set(stp) then stop

end
