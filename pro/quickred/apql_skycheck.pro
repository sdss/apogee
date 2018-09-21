;+
;
; APQL_SKYCHECK
;
; This program checks the sky for the quicklook software
;
; INPUTS:
;  str      The quicklook structure.
;  allstr   The quicklook structure of all previous reads of this exposure.
;  /silent  Don't print anything to the screen.
;
; OUTPUTS:
;  The sky values are updated in the STR structure.
;  =error   The error message if one occurred.
;
; USAGE:
;  IDL>apl_skycheck,str,allstr
;
; By D.Nidever  2010
;-
pro apql_skycheck,str,allstr,silent=silent,error=error

; Error Handling
;------------------
; Establish error handler. When errors occur, the index of the
; error is returned in the variable Error_status:  
;CATCH, Error_status 

;This statement begins the error handler:  
;if (Error_status ne 0) then begin 
;   error = !ERROR_STATE.MSG  
;   if not keyword_set(silent) then print,error
;   CATCH, /CANCEL 
;   return
;endif

; Not enough inputs
if n_elements(str) eq 0 then begin
  error = 'Not enough inputs'
  if not keyword_set(silent) then print,'Syntax - apql_skycheck,str,allstr,silent=silent,error=error'
  return
endif

; make sure the system variable exists before we proceed
DEFSYSV,'!apql',exists=apql_exists
if not apql_exists then begin
  error = 'APQL_SKYCHECK: system variable !apql is not defined'
  if not keyword_set(silent) then print,error
  return
endif

; Check that plugmap exists
if not PTR_VALID(!apql.plugmap.datastr) then begin
  error = 'APQL_SKYCHECK: no plugmap data available'
  if not keyword_set(silent) then print,error
  return
endif

; Check that the extracted spectra are there
if not PTR_VALID(str.frame) then begin
  error = 'APQL_SKYCHECK: NO extracted spectra'
  if not keyword_set(silent) then print,error
  return
endif

; Plot up the sky fibers and only show the regions around where
; you expect the sky lines to be
; Also show the estimated positions of the sky lines


;spectro_dir = APGETDIR('APQLSPECTRO_DIR',/exists,error=direrr)
;if n_elements(direrr) gt 0 then return
;linelist_dir = spectro_dir+'lib/linelists/'

npix = 2048L

;skyind = where(plugmap.fiberdata.objtype eq 'SKY',nskyind)
skyind = where( (*!apql.plugmap.datastr).fiberdata.objtype eq 'SKY' and $
                (*!apql.plugmap.datastr).fiberdata.holetype eq 'OBJECT' AND $
                (*!apql.plugmap.datastr).fiberdata.spectrographid eq 2 and $
                (*!apql.plugmap.datastr).fiberdata.fiberid ne -1,nskyind)
if nskyind eq 0 then begin
  error = 'APQL_SKYCHECK: No Sky fibers'
  if not keyword_set(silent) then print,error
  return
endif
;skyfiberind = plugmap.fiberdata[skyind].fiberid-1
;skyfiberind = (*!apql.plugmap.datastr).fiberdata[skyind].fiberid-1
; fiberid=1 is at the top of the detector or index=299
; index = 300-fiberid
skyfiberind = 300-(*!apql.plugmap.datastr).fiberdata[skyind].fiberid

skyframe1 = {flux:(*str.frame)[0].flux[*,skyfiberind],err:(*str.frame)[0].err[*,skyfiberind],$
             mask:(*str.frame)[0].mask[*,skyfiberind]}
skyframe2 = {flux:(*str.frame)[1].flux[*,skyfiberind],err:(*str.frame)[1].err[*,skyfiberind],$
             mask:(*str.frame)[1].mask[*,skyfiberind]}
skyframe3 = {flux:(*str.frame)[2].flux[*,skyfiberind],err:(*str.frame)[2].err[*,skyfiberind],$
             mask:(*str.frame)[2].mask[*,skyfiberind]}


; Maybe just median sky spectrum
medflux1 = median(skyframe1.flux,dim=2)
mederr1 = median(skyframe1.err,dim=2)
;medframe1 = {flux:medflux1,err:mederr1,mask:lonarr(npix)}
;APPEAKFIT,medframe1,medlinestr1,/nogauss,/silent

medflux2 = median(skyframe2.flux,dim=2)
mederr2 = median(skyframe2.err,dim=2)
;medframe2 = {flux:medflux2,err:mederr2,mask:lonarr(npix)}
;APPEAKFIT,medframe2,medlinestr2,/nogauss,/silent

medflux3 = median(skyframe3.flux,dim=2)
mederr3 = median(skyframe3.err,dim=2)
;medframe3 = {flux:medflux3,err:mederr3,mask:lonarr(npix)}
;APPEAKFIT,medframe3,medlinestr3,/nogauss,/silent

medflux = [medflux1, medflux2, medflux3]

;str.medsky[0].flux = medflux1
;str.medsky[1].flux = medflux2
;str.medsky[2].flux = medflux3
str.medsky[0].err = mederr1
str.medsky[1].err = mederr2
str.medsky[2].err = mederr3


; Cross-correlate all of the fibers against each other to get
; zero-point shifts and a "global" spectrum
;------------------------------------------
xshift = fltarr(nskyind)

fiber_arr1 = skyframe1.flux*0
fiber_arr2 = skyframe2.flux*0
fiber_arr3 = skyframe3.flux*0

fiber_arr1[*,0] = skyframe1.flux[*,0]
fiber_arr2[*,0] = skyframe2.flux[*,0]
fiber_arr3[*,0] = skyframe3.flux[*,0]

refspec = medflux2
refspec -= MEDFILT1D(refspec,150,/edge)

x = lindgen(npix)
for i=1,nskyind-1 do begin
  fiber1 = skyframe2.flux[*,i]
  fiber1 -= MEDFILT1D(fiber1,150,/edge)
  XCORLB,medflux2,fiber1,20,xsh,error=xerror
  if n_elements(xerror) gt 0 then xsh=0
  ;fiber2 = spline(x,skyframe2.flux[*,i],x-xsh)
  
  ; Shift the spectra
  spec1 = shift(skyframe1.flux[*,i],xsh)
  spec2 = shift(skyframe2.flux[*,i],xsh)
  spec3 = shift(skyframe3.flux[*,i],xsh)

  if xsh gt 0 then begin
    spec1[0:ceil(xsh)] = 0    ; fix the ends
    spec2[0:ceil(xsh)] = 0    ; fix the ends
    spec3[0:ceil(xsh)] = 0    ; fix the ends
    ;fiber2[0:ceil(xsh)]=0    ; fix the ends
  endif
  if xsh lt 0 then begin
    spec1[npix-floor(abs(xsh))-1:npix-1] = 0
    spec2[npix-floor(abs(xsh))-1:npix-1] = 0
    spec3[npix-floor(abs(xsh))-1:npix-1] = 0
    ;fiber2[npix-floor(abs(xsh))-1:npix-1] = 0
  endif

  fiber_arr1[*,i] = spec1
  fiber_arr2[*,i] = spec2
  fiber_arr3[*,i] = spec3
  ;fiber_arr[*,i] = fiber2

  xshift[i] = xsh

  ;print,xsh
  ;stop
end

;stop


; Take median of shifted sky spectra
medflux1 = median(fiber_arr1,dim=2)
medflux2 = median(fiber_arr2,dim=2)
medflux3 = median(fiber_arr3,dim=2)

str.medsky[0].flux = medflux1
str.medsky[1].flux = medflux2
str.medsky[2].flux = medflux3


chipgap1 = 147
chipgap2 = 140

;med=median(skyframe2.flux,dim=2)
;plot,medflux2,xr=[500,1300]
;oplot,med,co=250

;stop

end
