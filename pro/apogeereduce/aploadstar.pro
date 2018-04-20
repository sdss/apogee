pro aploadstar,starfile,str,pl=pl,stp=stp,error=error

;+
;
; APLOADSTAR
;
; Load an apStar file into an IDL structure
;
; INPUTS:
;  starfile   The filename of the apStar file
;  /pl        Plot the spectrum.
;  /stp       Stop at the end of the program.
;  =error     The error message if one occurred.
;
; OUTPUTS:
;  str        An IDL structure containing the information
;               from the apStar file.
;
; USAGE:
;  IDL>aploadstar,starfile,str
;
; By D.Nidever  July 2010
;-

apgundef,str

; Do we have enough inputs
nstarfile = n_elements(starfile)
if nstarfile eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax - aploadstar,starfile,str,pl=pl,stp=stp,error=error'
  return
end

; More than one file input
if nstarfile gt 1 then begin
  error = 'MORE THAN ONE FILE INPUT'
  print,error
  return
endif

; Does the visit file exist
test = file_test(starfile)
if test eq 0 then begin
  error = 'FILE '+starfile+' NOT FOUND'
  print,error
  return
endif

; Load the data
head = headfits(starfile,exten=0,errmsg=message0)
FITS_READ,starfile,spec,spec_head,exten=1,/no_abort,message=message1
FITS_READ,starfile,err,err_head,exten=2,/no_abort,message=message2
FITS_READ,starfile,mask,mask_head,exten=3,/no_abort,message=message3
;FITS_READ,starfile,wave,wave_head,exten=4,/no_abort,message=message4
FITS_READ,starfile,sky,sky_head,exten=4,/no_abort,message=message4
FITS_READ,starfile,skyerr,skyerr_head,exten=5,/no_abort,message=message5
FITS_READ,starfile,telluric,telluric_head,exten=6,/no_abort,message=message6
FITS_READ,starfile,telerr,telerr_head,exten=7,/no_abort,message=message7
FITS_READ,starfile,lsfcoef,lsfcoef_head,exten=8,/no_abort,message=message8
rvstr = MRDFITS(starfile,9,status=rvstatus,/silent)

message = message0+message1+message2+message3+message4+message5+message6+message7+message8
if message ne '' then begin
  error = 'ERROR LOADING '+starfile+'  '+message
  print,error
endif

; Construct the wavelength array
nwcoef = sxpar(head,'NWCOEF')
;  old version
if nwcoef gt 0 then begin
  nwcoef = sxpar(head,'NWCOEF')
  wcoef = dblarr(nwcoef)
  for i=0,nwcoef-1 do wcoef[i] = sxpar(head,'WCOEF'+strtrim(i,2))
  nspec = n_elements(spec)
  wave = POLY(dindgen(nspec)+1,wcoef)  ; starts with 1
; new version
endif else begin
  crval1 = sxpar(head,'CRVAL1')
  cdelt1 = sxpar(head,'CDELT1')
  ctype1 = sxpar(head,'CTYPE1')
  crpix1 = sxpar(head,'CRPIX1',count=ncrpix1)
  if ncrpix1 eq 0 then crpix1=1
  sz = size(spec)
  npix = sz[1]
  logwave = double(crval1) + double(cdelt1)*(dindgen(npix) + 1.0d0 - double(crpix1))
  wave = 10^logwave
endelse

; Spectrum is in units of 1e-17
spec *= 1e-17
; Error also in units of 1e-17
err *= 1e-17
; sky/skyerr in units of 1e-17
sky *= 1e-17
skyerr *= 1e-17

; Make the output structure
str = {file:starfile,head:head,spec:spec,err:err,mask:mask,wavelength:wave,$
       sky:sky,skyerr:skyerr,telluric:telluric,telluricerr:telerr,$
       lsfcoef:lsfcoef,spec_head:spec_head,err_head:err_head,mask_head:mask_head,sky_head:sky_head,lsfcoef_head:lsfcoef_head}
if rvstatus eq 0 then begin
  str = create_struct(str,'rv',rvstr)
  ;tags = tag_names(rvstr)
  ;for i=0,n_elements(tags)-1 do str=create_struct(str,tags[i],rvstr.(i))
endif

; add info from header
telescop=sxpar(str.head,'TELESCOP',count=count)
if count eq 0 then telescop=''
add_tag,str,'TELESCOPE',telescop,str
add_tag,str,'FIELD',sxpar(str.head,'FIELD'),str
add_tag,str,'RA',sxpar(str.head,'RA'),str
add_tag,str,'DEC',sxpar(str.head,'DEC'),str
add_tag,str,'GLON',sxpar(str.head,'GLON'),str
add_tag,str,'GLAT',sxpar(str.head,'GLAT'),str
add_tag,str,'J',float(sxpar(str.head,'J')),str
add_tag,str,'J_ERR',float(sxpar(str.head,'J_ERR')),str
add_tag,str,'H',float(sxpar(str.head,'H')),str
add_tag,str,'H_ERR',float(sxpar(str.head,'H_ERR')),str
add_tag,str,'K',float(sxpar(str.head,'K')),str
add_tag,str,'K_ERR',float(sxpar(str.head,'K_ERR')),str
add_tag,str,'APOGEE_TARGET1',sxpar(str.head,'TARG1'),str
add_tag,str,'APOGEE_TARGET2',sxpar(str.head,'TARG2'),str
add_tag,str,'APOGEE_TARGET3',sxpar(str.head,'TARG3'),str
add_tag,str,'AK_TARG',sxpar(str.head,'AKTARG'),str
add_tag,str,'AK_WISE',sxpar(str.head,'AKWISE'),str
add_tag,str,'AK_TARG_METHOD',sxpar(str.head,'AKMETHOD'),str
add_tag,str,'SFD_EBV',sxpar(str.head,'SFD_EBV'),str
add_tag,str,'NVISITS',sxpar(str.head,'NVISITS'),str
add_tag,str,'COMBTYPE',sxpar(str.head,'COMBTYPE'),str
add_tag,str,'VHELIO',sxpar(str.head,'VHELIO'),str
add_tag,str,'VSCATTER',sxpar(str.head,'VSCATTER'),str
add_tag,str,'VERR',sxpar(str.head,'VERR'),str
add_tag,str,'VERR_MED',sxpar(str.head,'VERR_MED'),str
add_tag,str,'OBSVHELIO',sxpar(str.head,'OVHELIO'),str
add_tag,str,'OBSVSCATTER',sxpar(str.head,'OVSCAT'),str
add_tag,str,'OBSVERR',sxpar(str.head,'OVERR'),str
add_tag,str,'OBSVERR_MED',sxpar(str.head,'OVERR_ME'),str
add_tag,str,'SYNTHVHELIO',sxpar(str.head,'SVHELIO'),str
add_tag,str,'SYNTHVSCATTER',sxpar(str.head,'SVSCAT'),str
add_tag,str,'SYNTHVERR',sxpar(str.head,'SVERR'),str
add_tag,str,'SYNTHVERR_MED',sxpar(str.head,'SVERR_ME'),str
add_tag,str,'SNR',sxpar(str.head,'SNR'),str
add_tag,str,'CHISQ', 1.0*sxpar(str.head,'CHISQ'),str
add_tag,str,'RV_TEFF',sxpar(str.head,'RVTEFF'),str
add_tag,str,'RV_LOGG',sxpar(str.head,'RVLOGG'),str
add_tag,str,'RV_FEH',sxpar(str.head,'RVFEH'),str
add_tag,str,'RV_ALPHA',sxpar(str.head,'RVALPH'),str
add_tag,str,'RV_CARB',sxpar(str.head,'RVCARB'),str
add_tag,str,'CCPFWHM',sxpar(str.head,'CCPFWHM'),str
add_tag,str,'AUTOFWHM',sxpar(str.head,'AUTOFWHM'),str
add_tag,str,'SYNTHSCATTER',sxpar(str.head,'SYNTHSCA'),str
add_tag,str,'STARFLAG',sxpar(str.head,'STARFLAG'),str
add_tag,str,'ANDFLAG',sxpar(str.head,'ANDFLAG'),str

; Plot the spectrum
if keyword_set(pl) then begin
  objid = strtrim(sxpar(str.head,'OBJID'),2)
  vhelio = strtrim(string(sxpar(str.head,'VHELIO'),format='(F10.2)'),2)
  verr = strtrim(string(sxpar(str.head,'VERR'),format='(F10.2)'),2)
  gap1beg = sxpar(str.head,'GAP1BEG')
  gap1end = sxpar(str.head,'GAP1END')
  gap2beg = sxpar(str.head,'GAP2BEG')
  gap2end = sxpar(str.head,'GAP2END')
  xr = minmax(str.wavelength)
  yr = [0.0, max(str.spec/1e-17)*1.05]
  plot,[0],[0],/nodata,xr=xr,yr=yr,xs=1,ys=1,xtit='Wavelength (Ang)',ytit='Flux (10!u-17!n ergs/s/cm^2/Ang)',$
       tit='OBJID='+objid+' Vhelio='+vhelio+'+/-'+verr+' km/s'
  ; Chip a
  oplot,str.wavelength[0:gap1beg-2],str.spec[0:gap1beg-1]/1e-17
  oplot,str.wavelength[0:gap1beg-2],str.err[0:gap1beg-1]/1e-17,co=250
  ; Chip b
  oplot,str.wavelength[gap1end:gap2beg-2],str.spec[gap1end:gap2beg-2]/1e-17
  oplot,str.wavelength[gap1end:gap2beg-2],str.err[gap1end:gap2beg-2]/1e-17,co=250
  ; Chip c
  oplot,str.wavelength[gap2end:*],str.spec[gap2end:*]/1e-17
  oplot,str.wavelength[gap2end:*],str.err[gap2end:*]/1e-17,co=250
endif

if keyword_set(stp) then stop

end
