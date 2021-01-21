pro aploadvisit,visitfile,str,error=error,stp=stp

;+
;
; APLOADVISIT
;
; Load a visit file into an IDL structure
;
; INPUTS:
;  visitfile  The filename of the apVisit file
;  =error     The error message if one occurred.
;
; OUTPUTS:
;  str        An IDL structure containing the information
;               from the visit file.
;
; USAGE:
;  IDL>aploadvisit,visitfile,str
;
; By D.Nidever  July 2010
;-

; Do we have enough inputs
nvisitfile = n_elements(visitfile)
if nvisitfile eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax - aploadvisit,visitfile,str,error=error'
  return
end

; More than one file input
if nvisitfile gt 1 then begin
  error = 'MORE THAN ONE FILE INPUT'
  print,error
  return
endif

; Does the visit file exist
test = file_test(visitfile)
if test eq 0 then begin
  error = 'FILE '+visitfile+' NOT FOUND'
  print,error
  return
endif

; Load the data
head = headfits(visitfile,exten=0,errmsg=message0)
FITS_READ,visitfile,spec,spec_head,exten=1,/no_abort,message=message1
FITS_READ,visitfile,err,err_head,exten=2,/no_abort,message=message2
FITS_READ,visitfile,mask,mask_head,exten=3,/no_abort,message=message3
FITS_READ,visitfile,wave,wave_head,exten=4,/no_abort,message=message4
FITS_READ,visitfile,sky,sky_head,exten=5,/no_abort,message=message5
FITS_READ,visitfile,skyerr,skyerr_head,exten=6,/no_abort,message=message6
FITS_READ,visitfile,telluric,telluric_head,exten=7,/no_abort,message=message7
FITS_READ,visitfile,telerr,telerr_head,exten=8,/no_abort,message=message8
FITS_READ,visitfile,wcoef,wcoef_head,exten=9,/no_abort,message=message9
FITS_READ,visitfile,lcoef,lcoef_head,exten=10,/no_abort,message=message10
FITS_READ,visitfile,fluxflam,fluxflam_head,exten=11,/no_abort,message=message11

message = message0+message1+message2+message3+message4+message5+message6+message7+message8+message9+message10
if message ne '' then begin
  error = 'ERROR LOADING '+visitfile+'  '+message
  print,error
endif


; Spectrum is in units of 1e-17
spec *= 1e-17
; Error also in units of 1e-17
err *= 1e-17
; sky/skyerr in units of 1e-17
sky *= 1e-17
skyerr *= 1e-17

str = {file:visitfile,head:head,spec:spec,spec_head:spec_head,err:err,err_head:err_head,mask:mask,$
       mask_head:mask_head,wave:wave,wave_head:wave_head,sky:sky,sky_head:sky_head,skyerr:skyerr,$
       skyerr_head:skyerr_head,telluric:telluric,telluric_head:telluric_head,telerr:telerr,$
       telerr_head:telerr_head,wcoef:wcoef,wcoef_head:wcoef_head,lcoef:lcoef,lcoef_head:lcoef_head}
if message11 eq '' then str = create_struct(str,'FLUXFLAM',fluxflam,'FLUXFLAM_HEAD',fluxflam_head)


; Fix NANs/Inf
bdnan = where(finite(str.spec) eq 0 or finite(str.err) eq 0 or str.err le 0 ,nbdnan)
if nbdnan gt 0 then begin
  str.spec[bdnan] = 0
  str.err[bdnan] = baderr()
  str.mask[bdnan] = 1
endif

; Get info from the header for the structure
ra = sxpar(str.head,'RA')
dec = sxpar(str.head,'DEC')
dateobs = sxpar(str.head,'DATE-OBS')
add_tag,str,'DATEOBS',dateobs,str
jd = date2jd(dateobs)
add_tag,str,'RA',ra,str
add_tag,str,'DEC',dec,str
add_tag,str,'JD',jd,str
jdmid = sxpar(head,'JD-MID',count=njdmid)
if njdmid gt 0 then add_tag,str,'JDMID',jdmid,str
GLACTC,ra,dec,2000.0,glon,glat,1,/deg
add_tag,str,'GLON',glon,str
add_tag,str,'GLAT',glat,str
jmag = sxpar(str.head,'J')
jmagerr = sxpar(str.head,'J_ERR')
hmag = sxpar(str.head,'H')
hmagerr = sxpar(str.head,'H_ERR')
kmag = sxpar(str.head,'K')
kmagerr = sxpar(str.head,'K_ERR')
add_tag,str,'J',jmag,str
add_tag,str,'J_ERR',jmagerr,str
add_tag,str,'H',hmag,str
add_tag,str,'H_ERR',hmagerr,str
add_tag,str,'K',kmag,str
add_tag,str,'K_ERR',kmagerr,str
add_tag,str,'RA_TARG',sxpar(str.head,'RA_TARG'),str
add_tag,str,'DEC_TARG',sxpar(str.head,'DEC_TARG'),str
add_tag,str,'AK_TARG',sxpar(str.head,'AKTARG'),str
add_tag,str,'AK_WISE',sxpar(str.head,'AKWISE'),str
add_tag,str,'AK_TARG_METHOD',sxpar(str.head,'AKMETHOD'),str
add_tag,str,'SFD_EBV',sxpar(str.head,'SFD_EBV'),str

add_tag,str,'RELFLUX',sxpar(str.head,'RELFLUX'),str
add_tag,str,'MTPFLUX',sxpar(str.head,'MTPFLUX'),str
add_tag,str,'SNR',sxpar(str.head,'SNR'),str
add_tag,str,'STARFLAG',sxpar(str.head,'STARFLAG'),str

;stop

if keyword_set(stp) then stop

end
