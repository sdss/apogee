;+
;
; APEXTRACT:
;
; This program extracts the flux from a 2D image if given
; a trace structure.
;
; INPUTS:
;  str       A structure that contains the flux, variance,
;              and mask 2D arrays with tags of FLUX, ERR and MASK.
;  tracestr  The trace structure as output by APFINDTRACE.PRO
;             for a flat field at the same dither position
;  =fibers   Fibers to extract.  The default is to extract all fibers.
;  /silent   Don't print anything to the screen
;  /stp      Stop at the end of the program.
;
; OUTPUTS:
;  outstr    The same type of structure as STR, but now with
;              the extracted spectra.
;  =error    The error message if one occurred.
;
; USAGE:
;  IDL>apextract,str,tracestr,outstr
;
; By D.Nidever  March 2010
;-
pro apextract,str,tracestr,outstr,stp=stp,fibers=fibers,silent=silent,error=error

; Error Handling
;------------------
; Establish error handler. When errors occur, the index of the  
; error is returned in the variable Error_status:  
CATCH, Error_status 

;This statement begins the error handler:  
if (Error_status ne 0) then begin 
   error = !ERROR_STATE.MSG  
   if not keyword_set(silent) then print,error
   CATCH, /CANCEL 
   return
endif


nstr = n_elements(str)
npeak = n_elements(tracestr)

; Not enough inputs
if nstr eq 0 or npeak eq 0 then begin
  error = 'Not enough inputs'
  if not keyword_set(silent) then $
    print,'Syntax - apextract,str,tracestr,outstr,fibers=fibers,stp=stp,silent=silent,error=error'
  return
endif


; Check that we have all the tags we need
tags = TAG_NAMES(str)
dum = where(tags eq 'FLUX',nflux)
if nflux eq 0 then begin
  error = 'FLUX tag not found'
  if not keyword_set(silent) then print,error
  return
endif
dum = where(tags eq 'ERR',nvar)
if nvar eq 0 then begin
  error = 'ERR tag not found'
  if not keyword_set(silent) then print,error
  return
endif
dum = where(tags eq 'MASK',nmask)
if nmask eq 0 then begin
  error = 'MASK tag not found'
  if not keyword_set(silent) then print,error
  return
endif

; Fibers to extract
ntraces = n_elements(tracestr)
if n_elements(fibers) eq 0 then fibers=lindgen(ntraces)  ; the default is to do all

; Make the new structure
;  everything should be identical except
;  that the FLUX/ERR/MASK should be [Npix,Npeak]
for i=0,n_elements(tags)-1 do begin
  sz = size(str.(i))
  type = size(str.(i),/type)
  if tags[i] eq 'FLUX' or tags[i] eq 'ERR' or tags[i] eq 'MASK' then begin
    arr = make_array(sz[1],n_elements(fibers),type=type)
  endif else arr=str.(i)

  if i eq 0 then begin
    outstr = CREATE_STRUCT(tags[i],arr)
  endif else begin
    outstr = CREATE_STRUCT(outstr,tags[i],arr)
  endelse
end


; The trace runs along the X-axis
szflux = size(str.flux)
nx = szflux[1]
ny = szflux[2]
x = lindgen(nx)

; Loop through the Fibers
For i=0,n_elements(fibers)-1 do begin

  fwhm = tracestr[fibers[i]].fwhm
  coef = tracestr[fibers[i]].coef
  ymid = POLY(x,coef)
  ylo = MIN(floor(ymid-fwhm) > 0)
  yhi = MAX(ceil(ymid+fwhm) < (ny-1))
  num = yhi-ylo+1

  ; Make a MASK based on the trace and FWHM
  yy = replicate(1.0,nx)#(lindgen(num)+ylo)
  ymid2d = ymid#replicate(1.0,num)
  mask = long(yy ge floor(ymid2d-fwhm) and yy le ceil(ymid2d+fwhm))


  ; Flux
  ;------
  flux = TOTAL( str.flux[*,ylo:yhi]*mask, 2)
  outstr.flux[*,i] = flux

  ; Error
  ;---------
  ;  add in quadrature
  err = sqrt( TOTAL( (str.err[*,ylo:yhi]^2)*mask, 2) )
  outstr.err[*,i] = err > 1    ; make sure it's greater than 0

  ; Flag Mask
  ;-----------
  ; Combine flags bitwise with OR
  flags = reform(str.mask[*,ylo:yhi])*mask
  outflag = lonarr(nx)
  ;for j=0,num-1 do outflag=outflag OR reform(long(im[xlo+j,*,2]))
  for j=0,num-1 do outflag=outflag OR reform(long(flags[*,j]))
  outstr.mask[*,i] = outflag

  ;stop

End

if keyword_set(stp) then stop

end
