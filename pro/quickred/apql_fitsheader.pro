;+
;
; APQL_FITSHEADER
;
; This program checks that the FITS header has the required keywords.
;
; INPUTS:
;  str               The quicklook structure.
;  =fitskw_err       A structure giving the required FITS keyword
;                      error codes.  This is normally obtained from
;                      the apogeeql database.
;  =required_fitskw  A structure giving the required FITS keywords.
;                      This is normally obtained from the apogeeql database.
;  /silent  Don't print anything to the screen.
;
; OUTPUTS:
;  Any errors in the FITS keywords are stored in the STR structure. 
;  =error   The error message if one occurred.
;
; USAGE:
;  IDL>apl_fitsheader,str,required_fitskw=required_fitskw
;
; By D.Nidever   2010
; updates by S.Beland   April 2011
;-
pro apql_fitsheader,str,required_fitskw=required_fitskw,fitskw_err=fitskw_err,silent=silent,error=error


; Initialize all values to BAD
str.fitsheader_status = 0
; this is a pointer
;str.fitsheader_errors = 

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
if n_elements(str) eq 0 or n_elements(required_fitskw) eq 0 or n_elements(fitskw_err) eq 0 then begin
  error = 'Not enough inputs'
  if not keyword_set(silent) then $
    print,'Syntax - apql_fitsheader,str,required_fitskw=required_fitskw,fitskw_err=fitskw_err,silent=silent,error=error'
  return
endif


nkeyword = n_elements(required_fitskw.name)
keyword_bad = intarr(nkeyword)  ; 1 indicates bad
errortype_pk = intarr(nkeyword)

for i=0,nkeyword-1 do begin
  value = sxpar(str.header,required_fitskw.name[i],count=count)
  dtype = size(value,/tname)

  ; check for errors
  if count eq 0 then begin
     p=where(strpos(fitskw_err.name,'MISSING') ge 0,count)
     if count ge 0 then errortype_pk[i]=p[0]
     keyword_bad[i]=1
  endif else if count gt 1 then begin
     p=where(strpos(fitskw_err.name,'DUPLICATE') ge 0,count)
     if count ge 0 then errortype_pk[i]=p[0]
     keyword_bad[i]=1
  endif else if dtype ne required_fitskw.dtype[i] then begin
     p=where(strpos(fitskw_err.name,'DATATYPE') ge 0,count)
     if count ge 0 then errortype_pk[i]=p[0]
     keyword_bad[i]=1
  endif else if dtype eq 'FLOAT' or dtype eq 'DOUBLE' or dtype eq 'INT' or dtype eq 'LONG' or dtype eq 'BYTE' then begin
     if double(value) lt double(required_fitskw.lowval[i]) or double(value) gt double(required_fitskw.highval[i]) then begin
        p=where(strpos(fitskw_err.name,'INVALID') ge 0,count)
        if count ge 0 then errortype_pk[i]=p[0]
        keyword_bad[i]=1
     endif
  endif

endfor  ; loop through the required keywords


;  store the list of errors as a pointer in the returned structure
; 0-good, 1-bad
p=where(keyword_bad gt 0,count)
if count gt 0 then begin
   str.fitsheader_status=1 
   if PTR_VALID(str.fitsheader_errors) then PTR_FREE,str.fitsheader_errors
   str.fitsheader_errors = PTR_NEW({name:required_fitskw.name[p], errortype_pk:errortype_pk[p]})
endif else begin
   str.fitsheader_status=0
   if PTR_VALID(str.fitsheader_errors) then PTR_FREE,str.fitsheader_errors
endelse

end
