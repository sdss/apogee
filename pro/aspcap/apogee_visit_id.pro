;+
; NAME:
;   apogee_visit_id
; PURPOSE:
;   Return the visit ID for a reduction of a plate or spectrum
; CALLING SEQUENCE:
;   apogee_visit_id, plate=, mjd=, fiberid=, apred_version= [, telescope=, /commissioning ]
; INPUTS:
;   plate, mjd - which spectrum 
;   apred_version - version of visit reductions ("rN")
; OPTIONAL INPUTS:
;   fiberid - which fiber (if this is a spectrum visit id) 
;   telescope - which telescope (default 'n' for Northern)
; OPTIONAL KEYWORDS:
;   /commissioning - set if commissioning data
; COMMENTS:
;   Returns a string of form:
;     apogee.[telescope].[c|s].[apred_version].plate.mjd.fiberid
;   where 'telescope' is 'n' or 's' for North or South, and 
;   'c' or 's' means commissioning data or survey data.
;   If fiberid is omitted, then the id refers to a plate as a 
;   whole and is returned as:
;     apogee.[telescope].[c|s].[apred_version].plate.mjd
; REVISION HISTORY:
;   M. Blanton, NYU, February 2013
;-
function apogee_visit_id, plate=plate, mjd=mjd, fiberid=fiberid, file=file, $
                          apogeeid=apogeeid, $
                          apred_version=apred_version, telescope=telescope, $
                          commissioning=commissioning

if(NOT keyword_set(telescope)) then $
   telescope='apo25m'

if(keyword_set(commissioning)) then $
   csstr= 'c' $
else $
   csstr= 's'

; special handling for multiple observations on same date for apo1m
if telescope eq 'apo1m' and  keyword_set(file) then begin
   tmp=strsplit(file,'-',/extract)
   cmjd=tmp[2]
endif else cmjd=strtrim(string(mjd),2)

id= 'apogee.'+telescope+'.'+csstr+'.'+ $
    apred_version+'.'+ $
    strtrim(string(plate),2)+'.'+  $
    cmjd

; for apo1m, use apogeeid instead of fiberid
if telescope eq 'apo1m' and keyword_set(apogeeid) then $
   id=id+ '.'+ strtrim(apogeeid,2) else $
if (keyword_set(fiberid)) then $
   id=id+ '.'+ strtrim(string(fiberid),2)

return, id

end
