;+
; NAME:
;   apogee_target_id
; PURPOSE:
;   Return the apStar target ID for a reduction of a spectrum
; CALLING SEQUENCE:
;   apogee_target_id, locid=, star= 
; INPUTS:
;   locid - location id of field star is in
;   star - 2MASS-style string with name of star
; COMMENTS:
;   Returns a string of form:
;     apogee.[locid].[star]
; REVISION HISTORY:
;   J. Holtzman, June 2013
;-
function apogee_target_id, locid=locid, star=star, field=field

if locid eq 1 then $
id=strtrim(string(locid),2)+'.'+strtrim(field,2)+'.'+strtrim(star,2) else $
id=strtrim(string(locid),2)+'.'+strtrim(star,2)

return, id

end
