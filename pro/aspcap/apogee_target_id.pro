;+
; NAME:
;   apogee_target_id
; PURPOSE:
;   Return the apStar target ID for a reduction of a spectrum
; CALLING SEQUENCE:
;   apogee_target_id, telescope=, field=, star= 
; INPUTS:
;   locid - location id of field star is in
;   star - 2MASS-style string with name of star
; COMMENTS:
;   Returns a string of form:
;     apogee.[locid].[star]
; REVISION HISTORY:
;   J. Holtzman, June 2013
;                April 2019
;-
function apogee_target_id, telescope=telescope, star=star, field=field, locid=locid

id=strtrim(string(telescope),2)+'.'+strtrim(string(locid),2)+'.'+strtrim(field,2)+'.'+strtrim(star,2) 

return, id

end
