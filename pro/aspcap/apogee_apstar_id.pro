;+
; NAME:
;   apogee_apstar_id
; PURPOSE:
;   Return the apStar ID for a reduction of a spectrum
; CALLING SEQUENCE:
;   apogee_apstar_id, locid=, star=, apstar_version= [, telescope=, /commissioning ]
; INPUTS:
;   locid - location id of field star is in
;   star - 2MASS-style string with name of star
;   apstar_version - version of apStar combine reductions ("sN")
; OPTIONAL INPUTS:
;   telescope - which telescope (default 'n' for Northern)
; OPTIONAL KEYWORDS:
;   /commissioning - set if commissioning data 
; COMMENTS:
;   Returns a string of form:
;     apogee.[telescope].[c|s].[apstar_version].[locid].[star]
;   where 'telescope' is 'apo25m' or 'apo1m' or 'lco25m', and 
;   'c' or 's' means commissioning data or survey data.
; REVISION HISTORY:
;   M. Blanton, NYU, February 2013
;-
function apogee_apstar_id, locid=locid, star=star, $
                           apstar_version=apstar_version, $
                           telescope=telescope, $
                           commissioning=commissioning

if(NOT keyword_set(telescope)) then $
   telescope='apo25m'

if(keyword_set(commissioning)) then $
   csstr= 'c' $
else $
   csstr= 's'

id= 'apogee.'+telescope+'.'+csstr+'.'+ $
    apstar_version+'.'+ $
    strtrim(string(locid),2)+'.'+strtrim(star,2)

return, id

end
