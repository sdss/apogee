;+
; NAME:
;   apogee_aspcap_id
; PURPOSE:
;   Return the apStar ID for a reduction of a spectrum
; CALLING SEQUENCE:
;   apogee_aspcap_id, locid=, star=, results_version= [, telescope=, /commissioning ]
; INPUTS:
;   locid - location id of field star is in
;   star - 2MASS-style string with name of star
;   results_version - version of apStar combine reductions ("vN")
; OPTIONAL INPUTS:
;   telescope - which telescope (default 'n' for Northern)
; OPTIONAL KEYWORDS:
;   /commissioning - set if commissioning data 
; COMMENTS:
;   Returns a string of form:
;     apogee.[telescope].[c|s].[results_version].[locid].[star]
;   where 'telescope' is 'n' or 's' for North or South, and 
;   'c' or 's' means commissioning data or survey data.
; REVISION HISTORY:
;   M. Blanton, NYU, February 2013
;-
function apogee_aspcap_id, locid=locid, star=star, $
                           results_version=results_version, $
                           telescope=telescope, $
                           commissioning=commissioning

if(NOT keyword_set(telescope)) then $
   telescope='apo25m'

if(keyword_set(commissioning)) then $
   csstr= 'c' $
else $
   csstr= 's'

id= 'apogee.'+strtrim(telescope,2)+'.'+csstr+'.'+ $
    results_version+'.'+ $
    strtrim(string(locid),2)+'.'+strtrim(star,2)

return, id

end
