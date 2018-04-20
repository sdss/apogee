pro apsetver,vers=vers,telescope=telescope,instrument=instrument
common apver,ver,telescop,instrume
;
; main routine for defining key names for file construction:
;   version, telescope, instrument
; these are passed in common for use by routine getdir

; default to apo25m/apogee-n on first call
if n_elements(telescop) eq 0 then telescop='apo25m'
if n_elements(instrume) eq 0 then instrume='apogee-n'

; set specified names
if keyword_set(vers) gt 0 then ver=vers
if keyword_set(telescope) gt 0 then begin
  telescop=telescope
  if telescope eq 'apo25m' then instrume='apogee-n'
  if telescope eq 'apo1m' then instrume='apogee-n'
  if telescope eq 'lco25m' then instrume='apogee-s'
endif
if keyword_set(instrument) gt 0 then instrume=instrument

end
