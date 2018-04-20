function get_nightlog,cmjd
;
; routine to return web addess of night log for specified MJD
;
dirs=getdir()
if dirs.instrument eq 'apogee-s' then begin
  spawn,"grep -Ei '<TITLE> \[lco-operations ([0-9]+)\] lco/apogee-2s night log [ms]jd ([0-9]+)' "+getenv('SAS_ROOT')+"/lco_staging/reports/*.html | grep "+cmjd+" | awk -F'html:' '{print $1}'",out
  return,'https://data.sdss.org/sas/sdsswork/lco_staging/reports/'+file_basename(out)+'html'
endif else begin
  spawn,"grep -Ei 'subject: 2\.5m obslog ([0-9]+) \([ms]jd ([0-9]+)\)' "+getenv('SAS_ROOT')+"/apo_staging/reports/*.log | grep "+cmjd+" | awk -F'log:' '{print $1}'",out
  return,'https://data.sdss.org/sas/sdsswork/apo_staging/reports/'+file_basename(out)+'log'
endelse

end

