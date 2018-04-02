function apogee_filename, filetype, chip=chip, num=num, base=base, dir=dir, nochip=nochip, _EXTRA=keywords
;
; Function to return full APOGEE file name given input file types, using sdss_filename,
;    but allowing for the possibility of multiple chips, and adding derived keyword information
;    as needed for sdss_filename
;

common com_sdss_filename, config

if n_elements(num) gt 0 and ~keyword_set(mjd) then mjd=string(long(num)/10000 + 55562L)

dirs=getdir()
prefix=dirs.prefix
if filetype eq 'R' then sdsstype = dirs.prefix+filetype else sdsstype = 'ap'+filetype
; special filenames for 1m (which doesn't use plates)
if filetype eq 'R' and dirs.telescope eq 'apo1m' then sdsstype='apR-1m'
if filetype eq 'Plan' and dirs.telescope eq 'apo1m' then sdsstype='apPlan-1m'
if filetype eq 'PlateSum' and dirs.telescope eq 'apo1m' then sdsstype='apPlateSum-1m'
if filetype eq 'Visit' and dirs.telescope eq 'apo1m' then sdsstype='apVisit-1m'
if filetype eq 'VisitSum' and dirs.telescope eq 'apo1m' then sdsstype='apVisitSum-1m'
if filetype eq 'Tellstar' and dirs.telescope eq 'apo1m' then sdsstype='apTellstar-1m'
if filetype eq 'Cframe' and dirs.telescope eq 'apo1m' then sdsstype='apCframe-1m'
if filetype eq 'Plate' and dirs.telescope eq 'apo1m' then sdsstype='apPlate-1m'

if n_elements(keywords) gt 0 then if tag_exist(keywords,'plate') then keywords=create_struct(keywords,'FIELD',apogee_field(0,keywords.plate))

if keyword_set(nochip) then chip='-asdf'

if n_elements(chip) le 1 then begin
  files=sdss_filename(sdsstype,prefix=prefix,chip=chip,num=num,mjd=mjd,apred=dirs.apred,telescope=dirs.telescope,instrument=dirs.instrument,_EXTRA=keywords)
  if keyword_set(nochip) then begin
    j=strpos(files,'-asdf')
    files=strmid(files,0,j-1)+strmid(files,j+5,strlen(files))
  endif
  if keyword_set(base) then files=file_basename(files) else if keyword_set(dir) then files=file_dirname(files)+'/' else files=files
endif else begin
  files=[]
  for i=0,n_elements(chip)-1 do begin
    file=sdss_filename(sdsstype,prefix=prefix,chip=chip[i],num=num,mjd=mjd,apred=dirs.apred,telescope=dirs.telescope,instrument=dirs.instrument,_EXTRA=keywords)
    if keyword_set(base) then files=[files,file_basename(file)] else if keyword_set(dir) then files=[files,file_dirname(file)+'/'] else files=[files,file]
  endfor
endelse

return,files
end


