pro apdailycals,waves=waves,psf=psf,lsfs=lsfs

if keyword_set(waves) and not keyword_set(psf) then begin
  print,'psf keyword must be given with waves keyword'
  return
endif

dirs=getdir(apogee_dir,cal_dir,spectro_dir,apogee_vers,lib_dir)
if not keyword_set(vers) then vers=dirs.apred
if not file_test(cal_dir,/dir) then file_mkdir,cal_dir
file=cal_dir+'/dailycal.par'
openw,lun,/get_lun,file,/append

if keyword_set(waves) then begin
  sz=size(waves,/dim)
  if n_elements(sz) eq 1 then i2=0 else i2=sz[1]-1
  for i=0,i2 do begin
    printf,lun,'wave     99999 99999','   ',string(format='(i8.8)',waves[0,i]),$
       '    ',strtrim(waves[0,i],2),',',strtrim(waves[1,i],2),'   ',psf
    printf,lun,'lsf     99999 99999','   ',string(format='(i8.8)',waves[0,i]),$
       '    ',strtrim(waves[0,i],2),'   ',psf
    printf,lun,'lsf     99999 99999','   ',string(format='(i8.8)',waves[1,i]),$
       '    ',strtrim(waves[1,i],2),'   ',psf
  endfor
endif
if keyword_set(lsfs) then begin
  sz=size(lsfs,/dim)
  for i=0,sz[0]-1 do begin
    printf,lun,'lsf      99999 99999','   ',string(format='(i8.8)',lsfs[i]),$
       '    ',strtrim(lsfs[i],2),'   ',psf
  endfor
endif
free_lun,lun

; get unique entries
;spawn,'sort --key=1 --key=4 '+file+' | uniq ',outsort
;openw,lun,/get_lun,file
;for i=0,n_elements(outsort)-1 do printf,lun,outsort[i]
;free_lun,lun


end
