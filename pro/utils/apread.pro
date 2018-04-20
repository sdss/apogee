function apread,type,_extra=keywords

; general read routine for APOGEE files, packs chips and extensions into structure

print,'apread: ', type
if type eq 'Star' then files=apogee_filename(type,_extra=keywords) else $
files=apogee_filename(type,_extra=keywords,chip=['a','b','c'])
print,files
ext=1

for ichip=0,n_elements(files)-1 do begin
  file=files[ichip]

  if not file_test(file) then file=file+'.fz'
  if not file_test(file) then begin
    print,'halt: Cant open file: ', file
    stop
    return,0
  endif

  d=mrdfits(file,0,head,/silent)

  if type eq 'Dark' then begin
    dark=mrdfits(file,ext,/silent)
    err=0.
    mask=mrdfits(file,ext+2,/silent)
    sz=size(dark) & nread=sz[3]
    d=dark[*,*,nread-1]-dark[*,*,1]
    if keyword_set(domask) then begin
       bad=where(mask and badmask())
       d[bad]=!values.f_nan
    endif 
    temp={hdr: head, flux: d, err: err, mask: mask} 
    if ichip eq 0 then out=temp else out=struct_append(out,temp)
  endif else if type eq 'Raw' then begin
    if not file_test(file) then apunzip,file+'.apz'
    first=mrdfits(file,2,head,/silent)
    nr=sxpar(head,'NFRAMES') 
    d=mrdfits(file,nr,/silent)
    d-=first
    err=0.
    mask=0
    temp={hdr: head, flux: d, err: err, mask: mask} 
    if ichip eq 0 then out=temp else out=struct_append(out,temp)
  endif else begin
    d=mrdfits(file,ext,/silent)
    if type eq '2Dmodel' then begin
      d=mrdfits(file,0,/silent)
      err=0.
      mask=0
    endif else  begin
      err=mrdfits(file,ext+1,/silent)
      mask=mrdfits(file,ext+2,/silent)
      if keyword_set(domask) then begin
       bad=where(mask and badmask())
       d[bad]=!values.f_nan
       err[bad]=!values.f_nan
      endif 
    endelse
    ;out[*,*,ichip]=d
    if type eq '1D' or type eq 'Cframe' or type eq 'Plate' or type eq 'Star' then begin
      wave=mrdfits(file,ext+3,/silent)
      wcoef=mrdfits(file,ext+4,/silent)
      if type eq 'Cframe' or type eq 'Plate' or type eq 'Star' then begin
        sky=mrdfits(file,ext+4,/silent)
        skyerr=mrdfits(file,ext+5,/silent)
        tell=mrdfits(file,ext+6,/silent)
        tellerr=mrdfits(file,ext+7,/silent)
        temp={hdr: head, flux: d, err: err, mask: mask, wave: wave, sky: sky, skyerr:skyerr, tell:tell, tellerr:tellerr} 
      endif else temp={hdr: head, flux: d, err: err, mask: mask, wave: wave, wcoef: wcoef} 
    endif else begin
      temp={hdr: head, flux: d, err: err, mask: mask} 
    endelse
    if ichip eq 0 then out=temp else out=struct_append(out,temp)
     
  endelse
endfor
return,out
end
