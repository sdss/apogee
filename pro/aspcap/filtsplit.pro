function filtsplit,el,maxwind=maxwind,outdir=outdir

; Read in filter (mask) file for given element
; Split it into discrete non-zero chunks
; output the chunks into {el}.wind file
; output a filter file for each chunk into {el}_{iwind}.filt files

if ~keyword_set(maxwind) then maxwind=0
if ~keyword_set(outdir) then outdir='./'

; read in master filter and wavelengths
readcol,getenv('APOGEE_DIR')+'/data/windows/filters_26042016/'+el+'.filt',mask,/silent
readcol,getenv('APOGEE_DIR')+'/data/windows/filters_26042016/'+'wave.dat',wave,/silent

; loop through mask finding non-zero windows
nfilt=0
wall=[]
weight=[]
wtot=0
if mask[0] gt 0 then w=[0]
for i=1,n_elements(mask)-1 do begin
  if mask[i] gt 0 and mask[i-1] eq 0 then begin
    nfilt=nfilt+1
    w=[i]
    wtot=mask[i]
  endif else if mask[i] eq 0 and mask[i-1] ne 0 then begin
    weight=[weight,wtot]
    w=[w,i]
    wall=[[wall],[w]]
  endif
  wtot+=mask[i]
endfor

; output the windows and corresponding filter files
sz=size(wall,/dim)
if n_elements(sz) eq 2 then nwind=sz[1] else if sz eq 0 then nwind=0 else nwind=1
print,'nwind:',nwind,maxwind
if nwind gt 1 and (nwind le maxwind or maxwind eq 0) then begin
  ; open output file
  openw,sum,outdir+el+'.wind',/get_lun

  for iwind=0,nwind-1 do begin
    openw,out,outdir+el+'_'+string(format='(i2.2)',iwind+1)+'.mask',/get_lun
    printf,sum,wave[wall[0,iwind]],wave[wall[1,iwind]],weight[iwind]
    for i=0,n_elements(mask)-1 do begin
      if i lt wall[0,iwind] or i gt wall[1,iwind] then val=0. else val=mask[i]
      printf,out,val
    endfor
    free_lun,out
  endfor
  free_lun,sum
endif

return,nwind

end


