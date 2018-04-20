
pro getflux,d,skyline,rows

if skyline.w1 gt d[0].wave[2047,150] then ichip=0 $
else if skyline.w1 gt d[1].wave[2047,150] then ichip=1 $
else ichip=2

cont=fltarr(n_elements(rows))
line=fltarr(n_elements(rows))
nline=fltarr(n_elements(rows))
for i=0,n_elements(rows)-1 do begin
  wave=d[ichip].wave[*,rows[i]]
  data=d[ichip].flux[*,rows[i]]
  icont=where((wave gt skyline.c1 and wave lt skyline.c2) or $
              (wave lt skyline.c3 and wave lt skyline.c4) )
  if icont[0] ge 0 then $
    cont[i]=median(d[ichip].flux[icont,rows[i]])
  iline=where(wave gt skyline.w1 and wave lt skyline.w2)
  if iline[0] ge 0 then begin
    line[i]=TOTAL(d[ichip].flux[iline,rows[i]],/nan)
    nline[i]=TOTAL(d[ichip].flux[iline,rows[i]]/d[ichip].flux[iline,rows[i]],/nan)
  endif
endfor
skyline.flux = line - nline*cont
if skyline.type eq 0 then skyline.flux /= cont

end
                                                                                                                                        776,1         Bot

