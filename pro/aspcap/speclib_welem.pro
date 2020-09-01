function speclib_welem,wave,spec,welem
;
; Given input wavelength/spectrum and welem windows, returns spectrum with only window regions
;
out=[]
wout=[]
sz=size(welem,/dim)
help,spec
stop
for i=0,sz[1]-1 do begin
  j=where(wave ge welem[0,i] and wave le welem[1,i],nj)
  if nj gt 0 then begin
    print,nj
    print,wave[j[0]],welem[0,i],wave[j[0]]-welem[0,i]
    print,wave[j[-1]],welem[1,i],wave[j[-1]]-welem[1,i]
    out=[out,spec[j]]
    wout=[out,wave[j]]
  endif
endfor

return,out
end
