function apfixmasked,spec,mask, setnan = setnan

IF keyword_set(setnan) THEN BEGIN
 new = spec
 new[where(mask)] = !Values.F_Nan
 return,new
ENDIF

istart=intarr(n_elements(mask))
iend=intarr(n_elements(mask))

nregion=-1
last=-10
for i=0,n_elements(mask)-1 do begin
  if mask[i] then begin
   if i gt last+1 then begin
     nregion+=1
     istart[nregion]=i
     iend[nregion]=i
   endif else begin
     iend[nregion]=i
   endelse
   last=i
  endif
endfor

new=spec
for i=0,nregion do begin
 nmask=iend[i]-istart[i]+1
 x0=0
 x1=nmask+1
 y0=0 & y1=0
 ; return from median to immediate pixels based on poorer proto-DR13 results for vscatter c.f. DR12
 if istart[i]-1 ge 0 then y0=spec[istart[i]-1]                                   ;These lines commented out
 if iend[i]+1 lt n_elements(spec)-1 then y1=spec[iend[i]+1]                      ;in proto-dr13
 ;if istart[i] eq 0 then new[istart[i]:iend[i]] = y1 else $                        ;Never
 ;if iend[i] eq n_elements(spec)-1 then new[istart[i]:iend[i]] = y0 else begin     ;Used
 ;if istart[i]-10 ge 0 then y0=median(spec[istart[i]-10:istart[i]-1])              ; These lines
 ;if iend[i]+10 lt n_elements(spec)-1 then y1=median(spec[iend[i]+1:iend[i]+10])   ; used for proto-dr13 
 if istart[i] eq 0 then new[istart[i]:iend[i]] = y1 else $                         ; These lines used 
 if iend[i] eq n_elements(spec)-1 then new[istart[i]:iend[i]] = y0 else begin      ; in both
   x=indgen(nmask)+1
   slope=(y1-y0)/(x1-x0)
   new[istart[i]:iend[i]] = y0+x*slope
 endelse
;print,i,nmask,istart[i],iend[i]
endfor
return,new
end
