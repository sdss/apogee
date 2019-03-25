function aspcap_clip,val,lo,hi,flag

tmp=val
bd=where(tmp lt lo,nbd)
if nbd gt 0 then begin
  tmp[bd]=lo
  for ibd=0,nbd-1 do flag[bd[ibd]]=flag[bd[ibd]] or paramflagval('CALRANGE_WARN')
endif
bd=where(val gt hi,nbd)
if nbd gt 0 then begin
  tmp[bd]=hi
  for ibd=0,nbd-1 do flag[bd[ibd]]=flag[bd[ibd]] or paramflagval('CALRANGE_WARN')
endif

return,tmp
end

