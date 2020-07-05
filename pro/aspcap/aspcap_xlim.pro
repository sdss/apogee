function aspcap_xlim,x,xmin,xmax,flag,flagval

  ret=x
  j=where(ret lt xmin, nbd)
  if nbd gt 0 then begin
    ret[j]=xmin
    for jj=0,n_elements(j)-1 do flag[j[jj]]=flag[j[jj]] or flagval
  endif
  j=where(ret gt xmax, nbd)
  if nbd gt 0 then begin
    ret[j]=xmax
    for jj=0,n_elements(j)-1 do flag[j[jj]]=flag[j[jj]] or flagval
  endif

  return,ret
end

