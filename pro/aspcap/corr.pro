function corr,cov
corr=fltarr(7,7)
for i=0,6 do begin
  for j=0,6 do begin
   corr[i,j]=cov[i,j]/sqrt(cov[i,i])/sqrt(cov[j,j])
  endfor
endfor
return,corr
end
