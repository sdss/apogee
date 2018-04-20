function zap,image,width,sig=sig

if n_elements(width) ne 2 then begin
  print,'must specify two elements of width...'
  return,-1
endif

s=image
sz=size(image)
nx=sz[1] & ny=sz[2]
for i=0,nx-1 do begin
  i1=i-width[0]/2 gt 0?i-width[0]/2:0
  i2=i+width[0]/2 gt nx-1?nx-1:i+width[0]/2
;  i2=i1+width[0]-1
;  if i2 gt nx-1 then begin
;    i2=nx-1 & i1=i2-width[0]+1
;  endif
  for j=0,ny-1 do begin
    j1=j-width[1]/2 gt 0?j-width[1]/2:0
    j2=j+width[1]/2 gt ny-1?ny-1:j+width[1]/2
;    j2=j1+width[1]-1
;    if j2 gt ny-1 then begin
;      j2=ny-1 & j1=j2-width[1]+1
;    endif
    region=image[i1:i2,j1:j2]
    amed = median(region)
    if keyword_set(sig) then begin
      region[i-i1,j-j1] = !values.f_nan
      sigma = stddev(region,/nan)
      if s[i,j]-amed gt sig*sigma then s[i,j] = amed
    endif else s[i,j] = amed
  endfor
endfor
return,s
end
