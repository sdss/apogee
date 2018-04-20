pro scat_remove,a,scat=scat,mask=mask

; remove scattered light

if scat eq 1 then begin
  ; simple stupid single level removal!
  if keyword_set(mask) then begin
    flux=a
    j=where((mask and badmask()) gt 0)
    flux[j]=!values.f_nan
  endif else flux=a
  bot=median(flux[100:1947,5:10])
  top=median(flux[100:1947,2038:2042])
  scatlevel=(bot+top)/2.
  print,'scatlevel: ',scatlevel
  a-=scatlevel

endif else begin

  ; variable scattered light, but only works for sparse exposures
  sz=size(a,/dim)
  t=a
  bad=where(finite(t) eq 0 or t lt -10)
  t[bad]=1e10
  nbox=51
  grid=fltarr(41,41)
  ii=0
  for i=4,2044,nbox do begin
  print,i
    jj=0
    for j=4,2044,nbox do begin
      i1=i-nbox/2 & i2=i+nbox/2
      j1=j-nbox/2 & j2=j+nbox/2
      i1=max([4,i1]) & i2=min([2044,i2])
      j1=max([4,j1]) & j2=min([2044,j2])
      sky=t[i1:i2,j1:j2]
      mmm,sky,val,sig,skew,highbad=1e5
      if sig gt 0 then grid[ii,jj]=val
      jj+=1
    endfor
    ii+=1
  endfor
  
  vec1=indgen(nbox)
  vec2=replicate(1.,nbox)
  xramp=vec1#vec2
  yramp=vec1##vec2
  
  w1=(nbox-xramp)/nbox*(nbox-yramp)/nbox
  w2=xramp/nbox*(nbox-yramp)/nbox
  w3=(nbox-xramp)/nbox*yramp/nbox
  w4=xramp/nbox*yramp/nbox
  
  out=fltarr(2048,2048)
  ii=0
  for i=4+nbox/2,2044-nbox/2,nbox do begin
   jj=0
   for j=4+nbox/2,2044-nbox/2,nbox do begin
     v1=grid[ii,jj]
     v2=grid[ii+1,jj]
     v3=grid[ii,jj+1]
     v4=grid[ii+1,jj+1]
     if v1 gt 1e9 then v1=v2
     if v2 gt 1e9 then v2=v1
     out[i-nbox/2:i+nbox/2,j-nbox/2:j+nbox/2]=v1*w1+v2*w2+v3*w3+v4*w4
     jj+=1
   endfor
   ii+=1
  endfor
  ;stop
  a-=out

endelse

end
