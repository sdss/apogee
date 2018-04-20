function normalize_sincos,wave,spec,nwave=nwave,L=l,mask=mask,design=design

if ~keyword_set(nwave) then nwave=3
if ~keyword_set(L) then L=1400

if ~keyword_set(design) then begin
  des=fltarr(n_elements(spec),2*nwave+1)
  help,des
  for i=0,nwave do begin
    if i eq 0 then begin
      des[*,i] = cos(2*!pi*i*wave/L)
    endif else begin
      des[*,i*2-1] = sin(2*!pi*i*wave/L)
      des[*,i*2] = cos(2*!pi*i*wave/L)
    endelse
  endfor
endif

if ~keyword_set(mask) then mask=replicate(1,n_elements(spec))
gd=where(mask gt 0)
print,n_elements(gd)
alpha=matrix_multiply(des[gd,*],des[gd,*],/atranspose)
beta=matrix_multiply(des[gd,*],spec[gd],/atranspose)
par = la_linear_equation(alpha,beta)

out=par[0]
for i=1,nwave do begin
  out += par[i*2-1]*sin(2*!pi*i*wave/L)
  out += par[i*2]*cos(2*!pi*i*wave/L)
endfor

return,out

end
