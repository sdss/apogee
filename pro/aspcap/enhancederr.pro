function enhancederr,flux,err,sky,sig=sig,skyfact=skyfact

; enhance error around sky lines into skyerr

if not keyword_set(sig) then sigfilt=5 else sigfilt=sig
if not keyword_set(skyfact) then fact=10 else fact=skyfact

sz=size(flux)
if sz[0] gt 1 then ndim=sz[2] else ndim=1
sky=reform(sky,sz[1],ndim)
err=reform(err,sz[1],ndim)
skyerr=reform(fltarr(sz[1],ndim),sz[1],ndim)

for idim=0,ndim-1 do begin

  bderr=where(err[*,idim] eq baderr(),nbderr)

  tsky=sky[*,idim]
  terr=err[*,idim]
  ; find sky lines from peaks in sky spectrum
  lp=findpeak(tsky,level=median(tsky)+5*mad(tsky),nlp)
  ; find sky lines that are a significant fraction of stellar flux
  if nlp gt 0 then lb=where(tsky[lp]/flux[lp] gt 1.0,nlb) else nlb=0
  if nlb gt 0 then begin
    l=lp[lb]
    ;l=findpeak(sky/flux,level=0.5)
  
    ; create a gaussian around the significant sky lines
    s=fltarr(n_elements(tsky))
    i=indgen(21)-10
    s[l]=1.
    skyfilt=convol(s,exp(-0.5*i^2/sigfilt^2))
    ; enhance the error around the lines
    tskyerr=terr*(1+skyfilt*fact)
    if nbderr gt 0 then tskyerr[bderr]=baderr()
  endif else tskyerr=terr
  skyerr[*,idim]=tskyerr
endfor
return,skyerr
end

