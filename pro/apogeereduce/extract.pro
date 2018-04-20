;    extract    : 2D -> 1D extraction routine
;      pmul

; pmul is an auxiliary function for extraction that sums the
;   product of two PSF basis functions
function pmul,p1lo,p1hi,img,p2

lo=MAX([p1lo,p2.lo])
k1=lo-p1lo
l1=lo-p2.lo
hi=MIN([p1hi,p2.hi])
k2=hi-p1lo
l2=hi-p2.lo
if lo gt hi then return,fltarr(2048)
img2=*p2.img
if lo eq hi then $
return,TOTAL(img[*,k1:k2]*img2[*,l1:l2],1,/nan) $
else $
return,TOTAL(img[*,k1:k2]*img2[*,l1:l2],2,/nan)
end

; extract is the main empirical PSF-based extraction
function extract,red,ntrace,psf,var,back,model=model,doback=doback,skip=skip

; the basic implementation here assumes that the extraction matrix
;   is tridiagonal, i.e. a given trace is only affected by its
;   two adjacent neighbors

; back not properly implemented -- it breaks tridiagonal!
if keyword_set(doback) then begin
  print,' background option not implemented correctly yet!'
  stop
endif
time0=systime(/seconds)
spec=fltarr(2048,ntrace)
err=fltarr(2048,ntrace)

; calculate extraction matrix
;alpha=fltarr(ntrace,ntrace,2048)
if keyword_set(doback) then begin
  nback=1 
  back=fltarr(2048)
endif else nback=0
beta=fltarr(ntrace+nback,2048)
betavar=fltarr(ntrace+nback,2048)
psftot=fltarr(ntrace+nback,2048)
tridiag=fltarr(3,ntrace+nback,2048)
; loop over all traces and load least squares matrices
;   beta[k] = sum_i (y_i * PSF_k)
;   alpha[k,l] = sum_i (PSF_k * PSF_l)  but stored as 3 vectors for tridiagonal case
; we are getting matrices for all columns simultaneously
for k=0,ntrace-1+nback do begin
  if k gt ntrace-1 then begin
    beta[k,*]=TOTAL(red[*,lo:hi],2,/nan)
    betavar[k,*]=TOTAL(var[*,lo:hi],2,/nan)
    psftot[k,*]=1.
  endif else begin
    p1=psf[k]
    lo=psf[k].lo & hi=psf[k].hi
    bad=where(finite(red[*,lo:hi]) eq 0)
    img=*p1.img
    if bad[0] ge 0 then begin
      img[bad]=!values.f_nan
    endif
    psftot[k,*]=TOTAL(img,2,/nan)
    beta[k,*]=TOTAL(red[*,lo:hi]*img,2,/nan)
    betavar[k,*]=TOTAL(var[*,lo:hi]*img^2,2,/nan)
  endelse
  if k eq 0 then begin
    ; first trace
    ll=1
    for l=k,k+1 do begin
      tridiag[ll,k,*]=pmul(p1.lo,p1.hi,img,psf[l])
      ll+=1
;      alpha[k,l,*]=pmul(p1.lo,p1.hi,img,psf[l])
    endfor
  endif else if k eq ntrace-1 then begin
    ; last trace
    ll=0
    for l=k-1,k do begin
      tridiag[ll,k,*]=pmul(p1.lo,p1.hi,img,psf[l])
      ll+=1
;      alpha[k,l,*]=pmul(p1.lo,p1.hi,img,psf[l])
    endfor
  endif else if k gt ntrace-1 then begin
    ; background
    tridiag[1,k,*]=hi-lo+1
  endif else begin
    ; all other traces
    ll=0
    for l=k-1,k+1 do begin
      tridiag[ll,k,*]=pmul(p1.lo,p1.hi,img,psf[l])
      ll+=1
;      alpha[k,l,*]=pmul(p1.lo,p1.hi,img,psf[l])
    endfor
  endelse
endfor
print,systime(/seconds)-time0

; now get the solution for each column
print,'invert'
if keyword_set(model) then model=red*0
for i=4,2043 do begin
;  if i mod 10 eq 0 then print,'invert/extract ',i
  good=where(psftot[*,i] gt 0.1)
  if good[0] ge 0 then begin
    ngood=n_elements(good)
;    alphagood=alpha[good,good,i]
;    betagood=beta[good,i]
;    eps=invert(alphagood)
;    spec[i,good]=eps##betagood

;    a=reform(tridiag[0,good,i])
;    b=reform(tridiag[1,good,i])
;    c=reform(tridiag[2,good,i])
;    v=reform(beta[good,i])
;    spec[i,good]=trisol(a,b,c,v,/double)

    ; solve tridiagonal matrix
    a=tridiag[0,good,i]
    b=tridiag[1,good,i]
    c=tridiag[2,good,i]
    v=beta[good,i]
    vvar=betavar[good,i]
    for j=1,ngood-1 do begin
       m=a[j]/b[j-1]
       b[j]=b[j]-m*c[j-1]
       v[j]=v[j]-m*v[j-1]
       vvar[j]=vvar[j]+m^2*vvar[j-1]
    endfor
    x=fltarr(ngood)
    xvar=fltarr(ngood)
    x[ngood-1]=v[ngood-1]/b[ngood-1]
    xvar[ngood-1]=vvar[ngood-1]/b[ngood-1]^2
    for j=ngood-2,0,-1 do begin
      x[j]=(v[j]-c[j]*x[j+1])/b[j]
      xvar[j]=(vvar[j]+c[j]^2*xvar[j+1])/b[j]^2
    endfor
    spec[i,good]=x
    err[i,good]=sqrt(xvar)
  endif else begin
    spec[i,*]=!values.f_nan
  endelse

  if keyword_set(doback) then back[i]=x[ngood-1]

endfor

; make 2D model image if requested
if keyword_set(model) then begin
  if not keyword_set(skip) then skip=-1
  t=spec
  bad=where(t le 0)
  t[bad] = 0
  for k=0,ntrace-1 do begin
    s=where(skip eq k)
    if skip[0] eq -1 or s eq -1 then begin 
      p1=psf[k]
      lo=psf[k].lo & hi=psf[k].hi
      img=*p1.img
      rows=fltarr(hi-lo+1)+1
      model[*,lo:hi]+=img[*,*]*(rows##t[*,k])
    endif else print,'skipping ',k
  endfor
endif
print,systime(/seconds)-time0
var=err

return,spec

end

