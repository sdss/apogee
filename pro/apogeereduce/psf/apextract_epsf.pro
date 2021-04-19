;+
;
; APEXTRACT_EPSF
;
; This extracts spectra using an empirical PSF
;
; INPUTS:
;  frame   The 2D input structure with FLUX, VAR and MASK
;  epsf    A structure with the empirical PSF.
;  /doback Subtract the background
;
; OUTPUTS:
;  outstr  The 1D output structure with FLUX, VAR and MASK
;  back    The background
;  =model  The model 2D image
;
; USAGE:
;  IDL>apextract_epsf,frame,epsf,outstr,model=model
;
; By J. Holtzman  2011
; Incorporated into ap2dproc.pro  D.Nidever May 2011
;-

function apextract_pmul,p1lo,p1hi,img,p2

lo=MAX([p1lo,p2.lo])
k1=lo-p1lo
l1=lo-p2.lo
hi=MIN([p1hi,p2.hi])
k2=hi-p1lo
l2=hi-p2.lo
if lo gt hi then return,fltarr(2048)
img2=*p2.img
;if lo eq hi then return,TOTAL(img[*,k1:k2]*img2[*,l1:l2],1,/nan) $
if lo eq hi then return,img[*,k1:k2]*img2[*,l1:l2] $
else return,TOTAL(img[*,k1:k2]*img2[*,l1:l2],2,/nan)
end

pro apextract_epsf,frame,epsf,outstr,back,model=model,doback=doback,skip=skip,scat=scat,subonly=subonly

;
; extract spectrum under the assumption that a given pixel only contributes
;  to two neighboring traces, leading to a tridiagonal matrix inversion
;
nframe = n_elements(frame)

; Not enough inputs
if nframe eq 0 then begin
  print,'syntax - apextract_epsf,frame,epsf,outstr,model=model,doback=doback)'
  return
endif

; back not properly implemented -- it breaks tridiagonal!
;if keyword_set(doback) then begin
;  print,' background option not implemented correctly yet!'
;  stop
;endif
;time0 = systime(/seconds)
;sz=size(trace)
;ntrace=sz[2]
;ntrace = n_elements(epsf)
ntrace = n_elements(epsf)

fibers=epsf.fiber
;if n_elements(fibers) eq 0 then fibers=indgen(ntrace)
;nfibers = n_elements(fibers)
red = frame.flux
var = frame.err^2
inmask = frame.mask

if keyword_set(scat) then scat_remove,red,scat=scat,mask=inmask

; Initialize output arrays
spec = fltarr(2048,300)
err = fltarr(2048,300)+baderr()
outmask = intarr(2048,300)+1

; calculate extraction matrix
;alpha=fltarr(ntrace,ntrace,2048)
if keyword_set(doback) then begin
  nback=1 
  back=fltarr(2048)
endif else nback=0
beta = dblarr(ntrace+nback,2048)
betavar = dblarr(ntrace+nback,2048)
psftot = dblarr(ntrace+nback,2048)
tridiag = dblarr(3,ntrace+nback,2048)
warnmasked = intarr(ntrace+nback,2048)
badmasked = intarr(ntrace+nback,2048)
inmask_warn = inmask AND warnmask()
inmask_bad = inmask AND badmask()
for k=0,ntrace-1+nback do begin
  ; Background
  if k gt ntrace-1 then begin
    beta[k,*] = TOTAL(red[*,lo:hi],2,/nan)
    betavar[k,*] = TOTAL(var[*,lo:hi],2,/nan)
    psftot[k,*] = 1.

  ; Fibers
  endif else begin

    ; get EPSF and set bad pixels to NaN
    p1 = epsf[k]
    lo = epsf[k].lo & hi = epsf[k].hi
    bad = where(finite(red[*,lo:hi]) eq 0 or red[*,lo:hi] eq 0 or $
          (inmask[*,lo:hi] AND badmask()) gt 0 , nbad)
    img = *p1.img
    if nbad gt 0 then img[bad]=!values.f_nan

    ; are there any warning flags for this trace? If so, flag the output
    for i=0,2047 do begin
     warnmasked[k,i]=inmask_warn[i,lo]
     badmasked[k,i]=inmask_bad[i,lo]
     for j=lo+1,hi do begin
       warnmasked[k,i] = warnmasked[k,i] OR inmask_warn[i,j]
       badmasked[k,i] = badmasked[k,i] OR inmask_bad[i,j]
     endfor
    endfor
    psftot[k,*] = TOTAL(img,2,/nan)
    beta[k,*] = TOTAL(red[*,lo:hi]*img,2,/nan)
    betavar[k,*] = TOTAL(var[*,lo:hi]*img^2,2,/nan)
  endelse

  ; First fiber (on the bottom edge)
  if k eq 0 then begin
    ll=1
    for l=k,k+1 do begin
      tridiag[ll,k,*] = apextract_pmul(p1.lo,p1.hi,img,epsf[l])
      ll+=1
    endfor

  ; Last fiber (on top edge)
  endif else if k eq ntrace-1 then begin
    ll=0
    for l=k-1,k do begin
      tridiag[ll,k,*] = apextract_pmul(p1.lo,p1.hi,img,epsf[l])
      ll+=1
    endfor

  ; Background terms
  endif else if k gt ntrace-1 then begin
    tridiag[1,k,*] = hi-lo+1

  ; Middle fibers (not first or last)
  endif else begin
    ll=0
    for l=k-1,k+1 do begin
      tridiag[ll,k,*] = apextract_pmul(p1.lo,p1.hi,img,epsf[l])
      ll+=1
    endfor
  endelse
endfor

;print,systime(/seconds)-time0
;print,'invert'

;if keyword_set(model) then model=red*0
for i=4,2043 do begin

  ; Good fibers
  good = where(psftot[*,i] gt 0.5,ngood,complement=bad,ncomplement=nbad)
  if nbad gt 0 then begin
    bad0=where(bad gt 0,nbad0)
    if nbad0 gt 0 then  tridiag[2,bad[bad0]-1,i]=0 
    bad1=where(bad lt ntrace-1,nbad1)
    if nbad1 gt 0 then tridiag[0,bad[bad1]+1,i]=0 
  endif
  if ngood gt 0 then begin

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
    spec[i,fibers[good]] = x
    err[i,fibers[good]] = sqrt(xvar)
    ; mask the bad pixels
    outmask[i,fibers[good]] = 0
    if nbad gt 0 then outmask[i,fibers[bad]]=maskval('NOT_ENOUGH_PSF') OR badmasked[bad,i]
    ; put the warning bits into the mask
    outmask[i,fibers]=outmask[i,fibers] OR warnmasked[*,i]
    ;var[i,good]=sqrt(xvar)

  ; No good fibers for this column
  endif else begin
    spec[i,*] = 0
    err[i,*] = baderr()
    outmask[i,fibers] = maskval('NOT_ENOUGH_PSF') OR badmasked[*,i]
    ;spec[i,*]=!values.f_nan
  endelse

  if keyword_set(doback) then back[i]=x[ngood-1]

;  if keyword_set(model) then begin
;    for k=0,ntrace-1 do begin
;      p1=psf[k]
;      lo=psf[k].lo & hi=psf[k].hi
;      img=*p1.img
;      t=spec[i,k] gt 0 ? spec[i,k] : 0
;      model[i,lo:hi]+=img[i,*]*t
;;  if k eq 0 then plot,red[i,*],xrange=[0,500],yrange=[0,30] else oplot,red[i,*]
;    endfor
;  endif

endfor

; catch any NaNs (shouldn't be there, but ....)
bad=where(finite(spec) eq 0, nbad)
if nbad gt 0 then begin
  spec[bad] = 0.
  err[bad] = baderr()
  outmask[bad] = 1
endif

; Put together the output structure
outstr = {flux:spec, err:err, mask:outmask}

; Create the Model 2D image
if arg_present(model) then begin
  model = red*0
  t=spec
  bad=where(t le 0)
  t[bad] = 0
  for k=0,ntrace-1 do begin
    nf=1
    ns=0
    if keyword_set(subonly) then junk=where(subonly eq k, nf)
    if keyword_set(skip) then junk=where(skip eq k, ns)
    if nf gt 0  and ns eq 0 then begin
      p1 = epsf[k]
      lo = epsf[k].lo & hi=epsf[k].hi
      img=*p1.img
      ;img = p1.img
      ;t=spec[*,k] gt 0 ? spec[*,k] : 0
      rows = fltarr(hi-lo+1)+1
      ;model[*,lo:hi] += img[*,*]*(rows##t[*,k])
      fiber=epsf[k].fiber
      model[*,lo:hi] += img[*,*]*(rows##t[*,fiber])
;  if k eq 0 then plot,red[i,*],xrange=[0,500],yrange=[0,30] else oplot,red[i,*]
    endif
  endfor

endif

;print,systime(/seconds)-time0

;return,spec

end
