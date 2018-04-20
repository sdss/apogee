;+
;
; APMKPSF_EPSF
;
; This creates an empirical PSF for APMKPSF.
;
; INPUTS
;  frame    A structure containing the 2D chip image of
;            the flat field image from which to make the PSF.
;            This should include FLUX, ERR and MASK.
;  /silent  Don't print anything to the screen
;
; OUTPUTS
;  outpsf   An Ntrace structure giving the empirical PSF
;             for each trace.
;
; USAGE:
;  IDL>apmkpsf_epsf,frame,outpsf
;
; By J. Holtzman 2011
; Incorporated into APMKPSF by D. Nidever  May 2011
;-

pro apmkpsf_epsf,frame,caldir,im,ichip,silent=silent,sparseid=sparseid,dmax=dmax,sdmax=sdmax,thresh=thresh,smooth=smooth,fiberid=fiberid,scat=scat,average=average

; Not enough inputs
if n_elements(frame) eq 0 then begin
  if not keyword_set(silent) then $
    print,'Synatax - apmkpsf_epsf,frame,outfile'
  return
endif
if not keyword_set(fiberid) then fiberid=0

dirs=getdir()
chip=['a','b','c']

if not keyword_set(silent) then print,'Creating Empirical PSF'

if size(frame,/type) eq 8 then begin
  red = frame.flux 
  mask= frame.mask
  head=frame.header
endif else begin
  red=frame
  mask=fix(frame*0.)
  bd=where(finite(red) eq 0,nbd)
  if nbd gt 0 then mask[bd]=1
  mkhdr,head,0
endelse

if keyword_set(scat) then scat_remove,red,scat=scat

wind = 2
if not keyword_set(average) then average=50
bad=where((mask and badmask()) gt 0,nbad)
if nbad gt 0 then red[bad]=!values.f_nan
;psf = smooth(red,[average,1],/nan)
; use zap to filt out bad pixels and reject remnant CRs, use [50,1] median filter to construct PSF
psf = zap(red,[average,1])

time0 = systime(/seconds)
; first find the fibers using average in central region by looking for peaks
p=fltarr(2048)
for i=0,2047 do begin
  robust_mean,psf[900:1100,i],tmp
  p[i]=tmp*201
endfor
;p = mash(psf,900,1100,/horz)
if not keyword_set(smooth) then smooth=4
ps= smooth(double(p),smooth,/nan)
;high = MAX(p,/nan)
; shifted to peak in smoothed spectrum 1/13
high = MAX(ps,/nan)
if not keyword_set(thresh) then thresh=0.05
ploc = findpeak(ps,level=max([thresh*high,5000.]))
ntrace = n_elements(ploc)
print,'Found ',strtrim(ntrace,2),' traces'

; determine which fiber each trace corresponds to
if fiberid gt 0 then begin
  file=apogee_filename('ETrace',chip=chip[ichip],num=fiberid)
  ref=mrdfits(file)
  fibers=ref[1000,*]
  fiber=intarr(ntrace)
  dmin=fltarr(n_elements(fiber))
  nmatched=0
  pmid=[]
  for i=0,n_elements(fiber)-1 do begin
    ; do a 3-point quadratic fit to get center a bit better than peak
    nfit=3
    fit=poly_fit(-nfit/2+indgen(nfit),ps[ploc[i]-nfit/2:ploc[i]+nfit/2],2)
    pcen=ploc[i]-fit[1]/(2*fit[2])
    if pcen lt ploc[i]-nfit/2 or pcen gt ploc[i]+nfit/2 then pcen=ploc[i]
    dist=abs(pcen-fibers)
    ;dist=abs(ploc[i]-fibers)
    dmin[nmatched]=min(dist,imin)
    if dmin[nmatched] lt 2 then begin
      fiber[nmatched]=imin 
      pmid=[pmid,ploc[i]*1.]
      nmatched+=1
    endif else begin
      ;if can't find a match, call this trace bad, but move on
      print,'not halted: cant find corresponding fiber for trace', i, ploc[i]
      ;stop
    endelse
  endfor
  if nmatched lt ntrace then begin
    ntrace=nmatched
    dmin=dmin[0:nmatched-1]
    fiber=fiber[0:nmatched-1]
  endif
endif else begin
  fiber=indgen(ntrace)
  dmin=fltarr(n_elements(fiber))
  pmid=ploc
endelse
print,'Mean distance: ',MEAN(dmin)
if ntrace gt 300 then begin
  print, 'halted: chip ', ichip,' has ',ntrace, ' traces!'
  stop
endif

;if we have an entire bad row, set it to zero so it doesn't kill the nearby traces
for j=4,2043 do begin
  bad=where(mask[*,j] and badmask() gt 0,nbad)
  if nbad eq 2048 then psf[*,j]=0.
endfor

; now march to higher columns using locations from last column to get location in next column 
;  get centroids using pixels within cent-wind < y < cent+wind, and save in trace[ncol,ntrace]
y = indgen(2048)
trace = fltarr(2048,ntrace)
;ptmp = ploc[0:ntrace-1]*1.  ;float!  do above to only include matched traces
ptmp=pmid
for i=1024,2047-average do begin
  yp = y*psf[i,*]
  for j=0,ntrace-1 do begin
    pc = (nint(ptmp[j]) > 0) < 2047
    lo = (pc-wind) > 0
    hi = (pc+wind) < 2047
    ;trace[i,j] = TOTAL(yp[lo:hi],/nan)/TOTAL(psf[i,lo:hi],/nan)
    trace[i,j] = TOTAL(yp[lo:hi])/TOTAL(psf[i,lo:hi])
  endfor
  new = trace[i,*]
  good = where(finite(new))
  ptmp[good] = new[good]
endfor

; now march to lower columns using locations from last column to get location in next column 
;ptmp = ploc*1.  ;float!
ptmp=pmid
for i=1023,average,-1 do begin
  yp = y*psf[i,*]
  for j=0,ntrace-1 do begin
    pc = (nint(ptmp[j]) > 0) < 2047
    lo = (pc-wind) > 0
    hi = (pc+wind) < 2047
    ;trace[i,j] = TOTAL(yp[lo:hi],/nan)/TOTAL(psf[i,lo:hi],/nan)
    trace[i,j] = TOTAL(yp[lo:hi])/TOTAL(psf[i,lo:hi])
  endfor
  new = trace[i,*]
  good = where(finite(new))
  ptmp[good] = new[good]
endfor

; do polynomical fit to traces, and replace into trace
ftrace=trace*0.
for j=0,ntrace-1 do begin
  gd=where(y gt average and y lt 2047-average and finite(trace[*,j]) eq 1)
  coef = poly_fit(y[gd],trace[gd,j],2,yfit=yfit)
  ;  plot,trace[*,j],yrange=[trace[1000,j]-10,trace[1000,j]+10]
  ftrace[*,j] = coef[0]+coef[1]*y+coef[2]*y*y
  ;  oplot,trace[*,j],color=3
endfor
trace = ftrace
sxaddpar,head,'NTRACE',ntrace
sxaddpar,head,'AVGDIST',MEAN(dmin)
file=apogee_filename('ETrace',chip=chip[ichip],num=im)
mwrfits,trace,file,head,/create

; get sparse pack PSF if desired and available
if keyword_set(sparseid) then begin
  file=apogee_filename('EPSF',chip=chip[ichip],num=sparseid)
  tmp=mrdfits(file,0,phead)
  nsparse=sxpar(phead,'NTRACE')
  img=ptrarr(nsparse,/allocate_heap)
  for i=0,nsparse-1 do begin
    ptmp=mrdfits(file,i+1,/silent)
    *img[i] = ptmp.img
    s ={cent: ptmp.cent, lo: ptmp.lo, hi: ptmp.hi, img: img[i]}
    if i eq 0 then spsf=replicate(s,nsparse)
    spsf[i] = s
  endfor
endif

; assign observed PSF to separate fibers
; flux from each pixel will be split between the two neighboring traces
; with sparseid, use sparse pack image to determine the weights
; if no sparse id, use gaussian with sigma=2 to set weights
sig2 = 2*2^2
if not keyword_set(dmax) then dmax = 7
bpsf = fltarr(2048,2048,ntrace)
; Column loop
for i=0,2047 do begin
  if not keyword_set(silent) then $
    print,format='(%"PSF column: %5d\r",$)',i

  k2=0
  tmp = fltarr(2048,ntrace)
  ; Row loop
  ;for j=0,2047 do begin
  for j=4,2043 do begin
   ; if whole row is masked, skip it
   gd=where(finite(psf[i,*]),ngd)
   if ngd gt 0 then begin
   ;if finite(psf[i,j]) then begin
     ; find the two traces on either side of this pixel
     while trace[i,k2] lt j and k2 lt ntrace-1 do k2+=1
     k1 = k2-1 lt 0 ? 0 : k2-1
     ; find the distance from the surrounding traces
     ; if this is the first trace, only use one
     if k1 ne k2 then d1 = j-trace[i,k1] else d1=dmax+1
     d2 = j-trace[i,k2]
     ; determine the weight to be given to each of the neighboring traces
     w1=0 & w2=0
     if abs(d1) lt dmax or abs(d2) lt dmax then begin
      if keyword_set(sparseid) then begin
       ; find the nearest sparse trace
       d=abs(spsf.cent[i]-j)
       dmin=min(d,is)
       sd1=spsf[is].cent[i]+d1-spsf[is].lo
;       sd2=spsf[is].cent[i]-d2-spsf[is].lo
       sd2=spsf[is].cent[i]+d2-spsf[is].lo
       a=spsf[is].img
       sz=size(*a)
       if abs(d1) lt dmax and nint(sd1) ge 0 and nint(sd1) lt sz[2] then begin
         wlo=(*a)[i,fix(sd1)]
         whi=(*a)[i,fix(sd1)+1]
         w1=wlo+(sd1-fix(sd1))*(whi-wlo)
         ;w1=(*a)[i,nint(sd1)] 
       endif
       if abs(d2) lt dmax and nint(sd2) ge 0 and nint(sd2) lt sz[2] then begin
         wlo=(*a)[i,fix(sd2)]
         whi=(*a)[i,fix(sd2)+1]
         w2=wlo+(sd2-fix(sd2))*(whi-wlo)
         ;w2=(*a)[i,nint(sd2)]
       endif
      endif else begin
       if (abs(d1) lt dmax) then w1=exp(-d1^2/sig2)
       if (abs(d2) lt dmax) then w2=exp(-d2^2/sig2)
      endelse
     endif
     wtot = w1+w2
     ; add this pixel into the empirical PSFs for the neighboring traces with
     ;    appropriate weights
     if wtot gt 0 then begin
       tmp[j,k1] = psf[i,j]*w1/wtot
       tmp[j,k2] = psf[i,j]*w2/wtot
     endif
   endif
  endfor
  ;norm = TOTAL(tmp,1,/nan)
  norm = TOTAL(tmp,1)
  for k=0,ntrace-1 do begin
    bp = where(tmp[*,k] gt 0,nbp)
    if nbp gt 0 then bpsf[i,bp,k] = tmp[bp,k]/norm[k]
  endfor

endfor
;print,systime(/seconds)-time0
file=apogee_filename('EPSF',chip=chip[ichip],num=im)
mwrfits,0,file,head,/create

; Put the PSFs in the output structure
for k=0,ntrace-1 do begin
  m = TOTAL(bpsf[*,*,k],1,/nan)
  ind = where(finite(m) and m ne 0)
  i1=MIN(ind)
  i2=MAX(ind)
  if i1 ge 0 then begin
    outpsf={fiber: fiber[k], cent: trace[*,k], lo: i1, hi: i2, img: bpsf[*,i1:i2,k]}
    mwrfits,outpsf,file,/silent
  endif else print,'not halted, but bad PSF at: ',k
endfor

end

