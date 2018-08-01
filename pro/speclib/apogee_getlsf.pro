function apogee_getlsf,lsfid,waveid,fibers,vers=vers,wave=wave,highres=highres
;
; Routine to return 2D LSF array averaged over requested fibers for sampling highres samples per pixel
; Built-in apStar wavelength scale
;
if ~keyword_set(vers) then vers='r5'
nLSFpix=15
npix=8575
;nLSFhigh=nLSFpix*highres
nLSFhigh=(nLSFpix-1)*highres+1
lsf2d=fltarr(npix,nLSFhigh)
wave=4.179d0+indgen(npix)*6.d-6
wave=10.^wave

if size(lsfid,/type) ne 0 then lsfid=string(format='(i8.8)',lsfid)
if size(waveid,/type) ne 0 then waveid=string(format='(i8.8)',waveid)

for ichip=0,2 do begin
 gdfibers=where(fibers ge 0)
 for jfiber=0,n_elements(gdfibers)-1 do begin
  ifiber=fibers[gdfibers[jfiber]]
  if ichip eq 0 then chip='a'
  if ichip eq 1 then chip='b'
  if ichip eq 2 then chip='c'
  lsffile=apogee_filename('LSF',chip=chip,num=lsfid)
  lsfpars=mrdfits(lsffile,0)
  wavefile=apogee_filename('Wave',chip=chip,num=waveid)
  lsfwave=mrdfits(wavefile,2)
  ; use wavelength scale of first fiber to set wave->pixel conversion (to make sure all lsf2d have same pixels, tiny variation is negligible)
  pix=wave2pix(wave,lsfwave[*,ifiber[0]])
  dx=slope(pix)
  dx=[dx[0],dx]
  xlsf=replicate(1.d0,n_elements(pix))#(dindgen(nLSFhigh)-nLSFhigh/2)
  xlsf*=dx/highres#replicate(1,nLSFhigh)
  xlsf+=pix#replicate(1.d0,nLSFhigh)
  gd=where(finite(pix) eq 1,complement=bd,ncomp=nbd)
  if jfiber eq 0 then $
    lsf2d[gd,*]=lsf_gh(xlsf[gd,*],pix[gd],lsfpars[ifiber,*]) else $
    lsf2d[gd,*]+=lsf_gh(xlsf[gd,*],pix[gd],lsfpars[ifiber,*])
 endfor
endfor
lsf2d>=0
for i=0,npix-1 do lsf2d[i,*]/=total(lsf2d[i,*])

; find bad pixels and fill for convolution near the edges
tot=fltarr(npix)
for i=0,npix-1 do tot[i]=total(lsf2d[i,*])
bd=where(finite(tot) eq 0,nbd,comp=gd)
for i=0,n_elements(bd)-1 do begin
  junk=min(abs(bd[i]-gd),imin)
  lsf2d[bd[i],*]=lsf2d[gd[imin],*]
endfor

return,lsf2d

end

