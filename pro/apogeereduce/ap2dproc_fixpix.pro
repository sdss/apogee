pro ap2dproc_fixpix,chstr

; Fix the Bad pixels and "Unfixable" pixels

sz = size(chstr.flux)
bpmmask = long((long(chstr.mask) AND 1) eq 1)  ; bad pixels
crmask = long((long(chstr.mask) AND 2) eq 2)   ; CR pixels
satmask = long((long(chstr.mask) AND 4) eq 4)  ; sat pixels
unfmask = long((long(chstr.mask) AND 8) eq 8)  ; unfixable
bdpix = where(bpmmask eq 1 or unfmask eq 1,nbdpix)
flux0 = chstr.flux
print,'Fixing ',strtrim(nbdpix,2),' bad pixels'
; Loop over the bad pixels
for j=0L,nbdpix-1 do begin
  ;print,j+1
  bdpix2d = array_indices(bpmmask,bdpix)
  halfwid = 5  ; 3
  xlo = bdpix2d[0,j]-halfwid > 0
  xhi = bdpix2d[0,j]+halfwid < (sz[1]-1)
  ylo = bdpix2d[1,j]-halfwid > 0
  yhi = bdpix2d[1,j]+halfwid < (sz[2]-1)
  ;bdxcen = bdpix2d[0,j]-xlo
  ;bdycen = bdpix2d[1,j]-ylo

  gdpix = where(bpmmask[xlo:xhi,ylo:yhi] eq 0 and unfmask[xlo:xhi,ylo:yhi] eq 0 and $
                satmask[xlo:xhi,ylo:yhi] eq 0,ngdpix)
  if ngdpix gt 4 then begin
    z = (chstr.flux[xlo:xhi,ylo:yhi])[gdpix]
    x = ( (findgen(xhi-xlo+1)+xlo)#replicate(1,yhi-ylo+1) )[gdpix]
    y = ( replicate(1,xhi-xlo+1)#(findgen(yhi-ylo+1)+ylo) )[gdpix]
    if range(x) gt 0 and range(y) gt 0 then begin  ; can't be co-linear

      ; Do the interpolation using the good pixels
      result = tri_surf(z,x,y,gs=[1,1],bounds=[bdpix2d[0,j],bdpix2d[1,j],bdpix2d[0,j]+1,bdpix2d[1,j]+1])
      chstr.flux[bdpix[j]] = result[0,0]

      ;displayc,chstr.flux[xlo:xhi,ylo:yhi],/z,tit=strtrim(j+1,2)+' X='+strtrim(bdpix2d[0,j],2)+' '+strtrim(bdpix2d[1,j],2)
      ;stop
    endif
  endif

end ; bad pixel loop

;stop

end
