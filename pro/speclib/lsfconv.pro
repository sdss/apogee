pro lsfconv,wave,spec,wlsf,lsf,xout,yout,highres=highres

  if keyword_set(highres) then begin
    highlsf=lsf
    sz=size(lsf,/dim)
    npix=sz[0]
    nlsf=sz[1]/highres
  endif else begin
    highres=3
    sz=size(lsf,/dim)
    npix=sz[0]
    nlsf=sz[1]
    highlsf=fltarr(npix,nlsf*highres)
    for ipix=0,n_elements(wlsf)-1 do $
      highlsf[ipix,*]=interpol(lsf[ipix,*],findgen(nlsf)-nlsf/2,$
        (findgen(nlsf*highres)-(nlsf*highres)/2)/float(highres))/float(highres)
  endelse

  disp=6.d-6/highres
  dw=alog10(wave[n_elements(wave)-1])-alog10(wave[0])
  xout=alog10(wave[0])+findgen(dw/disp)*disp
  xout*=alog(10.)
  yout=xout*0.
  nlsf2=(nlsf*highres)/2
  tmp=interpol(spec,wave,exp(xout))
  for ipix=nlsf*highres-nlsf2,n_elements(yout)-(nlsf*highres-nlsf2) do begin
    junk=min(abs(exp(xout[ipix])-wlsf),imin)
    yout[ipix-nlsf2:ipix+nlsf2]+=tmp[ipix]*highlsf[imin,*] 
  endfor

end
