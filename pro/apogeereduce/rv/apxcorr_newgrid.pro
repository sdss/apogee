pro apxcorr_newgrid,grid,spec,err,rvout,bestgrid,mask=mask,sum=sum,pl=pl
;+
; APXCORR_NEWGRID
;
; Gets cross correlation of observed spectrum with a number of
;  template spectra, and determines best fitting one
;-

wave=grid.wave
pixlim = grid.pixlim

; Make median filtered error array
obserr = err*0
for j=0,2 do begin
  lo=pixlim[0,j] & hi=pixlim[1,j] 
  IF lo GE hi THEN continue ;if any pixlim has lower limit higher than upper limit, then the chip is bad. Skip to the next chip, and use default value for obserr
  obserr1=err[lo:hi]
  bderr = where(obserr1 gt 1 or obserr1 le 0.0,nbderr,comp=gderr,ncomp=ngderr)                                      
  if ngderr gt 0 then begin
    if nbderr gt 0 then obserr1[bderr]=median([obserr1[gderr]],/even)
    obserr[lo:hi] = medfilt1d(obserr1,51,/edge)
  endif else obserr[lo:hi]=0.01
endfor

;pl=1
sz=size(grid.data)
if sz[0] gt 1 then ntemp=sz[1] else ntemp=1
for i=0,ntemp-1 do begin
  ; chip loop
  for j=0,2 do begin
    lo=pixlim[0,j] & hi=pixlim[1,j]
    IF lo GE hi THEN continue ;if any pixlim has lower limit higher than upper limit, then the chip is bad. Skip to the next chip, and use default values for template
    template = reform(grid.data[i,lo:hi])
    ; make "seed" a constant so the random numbers are reproducible
    template_with_noise = template + randomn(0,n_elements(template))*obserr[lo:hi]
    tempstr = {spec:template_with_noise, err:obserr, wave: wave[lo:hi]}
    APNORMSPEC,tempstr,/noerrcorr
    grid.ndata[i,lo:hi] = template / tempstr.continuum
  endfor

  ; do the xcorrs fast to get chisq without fitting peak
;pl=1
  apxcorr,wave,grid.ndata[i,*],spec,err,rvout,sum=sum,pixlim=pixlim,/nofit;,pl=pl
  ;apxcorr,wave,temp[i,*],spec,err,rvout,sum=sum,pl=pl,/nofit
  ;print,'template: ',i,rvout.vrel,rvout.chisq,ntemp
  if i eq 0 then allrv=rvout else allrv=[allrv,rvout]
endfor

chisqmin=min(allrv.chisq,bestgrid)
rvout=allrv[bestgrid]
; now redo the best fit one to get the peak
;pl=1
apxcorr,wave,grid.ndata[bestgrid,*],spec,err,rvout,sum=sum,pl=pl,pixlim=pixlim
if keyword_set(pl) then begin
  !p.multi=[0,0,0]
  sz=size(spec)
  if sz[0] eq 1 then begin
    plot,wave,spec,color = cgcolor('black'), background = cgcolor('white')
  endif else begin
    plot,wave,spec[*,0], color = cgcolor('black'), background = cgcolor('white')
    oplot,wave,spec[*,1]
    oplot,wave,spec[*,2]
  endelse
  oplot,wave,shift(grid.data[bestgrid,*],rvout.xshift0),color=cgcolor('red')
  stop
endif

; get autocorrelation of template
err=fltarr(n_elements(grid.ndata[bestgrid,*]))+0.001
apxcorr,wave,grid.ndata[bestgrid,*],reform(grid.ndata[bestgrid,*]),err,autoout,sum=sum,pl=pl,pixlim=pixlim
ADD_TAG,rvout,'AUTOFWHM',autoout.ccpfwhm,rvout

print,' Best Xshift = ',stringize(rvout.xshift0,ndec=2),' pixels'
print,' Best Gaussian fit Xshift = ',stringize(rvout.xshift,ndec=2),' +/- ',stringize(rvout.xshifterr,ndec=2),' pixels'
print,' Best synth parameters: [Fe/H]=',stringize(grid.metals[bestgrid],ndec=1),' Teff=',stringize(grid.teff[bestgrid],ndec=0), ' K logg=',stringize(grid.logg[bestgrid],ndec=1)
print,' Vrel = ',stringize(rvout.vrel,ndec=3),' +/- ',stringize(rvout.vrelerr,ndec=3),' km/s'


end
