pro aprvrefine_plot,allvisits,starstr,grid,bestgrid,vout0,outdir

; This makes the RV and CCF plot for aprvrefine.pro

file_mkdir,outdir
outfile = outdir+strtrim(allvisits[0].apogee_id,2)+'_rvccf'

nvisits = n_elements(allvisits)

; Spectrum and CCF plots
file = outfile
ps_open,file,/color,thick=3,/encap
cleanplot,/silent
loadct,39
!p.font = 0
device,/inches,xsize=14,ysize=8
nspec = starstr.spec[*,0]/starstr.continuum[*,0]
temp = reform(grid.ndata[bestgrid,*])
bdt = where(temp le 0.0,nbdt)
if nbdt gt 0 then temp[bdt]=!values.f_nan
xr = minmax(starstr.wave/1e4)
yr = [0.0>min([nspec,temp],/nan)*0.9<0.5,1.4>max([nspec,temp],/nan)*1.1<2.0]
plot,[0],[0],/nodata,xr=xr,xs=1,yr=yr,ys=1,xtit='Wavelength (microns)',$
     ytit='Normalized Flux',charsize=1.3,position=[0.08,0.57,0.98,0.98],xticklen=0.04,yticklen=0.015
oplot,starstr.wave/1e4,nspec,thick=2
oplot,grid.wave/1e4,temp,co=250,thick=2
vavg = mean(allvisits.vhelio)
snr = median(starstr.spec[*,0]/starstr.err[*,0])
;xyouts,1.52,yr[1]-0.09*range(yr),allvisits[0].apogee_id+'  V!dhelio!n='+stringize(vavg,ndec=2)+' km/s  chisq='+$
;       stringize(starstr.chisq,ndec=2)+'  S/N='+stringize(snr,ndec=1)+'  Nvisits='+$
;       strtrim(nvisits,2),align=0,charsize=1.3,charthick=1.3
;xyouts,1.52,yr[1]-0.17*range(yr),'Teff='+stringize(starstr.rv_teff,ndec=0)+'  logg='+stringize(starstr.rv_logg,ndec=1)+$
;       '  [Fe/H]='+stringize(starstr.rv_feh,ndec=2)+'  [alpha/Fe]='+stringize(starstr.rv_alpha,ndec=1)+$
;       '  [C/Fe]='+stringize(starstr.rv_carb,ndec=1),align=0,charsize=1.3,charthick=1.3
line1 = allvisits[0].apogee_id+'  V!dhelio!n='+stringize(vavg,ndec=2)+' km/s  chisq='+$
       stringize(starstr.chisq,ndec=2)+'  S/N='+stringize(snr,ndec=1)+'  Nvisits='+$
       strtrim(nvisits,2)
line2 = 'Teff='+stringize(starstr.rv_teff,ndec=0)+' logg='+stringize(starstr.rv_logg,ndec=1)+$
       ' [Fe/H]='+stringize(starstr.rv_feh,ndec=2)+' [alpha/Fe]='+stringize(starstr.rv_alpha,ndec=1)+$
       ' [C/Fe]='+stringize(starstr.rv_carb,ndec=1)+$
       ' FW='+stringize(starstr.ccpfwhm)+' AUTOFW='+stringize(starstr.autofwhm)
;dum = fsc_color('white',255)  ; this seems to be needed to make the legend box background white
al_legend,[line1,line2],position=[1.515,yr[1]-0.05*range(yr)],box=1,charsize=1.3,charthick=4,/clear
loadct,39,/silent
if nvisits gt 1 then nspec=nvisits+2 else nspec=1
co = scale_vector(findgen(nspec),50,250)

; apStar-relative CCFs
if n_elements(vout0) gt 0 then begin
  ccf0 = vout0.ccf
  mx = max(ccf0,dim=1)
  nlag = n_elements(ccf0[*,0])
   ccf0 /= replicate(1,nlag)#mx
  yr = minmax(ccf0)
  yr = [yr[0]-0.05*range(yr),yr[1]+0.05*range(yr)]
  plot,vout0.cclag,ccf0[*,0],xs=1,ys=1,yr=yr,position=[0.08,0.08,0.35,0.45],xtit='Lag (pixels)',$
       ytit='Normalized CCF',tit='apStar-relative CCFs',/noerase
  for j=0,nvisits-1 do oplot,vout0.cclag,ccf0[*,j],co=co[j],thick=2
endif else begin
  plot,[0],[0],/nodata,xs=1,ys=1,xr=[-100,100],yr=[0,1],position=[0.08,0.08,0.35,0.45],xtit='Lag (pixels)',$
       ytit='Normalized CCF',tit='apStar-relative CCFs',/noerase
  xyouts,0,0.5,'No apStar CCFs.  Only ONE Visit',align=0.5,charsize=1.3
endelse

; Synth CCFs
;  Nvisits+2
ccf1 = starstr.ccf
nlag = n_elements(ccf1[*,0])
if nvisits gt 1 then begin
  mx = max(ccf1,dim=1)
  ccf1 /= replicate(1,nlag)#mx
endif else ccf1/=max(ccf1)
yr = minmax(ccf1)
yr = [yr[0]-0.05*range(yr),yr[1]+0.05*range(yr)]
plot,starstr.ccflag,ccf1[*,0],xs=1,ys=1,yr=yr,position=[0.395,0.08,0.665,0.45],xtit='Lag (pixels)',$
     ytit=' ',tit='Synthetic CCFs',/noerase
oplot,[-500,500],[0,0],linestyle=2
for j=0,nspec-1 do oplot,starstr.ccflag,ccf1[*,j],co=co[j],thick=2

; Difference Synth CCFs
if nvisits gt 1 then begin
  ccf2 = ccf1
  ccf2 -= ccf1[*,0]#replicate(1,nvisits+2)
  yr = minmax(ccf2)
  yr = [yr[0]-0.05*range(yr),yr[1]+0.05*range(yr)]
endif else begin
  ccf2 = ccf1*0
  yr = [0,1]
endelse
plot,starstr.ccflag,ccf2[*,0],xs=1,ys=1,yr=yr,position=[0.71,0.08,0.98,0.45],xtit='Lag (pixels)',$
     ytit=' ',tit='Difference Synthetic CCFs',/noerase
oplot,[-500,500],[0,0],linestyle=2
for j=0,nspec-1 do oplot,starstr.ccflag,ccf2[*,j],co=co[j],thick=2

;stop,'are we looking over all CCFs here???'

ps_close
;ps2jpg,file+'.eps',/eps
ps2gif,file+'.eps',/eps
file_delete,file+'.eps',/allow

end
