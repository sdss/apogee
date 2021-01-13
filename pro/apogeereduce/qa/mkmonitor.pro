pro mkmonitor, term=term, read=read, nofiber=nofiber, vers=vers

apsetver,vers=vers

; Get long term trends from internal cals, using compiled apQAcal.fits files, make plots
;  and web page

if keyword_set(term) then hard=0 else hard=1
if keyword_set(read) then read=1 else read=0
if hard eq 1 then set_plot,'PS' else set_plot,'X'
if hard eq 1 then smcolor,/ps else smcolor

dirs=getdir(a,c,specdir)
sdir=specdir+'/monitor/'+dirs.instrument+'/'
mdir='monitor/'+dirs.instrument+'/'
file_mkdir,sdir
if read eq 1 then begin
  allcal=mrdfits(specdir+'/'+dirs.instrument+'Cal.fits',1)
  alldark=mrdfits(specdir+'/'+dirs.instrument+'Cal.fits',2)
  allexp=mrdfits(specdir+'/'+dirs.instrument+'Exp.fits',1)
  allsci=mrdfits(specdir+'/'+dirs.instrument+'Sci.fits',1)
endif else begin

  ; Append together the individual summary files
  files=file_search(specdir+'/cal/'+dirs.instrument+'/*/*QAcal*.fits')
  if files[0] eq '' then stop,'No files! do you have correct version set?' else begin
    for i=0,n_elements(files)-1 do begin
    print,files[i]
      a=mrdfits(files[i],1)
      if i eq 0 then allcal=a else allcal=[allcal,a]
    endfor
    mwrfits,allcal,specdir+'/'+dirs.instrument+'Cal.fits',/create
  endelse
  files=file_search(specdir+'/cal/'+dirs.instrument+'/*/*QAdarkflat*.fits')
  if files[0] eq '' then stop,'No files! do you have correct version set?' else begin
    for i=0,n_elements(files)-1 do begin
    print,files[i]
      a=mrdfits(files[i],1)
      if i eq 0 then alldark=a else alldark=[alldark,a]
    endfor
    mwrfits,alldark,specdir+'/'+dirs.instrument+'Cal.fits'
  endelse

  ; Get long term trends from dome flats
  ; Append together the individual summary files
  files=file_search(specdir+'/exposures/'+dirs.instrument+'/*/5*exp.fits')
  for i=0,n_elements(files)-1 do begin
    a=mrdfits(files[i],1)
    if i eq 0 then allexp=a else allexp=[allexp,a]
  endfor
  mwrfits,allexp,specdir+'/'+dirs.instrument+'Exp.fits',/create
  ; zeropoints
  files=file_search(specdir+'/visit/'+dirs.telescope+'/*/[1-9]*/*/'+dirs.prefix+'PlateSum*.fits')
  if files[0] eq '' then stop,'No files! do you have correct version set?' else begin
    for i=0,n_elements(files)-1 do begin
      a=mrdfits(files[i],1)
      if i eq 0 then allsci=a else allsci=[allsci,a]
    endfor
    mwrfits,allsci,specdir+'/'+dirs.instrument+'Sci.fits',/create
  endelse
  
endelse

openw,html,specdir+'/'+dirs.instrument+'-monitor.html',/get_lun
printf,html,'<HTML><BODY>'
printf,html,'<ul>'
printf,html,'<li> Throughput / lamp monitors'
printf,html,'<ul>'
printf,html,'<li> <a href=#quartz> Cal channel quartz</a>'
printf,html,'<li> <a href=monitor/fiber/fiber.html>Individual fiber throughputs</a>'
printf,html,'<li> <a href=#tharflux> Cal channel ThAr</a>'
printf,html,'<li> <a href=#uneflux> Cal channel UNe</a>'
printf,html,'<li> <a href=#dome>Dome flats</a>'
printf,html,'<li> <a href=#zero>Plate zeropoints</a>'
printf,html,'</ul>'
printf,html,'<li> Positions'
printf,html,'<ul>'
printf,html,'<li> <a href=#tpos>ThAr line position</a>'
printf,html,'</ul>'
printf,html,'<li> Line widths'
printf,html,'<ul>'
printf,html,'<li> <a href=#tfwhm>ThAr line FWHM</a>'
printf,html,'</ul>'
printf,html,'<li> <a href=#trace>Trace locations</a>'
printf,html,'<li> <a href=#detectors> Detectors'
printf,html,'<li> <a href=#sky> Sky brightness'
printf,html,'</ul>'

; find the different lamp types
t=where(allcal.thar eq 1)
u=where(allcal.une eq 1)
q=where(allcal.qrtz eq 1)

; Median quartz brightness
if hard eq 1 then begin
  device,file=sdir+'qflux.eps',/encap,/color,xsize=20,ysize=8,/in
  !p.multi=[0,1,3]
  printf,html,'<h3> <a name=qflux></a> Quartz lamp median brightness (per 10 reads) in extracted frame </h3>'
endif
for ichip=0,2 do begin
  if dirs.instrument eq 'apogee-s' then yr=[-1000,100000] else yr=[-1000,25000]
  plot,allcal[q].jd-2400000,allcal[q].flux[150,ichip]/allcal[q].nread*10.,psym=6,xtickformat='(i5)',xtitle='MJD',ytitle='Median flux',charsize=2,symsize=0.3,yr=yr,ys=1,thick=2,charthick=2
  oplot,allcal[q].jd-2400000,allcal[q].flux[10,ichip]/allcal[q].nread*10.,psym=6,color=2,symsize=0.3,thick=2
  oplot,allcal[q].jd-2400000,allcal[q].flux[150,ichip]/allcal[q].nread*10.,psym=6,color=3,symsize=0.3,thick=2
  oplot,allcal[q].jd-2400000,allcal[q].flux[290,ichip]/allcal[q].nread*10.,psym=6,color=4,symsize=0.3,thick=2
  al_legend,['Fiber 290','Fiber 150','Fiber 10'],textcolors=[2,3,4],/right
  al_legend,['Chip '+string(ichip)],/bottom
endfor
if hard eq 1 then begin
  device,/close
  ps2gif,sdir+'qflux.eps',/eps,/delete
  printf,html,'<img src='+mdir+'/qflux.gif>'
endif else stop

; Individual fiber throughputs
j=where(allcal.qrtz gt 0)                                                       
bd=where(allcal.qrtz gt 0 and max(allcal.flux[*,1],dim=1) lt 500)
file_mkdir,sdir+'/fiber'
printf,html,'<h3> <a href='+mdir+'/fiber/fiber.html> Individual fiber throughputs from quartz </h3>'
openw,fhtml,sdir+'/fiber/fiber.html',/get_lun
printf,fhtml,'<HTML><BODY><TABLE BORDER=2>'
openw,convert,sdir+'/fiber/fiber.csh',/get_lun
printf,convert,'#!/bin/csh'
printf,convert,'cd '+sdir+'/fiber'
for ifiber=1,300 do begin
  plotfile='fiber'+string(format='(i3.3)',ifiber)
  device,file=sdir+'/fiber/'+plotfile+'.eps',/encap,/color,xsize=20,ysize=4,/in
  !p.multi=[0,0,0]
  plot,allcal[j].jd-2400000,allcal[j].flux[300-ifiber,1],psym=6,xtickformat='(i5)',xtitle='MJD',ytitle='flux',charsize=2,symsize=0.3,thick=2,charthick=2
  device,/close
  printf,convert,'convert '+plotfile+'.eps '+plotfile+'.jpg'
  printf,convert,'rm '+plotfile+'.eps'
  ;ps2gif,sdir+'/fiber/'+plotfile+'.eps',/eps,/delete
  printf,fhtml,'<TR><TD>'+string(ifiber)+'<TD><img src='+plotfile+'.jpg>'
endfor
printf,fhtml,'</TABLE></BODY></HTML>'
free_lun,fhtml
free_lun,convert
if ~keyword_set(nofiber) then spawn,'csh '+sdir+'/fiber/fiber.csh'

; ThAr line brightness
if hard eq 1 then begin
  device,file=sdir+'tharflux.eps',/encap,/color,xsize=20,ysize=8,/in
  !p.multi=[0,1,3]
  printf,html,'<h3> <a name=tharflux></a>ThAr line brightness (per 10 reads) in extracted frame </h3>'
endif
flux=reform(allcal[t].gauss[0,*,*]*allcal[t].gauss[2,*,*]^2)
for ichip=0,2 do begin
  plot,allcal[t].jd-2400000,flux[0,ichip,*]/allcal[t].nread*10.,psym=6,xtickformat='(i5)',xtitle='MJD',ytitle='line flux',charsize=2,symsize=0.3,thick=2,charthick=2
  oplot,allcal[t].jd-2400000,flux[0,ichip,*]/allcal[t].nread*10.,psym=6,color=2,symsize=0.3,thick=2
  oplot,allcal[t].jd-2400000,flux[1,ichip,*]/allcal[t].nread*10.,psym=6,color=3,symsize=0.3,thick=2
  oplot,allcal[t].jd-2400000,flux[2,ichip,*]/allcal[t].nread*10.,psym=6,color=4,symsize=0.3,thick=2
  oplot,allcal[t].jd-2400000,flux[3,ichip,*]/allcal[t].nread*10.,psym=6,color=5,symsize=0.3,thick=2
  oplot,allcal[t].jd-2400000,flux[4,ichip,*]/allcal[t].nread*10.,psym=6,color=6,symsize=0.3,thick=2
  al_legend,['Fiber 290','Fiber 150','Fiber 10'],textcolors=[2,3,4],/right
  al_legend,['Chip '+string(ichip)],/bottom
endfor
if hard eq 1 then begin
  device,/close
  ps2gif,sdir+'tharflux.eps',/eps,/delete
  printf,html,'<img src='+mdir+'/tharflux.gif>'
endif else stop

; UNe line brightness
if hard eq 1 then begin
  device,file=sdir+'uneflux.eps',/encap,/color,xsize=20,ysize=8,/in
  !p.multi=[0,1,3]
  printf,html,'<h3> <a name=uneflux></a>UNe line brightness (per 10 reads) in extracted frame </h3>'
endif
flux=reform(allcal[u].gauss[0,*,*]*allcal[u].gauss[2,*,*]^2)
for ichip=0,2 do begin
  plot,allcal[u].jd-2400000,flux[0,ichip,*]/allcal[u].nread*10.,psym=6,xtickformat='(i5)',xtitle='MJD',ytitle='line flux',charsize=2,symsize=0.3,thick=2,charthick=2
  oplot,allcal[u].jd-2400000,flux[0,ichip,*]/allcal[u].nread*10.,psym=6,color=2,symsize=0.3,thick=2
  oplot,allcal[u].jd-2400000,flux[2,ichip,*]/allcal[u].nread*10.,psym=6,color=4,symsize=0.3,thick=2
  oplot,allcal[u].jd-2400000,flux[4,ichip,*]/allcal[u].nread*10.,psym=6,color=6,symsize=0.3,thick=2
  al_legend,['Fiber 290','Fiber 150','Fiber 10'],textcolors=[2,4,6],/right
  al_legend,['Chip '+string(ichip)],/bottom
endfor
if hard eq 1 then begin
  device,/close
  ps2gif,sdir+'uneflux.eps',/eps,/delete
  printf,html,'<img src='+mdir+'/uneflux.gif>'
endif else stop

; Get long term trends from dome flats
; Append together the individual summary files
; find the domeflats
d=where(strtrim(allexp.imagetyp,2) eq 'DomeFlat')

if hard eq 1 then begin
  printf,html,'<h3> <a name=dome></a>Dome flat median brightness</h3>'
  device,file=sdir+'dome.eps',/encap,/color,xsize=20,ysize=8,/in
  !p.multi=[0,1,3]
endif 
for ichip=0,2 do begin
  w=median(allexp[d].med[*,ichip])
  plot,allexp[d].jd-2400000,allexp[d].med[150,ichip],psym=6,xtickformat='(i5)',charsize=2,yrange=[0,2*w],symsize=0.3,thick=2,charthick=2
  oplot,!x.crange,[w,w],symsize=0.3,thick=2
  oplot,allexp[d].jd-2400000,allexp[d].med[10,ichip],psym=6,color=2,symsize=0.3,thick=2
  oplot,allexp[d].jd-2400000,allexp[d].med[150,ichip],psym=6,color=4,symsize=0.3,thick=2
  oplot,allexp[d].jd-2400000,allexp[d].med[290,ichip],psym=6,color=6,symsize=0.3,thick=2
  al_legend,['Fiber 290','Fiber 150','Fiber 10'],textcolors=[2,4,6]
  al_legend,['Chip '+string(ichip)],/bottom
endfor
if hard eq 1 then begin
  device,/close
  ps2gif,sdir+'dome.eps',/eps,/delete
  printf,html,'<img src='+mdir+'/dome.gif>'
endif else stop

; zeropoints
if hard eq 1 then begin
  printf,html,'<h3> <a name=zero></a>Science frame zero point</h3>'
  device,file=sdir+'zero.eps',/encap,/color,xsize=20,ysize=8,/in
  !p.multi=[0,0,0]
endif 
plot,allsci.mjd,allsci.zero,psym=6,xtickformat='(i5)',charsize=2,yrange=[13,20],symsize=0.3,thick=2,charthick=2
if hard eq 1 then begin
  device,/close
  ps2gif,sdir+'zero.eps',/eps,/delete
  printf,html,'<img src='+mdir+'/zero.gif>'
endif else stop

; Position of a bright line in each chip
if hard eq 1 then begin
  printf,html,'<h3> <a name=tpos></a>ThArNe lamp line position</h3>'
  device,file=sdir+'tpos.eps',/encap,/color,xsize=20,ysize=16,/in
  !p.multi=[0,1,3]
endif 
for ichip=0,2 do begin
  w=median(allcal[t].gauss[1,*,ichip])
  plot,allcal[t].jd-2400000,allcal[t].gauss[1,0,ichip],psym=6,xtickformat='(i5)',yrange=[w-10,w+10],xtitle='MJD',ytitle='Position',ystyle=1,charsize=2,symsize=0.3,thick=2,charthick=2
  oplot,allcal[t].jd-2400000,allcal[t].gauss[1,0,ichip],psym=6,color=2,symsize=0.3,thick=2
  oplot,allcal[t].jd-2400000,allcal[t].gauss[1,2,ichip],psym=6,color=4,symsize=0.3,thick=2
  oplot,allcal[t].jd-2400000,allcal[t].gauss[1,4,ichip],psym=6,color=6,symsize=0.3,thick=2
  al_legend,['Fiber 290','Fiber 150','Fiber 10'],textcolors=[2,4,6]
  al_legend,['Chip '+string(ichip)],/bottom
endfor
if hard eq 1 then begin
  device,/close
  ps2gif,sdir+'tpos.eps',/eps,/delete
  printf,html,'<img src='+mdir+'/tpos.gif>'
endif else stop

; FWHM of a bright line in each chip
for iline=0,1 do begin
  plotfile='tfwhm'+string(format='(i1)',iline)
  if hard eq 1 then begin
    printf,html,'<h3> <a name=tfwhm></a> ThArNe lamp line FWHM, line position (x pixel): '+string(format='(3f8.1)',allcal[t[0]].lines[*,iline])+'</h3>'
    device,file=sdir+plotfile+'.eps',/encap,/color,xsize=20,ysize=8,/in
    !p.multi=[0,1,3]
  endif
  for ichip=0,2 do begin
    w=median(2.354*allcal[t].gauss[2,*,ichip,iline])
    plot,allcal[t].jd-2400000,2.354*allcal[t].gauss[2,0,ichip,iline],psym=6,xtickformat='(i5)',yrange=[0.7*w,1.3*w],xtitle='MJD',ytitle='FWHM',ystyle=1,charsize=2,symsize=0.3,charthick=2,thick=2
    oplot,allcal[t].jd-2400000,2.354*allcal[t].gauss[2,0,ichip,iline],psym=6,color=2,symsize=0.3,thick=2
    w=median(2.354*allcal[t].gauss[2,0,ichip,iline])
    oplot,!x.crange,[w,w],symsize=0.3,thick=2
    oplot,allcal[t].jd-2400000,2.354*allcal[t].gauss[2,1,ichip,iline],psym=6,color=3,symsize=0.3,thick=2
    w=median(2.354*allcal[t].gauss[2,1,ichip,iline])
    oplot,!x.crange,[w,w],color=2,thick=2
    oplot,allcal[t].jd-2400000,2.354*allcal[t].gauss[2,2,ichip,iline],psym=6,color=4,symsize=0.3,thick=2
    w=median(2.354*allcal[t].gauss[2,2,ichip,iline])
    oplot,!x.crange,[w,w],color=4,thick=2
    oplot,allcal[t].jd-2400000,2.354*allcal[t].gauss[2,3,ichip,iline],psym=6,color=5,symsize=0.3,thick=2
    w=median(2.354*allcal[t].gauss[2,3,ichip,iline])
    oplot,!x.crange,[w,w],color=4,thick=2
    oplot,allcal[t].jd-2400000,2.354*allcal[t].gauss[2,4,ichip,iline],psym=6,color=6,symsize=0.3,thick=2
    w=median(2.354*allcal[t].gauss[2,4,ichip,iline])
    oplot,!x.crange,[w,w],color=4,thick=2
    al_legend,['Fiber 290','Fiber 220','Fiber 150','Fiber 80','Fiber 10'],textcolors=[2,3,4,5,6]
    al_legend,['Chip '+string(ichip)],/bottom
  endfor
  if hard eq 1 then begin
    device,/close
    ps2gif,sdir+plotfile+'.eps',/eps,/delete
    printf,html,'<img src='+mdir+'/'+plotfile+'.gif>'
  endif else stop
endfor

; trace locations
files=file_search(specdir+'/cal/psf/'+dirs.prefix+'EPSF-b-*.fits')
cent=[]
mjd=[]
oldnum=0
tmp={num: 0L, mjd: 0.d0, cent: 0., ln2level: 0.}
allepsf=[]
for i=0,n_elements(files)-1 do begin
  print,files[i]
  num=long(strmid(file_basename(files[i]),9,8))/10000
  if num gt 1000 then begin
    a=mrdfits(files[i],0,hdr,status=status,/silent)
    for j=147,155 do begin
      a=mrdfits(files[i],j,status=status,/silent)
      if status ge 0 then begin
        if a.fiber eq 150 then begin
          tmp.num=long(strmid(file_basename(files[i]),9,8))
          tmp.cent=a.cent[1000]
          tmp.mjd=sxpar(hdr,'JD-MID')-2400000.5
          tmp.ln2level=sxpar(hdr,'LN2LEVEL')
          oldnum=num
          break
        endif
      endif
    endfor
    allepsf=struct_append(allepsf,tmp)
  endif
endfor
mwrfits,allepsf,specdir+'/'+dirs.instrument+'Trace.fits',/create

if hard eq 1 then begin
  plotfile='trace'
  printf,html,'<h3> <a name=trace></a> Trace position, fiber 150, column 1000</h3>'
  device,file=sdir+plotfile+'.eps',/encap,/color,xsize=20,ysize=4,/in
  !p.multi=[0,1,2]
  plot,allepsf.mjd,allepsf.cent,psym=6,yr=[median(allepsf.cent)-1,median(allepsf.cent)+1],xtickformat='(i5)',xtitle='MJD',ytitle='Trace center',thick=2,charthick=2,ys=1
  plot,allepsf.ln2level,allepsf.cent,psym=6,yr=[median(allepsf.cent)-1,median(allepsf.cent)+1],xtitle='LN2 Level',ytitle='Trace center',thick=2,charthick=2,ys=1
  device,/close
  ps2gif,sdir+plotfile+'.eps',/eps,/delete
  printf,html,'<img src='+mdir+'/'+plotfile+'.gif>'
endif                                                                                                                                                                                                                                    
; detectors
if hard eq 1 then begin
  printf,html,'<h3> <a name=detectors></a>Detectors</h3>'
  printf,html,'<h4> Dark mean </h4>'
  device,file=sdir+'biasmean.eps',/encap,/color,xsize=20,ysize=8,/in
endif
dark=where(strpos(strupcase(alldark.exptype),'DARK') ge 0)
!p.multi=[0,1,3]
for ichip=0,2 do begin
 plot,alldark[dark].jd-2400000.5,alldark[dark].mean[ichip,2],psym=6,xtit='MJD',ytit=textoidl('mean(column median)'),/ylog,charsize=2,symsize=0.3,yr=[0.1,1000.],xtickformat='(i5)',thick=2,charthick=2
 for ibias=0,3 do begin
   oplot,alldark[dark].jd-2400000.5,alldark[dark].mean[ichip,ibias],psym=6,color=ibias+1,symsize=0.3,thick=2
 endfor
 al_legend,['Chip '+string(ichip)],/bottom
endfor
if hard eq 1 then begin
  device,/close
  ps2gif,sdir+'biasmean.eps',/eps,/del
  printf,html,'<img src='+mdir+'/biasmean.gif>'
endif

if hard eq 1 then begin
  printf,html,'<h4> Dark sigma </h4>'
  device,file=sdir+'biassig.eps',/encap,/color,xsize=20,ysize=8,/in
endif
dark=where(strpos(strupcase(alldark.exptype),'DARK') ge 0)
!p.multi=[0,1,3]
for ichip=0,2 do begin
 plot,alldark[dark].jd-2400000.5,alldark[dark].sig[ichip,2],psym=6,xtit='MJD',ytit=textoidl('\sigma(column median)'),/ylog,charsize=2,symsize=0.3,yr=[0.1,1000.],xtickformat='(i5)',thick=2,charthick=2
 for ibias=0,3 do begin
   oplot,alldark[dark].jd-2400000.5,alldark[dark].sig[ichip,ibias],psym=6,color=ibias+1,symsize=0.3,thick=2
 endfor
 al_legend,['Chip '+string(ichip)],/bottom
endfor
if hard eq 1 then begin
  device,/close
  ps2gif,sdir+'biassig.eps',/eps,/del
  printf,html,'<img src='+mdir+'/biassig.gif>'
endif

if hard eq 1 then begin
  printf,html,'<h3> <a name=sky></a>Sky brightness</h3>'
  device,file=sdir+'moonsky.eps',/encap,/color,xsize=20,ysize=8,/in
endif
!p.multi=[0,1,2]
loadct,39
plotc,allsci.moonphase,allsci.sky[1],allsci.moondist,psym=6,yr=[16,12],xtitle='moon phase',ytitle='Sky brightness',title='moon distance',min=0,max=90,symsize=0.3
gd=where(finite(allsci.zero) eq 1 and allsci.zero lt 20 and allsci.zero gt 0)
plotc,allsci[gd].moonphase,allsci[gd].sky[1],allsci[gd].zero,psym=6,yr=[16,12],min=17,max=19,xtitle='moon phase', ytitle='Sky brightness',title='zeropoint (cloudiness)',symsize=0.3
if hard eq 1 then begin
  device,/close
  ps2gif,sdir+'moonsky.eps',/eps,/delete
  printf,html,'<img src='+mdir+'/moonsky.gif>'
endif

printf,html,'</BODY></HTML>'
free_lun,html
end
