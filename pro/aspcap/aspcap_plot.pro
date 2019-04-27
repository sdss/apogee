;+
; aspcap_plot
;
; makes a bunch of summary plots/html page for an allstar structure: params vs params, etc.
;-

function countflag,flag,name
  gd=where((flag and aspcapflagval('NO_ASPCAP_RESULT')) eq 0)
  j=where((flag[gd] and (aspcapflagval(name+'_BAD') or aspcapflagval(name+'_WARN'))) eq 0,ngood)
  j=where((flag[gd] and aspcapflagval(name+'_WARN')) gt 0,nwarn)
  j=where((flag[gd] and aspcapflagval(name+'_BAD')) gt 0,nbad)
  return,[ngood,nwarn,nbad]
end

pro aspcap_plot,str,file,hard=hard,symsize=symsize,charsize=charsize,results_vers=results_vers,corronly=corronly,fit=fit,errfact=errfact,elems=elem_symbol

if n_elements(str) eq 0 then begin
  print,'no allStar structure specified!'
  print,'making master index.html only'
  goto,index
endif

; makes summary ASPCAP web page and plots for all objects passed in input structure
if ~keyword_set(errfact) then errfact=1

file_mkdir,file_dirname(file)

; with fit keyword, plot using fit (not calibrated) quantities
if keyword_set(fit) then begin
  teff=str.fparam[0]
  teff_err=sqrt(str.fparam_cov[0,0])*errfact
  logg=str.fparam[1]
  logg_err=sqrt(str.fparam_cov[1,1])*errfact
  feh=str.fparam[3]
  feh_err=sqrt(str.fparam_cov[3,3])*errfact
  cfe=str.fparam[4]
  cfe_err=sqrt(str.fparam_cov[4,4])*errfact
  nfe=str.fparam[5]
  nfe_err=sqrt(str.fparam_cov[5,5])*errfact
  afe=str.fparam[6]
  afe_err=sqrt(str.fparam_cov[6,6])*errfact
  file=file+'fit'
  overplot=0
endif else begin
  teff=str.teff
  teff_err=str.teff_err
  logg=str.logg
  logg_err=str.logg_err
  feh=str.m_h
  feh_err=str.m_h_err
  afe=str.alpha_m
  afe_err=str.alpha_m_err
  cfe=str.c_m
  cfe_err=str.c_m_err
  nfe=str.n_m
  nfe_err=str.n_m_err
  overplot=1
endelse
teffflag=str.aspcapflag and (aspcapflagval('TEFF_WARN') or aspcapflagval('TEFF_BAD'))
loggflag=str.aspcapflag and (aspcapflagval('LOGG_WARN') or aspcapflagval('LOGG_BAD'))
fehflag=str.aspcapflag and (aspcapflagval('M_H_WARN') or aspcapflagval('M_H_BAD'))
afeflag=str.aspcapflag and (aspcapflagval('ALPHA_M_WARN') or aspcapflagval('ALPHA_M_BAD'))
cfeflag=str.aspcapflag and (aspcapflagval('C_M_WARN') or aspcapflagval('C_M_BAD'))
nfeflag=str.aspcapflag and (aspcapflagval('N_M_WARN') or aspcapflagval('N_M_BAD'))
chi2=str.aspcap_chi2
sn=str.snr
jk0=(str.j-str.k)-1.5*str.ak_targ
tmin=3500
tmax=10000

if not keyword_set(symsize) then symsize=0.2
if not keyword_set(charsize) then charsize=0.8

gd=where((str.aspcapflag and (warnaspcapflag() or badaspcapflag())) eq 0)
if keyword_set(hard) then set_plot,'PS'
if keyword_set(corronly) then goto,docorr

cleanplot,/silent
if keyword_set(hard) then begin
  openw,html,file+'.html',/get_lun
  printf,html,'<HTML><BODY>'
  openw,htmlerr,file+'_err.html',/get_lun
  printf,htmlerr,'<HTML><BODY>'

  printf,html,'<a href=../allStar-'+results_vers+'.fits> allStar-'+results_vers+'.fits </a><br>'
  printf,html,'<a href=../allVisit-'+results_vers+'.fits> allVisit-'+results_vers+'.fits </a><br>'
  ntot=n_elements(str)
  j=where((str.aspcapflag and aspcapflagval('NO_ASPCAP_RESULT')) eq 0,naspcap)
  main=where((str.aspcapflag and aspcapflagval('NO_ASPCAP_RESULT')) eq 0 and (strpos(str.targflags,'APOGEE_SHORT') ge 0) or (strpos(str.targflags,'APOGEE_INTERMEDIATE') ge 0) or (strpos(str.targflags,'APOGEE_LONG') ge 0),nmain)
  gmain=where((str.aspcapflag and aspcapflagval('NO_ASPCAP_RESULT')) eq 0 and (strpos(str.targflags,'APOGEE_SHORT') ge 0 or strpos(str.targflags,'APOGEE_INTERMEDIATE') ge 0 or strpos(str.targflags,'APOGEE_LONG') ge 0) and str.fparam[1] lt 3.8,ngmain)
  giant=where((str.aspcapflag and aspcapflagval('NO_ASPCAP_RESULT')) eq 0 and str.fparam[1] lt 3.8,ngiant)
  printf,html,'<TABLE BORDER=2>'
  printf,html,'<TR><TD>N(TOTAL)<TD>', ntot
  printf,html,'<TR><TD>N(ASPCAP)<TD>', naspcap
  printf,html,'<TR><TD>N(ASPCAP+MAIN)<TD>', nmain
  printf,html,'<TR><TD>N(ASPCAP+MAIN+GIANT)<TD>', ngmain
  printf,html,'<TR><TD>N(ASPCAP+GIANT)<TD>', ngiant
  j=where((str.starflag and starflagval('SUSPECT_RV_COMBINATION')) gt 0,n)
  printf,html,'<TR><TD>SUSPECT_RV_COMBINATION <TD>',n
  j=where((str.starflag and starflagval('PERSIST_HIGH')) gt 0,n)
  printf,html,'<TR><TD>PERSIST_HIGH <TD>',n
  printf,html,'</TABLE>'

  printf,html,'ASPCAP flags'
  printf,html,'<TABLE BORDER=2>'
  printf,html,'<TR><TD><TD colspan=3>All ASPCAP'
  printf,html,'<TD colspan=3>Main sample'
  printf,html,'<TD colspan=3>Main AND giant sample'
  printf,html,'<TD colspan=3>Giant sample'
  printf,html,'<TR><TD>Condition <TD> Good <TD> WARN set <TD> BAD set'
  printf,html,'<TD> Good <TD> WARN set <TD> BAD set'
  printf,html,'<TD> Good <TD> WARN set <TD> BAD set'
  printf,html,'<TD> Good <TD> WARN set <TD> BAD set'

  flags=['STAR','ROTATION','CHI2','SN','COLORTE','TEFF','LOGG','M_H','ALPHA_M','C_M','N_M']
  for iflag=0,n_elements(flags)-1 do begin
    n=countflag(str.aspcapflag,flags[iflag])
    printf,html,'<TR><TD>'+flags[iflag]+'<TD>',n[0],'<TD>',n[1],'<TD>',n[2]
    n=countflag(str[main].aspcapflag,flags[iflag])
    printf,html,'<TD>',n[0],'<TD>',n[1],'<TD>',n[2]
    n=countflag(str[gmain].aspcapflag,flags[iflag])
    printf,html,'<TD>',n[0],'<TD>',n[1],'<TD>',n[2]
    n=countflag(str[giant].aspcapflag,flags[iflag])
    printf,html,'<TD>',n[0],'<TD>',n[1],'<TD>',n[2]
  endfor
  printf,html,'</TABLE>'

  printf,html,'<TABLE>'
  ; CHI2 histogram
  printf,html,'<TR><TD>Chi^2 distribution (ASPCAP_FLAG not bad in green)'
  printf,html,'<TD>S/N distribution (ASPCAP_FLAG not bad in green)'
  set_plot,'PS'
  device,file=file+'chi2hist.eps',/encap,xsize=8,ysize=8,/inches
  smcolor,/ps
  plothist,fix(nint(chi2)),bin=1,xrange=[1,100],xstyle=1,thick=2,charthick=2
  plothist,fix(nint(chi2[gd])),bin=1,/overplot,color=3
  device,/close
  ps2gif,file+'chi2hist.eps',/eps,/delete
  printf,html,'<TR><TD><a href='+file_basename(file)+'.chi2hist.gif><IMG src='+file_basename(file)+'chi2hist.gif width=500></a>'

  ; S/N histogram
  device,file=file+'sn.eps',/encap,xsize=8,ysize=8,/inches
  smcolor,/ps
  plothist,fix(nint(sn)),bin=1,xrange=[1,500],xstyle=1,thick=2,charthick=2
  plothist,fix(nint(sn[gd])),bin=1,/overplot,color=3
  device,/close
  ps2gif,file+'sn.eps',/eps,/delete
  printf,html,'<TD><a href='+file_basename(file)+'.sn.gif><IMG src='+file_basename(file)+'sn.gif width=500></a>'

  ; CHI2 with temperature
  device,file=file+'chi2.eps',/encap,xsize=8,ysize=8,/inches
  smcolor,/ps
  dwarfs=where(str.fparam[1] gt 3.8)
  plot,str.teff,str.aspcap_chi2,psym=3,xrange=[7500,3300],yrange=[0,100],xtitle=textoidl('T_{\rm{eff}}'),ytitle=textoidl('\Chi^2'),xstyle=1,thick=2
  oplot,str[dwarfs].teff,str[dwarfs].aspcap_chi2,psym=3,color=3
  ;oplot,[7500,3300],[10,10],color=2,thick=3,linestyle=2
  ;oplot,[7500,3300],[30,30],color=2,thick=3
  xyouts,7000,80,'Giants (log g < 3.8): black'
  xyouts,7000,70,'Dwarfs (log g > 3.8): green',color=3
  dt=100 & medt=3550+indgen(70)*dt
  medchi=fltarr(n_elements(medt)) & medn=intarr(n_elements(medt))
  for i=0,n_elements(medt)-1 do begin
    j=where(str.teff ge medt[i]-dt/2. and str.teff lt medt[i]+dt/2.,nj)
    medn[i]=nj
    if nj gt 1 then medchi[i]=median(str[j].aspcap_chi2)
  endfor
  chicrit=10+fltarr(n_elements(medt))
  oplot,medt,max([[medchi*2],[chicrit]],dim=2),psym=10,thick=3,color=2,linestyle=2
  chicrit=30+fltarr(n_elements(medt))
  oplot,medt,max([[medchi*3],[chicrit]],dim=2),psym=10,thick=3,color=2
  device,/close
  ps2gif,file+'chi2.eps',/eps,/delete
  printf,html,'<TD><a href='+file_basename(file)+'chi2.gif><IMG src='+file_basename(file)+'chi2.gif width=500></a>'

  ; color-temp relations
  printf,html,'<TR><TD> Color-temp relation, color-coded by chi^2'
  printf,html,'<TD> Color-temp relation, color-coded by [Fe/H] (with GHB IRFM relations)'
  printf,html,'<TD> Teff-Teff(GHB)'
  printf,html,'<TD> Color-temp relation, color-coded by logg'
  printf,html,'<TR>'

  device,file=file+'temp.eps',/encap,xsize=8,ysize=8,/inches
  loadct,39
  x=jk0 & xlab='(J-K)_0'  & xrange=[0,1.5]
  y=teff & ylab='Teff' & yrange=[tmin,tmax]
  z=chi2 & zmin=0 & zmax=30
  plotc,x,y,z,psym=6,xtitle=xlab,ytitle=ylab,min=zmin,max=zmax,$
       thick=2,charthick=2,xrange=xrange,yrange=yrange,symsize=symsize,charsize=charsize
  device,/close
  ps2gif,file+'temp.eps',/eps,/delete
  printf,html,'<TD><a href='+file_basename(file)+'.temp.gif><IMG src='+file_basename(file)+'temp.gif width=500></a>'

  device,file=file+'temp2.eps',/encap,xsize=8,ysize=8,/inches
  loadct,39
  x=jk0 & xlab='(J-K)_0'  & xrange=[0,1.5]
  y=teff & ylab='Teff' & yrange=[tmin,tmax]
  z=feh & zmin=-2 & zmax=0.5
  plotc,x,y,z,psym=6,xtitle=xlab,ytitle=ylab,min=zmin,max=zmax,$
       thick=2,charthick=2,xrange=xrange,yrange=yrange,symsize=symsize,charsize=charsize
  ; overplot Gonzalez Hernandez & Bonifacio relation for giants with [Fe/H]=0
  smcolor
  x=0.+indgen(100)*0.01
  fe=[0.,-1.,-2.] & fecolor=[2,3,4]
  for ife=0,n_elements(fe)-1 do begin
    ;oplot,x,aspcap_colorte(x,fe[ife]),color=fecolor[ife],thick=2,linestyle=1
    oplot,x,aspcap_colorte(x,fe[ife])-500,color=2,thick=2,linestyle=2
    oplot,x,aspcap_colorte(x,fe[ife])+500,color=2,thick=2,linestyle=2
    oplot,x,aspcap_colorte(x,fe[ife])-1000,color=2,thick=2
    oplot,x,aspcap_colorte(x,fe[ife])+1000,color=2,thick=2
  endfor
  device,/close
  ps2gif,file+'temp2.eps',/eps,/delete
  printf,html,'<TD><a href='+file_basename(file)+'.temp2.gif><IMG src='+file_basename(file)+'temp2.gif width=500></a>'

  device,file=file+'temp2diff.eps',/encap,xsize=8,ysize=8,/inches
  !p.multi=[0,1,2]
  x=teff & xrange=[7000,3500] & xlab='Teff'
  y=teff-aspcap_colorte(jk0,feh) & yrange=[-1000,1000] & ylab='Teff-Teff(GH&B)'
  z=feh
  plt=where(str.ak_targ lt 0.05 and str.fparam[1] lt 3.8)
  plotc,x[plt],y[plt],z[plt],psym=6,xtitle=xlab,ytitle=ylab,xstyle=1,min=zmin,max=zmax,$
       thick=2,charthick=2,xrange=xrange,yrange=yrange,symsize=symsize,charsize=charsize
  plt=where(str.ak_targ lt 0.05 and str.fparam[1] gt 3.8)
  plotc,x[plt],y[plt],z[plt],psym=6,xtitle=xlab,ytitle=ylab,xstyle=1,min=zmin,max=zmax,$
       thick=2,charthick=2,xrange=xrange,yrange=yrange,symsize=symsize,charsize=charsize
  device,/close
  !p.multi=[0,0,0]
  ps2gif,file+'temp2diff.eps',/eps,/delete
  printf,html,'<TD><a href='+file_basename(file)+'.temp2diff.gif><IMG src='+file_basename(file)+'temp2diff.gif width=500></a>'

  device,file=file+'temp3.eps',/encap,xsize=8,ysize=8,/inches
  loadct,39
  x=jk0 & xlab='(J-K)_0'  & xrange=[0,1.5]
  y=teff & ylab='Teff' & yrange=[tmin,tmax]
  z=logg & zmin=0 & zmax=5
  plotc,x,y,z,psym=6,xtitle=xlab,ytitle=ylab,min=zmin,max=zmax,$
       thick=2,charthick=2,xrange=xrange,yrange=yrange,symsize=symsize,charsize=charsize
  device,/close
  ps2gif,file+'temp3.eps',/eps,/delete
  printf,html,'<TD><a href='+file_basename(file)+'.temp3.gif><IMG src='+file_basename(file)+'temp3.gif width=500></a>'

  printf,html,'</TABLE>'

  ; microturbulence plots
  printf,html,'<br> Microturbulence relation'
  printf,html,'<TABLE><TR>'
  plotfile=file+'vmicro1'
  device,file=plotfile+'.eps',/color,/encap
  gd=where((str.aspcapflag and badaspcapflag()) eq 0 and str.snr gt 100 and str.fparam[0] lt 5500)
  plotc,str[gd].fparam[1],10.^(str[gd].fparam[2]),str[gd].fparam[3],psym=6,min=-2,max=0.5,xr=[0,5],yr=[0,5],xtitle='log g',ytitle='vmicro (km/s)',symsize=0.2
  device,/close
  ps2jpg,plotfile+'.eps',/eps,/delete
  printf,html,'<TD><a href='+file_basename(plotfile)+'.jpg><IMG SRC='+file_basename(plotfile)+'.jpg></a>'
  plotfile=file+'vmicro2'
  device,file=plotfile+'.eps',/color,/encap
  gd=where((str.aspcapflag and badaspcapflag()) eq 0 and str.snr gt 100 and str.fparam[0] lt 5500)
  plotc,str[gd].fparam[1],10.^(str[gd].fparam[2]),str[gd].fparam[0],psym=6,min=3500,max=5500,xr=[0,5],yr=[0,5],xtitle='log g',ytitle='vmicro (km/s)',symsize=0.2
  device,/close
  ps2jpg,plotfile+'.eps',/eps,/delete
  printf,html,'<TD><a href='+file_basename(plotfile)+'.jpg><IMG SRC='+file_basename(plotfile)+'.jpg></a>'
  printf,html,'</TABLE>'

  printf,html,'<br> Following plots only have stars with ASPCAP_FLAG=0<br>'
  printf,html,'<TABLE BORDER=2>'
  printf,html,'<TR><TD>x-axis<TD>Teff<TD>log g<TD>[Fe/H]<TD>[alpha/Fe]<TD>[C/Fe]<TD>[N/Fe]'
  printf,htmlerr,'<TABLE BORDER=2>'
  printf,htmlerr,'<TR><TD>x-axis<TD>Teff<TD>log g<TD>[Fe/H]<TD>[alpha/Fe]<TD>[C/Fe]<TD>[N/Fe]'
endif


; multiplots of parameters and parameter errors vs lots of things
for ix=0,8 do begin
  if ix eq 0 then begin
    x=teff & xlab='Teff' & xrange=[tmax,tmin] & xflag=teffflag
    z=chi2 & zmin=0 & zmax=20 & zlab='chi^2'
  endif
  if ix eq 1 then begin
    x=teff & xlab='Teff' & xrange=[tmax,tmin] & xflag=teffflag
    z=feh & zmin=-2 & zmax=0.5 & zlab='[Fe/H]'
  endif
  if ix eq 2 then begin
    x=teff & xlab='Teff' & xrange=[tmax,tmin] & xflag=teffflag
    z=logg & zmin=0 & zmax=5 & zlab='log g'
  endif
  if ix eq 3 then begin
    x=teff & xlab='Teff' & xrange=[tmax,tmin] & xflag=teffflag
    z=sn & zmin=50 & zmax=150 & zlab='S/N'
  endif
  if ix eq 4 then begin
    x=logg & xlab='log g' & xrange=[5,0] & xflag=loggflag
    z=teff & zmin=3500 & zmax=6500 & zlab='Teff'
  endif
  if ix eq 5 then begin
    x=feh & xlab='[Fe/H]' & xrange=[-2.5,0.5] & xflag=fehflag
    z=teff & zmin=3500 & zmax=6500 & zlab='Teff'
  endif
  if ix eq 6 then begin
    x=feh & xlab='[Fe/H]' & xrange=[-2.5,0.5] & xflag=fehflag
    z=logg & zmin=0 & zmax=5 & zlab='log g'
  endif
  if ix eq 7 then begin
    x=afe & xlab='[alpha/Fe]' & xrange=[-1.0,1.0] & xflag=afeflag
    z=teff & zmin=3500 & zmax=6500 & zlab='Teff'
  endif
  if ix eq 8 then begin
    x=afe & xlab='[alpha/Fe]' & xrange=[-1.0,1.0] & xflag=afeflag
    z=logg & zmin=0 & zmax=5 & zlab='log g'
  endif
  if keyword_set(hard) then begin
     printf,html,'<TR><TD>'+xlab+' coded by '+zlab
     printf,htmlerr,'<TR><TD>'+xlab+' coded by '+zlab
  endif
  for iy=0,5 do begin
    if iy eq 0 then begin
      y=teff & ylab='Teff' & yrange=[tmin,tmax] & yflag=teffflag & yerr=teff_err & yerrrange=[0,500]
    endif
    if iy eq 1 then begin
      y=logg & ylab='log g' & yrange=[5,0] & yflag=loggflag & yerr=logg_err & yerrrange=[0,1]
    endif
    if iy eq 2 then begin
      y=feh & ylab='[Fe/H]' & yrange=[-2.5,0.5] & yflag=fehflag & yerr=feh_err & yerrrange=[0,0.5]
    endif
    if iy eq 3 then begin
      y=afe & ylab='[alpha/Fe]' & yrange=[-1,1] & yflag=afeflag & yerr=afe_err & yerrrange=[0,0.5]
    endif
    if iy eq 4 then begin
      y=cfe & ylab='[C/Fe]' & yrange=[-1,1] & yflag=cfeflag & yerr=cfe_err & yerrrange=[0,0.5]
    endif
    if iy eq 5 then begin
      y=nfe & ylab='[N/Fe]' & yrange=[-1,1] & yflag=nfeflag & yerr=nfe_err & yerrrange=[0,0.5]
    endif
;    if xlab ne ylab then begin

     if keyword_set(hard) then begin
       set_plot,'PS' 
       plotfile=file+string(format='(i1)',ix)+string(format='(i1)',iy)
       device,file=plotfile+'.eps',/encap,/port,xsize=8,ysize=8,/inches,/color
     endif else begin
       set_plot,'X'
     endelse

     gdplot=where(xflag eq 0 and yflag eq 0 and $
       (str.aspcapflag and $
       (aspcapflagval('CHI2_WARN') or aspcapflagval('SN_WARN') or aspcapflagval('ROTATION_WARN'))) eq 0,ngd)
     if keyword_set(overplot) then begin
       loadct,0 
       plotc,x,y,z,psym=6,xtitle=xlab,ytitle=ylab,min=zmin,max=1.e30,$
         thick=2,charthick=2,xrange=xrange,yrange=yrange,symsize=symsize,charsize=charsize,/nocolorbar
       loadct,39
       if ngd gt 0 then $
        plotc,x[gdplot],y[gdplot],z[gdplot],psym=6,xtitle=xlab,ytitle=ylab,min=zmin,max=zmax,$
         thick=2,charthick=2,xrange=xrange,yrange=yrange,symsize=symsize,charsize=charsize,/overplot
     endif else begin
       loadct,39
       plotc,x,y,z,psym=6,xtitle=xlab,ytitle=ylab,min=zmin,max=zmax,$
         thick=2,charthick=2,xrange=xrange,yrange=yrange,symsize=symsize,charsize=charsize
     endelse
     if keyword_set(hard) then begin
       device,/close
       ps2gif,plotfile+'.eps',/eps,/delete
       printf,html,'<TD><a href='+file_basename(plotfile)+'.gif><IMG SRC='+file_basename(plotfile)+'.gif width=250></a>'

       plotfile=file+'_err'+string(format='(i1)',ix)+string(format='(i1)',iy)
       device,file=plotfile+'.eps',/encap,/port,xsize=8,ysize=8,/inches,/color
       if keyword_set(overplot) then begin
         loadct,0
         plotc,x,yerr,z,psym=6,xtitle=xlab,ytitle=ylab,min=zmin,max=1.e30,$
           thick=2,charthick=2,xrange=xrange,yrange=yerrrange,symsize=symsize,charsize=charsize,/nocolorbar
         loadct,39
         if ngd gt 0 then $
          plotc,x[gdplot],yerr[gdplot],z[gdplot],psym=6,xtitle=xlab,ytitle=ylab,min=zmin,max=zmax,$    
           thick=2,charthick=2,xrange=xrange,yrange=yerrrange,symsize=symsize,charsize=charsize,/overplot
       endif else begin
         loadct,39
          plotc,x,yerr,z,psym=6,xtitle=xlab,ytitle=ylab,min=zmin,max=zmax,$    
           thick=2,charthick=2,xrange=xrange,yrange=yerrrange,symsize=symsize,charsize=charsize
       endelse
       device,/close
       ps2gif,plotfile+'.eps',/eps,/delete
       printf,htmlerr,'<TD><a href='+file_basename(plotfile)+'.gif><IMG SRC='+file_basename(plotfile)+'.gif width=250></a>'
     endif else stop
;    endif else begin
;     if keyword_set(hard) then begin
;       printf,html,'<TD>'
;       printf,htmlerr,'<TD>'
;     endif
;    endelse
  endfor
endfor

if keyword_set(hard) then begin
  printf,html,'</TABLE></BODY></HTML>'
  printf,htmlerr,'</TABLE></BODY></HTML>'
  free_lun,html
  free_lun,htmlerr
  set_plot,'X'
endif

docorr:
if keyword_set(hard) then begin
  openw,html,file+'corr.html',/get_lun
  printf,html,'<HTML><BODY>'
  printf,html,'<table border=2>'
  set_plot,'PS'
  loadct,39
endif

gd=where((str.aspcapflag and (warnaspcapflag() or badaspcapflag())) eq 0)
corrcoef=fltarr(n_elements(gd))
corrcoef2=fltarr(n_elements(gd))
fid=[100.,0.1,0.,0.01,0.01,0.01,0.01]
c=fltarr(7,7)
cparam=['Teff','log g','vmicro','[Z/H]','[alpha/Fe]','[C/Fe]','[N/Fe]']
if keyword_set(hard) then begin
  printf,html,'<TR><TD>'
  for i=0,6 do printf,html,'<TD>'+cparam[i]
endif
for i=0,6 do begin
 if keyword_set(hard) then printf,html,'<TR><TD>'+cparam[i]
 if i ne 2 then begin
  for j=0,6 do begin
   if keyword_set(hard) then printf,html,'<TD>'
   if j ne 2 then begin
    plotfile=file+'corr'+string(format='(i1)',i)+string(format='(i1)',j)
    for k=0L,n_elements(gd)-1 do begin
      c=corr(str[gd[k]].fparam_cov)
      corrcoef[k]=c[i,j]
      ;corrcoef2[k]=str[gd[k]].param_cov[i,j]/fid[i]/fid[j]
    endfor
    if keyword_set(hard) then device,file=plotfile+'.eps',/encap,/color
    if keyword_set(fit) then begin
      x=str[gd].fparam[0]
      z=str[gd].fparam[3]
    endif else begin
      x=str[gd].param[0]
      z=str[gd].param[3]
    endelse
    plotc,x,corrcoef,z,psym=6,zmin=-2,zmax=0.5,thick=3,xtitle='Teff',ytitle='Corr coeff '+cparam[i]+'-'+cparam[j],xrange=[3500,5500],xstyle=1,yrange=[-1,1],ystyle=1
    ;plotc,str[gd].param[0],corrcoef2,str[gd].param[3],psym=6,zmin=-2,zmax=0.5,thick=3,xtitle='Teff',ytitle='Corr coeff '+cparam[i]+'-'+cparam[j]
    if keyword_set(hard) then begin
      device,/close
      ps2gif,plotfile+'.eps',/eps,/delete
      printf,html,'<a href='+file_basename(plotfile)+'.gif><IMG SRC='+file_basename(plotfile)+'.gif width=250></a>'
    endif
   endif
  endfor
 endif
endfor

if keyword_set(hard) then begin
  printf,html,'</TABLE></BODY></HTML>'
  free_lun,html
  set_plot,'X'
endif

; make individual element plots
if keyword_set(elem_symbol) then elem,str,elem_symbol,/hard

; make residual plots
aspcap_sum,file_search('*/*/aspcapField*'),dir='res'

index:

openw,html,'000index.html',/get_lun
printf,html,'<HTML><BODY>'
printf,html,'<H2> Full sample (if available)</h2>'
printf,html,'<UL>'
printf,html,'<li>allStar file: <A HREF=allStar-'+results_vers+'.fits> allStar-'+results_vers+'.fits</A>'
printf,html,'<li><a href=html/all.html> full sample parameter plots</A>'
printf,html,'<li><a href=html/elem.html> full sample individual element plots</A>'
printf,html,'<li><a href=res/res.html> full sample mean spectra and residuals in Teff/logg grid</A>'
printf,html,'</UL>'
printf,html,'<H2> Calibration subsample (if available)</h2>'
printf,html,'<UL>'
printf,html,'<li>calibration summary file: <A HREF=cal.fits> cal.fits</A>'
printf,html,'<li><a href=html/calplots.html> calibration sample plot index</A>'
printf,html,'</UL>'
printf,html,'</BODY></HTML>'

free_lun,html

end
