; procedure to make plots and web pages for a single plate
;   showing a table with individual frames, zeropoints, S/N
;   plots, variations across frame etc
; also works for a combined image by specifying ims=[0]
; also makes a web page/plots for each individual exposure,
;   show each individual spectrum (unless /noplot is specified
; also makes a summary web page for the plate

pro plotmag,ims,plate,cmjd=cmjd,flat=flat,clobber=clobber,starfiber=starfiber,starnames=starnames,noplot=noplot,mapname=mapname,starmag=starmag,onem=onem,fixfiberid=fixfiberid,badfiberid=badfiberid,survey=survey,mapper_data=mapper_data

clobber=1

; set up directory names
if not keyword_set(cmjd) then cmjd=getcmjd(ims[0])
if size(plate,/type) eq 7 then cplate=plate else cplate=strtrim(string(format='(i)',plate),2)
dirs=getdir(apodir,caldir,spectrodir,vers,apred_vers=apred_vers)
reddir=spectrodir+'red/'+cmjd
telescope=dirs.telescope
platedir=apogee_filename('Plate',plate=plate,mjd=cmjd,chip='a',/dir)
outdir=platedir+'/plots/'
if file_test(outdir,/directory) eq 0 then file_mkdir,outdir
htmldir=platedir+'/html/'
if file_test(htmldir,/directory) eq 0 then file_mkdir,htmldir

; open the output HTML file for this plate
if keyword_set(flat) then $
  file=cplate+'-'+cmjd+'flat' $
;else if keyword_set(starfiber) then $
;  file=cplate+'-'+cmjd+'fiber' $
else if keyword_set(onem) then $
  file=cmjd+'-'+starnames[0] $
else $
  file=cplate+'-'+cmjd
platefile=file
if ims[0] eq 0 then file='sum'+file
openw,html,/get_lun,htmldir+file+'.html'
openw,htmlsum,/get_lun,htmldir+file+'sum.html'
printf,html,'<HTML><BODY>'
printf,htmlsum,'<HTML><BODY>'
if not keyword_set(starfiber) then begin
  printf,html,'Left plots: red are targets, blue are telluric. Observed mags are calculated from median value of green chip. Zeropoint gives overall throughput: bigger number is more throughput.'
  printf,html,'<br>First spatial plots: circles are objects, squares are tellurics, crosses are sky fibers. Colors give deviation of observed mag from expected 2MASS mag using the median zeropoint; red is brighter'
  printf,html,'<br>Second spatial plots: circles are sky fibers. Colors give sky line brightness relative to plate median sky line brightness'
endif
if not keyword_set(starfiber) then begin
  printf,html,'<TABLE BORDER=2>'
  printf,html,'<TR><TD>Frame<TD>Nreads<TD>Zeropoints<TD>Mag plots'
  printf,html,'<TD>Spatial mag deviation'
  printf,html,'<TD>Spatial sky 16325A emission deviations (filled: sky, open: star)'
  printf,html,'<TD>Spatial sky continuum emission '
  printf,html,'<TD>Spatial sky telluric CO2 absorption deviations (filled: H &lt 10) '
endif else begin
  printf,html,'<TABLE BORDER=2>'
  printf,html,'<TR><TD>Frame<TD>Fiber<TD>Star'
endelse
printf,htmlsum,'<TABLE BORDER=2>'
printf,htmlsum,'<TR bgcolor=lightgreen><TD>Frame<TD>Plate<TD>Cart<TD>sec z<TD>HA<TD>DESIGN HA<TD>seeing<TD>FWHM<TD>GDRMS<TD>Nreads<TD>Dither<TD>Zero<TD>Zerorms<TD>Zeronorm<TD>sky continuum<TD>S/N<TD>S/N(c)<TD>unplugged<TD>faint'

; define a filled circle point type
A = FINDGEN(17) * (!PI*2/16.)  

; get the fiber association for this plate
if ims[0] eq 0 then tot=apread('Plate',mjd=cmjd,plate=plate) else tot=apread('1D',mjd=cmjd,num=ims[0])
if size(tot,/type) ne 8 then begin
  printf,html,'<FONT COLOR=red> PROBLEM/FAILURE WITH: ',ims[0]
  printf,htmlsum,'<FONT COLOR=red> PROBLEM/FAILURE WITH: ',ims[0]
  free_lun,html
  free_lun,htmlsum
  stop,'Error in plotmag'
  return
endif
if keyword_set(mapname) then $
  if mapname[0] eq 'header' then plugid=sxpar(tot[0].hdr,'NAME') else plugid=mapname[0] $
else plugid=sxpar(tot[0].hdr,'NAME')
;fiber=getfiber(plate,cmjd,plugid=plugid,/mapa,starfiber=starfiber)

if keyword_set(onem) then begin
  telescope='apo1m'
  reduction_id=starnames[0]
  platedata=getplatedata(cplate,cmjd,plugid=plugid,obj1m=starnames[0],starfiber=starfiber,fixfiberid=fixfiberid) 
endif else $
  platedata=getplatedata(cplate,cmjd,plugid=plugid,fixfiberid=fixfiberid,badfiberid=badfiberid,mapper_data=mapper_data) 
gd=where(platedata.fiberdata.fiberid gt 0)
fiber=platedata.fiberdata[gd]
nfiber=n_elements(fiber)
rows=300-fiber.fiberid
guide=platedata.guidedata
add_tag,fiber,'sn',fltarr(n_elements(ims),3),fiber
add_tag,fiber,'obsmag',fltarr(n_elements(ims),3),fiber

unplugged=where(fiber.fiberid lt 0,nunplugged)
if keyword_set(flat) then begin
  fiber.hmag=12
  fiber.object='FLAT'
endif
;if keyword_set(starfiber) then begin
;  fiber[300-starfiber].hmag=starmag[0]
;  fiber.object='STAR'
;  fiber.objtype='SKY'
;  fiber[300-starfiber].objtype='STAR_BHB'
;  if keyword_set(starnames) then fiber[300-starfiber].object=starnames
;endif
fibertelluric=where(fiber.objtype eq 'SPECTROPHOTO_STD' or fiber.objtype eq 'HOT_STD',ntelluric)
telluric=rows[fibertelluric]
fiberobj=where(fiber.objtype eq 'STAR_BHB' or fiber.objtype eq 'STAR' or fiber.objtype eq 'EXTOBJ',nobj)
;obj=300-fiber[obj].fiberid
obj=rows[fiberobj]
fibersky=where(fiber.objtype eq 'SKY')
;sky=300-fiber[sky].fiberid
sky=rows[fibersky]

; define skylines structure which we will use to get crude sky levels in lines
temp={w1: 16230., w2:16240., c1: 16215., c2: 16225., c3: 16245., c4: 16255.,$
      flux: fltarr(nfiber), type: 1}
skylines=replicate(temp,2)
skylines[1].w1=15990. & skylines[1].w2=16028.
skylines[1].c1=15980. & skylines[1].c2=15990.
skylines[1].c3=0. & skylines[1].c4=0.
skylines[1].type=0

; loop through all the images for this plate, and make the plots
; load up and save information for this plate in a FITS table
allsky=fltarr(n_elements(ims),3)
allzero=fltarr(n_elements(ims),3)
allzerorms=fltarr(n_elements(ims),3)
ra=sxpar(tot[0].hdr,'RADEG')
dec=sxpar(tot[0].hdr,'DECDEG')
mjd=0L
reads,cmjd,mjd
; get moon information for this observation
moonpos,2400000+mjd,ramoon,decmoon
gcirc,2,ra,dec,ramoon,decmoon,moondist
moondist/=3600.
mphase,2400000+mjd,moonphase
; get guider information
if ~keyword_set(onem) then gcam=get_gcam(cmjd)
mjd0=99999
mjd1=0.
; FITS table structure
tab = {telescope: telescope, plate: plate, nreads: 0, dateobs: '', secz: 0., ha: -99., design_ha: [-99.,-99.,-99.], seeing: 0., fwhm: 0., gdrms: 0., cart: 0, plugid: plugid, dither: 0., mjd: mjd, im: 0L, zero: 0., zerorms: 0., zeronorm: 0., sky: [0.,0.,0.], $
  sn: [0.,0.,0.], snc: [0.,0.,0.], altsn: [0.,0.,0.], nsn: 0, snratio: 0., moondist: moondist, moonphase: moonphase, tellfit: fltarr(6,3)} 
platetab=replicate(tab,n_elements(ims))
for i=0,n_elements(ims)-1 do begin
 ; read image
 if ims[i] eq 0 then file=file_basename(apogee_filename('Plate',plate=cplate,mjd=cmjd,/nochip,/base),'.fits') else $
   file=file_basename(apogee_filename('1D',num=ims[i],/nochip,/base),'.fits')
 if keyword_set(clobber) or not file_test(outdir+file+'.tab') then begin
  if ims[0] eq 0 then d=apread('Plate',mjd=cmjd,plate=plate) else d=apread('1D',mjd=cmjd,num=ims[i])
  if size(d,/type) ne 8 then goto,badim
  if ims[0] eq 0 then cframe=apread('Plate',mjd=cmjd,plate=plate) else begin
    cframefile=apogee_filename('Cframe',plate=plate,mjd=cmjd,num=ims[i],chip='c')
    if file_test(cframefile[0]) then cframe=apread('Cframe',plate=plate,mjd=cmjd,num=ims[i]) else cframe=0
  endelse
  obs=fltarr(nfiber,3)
  sn=fltarr(nfiber,3)
  snc=fltarr(nfiber,3)
  snt=fltarr(nfiber,3)
  ;openw,out,/get_lun,outdir+file+'.dat'
  ; create the html file for plotting spectra for each frame
  ;if keyword_set(starfiber) then begin
  ;  objhtml=html 
  ;endif else begin
    openw,objhtml,/get_lun,htmldir+file+'.html'
    printf,objhtml,'<HTML>'
    printf,objhtml,'<HEAD><script type=text/javascript src=../../../../html/sorttable.js></script></head>'
    printf,objhtml,'<BODY>'
    if ims[i] eq 0 then begin
      printf,objhtml,'<H2>'+file+'</H2>'
      platefile=apogee_filename('Plate',plate=cplate,mjd=cmjd,chip=['a','b','c'],/file)
      for ichip=0,2 do printf,objhtml,'<A HREF=../'+platefile[ichip]+'>'+platefile[ichip]+'</A>'
    endif else begin
      printf,objhtml,'<H2>'+string(ims[i])+'</H2>'
      if ~keyword_set(noplot) then begin
        printf,objhtml,'<A HREF=../../../../red/'+cmjd+'/html/'+file+'.html> 1D frames </A>'
        printf,objhtml,'<BR><A HREF=../../../../red/'+cmjd+'/html/ap2D-'+string(format='(i8.8)',ims[i])+'.html> 2D frames </A>'
      endif
    endelse
    printf,objhtml,'<TABLE BORDER=2 CLASS=sortable>'
    printf, objhtml,'<TR><TD>Fiber<TD>Star<TD>H mag<TD>Diff<TD>S/N<TD>S/N (cframe)<TD>Target flags'
  ;endelse
  stars=reverse(indgen(300))
  ;if not keyword_set(starfiber) then $
  ;  stars=reverse(indgen(300)) $
  ;else stars=300-starfiber[i]
  for j=0,nfiber-1 do begin
    ; for each fiber, get an observed mag from a median value
    for ichip=0,2 do obs[j,ichip]=median(d[ichip].flux[*,rows[j]])
  endfor
  if not keyword_set(flat) then begin
    for iline=0,n_elements(skylines)-1 do begin
      skyline=skylines[iline]
      getflux,d,skyline,rows
      skylines[iline]=skyline
    endfor
  endif
  ; get a "magnitude" for each fiber from a median on each chip
  ; do a crude sky subtraction, calculate S/N
  for ichip=0,2 do begin
    if ims[0] eq 0 then medsky=0. else medsky=median(obs[fibersky,ichip])
    if nobj gt 0 then obs[fiberobj,ichip]=median(d[ichip].flux[*,obj],dim=1)-medsky
    if ntelluric gt 0 then obs[fibertelluric,ichip]=median(d[ichip].flux[*,telluric],dim=1)-medsky
    if nobj gt 0 then begin
     sn[fiberobj,ichip]=median((d[ichip].flux[*,obj]-medsky)/d[ichip].err[*,obj],dim=1)
     if n_elements(cframe) gt 1 then $
       snc[fiberobj,ichip]=median(cframe[ichip].flux[*,obj]/cframe[ichip].err[*,obj],dim=1)
    endif
    if ntelluric gt 0 then begin
      sn[fibertelluric,ichip]=median((d[ichip].flux[*,telluric]-medsky)/d[ichip].err[*,telluric],dim=1)
      if n_elements(cframe) gt 1 then begin
         snc[fibertelluric,ichip]=median(cframe[ichip].flux[*,telluric]/cframe[ichip].err[*,telluric],dim=1)
         medfilt=medfilt2d(cframe[ichip].flux[*,telluric],50,dim=1)
         sz=size(cframe[ichip].flux)
         i1=900*sz[1]/2048 & i2=1000*sz[1]/2048
         for itell=0,n_elements(telluric)-1 do $
           snt[fibertelluric[itell],ichip]=$
          MEAN(cframe[ichip].flux[i1:i2,telluric[itell]])/$
          STDDEV(cframe[ichip].flux[i1:i2,telluric[itell]]-medfilt[i1:i2,itell])
       endif else begin
         snc[fibertelluric,ichip]=sn[fibertelluric,ichip]
         medfilt=medfilt2d(d[ichip].flux[*,telluric],50,dim=1)
         sz=size(d[ichip].flux)
         i1=900*sz[1]/2048 & i2=1000*sz[1]/2048
         for itell=0,n_elements(telluric)-1 do $
           snt[fibertelluric[itell],ichip]=$
             MEAN(d[ichip].flux[i1:i2*sz[1]/2048,telluric[itell]])/$
             STDDEV(d[ichip].flux[i1:i2,telluric[itell]]-medfilt[i1:i2,itell])
       endelse
    endif
  endfor

  ; calculate zeropoints from known H band mags
  ; use a static zeropoint to calculate sky brightness
  nreads=sxpar(d[0].hdr,'NFRAMES')
  exptime=sxpar(d[0].hdr,'EXPTIME')
  skyzero=14.75+2.5*alog10(nreads)
  zero=0 & zerorms=0. & faint=-1 & nfaint=0 & achievedsn=[0.,0.,0.] & achievedsnc=[0.,0.,0.] & altsn=[0.,0.,0.] & nsn=0
  ;if not keyword_set(starfiber) then begin
    ;good=where(fiber.hmag ne 0 and fiber.hmag lt 20)
    ;igood=300-fiber[good].fiberid
    ;igood=rows[good]
    zero=median(fiber[[fiberobj,fibertelluric]].hmag+2.5*alog10(obs[[fiberobj,fibertelluric],1]))
    zerorms=robust_sigma(fiber[[fiberobj,fibertelluric]].hmag+2.5*alog10(obs[[fiberobj,fibertelluric],1]))
    faint=where((fiber[[fiberobj,fibertelluric]].hmag+2.5*alog10(obs[[fiberobj,fibertelluric],1])-zero) lt -0.5,nfaint)
  ;endif else begin
  ;  good=where(fiber.hmag ne 0 and fiber.hmag lt 20,ngood)
  ;  if ngood gt 0 then begin
  ;    zero=median(fiber[fiberobj[good]].hmag+2.5*alog10(obs[fiberobj[good],1]))
  ;    zerorms=robust_sigma(fiber[fiberobj[good]].hmag+2.5*alog10(obs[fiberobj[good],1]))
  ;    faint=where((fiber[fiberobj[good]].hmag+2.5*alog10(obs[fiberobj[good],1])-zero) lt -0.5,nfaint)
  ;  endif
  ;endelse
  zeronorm=zero-2.5*alog10(nreads)

  ; for each star, create the exposure entry on the web page and set up the plot of the spectrum
  openw,cfile,/get_lun,outdir+file+'.csh'
  ;for jj=0,n_elements(stars)-1 do begin
  ;  j=stars[jj]
  jsort=sort(fiber.fiberid)
  for jj=0,n_elements(fiber)-1 do begin
    j=jsort[jj]
    ;printf,out,format='(2i4,4f12.6,2f8.3,2x,a,f10.1)',j,fiber[j].fiberid,fiber[j].ra,fiber[j].dec,fiber[j].eta,fiber[j].zeta,fiber[j].hmag,-2.5*alog10(obs[j,1]),fiber[j].objtype

    ; html page entry
    printf,objhtml,'<TR>'
    ;if keyword_set(starfiber) then begin
    ;  printf,objhtml,'<TD><A HREF=../../../../red/'+cmjd+'/html/'+file+'.html>'+string(ims[i])+'</A>'
    ;  printf,objhtml,'<A HREF=../../../../red/'+cmjd+'/html/ap2D-'+string(format='(i8.8)',ims[i])+'.html> (2D frames) </A>'
    ;endif
    color='white'
    if fiber[j].objtype eq 'SPECTROPHOTO_STD' or fiber[j].objtype eq 'HOT_STD' then color='lightblue'
    if fiber[j].objtype eq 'SKY' then color='lightgreen'
    visitfile=apogee_filename('Visit',plate=cplate,mjd=cmjd,fiber=fiber[j].fiberid,reduction=reduction_id,/base)
    if ims[0] eq 0 then $
    printf,objhtml,'<TD><A HREF=../'+visitfile+'>'+string(format='(i3)',fiber[j].fiberid)+'</A>' $
    else $
    printf,objhtml,'<TD>'+string(format='(i3)',fiber[j].fiberid)
    if ims[0] eq 0 then $
      printf,objhtml,'<TD BGCOLOR='+color+'><a href=../plots/'+file_basename(visitfile,'.fits')+'.gif>'+fiber[j].object+'</A>' $
    else printf,objhtml,'<TD BGCOLOR='+color+'>'+fiber[j].object
    rastring=stringize(fiber[j].ra,ndec=5)
    decstring=stringize(fiber[j].dec,ndec=5)
    if fiber[j].object ne 'sky' and fiber[j].fiberid ge 0 then $
    printf,objhtml,'<BR><A HREF="http://simbad.cfa.harvard.edu/simbad/sim-basic?Ident='+RASTRING+'+%09'+DECSTRING+'++&submit=SIMBAD+search"> (SIMBAD) </A>'
    ;if not keyword_set(starfiber) then begin
      printf,objhtml,'<TD>'+string(format='(f8.2)',fiber[j].hmag)
      printf,objhtml,'<TD>'+string(format='(f8.2)',fiber[j].hmag+2.5*alog10(obs[j,1])-zero)
      printf,objhtml,'<TD>'+string(format='(f8.2)',sn[j,1])
      printf,objhtml,'<TD>'+string(format='(f8.2)',snc[j,1])
      printf,objhtml,'<TD>'+targflag(fiber[j].target1,fiber[j].target2,fiber[j].target3,survey=survey)
      if ims[0] eq 0 and fiber[j].fiberid ge 0 then begin
        ;vfile=platedir+'/apVisit-'+apred_vers+'-'+cplate+'-'+cmjd+'-'+string(format='(i3.3)',fiber[j].fiberid)+'.fits'
        vfile=apogee_filename('Visit',plate=cplate,mjd=cmjd,fiber=fiber[j].fiberid,reduction=reduction_id)
        if file_test(vfile) then begin
          h=headfits(vfile)
          if size(h,/type) eq 7 then printf,objhtml,'<BR>'+starflag(sxpar(h,'STARFLAG'))
        endif
      endif
    ;endif
 
    ; plot of spectrum?
    ;if (j mod 300 gt -1 or keyword_set(starfiber)) then begin
    if j mod 300 gt -1  then begin
;      if fiber[j].objtype ne 'SKY' then begin
       if not keyword_set(noplot) then begin
        set_plot,'ps'
        !p.multi=[0,0,0]
        device,file=outdir+file+'-'+string(format='(i3.3)',rows[j])+'.eps',/encap,xsize=48,ysize=6,/color,/in
        smcolor
        if n_elements(d[0].wave) gt 10 then begin
          ymin=2*min([obs[j,0],obs[j,1],obs[j,2],0.])
          ymax=2*max([obs[j,0],obs[j,1],obs[j,2]])
          yr=[ymin,ymax]
          if fiber[j].objtype eq 'SKY' then yr=[-50,50]
          ;plot,d[0].wave[*,rows[j]],d[0].flux[*,rows[j]],xrange=[15000,17000],xstyle=1,yrange=[0,2*max([obs[j,0],obs[j,1],obs[j,2]])],color=2
          plot,d[0].wave[*,rows[j]],d[0].flux[*,rows[j]],xrange=[15000,17000],xstyle=1,yrange=yr,color=2
          oplot,d[1].wave[*,rows[j]],d[1].flux[*,rows[j]],color=3
          oplot,d[2].wave[*,rows[j]],d[2].flux[*,rows[j]],color=4
          oplot,d[0].wave[*,rows[j]],d[0].err[*,rows[j]]*10,color=0,linestyle=1
          oplot,d[1].wave[*,rows[j]],d[1].err[*,rows[j]]*10,color=0,linestyle=1
          oplot,d[2].wave[*,rows[j]],d[2].err[*,rows[j]]*10,color=0,linestyle=1
          if n_elements(cframe) gt 1 then begin
            oplot,cframe[0].wave[*,rows[j]],cframe[0].flux[*,rows[j]],color=5
            oplot,cframe[1].wave[*,rows[j]],cframe[1].flux[*,rows[j]],color=5
            oplot,cframe[2].wave[*,rows[j]],cframe[2].flux[*,rows[j]],color=5
            if fiber[j].objtype eq 'SKY' then begin
              oplot,cframe[0].wave[*,rows[j]],cframe[0].flux[*,rows[j]]*0.
              oplot,cframe[1].wave[*,rows[j]],cframe[1].flux[*,rows[j]]*0.
              oplot,cframe[2].wave[*,rows[j]],cframe[2].flux[*,rows[j]]*0.
            endif else begin
              for ichip=0,2 do begin
                tmp=cframe[ichip].flux[*,rows[j]]
                bd=where(cframe[ichip].mask[*,rows[j]] and badmask(),nbd)
                if nbd gt 0 then tmp[bd] = !values.f_nan
                oplot,cframe[ichip].wave[*,rows[j]], tmp, color=3
                bd=where(cframe[ichip].mask[*,rows[j]] and maskval('SIG_SKYLINE'),nbd)
                if nbd gt 0 then tmp[bd] = !values.f_nan
                oplot,cframe[ichip].wave[*,rows[j]], tmp, color=4
              endfor
            endelse
          endif
        endif else begin
          x=indgen(2048)
          plot,x,d[0].flux[*,rows[j]],xrange=[0,2048*3+20],xstyle=1,yrange=[0,2*obs[j,1]],color=2
          x+=2058
          oplot,x,d[1].flux[*,rows[j]],color=3
          x+=2058
          oplot,x,d[2].flux[*,rows[j]],color=4
        endelse
        device,/close
        set_plot,'X'
;      spawn,'convert '+outdir+file+'-'+string(format='(i3.3)',j)+'.eps '+outdir+file+'-'+string(format='(i3.3)',j)+'.gif &'
        infile=outdir+file+'-'+string(format='(i3.3)',rows[j])+'.eps'
        outfile=outdir+file+'-'+string(format='(i3.3)',rows[j])+'.jpg'
        printf,cfile,'echo '+infile
        printf,cfile,'gs -sDEVICE=jpeg -sOutputFile='+outfile+' -r50 -dBATCH -dNOPAUSE -dDEVICEWIDTHPOINTS=3450 -dDEVICEHEIGHTPOINTS=431 '+infile+'  >& /dev/null'
        printf,cfile,'''rm'' '+infile
        printf,cfile,'chmod 664 '+outfile
        jfile='../plots/'+file+'-'+string(format='(i3.3)',rows[j])+'.jpg'
        printf,objhtml,'<TD><A HREF='+jfile+'><IMG SRC='+jfile+' WIDTH=1000></A>'
       endif else $
        printf,objhtml,'<TD>No plots for individual exposures, see plate plots'
;      endif 
    endif
  endfor
  ;if not keyword_set(starfiber) then free_lun,objhtml
  free_lun,objhtml
  ;free_lun,out
  free_lun,cfile
  ; do the JPG conversions in one batch (still takes a little while)
  if not keyword_set(noplot) then spawn,'csh '+outdir+file+'.csh'
  file_delete,outdir+file+'.csh',/allow_nonexistent

  ; create a series of plots for this exposure
  set_plot,'ps'
  device,file=outdir+file+'.eps',/encap,xsize=16,ysize=16,/color
  thick=3
  !p.thick=thick
  !P.CharThick =thick
  !X.Thick =thick
  !Y.Thick =thick
  !Z.Thick =thick
  smcolor
  usersym,cos(A),sin(A),/fill  
  if not keyword_set(flat) and not keyword_set(onem) then begin
    ; 5 panel plot:
    ;   observed mag vs H mag
    !p.multi=[0,1,5]
    plot,fiber.hmag,-2.5*alog10(obs[*,1]),psym=8,xrange=[7,16],xtitle='2MASS H',ytitle='mag=-2.5 log(cnts)'
    oplot,fiber.hmag,-2.5*alog10(obs[*,1]),psym=8,color=2
    if ntelluric gt 0 then oplot,fiber[telluric].hmag,-2.5*alog10(obs[telluric,1]),psym=8,color=4

    ;   observed mag - fit mag vs H mag
    plot,fiber.hmag,fiber.hmag+2.5*alog10(obs[*,1])-zero,psym=8,xrange=[7,16],ytitle='H-(mag+zero)',yrange=[-1.5,1.5]
    oplot,fiber.hmag,fiber.hmag+2.5*alog10(obs[*,1])-zero,psym=8,color=2
    if ntelluric gt 0 then oplot,fiber[telluric].hmag,fiber[telluric].hmag+2.5*alog10(obs[telluric,1])-zero,psym=8,color=4
    ;xyouts,0.2,0.4,'Z='+string(format='(f5.2)',zero),/normal
    ;xyouts,0.2,0.35,'Znorm='+string(format='(f5.2)',zeronorm),/normal

    ;   S/N as calculated from ap1D frame
    plot,fiber.hmag,sn,psym=8,xrange=[7,16],ytitle='S/N',yrange=[1,1000],/ylog
    oplot,fiber.hmag,sn[*,0],psym=8,color=2
    oplot,fiber.hmag,sn[*,1],psym=8,color=3
    oplot,fiber.hmag,sn[*,2],psym=8,color=4
    ; plot target line that has S/N=100 for 3 hour exposure at H=12.2
    sntarget=100*sqrt(exptime/(3.*3600))
    sntargetmag=12.2
    ;if ims[0] eq 0 then sntarget=100*sqrt(10*nreads/180) $
    ;else sntarget=100*sqrt(10./180) 
    oplot,[sntargetmag-10,sntargetmag+2.5],[sntarget*100,sntarget/sqrt(10)]
    ; get typical S/N for this plate
    snstars=where(fiber.hmag gt 12 and fiber.hmag lt 12.2,nsn)
    scale=1
    if nsn lt 3 then begin
      bright=where(fiber.hmag lt 12)
      hmax=max(fiber[bright].hmag)
      snstars=where(fiber.hmag gt hmax-0.2 and fiber.hmag le hmax,nsn)
      scale=sqrt(10^(0.4*(hmax-12.2)))
    endif
    achievedsn=median(sn[snstars,*],dimension=1)*scale
    ; alternative S/N as computed from median of all stars with H<12.2, scaled
    snstars=where(fiber.hmag lt 12.2)
    scale=sqrt(10^(0.4*(fiber[snstars].hmag-12.2)))
    altsn=achievedsn*0.
    for ichip=0,2 do altsn[ichip]=median(sn[snstars,ichip]*scale)

    ;   S/N as calculated from apCframe
    plot,fiber.hmag,snc,psym=8,xrange=[7,16],ytitle='S/N',yrange=[1,1000],/ylog
    oplot,fiber.hmag,snc[*,0],psym=8,color=2
    oplot,fiber.hmag,snc[*,1],psym=8,color=3
    oplot,fiber.hmag,snc[*,2],psym=8,color=4
    ; plot target line that has S/N=100 for 3 hour exposure at H=12.2
    ;if ims[0] eq 0 then sntarget=100*sqrt(10*nreads/180) $
    ;else sntarget=100*sqrt(10./180) 
    sntargetmag=12.2
    oplot,[sntargetmag-10,sntargetmag+2.5],[sntarget*100,sntarget/sqrt(10)]
    achievedsnc=median(snc[snstars,*],dimension=1)*scale

    ; S/N from scatter relative to calculated S/N for tellurics
    if ntelluric gt 0 then begin
      plot,fiber[telluric].hmag,snt[telluric,1]/snc[telluric,1],xrange=[7,16],psym=6,yrange=[0,1.5]
      oplot,[7,16],[1,1],linestyle=1
    endif
    device,/close
    ps2gif,outdir+file+'.eps',/eps,chmod='664'o,/delete
  endif else if onem then begin
    achievedsn=median([sn[obj,*]],dim=1)
  endif
   ;else begin
    ;ii=where(fiber.objtype eq 'STAR')
    ;plot,fiber.fiberid,fiber.hmag+2.5*alog10(obs[300-fiber[ii].fiberid,1])-zero,psym=8,xrange=[-5,305],ytitle='H-(mag+zero)',yrange=[-1.5,1.5]
  ;endelse

  ; more plots!
  !p.multi=[0,0,0]
  if not keyword_set(starfiber) and not keyword_set(onem) then begin
    ; spatial plot of residuals
    good=where(fiber.hmag gt 0)
    device,file=outdir+file+'.eps',/encap,xsize=16,ysize=16,/color
    loadct,39
    usersym,cos(A),sin(A),/fill  
    overplot=0
    if dirs.telescope eq 'lco25m' then lim=[-1.2,1.2] else lim=[-1.6,1.6]
    if nobj gt 0 then begin
      plotc,fiber[fiberobj].zeta,fiber[fiberobj].eta,fiber[fiberobj].hmag+2.5*alog10(obs[fiberobj,1])-zero,min=-0.5,max=0.5,$
       xr=lim,yr=lim,ps=8,xtit='Zeta(deg)',ytit='Eta(deg)',xstyle=1,ystyle=1
      overplot=1
    endif
    usersym,[-1,1,1,-1],[-1,-1,1,1],/fill  
    if ntelluric gt 0 then $
    plotc,fiber[fibertelluric].zeta,fiber[fibertelluric].eta,fiber[fibertelluric].hmag+2.5*alog10(obs[fibertelluric,1])-zero,min=-0.5,max=0.5,$
       xr=lim,yr=lim,ps=8,overplot=overplot,xstyle=1,ystyle=1
    if sky[0] ge 0 then oplot,fiber[fibersky].zeta,fiber[fibersky].eta,ps=1
    usersym,cos(A),sin(A)
    oplot,guide.zeta,guide.eta,ps=8,symsize=2
    device,/close
    ps2jpg,outdir+file+'.eps',/eps,chmod='664'o,/delete
  
    ; spatial plot of sky line emission
    device,file=outdir+file+'sky.eps',/encap,xsize=16,ysize=16,/color
    loadct,39
    usersym,cos(A),sin(A),/fill  
    medsky=median(skylines[0].flux[sky])
    plotc,fiber[fibersky].zeta,fiber[fibersky].eta,skylines[0].flux[fibersky]/medsky,min=0.9,max=1.1,$
     xr=lim,yr=lim,ps=8,xstyle=1,ystyle=1
    if nobj gt 0 then $
    plotc,fiber[fiberobj].zeta,fiber[fiberobj].eta,skylines[0].flux[fiberobj]/medsky,min=0.9,max=1.1,$
     xr=lim,yr=lim,ps=6,overplot=1,xstyle=1,ystyle=1
    if ntelluric gt 0 then $
    plotc,fiber[fibertelluric].zeta,fiber[fibertelluric].eta,skylines[0].flux[fibertelluric]/medsky,min=0.9,max=1.1,$
     xr=lim,yr=lim,ps=4,overplot=1,xstyle=1,ystyle=1
    device,/close
    ps2jpg,outdir+file+'sky.eps',/eps,chmod='664'o,/delete

    ; spatial plot of sky continuum emission
    device,file=outdir+file+'skycont.eps',/encap,xsize=16,ysize=16,/color
    loadct,39
    usersym,cos(A),sin(A),/fill  
    plotc,fiber[fibersky].zeta,fiber[fibersky].eta,-2.5*alog10(obs[fibersky,1])+skyzero,min=13,max=15,$
     xr=lim,yr=lim,ps=8,xstyle=1,ystyle=1
    device,/close
    ps2jpg,outdir+file+'skycont.eps',/eps,chmod='664'o,/delete

    ; spatial plot of telluric absorption
;    if ntelluric gt 0 then begin
;      medsky=median(skylines[1].flux[telluric])
;      device,file=outdir+file+'telluric.eps',/encap,xsize=16,ysize=16,/color
;      loadct,39
;      usersym,cos(A),sin(A),/fill  
;      bright=where(fiber[telluric].hmag lt 10)
;      plotc,fiber[telluric].zeta,fiber[telluric].eta,$
;        skylines[1].flux[telluric]/medsky,min=0.5,max=2,$
;        xr=[-1.6,1.6],yr=[-1.6,1.6],ps=4
;      if bright[0] ge 0 then $
;        plotc,fiber[telluric[bright]].zeta,fiber[telluric[bright]].eta,$
;        skylines[1].flux[telluric[bright]]/medsky,min=0.5,max=2,$
;        xr=[-1.6,1.6],yr=[-1.6,1.6],ps=8,overplot=1
;      device,/close
;      ps2jpg,outdir+file+'telluric.eps',/eps,chmod='664'o,/delete
;    endif
  endif
  
 endif

 ; put all of the info and plots on the plate web page
 medsky=fltarr(3)
 for ichip=0,2 do $
    if median(obs[fibersky,ichip]) gt 0 then medsky[ichip]=-2.5*alog10(median(obs[fibersky,ichip]))+skyzero else medsky[ichip]=99.999
 ;if not keyword_set(starfiber) then begin
   printf,html,'<TR><TD><A HREF=../html/'+file+'.html>',ims[i],'</A>'
   printf,html,'<TD>'+string(nreads)
   printf,html,'<TD><TABLE BORDER=1><TD><TD>Red<TD>Green<TD>Blue'
   printf,html,'<TR><TD>z<TD><TD>'+string(format='(f5.2)',zero)
   printf,html,'<TR><TD>znorm<TD><TD>'+string(format='(f5.2)',zeronorm)
   printf,html,'<TR><TD>sky'+$
             string(format='("<TD>",f5.1,"<TD>",f5.1,"<TD>",f5.1)',medsky)
   printf,html,'<TR><TD>S/N'+$
             string(format='("<TD>",f5.1,"<TD>",f5.1,"<TD>",f5.1)',achievedsn)
   printf,html,'<TR><TD>S/N(c)'+$
             string(format='("<TD>",f5.1,"<TD>",f5.1,"<TD>",f5.1)',achievedsnc)
   if ntelluric gt 0 then $
   printf,html,'<TR><TD>SN(E/C)<TD<TD>'+string(format='(f5.2)',median(snt[telluric,1]/snc[telluric,1])) $
   else printf,html,'<TR><TD>SN(E/C)<TD<TD>'
   printf,html,'</TABLE>'
 ;  printf,html,'<TD>z='+string(format='(f5.2)',zero)+$
 ;      '<BR>znorm='+string(format='(f5.2)',zeronorm)+$
 ;      '<BR>medsky='+string(format='("[",f5.2,",",f5.2,",",f5.2,"]")',medsky)+$
 ;      '<BR>S/N='+string(format='("[",f5.1,",",f5.1,",",f5.1,"]")',achievedsn)+$
 ;      '<BR>S/N(c)='+string(format='("[",f5.1,",",f5.1,",",f5.1,"]")',achievedsnc)
   printf,html,'<TD><IMG SRC=../plots/'+file+'.gif>'
   printf,html,'<TD> <IMG SRC=../plots/'+file+'.jpg>'
   printf,html,'<TD> <IMG SRC=../plots/'+file+'sky.jpg>'
   printf,html,'<TD> <IMG SRC=../plots/'+file+'skycont.jpg>'
   printf,html,'<TD> <IMG SRC=../plots/'+file+'telluric.jpg>'
 ;endif

 ; get guider infor
 if ~keyword_set(onem) then begin
  dateobs=sxpar(d[0].hdr,'DATE-OBS')
  mjdstart=date_conv(dateobs,'MODIFIED')
  mjdend=mjdstart+sxpar(d[0].hdr,'EXPTIME')/86400.
  mjd0=min([mjd0,mjdstart])
  mjd1=max([mjd1,mjdend])
  if size(gcam,/type) eq 8 then jcam=where(gcam.mjd gt mjdstart and gcam.mjd lt mjdend,nj) else nj=0
  if nj gt 1 then begin
   fwhm=median(gcam[jcam].fwhm_median) 
   gdrms=median(gcam[jcam].gdrms)
  endif else begin
   fwhm=-1.
   gdrms=-1.
   print,'not halted: no matching mjd range in gcam...'
  endelse
 endif else begin
   fwhm=-1
   gdrms=-1
 endelse

 ; summary plate web page
 printf,htmlsum,'<TR><TD><A HREF=../html/'+file+'.html>',ims[i],'</A>'
 printf,htmlsum,'<TD><A HREF=../../../../plates/'+cplate+'/'+cmjd+'/html/'+cplate+'-'+cmjd+'.html>',sxpar(d[0].hdr,'PLATEID'),'</A>'
 printf,htmlsum,'<TD>',sxpar(d[0].hdr,'CARTID')
 alt=sxpar(d[0].hdr,'ALT',count=count)
 if count gt 0 then secz=1./cos((90.-alt)*!pi/180.) else secz=sxpar(d[0].hdr,'ARMASS')
 seeing=sxpar(d[0].hdr,'SEEING')
 ha=sxpar(d[0].hdr,'HA')
 design_ha=platedata.ha
 dither=-99.
 if n_elements(cframe) gt 1 then dither = sxpar(cframe[0].hdr,'DITHSH') 
 printf,htmlsum,'<TD>',string(format='(f6.2)',secz)
 printf,htmlsum,'<TD>',string(format='(f6.2)',ha)
 printf,htmlsum,'<TD>',string(format='(f6.0,",",f6.0,",",f6.0)',design_ha)
 printf,htmlsum,'<TD>',string(format='(f6.2)',seeing)
 printf,htmlsum,'<TD>',string(format='(f6.2)',fwhm)
 printf,htmlsum,'<TD>',string(format='(f6.2)',gdrms)
 printf,htmlsum,'<TD>'+string(nreads)
 if n_elements(cframe) gt 1 then $
 printf,htmlsum,'<TD>'+string(format='(f8.2)',sxpar(cframe[0].hdr,'DITHSH')) else $
 printf,htmlsum,'<TD>'
 printf,htmlsum,'<TD>',string(format='(f5.2)',zero)
 printf,htmlsum,'<TD>',string(format='(f5.2)',zerorms)
 printf,htmlsum,'<TD>',string(format='(f5.2)',zeronorm)
 printf,htmlsum,'<TD>',string(format='("[",f5.2,",",f5.2,",",f5.2,"]")',medsky)
 printf,htmlsum,'<TD>',string(format='("[",f5.1,",",f5.1,",",f5.1,"]")',achievedsn)
 printf,htmlsum,'<TD>',string(format='("[",f5.1,",",f5.1,",",f5.1,"]")',achievedsnc)
 printf,htmlsum,'<TD>'
 for j=0,nunplugged-1 do printf,htmlsum,300-unplugged[j]
 printf,htmlsum,'<TD>'
 if faint[0] ge 0 then for j=0,nfaint-1 do printf,htmlsum,fiber[faint[j]].fiberid
 allsky[i,*]=medsky
 allzero[i,*]=zero
 allzerorms[i,*]=zerorms
 ; summary information in apPlateSum FITS file
 if ims[i] gt 0 then begin
   tellfile=apogee_filename('Tellstar',plate=cplate,mjd=cmjd,reduction=reduction_id)
   telstr=mrdfits(tellfile,1,status=status)
   if status eq 0 then begin
     jtell=where(telstr.im eq ims[i],ntell)
     if ntell gt 0 then platetab[i].tellfit = telstr[jtell].fitpars
   endif else stop,'Error reading Tellstar file: ',tellfile
 endif
 platetab[i].im = ims[i]
 platetab[i].nreads = nreads
 platetab[i].secz = secz
 platetab[i].ha = ha
 platetab[i].design_ha = design_ha
 platetab[i].seeing = seeing
 platetab[i].fwhm = fwhm
 platetab[i].gdrms = gdrms
 platetab[i].cart = sxpar(d[0].hdr,'CARTID')
 platetab[i].dateobs = sxpar(d[0].hdr,'DATE-OBS')
 platetab[i].dither = dither
 platetab[i].zero = zero
 platetab[i].zerorms = zerorms
 platetab[i].zeronorm = zeronorm
 platetab[i].sky = medsky
 platetab[i].sn = achievedsn
 platetab[i].altsn = altsn
 platetab[i].nsn = nsn
 platetab[i].snc = achievedsnc
 if ntelluric gt 0 then platetab[i].snratio = median(snt[telluric,1]/snc[telluric,1])

 for j=0,n_elements(fiber)-1 do begin
   fiber[j].sn[i,*]=sn[j,*]
   fiber[j].obsmag[i,*]=-2.5*alog10(obs[j,*])+zero
 endfor

 badim:
endfor


; write out the FITS table
platefile=apogee_filename('PlateSum',plate=plate,mjd=cmjd,reduction=reduction_id) 
;if ims[0] gt 0 and not keyword_set(starfiber) then begin
if ims[0] gt 0 then begin
  mwrfits,platetab,platefile,/create
  mwrfits,fiber,platefile
endif
if ims[0] eq 0 then begin
  mwrfits,platetab,platefile
  mwrfits,fiber,platefile
endif

printf,html,'</TABLE>'

; for individual frames, make plots of variation of sky and zeropoint
; for combined frames, make table of combination parameters
!p.multi=[0,0,0]
if keyword_set(onem) then name=starnames[0]+'-'+cmjd else name=cplate+'-'+cmjd
if ims[0] gt 0 then begin
  ; guider rms plot
  if ~keyword_set(onem) then begin
    if size(gcam,/type) eq 8 then begin
      jcam=where(gcam.mjd gt mjd0 and gcam.mjd lt mjd1,nj) 
      device,file=outdir+file+'guider.eps',/encap
      plot,gcam[jcam].mjd,gcam[jcam].gdrms
      device,/close
      file='guider-'+name
      ps2gif,outdir+file+'guider.eps',/eps,chmod='664'o,/delete
    endif
  endif
  ; make plot of sky levels for this plate
  printf,html,'<TABLE BORDER=2><TR>'
  file='sky-'+name
  set_plot,'ps'
  device,file=outdir+file+'.eps',/encap,ysize=8,/color
  plot,ims,allsky,yrange=[max(allsky)+0.3,min(allsky)-0.3],psym=8,xrange=[(ims[0] mod 10000)-1,(ims[n_elements(ims)-1] mod 10000)+1],xtitle='Image number',ytitl='Continuum sky per pixel'
  smcolor
  for i=0,2 do begin
   oplot,ims mod 10000,allsky[*,i],color=i+2,psym=8
  endfor
  device,/close
  ps2gif,outdir+file+'.eps',/eps,chmod='664'o,/delete
  printf,html,'<TD><IMG SRC=../plots/'+file+'.gif>'
  ; make plot of zeropoints for this plate
  file='zero-'+name
  device,file=outdir+file+'.eps',/encap,ysize=8,/color
  plot,ims,allzero,yrange=[max(allzero)+0.3,min(allzero)-0.3],psym=8,xrange=[(ims[0] mod 10000)-1,(ims[n_elements(ims)-1] mod 10000)+1],xtitle='Image number',ytitle='Zeropoint per pixel'
  smcolor
  for i=0,2 do begin
   oplot,ims mod 10000,allzero[*,i],color=i+2,psym=8
  endfor
  device,/close
  ps2gif,outdir+file+'.eps',/eps,chmod='664'o,/delete
  printf,html,'<TD><IMG SRC=../plots/'+file+'.gif>'
  printf,html,'</TABLE>'
endif else begin
  file=apogee_filename('Plate',plate=cplate,mjd=cmjd,chip='a')
  shiftstr=mrdfits(file,13)
  pairstr=mrdfits(file,14)
  npairs=n_elements(pairstr)
  if size(pairstr,/type) eq 8 and npairs gt 0 then begin
    ; Pair table
    printf,html,'<BR><TABLE BORDER=2>'
    printf,html,'<TR><TD>IPAIR<TD>NAME<TD>SHIFT<TD>NEWSHIFT<TD>S/N'
    printf,html,'<TD>NAME<TD>SHIFT<TD>NEWSHIFT<TD>S/N'
    for ipair=0,npairs-1 do begin
      printf,html,'<TR><TD>',ipair
      for j=0,1 do begin
        printf,html,'<TD>',pairstr[ipair].framename[j]
        printf,html,'<TD>',pairstr[ipair].oldshift[j]
        printf,html,'<TD>',pairstr[ipair].shift[j]
        printf,html,'<TD>',pairstr[ipair].sn[j]
      endfor
    endfor
  endif else begin
    ; table of combination parameters
    printf,html,'<BR><TABLE BORDER=2>'
    for iframe=0,n_elements(shiftstr)-1 do begin
      printf,html,'<TR><TD>',shiftstr[iframe].framenum
      printf,html,'<TD>',shiftstr[iframe].shift
      printf,html,'<TD>',shiftstr[iframe].sn
    endfor
  endelse
  printf,html,'</TABLE>'

endelse

printf,html,'</BODY></HTML>'
printf,htmlsum,'</TABLE>'

if keyword_set(onem) then begin
  file=cmjd+'-'+starnames[0]
  printf,htmlsum,'<a href=../plots/apVisit-'+apred_vers+'-'+file+'.gif><IMG src='+'../plots/apVisit-'+apred_vers+'-'+file+'.gif></A>'
endif
printf,htmlsum,'</BODY></HTML>'
free_lun,html
free_lun,htmlsum
set_plot,'X'

end
