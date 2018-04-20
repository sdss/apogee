;+
; APQA
;
;  call routines to make "QA" plots and web pages for a plate/MJD
;  for calibration frames, measures some features and makes a apQAcal file
;    with some summary information about the calibration data
;-
pro apqa,planfile,noplot=noplot,mapper_data=mapper_data

; load the planfile structure
aploadplan,planfile,planstr
dirs=getdir(apodir,caldir,spectrodir,vers)

; find all of the object frames and load IDs into ims array
objs=where(planstr.apexp.flavor eq 'object')
ims=lonarr(n_elements(objs))
for iframe=0,n_elements(objs)-1 do begin
;  reads,strmid(planstr.apexp[objs[iframe]].name[0],6,8),im
  im=0L
  reads,planstr.apexp[objs[iframe]].name[0],im
  ims[iframe]=im
endfor
;plateid=0
;reads,planstr.plateid,plateid
plateid=planstr.plateid
if size(plateid,/type) ne 7  then cplate=strtrim(string(format='(i6.4)',plateid),2) else cplate=strtrim(plateid)

cmjd=getcmjd(ims[0],mjd=mjd)

; check for tag for fixing fiberid in plugmap files
if tag_exist(planstr,'fixfiberid') then fixfiberid=planstr.fixfiberid
if tag_exist(planstr,'badfiberid') then badfiberid=planstr.badfiberid
if tag_exist(planstr,'survey') then survey=planstr.survey else survey='apogee'

; make plots for normal and single plates
if planstr.platetype eq 'normal' then begin
  plotmag,ims,plateid,/clobber,mapname=planstr.plugmap,/noplot,fixfiberid=fixfiberid,badfiberid=badfiberid,survey=survey,mapper_data=mapper_data
  plotmag,0,plateid,cmjd=cmjd,/clobber,mapname=planstr.plugmap,noplot=noplot,fixfiberid=fixfiberid,badfiberid=badfiberid,survey=survey,mapper_data=mapper_data
  plotflux,planfile
  mkhtmlplate,plateid,mjd,fluxid=planstr.fluxid
  platefile=apogee_filename('PlateSum',plate=plateid,mjd=cmjd)
  sntab,tabs=platefile,outfile=platefile+'.dat'
endif else if planstr.platetype eq 'single' then begin
  ;plotmag,ims,plateid,starfiber=planstr.apexp[objs].single,starnames=planstr.apexp[objs].singlename,/clobber,mapname=planstr.plugmap,noplot=noplot,/onem,starmag=planstr.hmag
  plotmag,ims,plateid,starnames=planstr.apexp[objs].singlename,starfiber=planstr.apexp[objs[0]].single,fixfiberid=fixfiberid,/clobber,mapname=planstr.plugmap,noplot=noplot,/onem,starmag=planstr.hmag,survey=survey
endif

; for calibration plates, measure lamp brightesses and/or line widths, etc.
if planstr.platetype eq 'cal' then begin
 nlines=2
 if dirs.instrument eq 'apogee-s' then begin
   tharline=[[944.,1112.,1102.],[1726.,608.,1745.]]
   uneline=[[607.,1229.,1088.],[1765.,620.,1860.]]
 endif else begin
   tharline=[[940.,1128.,1130.],[1724.,623.,1778.]]
   uneline=[[603.,1213.,1116.],[1763.,605.,1893.]]
 endelse
 fibers=[10,80,150,220,290]
 nfibers=n_elements(fibers)
 nchips=3
 tmp={name: ' ', mjd: ' ', jd: 0.d0, nframes: 0, nread: 0, exptime: 0., qrtz: 0, une: 0, thar: 0, flux: fltarr(300,nchips), gauss: fltarr(4,nfibers,nchips,nlines), wave: dblarr(nfibers,nchips,nlines), fibers: fibers, lines: fltarr(nchips,nlines)}
 str=replicate(tmp,n_elements(planstr.apexp))
 for i=0,n_elements(planstr.apexp)-1 do begin
  a=apread('1D',num=long(planstr.apexp[i].name))
  if size(a,/type) eq 8 then begin
   if ~tag_exist(a,'flux',index=fluxid) then junk =tag_exist(a,'data',index=fluxid)
   str[i].name=planstr.apexp[i].name
   str[i].mjd=planstr.mjd
   str[i].jd=sxpar(a[0].hdr,'JD-MID')
   str[i].nframes=sxpar(a[0].hdr,'NFRAMES')
   str[i].nread=sxpar(a[0].hdr,'NREAD')
   str[i].exptime=sxpar(a[0].hdr,'EXPTIME')
   str[i].qrtz=sxpar(a[0].hdr,'LAMPQRTZ')
   str[i].thar=sxpar(a[0].hdr,'LAMPTHAR')
   str[i].une=sxpar(a[0].hdr,'LAMPUNE')

   ; quartz exposures
   if str[i].qrtz eq 1 then begin
     str[i].flux=median(a.(fluxid),dim=1)
   endif

   ; arc lamp exposures
   if str[i].thar eq 1  or str[i].une eq 1 then begin
    if str[i].thar eq 1 then line=tharline else line=uneline
    str[i].lines=line
    sz=size(line)
    if sz[0] eq 1 then nlines=1 else nlines=sz[2]
    for iline=0,nlines-1 do begin
     for ichip=0,nchips-1 do begin
       print,'calling appeakfit'
       appeakfit,a[ichip],linestr,fibers=fibers,nsigthresh=10
       for ifiber=0,nfibers-1 do begin
         fiber=fibers[ifiber]
         j=where(linestr.fiber eq fiber,nj)
         if nj gt 0 then begin
           junk=min(abs(linestr[j].gaussx-line[ichip,iline]),jline)
           str[i].gauss[*,ifiber,ichip,iline] = linestr[j[jline]].gpar
           sz=size(a[ichip].wcoef,/dim)
           if sz[0] eq 2 then str[i].wave[ifiber,ichip,iline] = pix2wave(linestr[j[jline]].gaussx,a[ichip].wcoef[fiber,*])
           str[i].flux[fiber,ichip] = linestr[j[jline]].sumflux
         endif
       endfor
     endfor
    endfor
   endif
  endif
 endfor
 outfile=apogee_filename('QAcal',mjd=planstr.mjd)
 mwrfits,str,outfile,/create
endif

; dark and flats
if planstr.platetype eq 'dark' then begin
 nchips=3
 nquad=4
 tmp={name: ' ', mjd: ' ', jd: 0.d0, nframes: 0, nread: 0, exptime: 0., qrtz: 0, une: 0, thar: 0, exptype: '', mean: fltarr(nchips,nquad), sig: fltarr(nchips,nquad)}
 str=replicate(tmp,n_elements(planstr.apexp))
 for i=0,n_elements(planstr.apexp)-1 do begin
  a=apread('2D',num=long(planstr.apexp[i].name))
  if size(a,/type) eq 8 then begin
   if ~tag_exist(a,'flux',index=fluxid) then junk =tag_exist(a,'data',index=fluxid)
   str[i].name=planstr.apexp[i].name
   str[i].mjd=planstr.mjd
   str[i].jd=sxpar(a[0].hdr,'JD-MID')
   str[i].nframes=sxpar(a[0].hdr,'NFRAMES')
   str[i].nread=sxpar(a[0].hdr,'NREAD')
   str[i].exptime=sxpar(a[0].hdr,'EXPTIME')
   str[i].qrtz=sxpar(a[0].hdr,'LAMPQRTZ')
   str[i].thar=sxpar(a[0].hdr,'LAMPTHAR')
   str[i].une=sxpar(a[0].hdr,'LAMPUNE')
   str[i].exptype=sxpar(a[0].hdr,'EXPTYPE')
   ; for darks/flats, get mean and stdev of column-medianed quadrants
   for ichip=0,nchips-1 do begin
    i1=10
    i2=500
    for iquad=0,3 do begin
      sm=median(a[ichip].flux[i1:i2,10:2000],dim=1)
      str[i].mean[ichip,iquad] = mean(sm)
      str[i].sig[ichip,iquad] = stddev(sm)
      i1+=512
      i2+=512
    endfor
   endfor
  endif
 endfor
 outfile=file_dirname(planfile)+'/apQAdarkflat-'+string(planstr.mjd,format='(i5.5)')+'.fits'
 mwrfits,str,outfile,/create
endif

;mkhtml,mjd

;mkhtmlsum

end
