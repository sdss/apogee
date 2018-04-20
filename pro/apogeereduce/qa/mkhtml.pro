; procedure to make summary HTML files for a given night
pro mkhtml,mjd,vers=vers,nocheck=nocheck

nocheck=1

print,'mkhtml for: ',mjd
dirs=getdir(apodir,caldir,spectrodir,vers)

cmjd=string(format='(i5.5)',mjd)

platedir=dirs.spectrodir+'/visit/'+dirs.telescope+'/*/*/'+cmjd+'/'
reddir=dirs.expdir+'/'+cmjd+'/'
outdir=dirs.expdir+'/'+cmjd+'/html/'
datadir=dirs.datadir

if file_test(outdir,/directory) eq 0 then file_mkdir,outdir

; open plate summary page
while file_test(outdir+cmjd+'.html.lock') do apwait,outdir+cmjd+'.html.lock',10

openw,lock,/get_lun,outdir+cmjd+'.html.lock'
free_lun,lock

; get all apR file numbers for the night
print,'looking for apR files and missing frames ....'
files=file_search(datadir+cmjd+'/'+dirs.prefix+'R-*.apz')
nums=lonarr(n_elements(files))
checksums=intarr(n_elements(files))
if strlen(files[0]) eq 0 then begin
  file_delete,outdir+cmjd+'.html.lock',/allow_nonexistent
  return
endif

openw,html,/get_lun,outdir+cmjd+'.html'
printf,html,'<HTML><BODY><CENTER><H2>'+cmjd+'</H2></CENTER>'
for i=0,n_elements(files)-1 do begin
  file=file_basename(files[i])
  if not keyword_set(nocheck) then begin
   if not file_test(reddir+file+'.check') then begin
    checksum=apg_checksum(files[i],fitsdir=getlocaldir())
    openw,clun,/get_lun,reddir+file+'.check'
    printf,clun,checksum
    free_lun,clun
   endif
   readcol,reddir+file+'.check',format='(i)',check
   checksums[i]=check
  endif
  comp=strsplit(file,'-',/extract)
  name=strsplit(comp[2],'.',/extract)
  chip=comp[1]
  num=0L
  reads,name[0],num
  nums[i]=num
endfor
first=min(nums)
last=max(nums)
sortnums=sort(nums)
uniqnums=uniq(nums,sortnums)

if dirs.instrument eq 'apogee-s' then begin
  spawn,"grep -Ei '<TITLE> \[lco-operations ([0-9]+)\] lco/apogee-2s night log [ms]jd ([0-9]+)' "+getenv('SAS_ROOT')+"/lco_staging/reports/*.html | grep "+cmjd+" | awk -F'html:' '{print $1}'",out
  printf,html,' <a href=https://data.sdss.org/sas/sdsswork/lco_staging/reports/'+file_basename(out)+'html> <h3>LCO 2.5m Observing report </h3></a>'
endif else begin
  spawn,"grep -Ei 'subject: 2\.5m obslog ([0-9]+) \([ms]jd ([0-9]+)\)' "+getenv('SAS_ROOT')+"/apo_staging/reports/*.log | grep "+cmjd+" | awk -F'log:' '{print $1}'",out
  printf,html,' <a href=https://data.sdss.org/sas/sdsswork/apo_staging/reports/'+file_basename(out)+'log> <h3>APO 2.5m Observing report </h3></a>'
endelse

; look for missing raw frames (assuming contiguous sequence)
printf,html,'<h3>Raw frames:</h3> ',nums[sortnums[0]], ' to ', nums[sortnums[n_elements(nums)-1]]
printf,html,' (<a href=../../../../../../'+file_basename(dirs.datadir)+'/'+cmjd+'/'+cmjd+'.log.html> image log</a>)'
printf,html,'<BR>'


printf,html,'<h3>Missing raw data:</h3>'
chip=['a','b','c']
nmiss=0
for ichip=0,2 do begin
  printf,html,'<FONT color=red>'
  for j=first,last do begin
    if not file_test(datadir+cmjd+'/'+dirs.prefix+'R-'+chip[ichip]+'-'+string(format='(i8.8)',j)+'.apz') then begin
       printf,html,'',dirs.prefix+'R-'+chip[ichip]+'-'+string(format='(i8.8)',j)+'.apz'
      nmiss+=1
    endif
  endfor
  printf,html,'</font>'
endfor
if nmiss eq 0 then printf,html,'<font color=green> NONE</font>'
printf,html,'<BR>'

if not keyword_set(nocheck) then begin
 print,'looking for bad checksums...'
 bad=where(checksums ne 1, nbad)
 printf,html,'<h3> Bad CHECKSUMS:</h3>'
 if nbad gt 0 then begin
  printf,html,'<font color=red> '
  for i=0,n_elements(bad)-1 do printf,html,file_basename(files[bad[i]])
  printf,html,'</font><BR>'
 endif else printf,html,'<font color=green> NONE </font>'
endif

;  look for missing reduced frames
print,'looking for missing reduced data...'
printf,html,'<h3>Missing reduced data:</h3><BR><TABLE BORDER=2>'
printf,html,'<TR><TD>ID<TD>NFRAMES/NREAD<TD>TYPE<TD>PLATEID<TD>CARTID<TD>1D missing<TD>2D missing'
for i=1,n_elements(uniqnums)-1 do begin
  n=nums[uniqnums[i]]
  file1d=apogee_filename('1D',num=n,chip='c',/base)
  if not file_test(reddir+file1d) then begin 
   file2d=apogee_filename('2D',num=n,chip='c',/base)
   if ~file_test(reddir+file2d) and ~file_test(reddir+file2d+'.fz') then miss2d=1 else miss2d=0
    type='unknown'
    head=[' ',' ']
    rawfile=apogee_filename('R',num=n,chip='a')
    if file_test(rawfile) then begin
     ;a=mrdfits(datadir+'apR-a-'+string(format='(i8.8)',n)+'.apz',1,head,/silent)
     head=headfits(rawfile,exten=1)
     type=strtrim(sxpar(head,'IMAGETYP'),2)
    endif
    color='white'
    if type eq 'Object' then color='red'
    if type eq 'unknown' then color='magenta'
    if type eq 'Dark' and miss2d then color='yellow'
    if type ne 'Dark' or miss2d then begin
      printf,html,'<TR bgcolor='+color+'><TD> ',string(format='(i8.8)',n)
      printf,html,'<TD><CENTER>',sxpar(head,'NFRAMES'),'/',sxpar(head,'NREAD'),'</CENTER>'
      printf,html,'<TD><CENTER>',sxpar(head,'IMAGETYP'),'</CENTER>'
      printf,html,'<TD><CENTER>',sxpar(head,'PLATEID'),'</CENTER>'
      printf,html,'<TD><CENTER>',sxpar(head,'CARTID'),'</CENTER>'
      printf,html,'<TD> ',file1d
      if ~file_test(reddir+file2d) and ~file_test(reddir+file2d+'.fz') then printf,html,'<TD> ',file2d
    endif
  endif
endfor
printf,html,'</TABLE>'

; get all observed plates (from planfiles)
print,'getting observed plates ....'
planfiles=file_search(platedir+'*Plan*.par')
printf,html,'<TABLE BORDER=2>'
printf,html,'<TR><TD>Planfile<TD>Nframes<TD>Median zeropoint<TD>Median RMS zeropoint<TD>Cartridge<TD>Unmapped<TD>Missing'
for i=0,n_elements(planfiles)-1 do begin
 if planfiles[i] ne '' then begin
  printf,html,'<TR><TD>',file_basename(planfiles[i],'.par')
  aploadplan,planfiles[i],planstr
  cplate=strtrim(string(format='(i6.4)',planstr.plateid),2)
  platefile=apogee_filename('PlateSum',plate=cplate,mjd=planstr.mjd,/base)
  ;platefile='apPlateSum-'+cplate+'-'+string(format='(i5.5)',planstr.mjd)+'.fits'
  if planstr.platetype eq 'normal' and file_test(file_dirname(planfiles[i])+'/'+platefile) then begin
    platetab=mrdfits(platedir+platefile,1)
    platefiber=mrdfits(platedir+platefile,2)
    printf,html,'<TD>',n_elements(platetab)
    if n_elements(platetab.zero) gt 1 then $
      printf,html,'<TD>',string(format='(f8.2)',median(platetab.zero)) else $
      printf,html,'<TD>',string(format='(f8.2)',platetab.zero)
    if tag_exist(platetab,'zerorms') then begin
    if n_elements(platetab.zerorms) gt 1 then $
      printf,html,'<TD>',string(format='(f8.2)',median(platetab.zerorms)) else $
      printf,html,'<TD>',string(format='(f8.2)',platetab.zerorms)
    endif else printf,html,'<TD>'
    if tag_exist(platetab,'cart') then $
    printf,html,'<TD>',platetab.cart else printf,html,'<TD>'
    unplugged=where(platefiber.fiberid lt 0)
    printf,html,'<TD>'
    if unplugged[0] ge 0 then printf,html,300-unplugged
    printf,html,'<TD>'
    expfile=apogee_filename('1D',num=planstr.fluxid,chip='b')
    if file_test(expfile) then begin
      domeflat=mrdfits(expfile,1)
      level=median(domeflat,dim=1)
      bad=where(level eq 0)
      if bad[0] ge 0 then printf,html,300-bad
    endif
  endif
 endif
endfor
printf,html,'</TABLE>'

print,'wavehtml...'
file=caldir+'wave/html/wave'+string(format='(i5.5)',mjd)+'.html'
if file_test(file) then begin
  spawn,'cat '+file,wavehtml
  for i=1,n_elements(wavehtml)-2 do printf,html,wavehtml[i]
endif

; get all succesfully reduced plates
print,'getting successfully reduced plates...'
platefiles=file_search(platedir+'*PlateSum*.fits')
; make master plot of zeropoint and sky levels for the night
if n_elements(platefiles) ge 1 and platefiles[0] ne '' then begin
 for i=0,n_elements(platefiles)-1 do begin
  platetab=mrdfits(platefiles[i],1)
  cplate=strtrim(string(format='(i6.4)',platetab[0].plate),2)
  cmjd=string(format='(i5.5)',platetab[0].mjd)
  sntab,tabs=platefiles[i],outfile=file_dirname(platefiles[i])+'/sn-'+cplate+'-'+cmjd+'.dat'
  sntab,tabs=platefiles[i],outfile=file_dirname(platefiles[i])+'/altsn-'+cplate+'-'+cmjd+'.dat',/altsn
  if i eq 0 then begin
    zero=platetab.zero
    ims=platetab.im
    moondist=platetab.moondist
    skyr=reform(platetab.sky[0,*])
    skyg=reform(platetab.sky[1,*])
    skyb=reform(platetab.sky[2,*])
  endif else begin
    zero=[zero,platetab.zero]
    ims=[ims,platetab.im]
    skyr=[skyr,reform(platetab.sky[0,*])]
    skyg=[skyg,reform(platetab.sky[1,*])]
    skyb=[skyb,reform(platetab.sky[2,*])]
    moondist=[moondist,platetab.moondist]
  endelse
 endfor
 printf,html,'<h3>Zeropoints and sky levels: </h3><br>'
 printf,html,'<table border=2><tr><td>Zeropoints<td>Sky level<td> Sky level vs moon distance'

 print,'plots...'
 set_plot,'ps'
 smcolor
 if not file_test(reddir+'/plots',/dir) then file_mkdir,reddir+'/plots'
 device,file=reddir+'/plots/'+cmjd+'zero.eps',/encap,ysize=8,/color
 xmin=min(ims mod 10000)-1 & xmax=max(ims mod 10000)+1
 good=where(zero gt 0)
 ymin=min(zero(good)) & ymax=max(zero)
 if ymin gt 15 then ymin=15
 if ymax lt 20 then ymax=20
 plot,ims mod 10000,zero,psym=6,yrange=[ymin,ymax],xrange=[xmin,xmax],xtitle='Image number',ytitle='Zeropoint per pixel'
 device,/close
 ps2gif,reddir+'/plots/'+cmjd+'zero.eps',chmod='664'o,/delete,/eps
 printf,html,'<tr><td><img src=../plots/'+cmjd+'zero.gif>'

 device,file=reddir+'/plots/'+cmjd+'sky.eps',/encap,ysize=8,/color
 ymin=min(skyr) & ymax=max(skyr)
 if ymin gt 11 then ymin=11
 if ymax lt 16 then ymax=16
 plot,ims mod 10000,skyr,psym=6,yrange=[ymax,ymin],xrange=[xmin,xmax],xtitle='Image number',ytitle='Continuum sky per pixel '
 oplot,ims mod 10000,skyr,psym=6,color=2
 oplot,ims mod 10000,skyg,psym=6,color=3
 oplot,ims mod 10000,skyb,psym=6,color=4
 device,/close
 ps2gif,reddir+'/plots/'+cmjd+'sky.eps',chmod='664'o,/delete,/eps
 printf,html,'<td><img src=../plots/'+cmjd+'sky.gif>'

 device,file=reddir+'/plots/'+cmjd+'moonsky.eps',/encap,ysize=8,/color
 ymin=min(skyr) & ymax=max(skyr)
 if ymin gt 11 then ymin=11
 if ymax lt 16 then ymax=16
 plot,moondist,skyr,psym=6,yrange=[ymax,ymin],xtitle='Moon distance',ytitle='Continuum sky per pixel '
 oplot,moondist,skyr,psym=6,color=2
 oplot,moondist,skyg,psym=6,color=3
 oplot,moondist,skyb,psym=6,color=4
 device,/close
 ps2gif,reddir+'/plots/'+cmjd+'moonsky.eps',chmod='664'o,/delete,/eps
 printf,html,'<td><img src=../plots/'+cmjd+'moonsky.gif>'

 printf,html,'</table>'

 printf,html,'<br>Moon phase: ',platetab[0].moonphase,'<br>',format='(a,f8.2,a)'
endif

;spawn,'tail --lines=+2 '+platedir+'/html/*sum.html | head --lines=-1',result
printf,html,'<p><h3>Observed plates:</h3> '
printf,html,'<TABLE BORDER=2>' 
sumfiles=file_search(platedir+'/html/[0-9]*-?????sum.html')
for i=0,n_elements(sumfiles)-1 do begin
  if sumfiles[i] ne '' then begin
    spawn,'tail --lines=+3 '+sumfiles[i]+' | head --lines=-2',result
    printf,html,result
    printf,html,'<TR>'
  endif
endfor
printf,html,'</table>'


printf,html,'</body></html>'
free_lun,html

; web page with summary images/links
print,'web pages...'
files=file_search(platedir+'/apSum-a*.fits')
for i=0,n_elements(files)-1 do begin
 if files[i] ne '' then begin
  dirname=file_dirname(files[i]) 
  file=strmid(file_basename(files[i],'.fits'),8,10)
  openw,lhtml,/get_lun,dirname+'/html/apSum-'+file+'.html'
  printf,lhtml,'<HTML><BODY><H2>'+file+'</H2>'
  printf,lhtml,'<TABLE BORDER=2>'
  printf,lhtml,'<TR><TD><A HREF=../apSum-a-'+file+'.fits> apSum-a-'+file+'.fits</A>'
  printf,lhtml,'<TD><A HREF=../apSum-b-'+file+'.fits> apSum-b-'+file+'.fits</A>'
  printf,lhtml,'<TD><A HREF=../apSum-c-'+file+'.fits> apSum-c-'+file+'.fits</A>'
  printf,lhtml,'<TR><TD><A HREF=../apSum2D-a-'+file+'.fits> apSum2D-a-'+file+'.fits</A>'
  printf,lhtml,'<TD><A HREF=../apSum2D-b-'+file+'.fits> apSum2D-b-'+file+'.fits</A>'
  printf,lhtml,'<TD><A HREF=../apSum2D-c-'+file+'.fits> apSum2D-c-'+file+'.fits</A>'
  printf,lhtml,'<TR><TD><IMG SRC=../plots/apSum-a-'+file+'a.jpg>'
  printf,lhtml,'<TD><IMG SRC=../plots/apSum-b-'+file+'a.jpg>'
  printf,lhtml,'<TD><IMG SRC=../plots/apSum-c-'+file+'a.jpg>'
  printf,lhtml,'<TR><TD><IMG SRC=../plots/apSum-a-'+file+'b.jpg>'
  printf,lhtml,'<TD><IMG SRC=../plots/apSum-b-'+file+'b.jpg>'
  printf,lhtml,'<TD><IMG SRC=../plots/apSum-c-'+file+'b.jpg>'
  printf,lhtml,'</TABLE></BODY></HTML>'
  free_lun,lhtml
 endif
endfor

; exposure summary page: individual exposures MJD5exp.html web page
openw,html,/get_lun,outdir+cmjd+'exp.html'
printf,html,'<HTML><BODY><CENTER><H2>'+cmjd+'</H2></CENTER>'
printf,html,'<TABLE BORDER=2>'
printf,html,'<TR><TD>ID<TD>NREADS<TD>TYPE<TD>PLATEID<TD>CARTID<TD>2D images<TD>1D images'
for i=0,n_elements(uniqnums)-1 do begin
  num=nums[uniqnums[i]]
  str={mjd: 0L, dateobs: ' ', jd: 0.d0, num: 0L, nframes: 0, imagetyp: ' ', plateid: 0, cartid: 0, ra: 0.d0, dec: 0.d0, seeing: 0., alt: 0., qrtz: 0, thar: 0, une: 0, ffs: ' ', ln2level: 0., dithpix: 0., tracedist: 0., med: fltarr(300,3)}
  file=string(format='(i8.8)',num)
  printf,html,'<TR><TD>',num
  rawfile=apogee_filename('R',num=num,chip='a')
  if file_test(rawfile) then begin
    ;a=mrdfits(datadir+'apR-a-'+file+'.apz',1,head,/silent)
    head=headfits(rawfile,exten=1)
    printf,html,'<TD><CENTER>',sxpar(head,'NFRAMES'),'</CENTER>'
    printf,html,'<TD><CENTER>',sxpar(head,'IMAGETYP'),'</CENTER>'
    printf,html,'<TD><CENTER>',sxpar(head,'PLATEID'),'</CENTER>'
    printf,html,'<TD><CENTER>',sxpar(head,'CARTID'),'</CENTER>'
    str.mjd=mjd
    str.dateobs=sxpar(head,'DATE-OBS')
    str.jd=date_conv(str.dateobs,'J')
    str.num=num
    str.nframes=sxpar(head,'NFRAMES')
    str.imagetyp=sxpar(head,'IMAGETYP')
    str.plateid=sxpar(head,'PLATEID')
    str.cartid=sxpar(head,'CARTID')
    str.ra=sxpar(head,'RA')
    str.dec=sxpar(head,'DEC')
    str.seeing=sxpar(head,'SEEING')
    str.alt=sxpar(head,'ALT')
    str.qrtz=sxpar(head,'LAMPQRTZ')
    str.thar=sxpar(head,'LAMPTHAR')
    str.une=sxpar(head,'LAMPUNE')
    str.ffs=sxpar(head,'FFS')
    str.ln2level=sxpar(head,'LN2LEVEL')
    str.dithpix=sxpar(head,'DITHPIX')
  endif else printf,html,'<TD><TD><TD><TD>'
  ;a1=apread(num,/oned) 
  ;if size(a1,/type) eq 8 then str.med=median(a1.flux,dim=1) else str.med=fltarr(300,3)
  if strpos(str.imagetyp,'DomeFlat') ge 0 then begin
   file=apogee_filename('Flux',num=num,chip='c')
   if file_test(file) then begin
    a1=apread('Flux',num=num)
    if size(a1,/type) eq 8 then begin
      for ifiber=0,299 do a1.flux[*,ifiber,*]*=a1.mask
      str.med=median(a1.flux,dim=1) 
    endif else str.med=fltarr(300,3)
    ; get distance from reference trace
    file=apogee_filename('ETrace',num=num,chip='c')
    if file_test(file) then begin
     a1=apread('ETrace',num=num)
     if size(a1,/type) eq 8 then str.tracedist=sxpar(a1[0].hdr,'AVGDIST')
    endif
   endif
  endif 
  if i eq 0 then allstr=str else allstr=[allstr,str]
  if file_test(reddir+'/plots/ap2D-a-'+file+'a.jpg') then begin
    openw,lhtml,/get_lun,outdir+'ap2D-'+file+'.html'
    printf,lhtml,'<HTML><BODY><H2>'+file+'</H2>'

    printf,lhtml,'<TABLE BORDER=2>'
    printf,lhtml,'<TR><TD><A HREF=../ap2D-a-'+file+'.fits> ap2D-a-'+file+'.fits</A>'
    printf,lhtml,'<TD><A HREF=../ap2D-b-'+file+'.fits> ap2D-b-'+file+'.fits</A>'
    printf,lhtml,'<TD><A HREF=../ap2D-c-'+file+'.fits> ap2D-c-'+file+'.fits</A>'
    printf,lhtml,'<TR><TD><IMG SRC=../plots/ap2D-a-'+file+'a.jpg>'
    printf,lhtml,'<TD><IMG SRC=../plots/ap2D-b-'+file+'a.jpg>'
    printf,lhtml,'<TD><IMG SRC=../plots/ap2D-c-'+file+'a.jpg>'
    printf,lhtml,'<TR><TD><IMG SRC=../plots/ap2D-a-'+file+'b.jpg>'
    printf,lhtml,'<TD><IMG SRC=../plots/ap2D-b-'+file+'b.jpg>'
    printf,lhtml,'<TD><IMG SRC=../plots/ap2D-c-'+file+'b.jpg>'
    printf,lhtml,'</TABLE></BODY></HTML>'
    free_lun,lhtml
    printf,html,'<TD><A HREF=ap2D-'+file+'.html>2D images</a>'
  endif
  if file_test(reddir+'/plots/ap1D-a-'+file+'a.jpg') then begin
    openw,lhtml,/get_lun,outdir+'ap1D-'+file+'.html'
    printf,lhtml,'<HTML><BODY><H2>'+file+'</H2>'
    printf,lhtml,'<TABLE BORDER=2>'
    printf,lhtml,'<TR><TD><A HREF=../ap1D-a-'+file+'.fits> ap1D-a-'+file+'.fits</A>'
    printf,lhtml,'<TD><A HREF=../ap1D-b-'+file+'.fits> ap1D-b-'+file+'.fits</A>'
    printf,lhtml,'<TD><A HREF=../ap1D-c-'+file+'.fits> ap1D-c-'+file+'.fits</A>'
    printf,lhtml,'<TR><TD><IMG SRC=../plots/ap1D-a-'+file+'a.jpg>'
    printf,lhtml,'<TD><IMG SRC=../plots/ap1D-b-'+file+'a.jpg>'
    printf,lhtml,'<TD><IMG SRC=../plots/ap1D-c-'+file+'a.jpg>'
    printf,lhtml,'<TR><TD><IMG SRC=../plots/ap1D-a-'+file+'b.jpg>'
    printf,lhtml,'<TD><IMG SRC=../plots/ap1D-b-'+file+'b.jpg>'
    printf,lhtml,'<TD><IMG SRC=../plots/ap1D-c-'+file+'b.jpg>'
    printf,lhtml,'</TABLE></BODY></HTML>'
    free_lun,lhtml
    printf,html,'<TD><A HREF=ap1D-'+file+'.html>1D images</a>'
  endif
endfor
mwrfits,allstr,reddir+cmjd+'exp.fits',/create

printf,html,'</TABLE></BODY></HTML>'
free_lun,html

file_delete,outdir+cmjd+'.html.lock',/allow_nonexistent

end 

