pro mkhtmlplate,plate,mjd,fluxid=fluxid
; write master apQA web page

print,'mkhtmlplate for: ',plate,mjd

dirs=getdir(apodir,caldir,spectrodir,vers)
cmjd=string(format='(i5.5)',mjd)
cplate=strtrim(string(format='(i6.4)',plate),2)
platesum=apogee_filename('PlateSum',plate=plate,mjd=mjd)
platedir=file_dirname(platesum)+'/'

fieldname=apogee_field(0,plate)
print,'platedir: ',platedir

print,platedir+'apPlateSum-'+cplate+'-'+cmjd+'.fits'
if file_test(platesum) eq 0 then return

qafile=apogee_filename('QA',plate=plate,mjd=mjd)
print,'opening ',platedir+'/html/apQA'
if not file_test(file_dirname(qafile),/dir) then file_mkdir,file_dirname(qafile)
openw,html,/get_lun,qafile
printf,html,'<HTML><BODY>'
printf,html,'<H2> PLATE:',plate,' MJD: ', mjd, ' FIELD: ', fieldname,'</H2>'

tab1=mrdfits(platesum,1)
tab2=mrdfits(platesum,2)
tab3=mrdfits(platesum,3,status=status)
if status lt 0 then begin
  printf,html,'<FONT COLOR=red> ERROR READING TAB FILE'
  printf,html,'</BODY></HTML>'
  free_lun,html
  return
endif
platefile=apogee_filename('Plate',plate=plate,mjd=mjd,chip='a')
shiftstr=mrdfits(platefile,13)
pairstr=mrdfits(platefile,14,status=status)
if status lt 0 then begin
  printf,html,'<FONT COLOR=red> ERROR READING apPlate FILE'
  printf,html,'</BODY></HTML>'
  free_lun,html
  return
endif

; table of individual exposures
printf,html,'<TABLE BORDER=2>'
printf,html,'<TR bgcolor=lightgreen>'
printf,html,'<TD>Frame<TD>Cart<TD>sec z<TD>HA<TD>DESIGN HA<TD>seeing<TD>FWHM<TD>GDRMS<TD>Nreads<TD>Dither<TD>Pixshift<TD>Zero<TD>Zero rms<TD>sky continuum<TD>S/N<TD>S/N(cframe)'
for i=0,n_elements(tab1)-1 do begin
  printf,html,'<TR>'
  printf,html,'<TD>',tab1[i].im
  printf,html,'<TD>'+string(tab1[i].cart)
  printf,html,'<TD>',string(format='(f8.2)',tab1[i].secz)
  printf,html,'<TD>',string(format='(f8.2)',tab1[i].ha)
  printf,html,'<TD>',string(format='(f6.0,",",f6.0,",",f6.0)',tab1[i].design_ha)
  printf,html,'<TD>',string(format='(f8.2)',tab1[i].seeing)
  printf,html,'<TD>',string(format='(f8.2)',tab1[i].fwhm)
  printf,html,'<TD>',string(format='(f8.2)',tab1[i].gdrms)
  printf,html,'<TD>'+string(tab1[i].nreads)
  j=where(shiftstr.framenum eq tab1[i].im,nj)
  if nj gt 0 then begin
    printf,html,'<TD>'+string(format='(f5.2)',shiftstr[j].shift)
    printf,html,'<TD>'+string(format='(f5.2)',shiftstr[j].pixshift)
  endif else printf,html,'<TD><TD>'
  printf,html,'<TD>'+string(format='(f8.2)',tab1[i].zero)
  printf,html,'<TD>'+string(format='(f8.2)',tab1[i].zerorms)
  printf,html,'<TD>',string(format='("[",f5.2,",",f5.2,",",f5.2,"]")',tab1[i].sky)
  printf,html,'<TD>',string(format='("[",f5.1,",",f5.1,",",f5.1,"]")',tab1[i].sn)
  printf,html,'<TD>',string(format='("[",f5.1,",",f5.1,",",f5.1,"]")',tab1[i].snc)
endfor
printf,html,'</TABLE>'

; table of exposure pairs
npairs=n_elements(pairstr)
if size(pairstr,/type) eq 8 and npairs gt 0 then begin
  ; Pair table
  printf,html,'<BR><TABLE BORDER=2>'
  printf,html,'<TR bgcolor=lightgreen><TD>IPAIR<TD>NAME<TD>SHIFT<TD>NEWSHIFT<TD>S/N'
  printf,html,'<TD>NAME<TD>SHIFT<TD>NEWSHIFT<TD>S/N'
  for ipair=0,npairs-1 do begin
    printf,html,'<TR><TD>',ipair
    for j=0,1 do begin
      printf,html,'<TD>',pairstr[ipair].framename[j]
      printf,html,'<TD>',string(format='(f8.3)',pairstr[ipair].oldshift[j])
      printf,html,'<TD>',string(format='(f8.3)',pairstr[ipair].shift[j])
      printf,html,'<TD>',string(format='(f8.2)',pairstr[ipair].sn[j])
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

; link to combined spectra page
printf,html,'<P><A HREF='+dirs.prefix+'Plate-'+cplate+'-'+cmjd+'.html> Visit (multiple exposure and dither combined) spectra </a><p>'

; flat field plots
if keyword_set(fluxid) then begin
  file=file_basename(apogee_filename('Flux',num=fluxid,chip=['a','b','c'],/base),'.fits')
  printf,html,'<P> Flat field relative fluxes'
  printf,html,'<TABLE BORDER=2><TR>'
  for ichip=0,2 do printf,html,'<TD> <A HREF='+'../plots/'+file[ichip]+'.jpg><IMG SRC=../plots/'+file[ichip]+'.jpg WIDTH=400></A>'
  blockfile=file[0].replace('-a-','-block-')
  printf,html,'<TD> <A HREF='+'../plots/'+blockfile+'.jpg><IMG SRC=../plots/'+blockfile+'.jpg WIDTH=400></A>'
  printf,html,'</TABLE>'
endif
gfile='guider-'+cplate+'-'+cmjd+'.jpg'
printf,html,'<A HREF='+'../plots/'+gfile+'><IMG SRC=../plots/'+gfile+' WIDTH=400></A>'

; table of exposure plots
printf,html,'<TABLE BORDER=2>'

printf,html,'<TR><TD>Frame<TD>Zeropoints<TD>Mag plots'
printf,html,'<TD>Spatial mag deviation'
printf,html,'<TD>Spatial sky telluric CH4 '
printf,html,'<TD>Spatial sky telluric CO2 '
printf,html,'<TD>Spatial sky telluric H2O '
printf,html,'<TD>Spatial sky 16325A emission deviations (filled: sky, open: star)'
printf,html,'<TD>Spatial sky continuum emission '

for i=0,n_elements(tab1)-1 do begin
  ims=tab1[i].im
  file = file_basename(apogee_filename('1D',num=ims,/nochip),'.fits')
  printf,html,'<TR><TD><A HREF=../html/'+file+'.html>',ims,'</A>'
  printf,html,'<TD><TABLE BORDER=1><TD><TD>Red<TD>Green<TD>Blue'
  printf,html,'<TR><TD>z<TD><TD>'+string(format='(f5.2)',tab1[i].zero)
  printf,html,'<TR><TD>znorm<TD><TD>'+string(format='(f5.2)',tab1[i].zeronorm)
  printf,html,'<TR><TD>sky'+$
             string(format='("<TD>",f5.1,"<TD>",f5.1,"<TD>",f5.1)',tab1[i].sky)
  printf,html,'<TR><TD>S/N'+$
             string(format='("<TD>",f5.1,"<TD>",f5.1,"<TD>",f5.1)',tab1[i].sn)
  printf,html,'<TR><TD>S/N(c)'+$
             string(format='("<TD>",f5.1,"<TD>",f5.1,"<TD>",f5.1)',tab1[i].snc)
  if tag_exist(tab1[i],'snratio') then $
  printf,html,'<TR><TD>SN(E/C)<TD<TD>'+string(format='(f5.2)',tab1[i].snratio) $
  else printf,html
  printf,html,'</TABLE>'
  printf,html,'<TD><IMG SRC=../plots/'+file+'.gif>'
  printf,html,'<TD> <IMG SRC=../plots/'+file+'.jpg>'
  cim=string(format='(i8.8)',ims)
  printf,html,'<TD> <a href=../plots/'+dirs.prefix+'telluric_'+cim+'_skyfit_CH4.jpg> <IMG SRC=../plots/'+dirs.prefix+'telluric_'+cim+'_skyfit_CH4.jpg height=400></a>'
  printf,html,'<TD> <a href=../plots/'+dirs.prefix+'telluric_'+cim+'_skyfit_CO2.jpg> <IMG SRC=../plots/'+dirs.prefix+'telluric_'+cim+'_skyfit_CO2.jpg height=400></a>'
  printf,html,'<TD> <a href=../plots/'+dirs.prefix+'telluric_'+cim+'_skyfit_H2O.jpg> <IMG SRC=../plots/'+dirs.prefix+'telluric_'+cim+'_skyfit_H2O.jpg height=400></a>'
  printf,html,'<TD> <IMG SRC=../plots/'+file+'sky.jpg>'
  printf,html,'<TD> <IMG SRC=../plots/'+file+'skycont.jpg>'
endfor
printf,html,'</table>'

printf,html,'</BODY></HTML>'
free_lun,html
end
