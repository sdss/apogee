
apsetver,vers='t9'

apsetver,telescope='apo25m'
dirs=getdir()
a=mrdfits(dirs.spectrodir+'/apogee-nSci.fits',1)
j=where(a.nreads ge 47 and a.zero gt 18.5 and a.mjd gt 58000)
openw,lun,dirs.spectrodir+'/monitor/apogee-n/plateplots.html',/get_lun
printf,lun,'<HTML><BODY><TABLE BORDER=2>'
for ii=0,n_elements(j)-4,4 do begin
  printf,lun,'<TR>'
  for iii=0,3 do begin
    i=j[ii+iii]
    printf,lun,'<TD>'+apogee_field(0,a[i].plate)+'/'+string(a[i].plate)+'/'+string(a[i].mjd)
  endfor
  printf,lun,'<TR>'
  for iii=0,3 do begin
    i=j[ii+iii]
    printf,lun,'<TD><img src=../../visit/apo25m/'+apogee_field(0,a[i].plate)+'/'+strtrim(a[i].plate,2)+'/'+strtrim(a[i].mjd,2)+'/plots/ap1D-'+string(format='(i8.8)',a[i].im)+'.jpg width=300>'
  endfor
endfor
printf,lun,'</table></BODY></HTML>'
free_lun,lun

apsetver,telescope='lco25m'
a=mrdfits(dirs.spectrodir+'/apogee-sSci.fits',1)
j=where(a.nreads ge 47 and a.zero gt 18.0 and a.mjd gt 58000)
openw,lun,dirs.spectrodir+'/monitor/apogee-s/plateplots.html',/get_lun
printf,lun,'<HTML><BODY><TABLE BORDER=2>'
for ii=0,n_elements(j)-4,4 do begin
  printf,lun,'<TR>'
  for iii=0,3 do begin
    i=j[ii+iii]
    printf,lun,'<TD>'+apogee_field(0,a[i].plate)+'/'+string(a[i].plate)+'/'+string(a[i].mjd)
  endfor
  printf,lun,'<TR>'
  for iii=0,3 do begin
    i=j[ii+iii]
    printf,lun,'<TD><img src=../../visit/lco25m/'+apogee_field(0,a[i].plate)+'/'+strtrim(a[i].plate,2)+'/'+strtrim(a[i].mjd,2)+'/plots/as1D-'+string(format='(i8.8)',a[i].im)+'.jpg width=300>'
  endfor
endfor
printf,lun,'</table></BODY></HTML>'
free_lun,lun
end

