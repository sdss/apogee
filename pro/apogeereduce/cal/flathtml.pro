;======================================================================
pro flathtml,caldir,plots=plots

flats=file_search(caldir+'/flatcorr/*.tab')
if not file_test(caldir+'/flatcorr/html') then file_mkdir,caldir+'/flatcorr/html'
openw,lun,/get_lun,caldir+'/flatcorr/html/flats.html'
printf,lun,'<HTML><BODY><TABLE BORDER=1>'
printf,lun,'<TR><TD>ID'
printf,lun,'<TD>NFRAMES'
printf,lun,'<TD>A<TD> B<TD> C'
chips=['a','b','c']
for i=0,n_elements(flats)-1 do begin
  flatlog=mrdfits(flats[i],1)
  for ichip=0,2 do begin
    if keyword_set(plots) then begin
      file=caldir+'/flatcorr/'+flats[i].name+'.fits'
      flat=mrdfits(file,1)
      flatplot,a,file
    endif
    if ichip eq 0 then begin
      printf,lun,'<TR><TD>',flatlog[ichip].num 
      printf,lun,'<TD><center>',flatlog[ichip].nframes
    endif
    file=string(format='("apFlat-",a,"-",i8.8)',chips[ichip],flatlog[ichip].num) 
    printf,lun,'<TD><center><a href=../plots/'+file+'.jpg><img src=../plots/'+file+'.jpg width=100></a>'
  endfor
endfor
printf,lun,'</table></body></html>'
free_lun,lun

end
