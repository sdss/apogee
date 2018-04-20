;======================================================================
pro darkhtml,caldir

darks=file_search(caldir+'darkcorr/*.tab')
if not file_test(caldir+'/darkcorr/html',/dir) then file_mkdir,caldir+'/darkcorr/html'
openw,lun,/get_lun,caldir+'darkcorr/html/darks.html'
printf,lun,'<HTML><BODY><TABLE BORDER=1>'
printf,lun,'<TR><TD>ID'
printf,lun,'<TD>CHIP'
printf,lun,'<TD>NREADS'
printf,lun,'<TD>NFRAMES'
printf,lun,'<TD>MEDIAN RATE'
printf,lun,'<TD>NSAT'
printf,lun,'<TD>NHOT'
printf,lun,'<TD>NHOTNEIGH'
printf,lun,'<TD>NBAD'
printf,lun,'<TD>NNEG'
chips=['a','b','c']
for i=0,n_elements(darks)-1 do begin
  darklog=mrdfits(darks[i],1)
  for ichip=0,2 do begin
    if ichip eq 0 then printf,lun,'<TR><TD>',darklog[ichip].num else printf,lun,'<TR><TD>'
    printf,lun,'<TD><center>',chips[ichip]
    printf,lun,'<TD><center>',darklog[ichip].nreads
    printf,lun,'<TD><center>',darklog[ichip].nframes
    printf,lun,'<TD><center>',darklog[ichip].medrate
    printf,lun,'<TD><center>',darklog[ichip].nsat
    printf,lun,'<TD><center>',darklog[ichip].nhot
    printf,lun,'<TD><center>',darklog[ichip].nhotneigh
    printf,lun,'<TD><center>',darklog[ichip].nbad
    printf,lun,'<TD><center>',darklog[ichip].nneg
    file=string(format='("apDark-",a,"-",i8.8)',chips[ichip],darklog[ichip].num) 
    printf,lun,'<TD><center><a href=../plots/'+file+'.jpg>Image</a>'
    printf,lun,'<TD><center><a href=../plots/'+file+'.gif>Plots</a>'
  endfor
endfor
printf,lun,'</table></body></html>'
free_lun,lun

end
