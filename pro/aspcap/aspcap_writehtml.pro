pro aspcap_writehtml,finalstr,path,name,bestclass=bestclass,spectra=spectra,elem=elem,npar=npar

; Routine to make final HTML pages with summary ASPCAP results

;file_mkdir,path+'/html'
if keyword_set(spectra) then file=path+'/'+name+'_spec.html' else $
 if keyword_set(elem) then file=path+'/'+name+'_elem.html' else $
   file=path+'/'+name+'.html'
if tag_exist(finalstr.param[0],'field') then field=strtrim(finalstr.param[0].field,2) else field='temp'

; open output file
openw,html,file,/get_lun

; link for sortable table
if keyword_set(bestclass) then $
printf,html,'<HTML><HEAD><script type=text/javascript src=../../../html/sorttable.js></script></head>' else $
printf,html,'<HTML><HEAD><script type=text/javascript src=../../../../html/sorttable.js></script></head>' 
printf,html,'<BODY>'
printf,html,'<h1> ASPCAP results: ',name, ' (',field,')'
if keyword_set(bestclass) then printf,html,' Best class' else printf,html,' Class ',finalstr.param[0].class
printf,html,'</h1>'

; links to FITS files
if keyword_set(bestclass) then begin
  printf,html,'<br><a href=aspcapField-'+name+'.fits> aspcapField-'+name+'.fits</A> (FITS tables with results)'
  printf,html,'<br><a href=class_K> class_K files </a>'
  printf,html,'<br><a href=class_G> class_G files </a>'
  printf,html,'<br><a href=class_F> class_F files </a>'
  printf,html,'<br><a href=class_A> class_A files </a>'
  printf,html,'<br>'
endif else begin
  class=finalstr.param[0].class
  file=class+'-'+field
  printf,html,'<br><a href='+file+'.fits> FITS table for class ',class
  printf,html,'<br><a href='+file+'.ipf> ',file+'.ipf </a>'
  printf,html,'<a href='+file+'.frd> ',file+'.frd </a>'
  printf,html,'<a href='+file+'.err> ',file+'.err </a>'
  printf,html,'<a href='+file+'.con> ',file+'.con </a>'
  printf,html,'<a href='+file+'.spm> ',file+'.spm </a>'
  printf,html,'<a href='+file+'.mdl> ',file+'.mdl </a>'
endelse

; summary plots for this field
set_plot,'PS'
cleanplot,/silent
if tag_exist(finalstr.lib,'elem_symbol') then begin
  device,file=path+'/'+field+'.eps',/encap,/port,xsize=32,ysize=3,/inches 
  !p.multi=[0,5+n_elements(finalstr.lib.elem_symbol),1]
endif else begin
  device,file=path+'/'+field+'.eps',/encap,/port,xsize=16,ysize=3,/inches
  !p.multi=[0,5,1]
endelse
smcolor,/ps
bd=where(finalstr.param.aspcapflag,nbd)
plot,[finalstr.param.fparam[0]],[finalstr.param.fparam[1]],psym=6,xtitle='Teff',ytitle='log g',$
   thick=2,charthick=2,xrange=[6500,3500],yrange=[5,0],symsize=0.3
if nbd gt 0 then oplot,[finalstr.param[bd].fparam[0]],[finalstr.param[bd].fparam[1]],$
   psym=6,color=2, thick=2,symsize=0.3
plot,[finalstr.param.fparam[0]],[finalstr.param.fparam[3]],psym=6,xtitle='Teff',ytitle='[Fe/H]',$
  thick=2,charthick=2,xrange=[6500,3500],yrange=[-2,0.5],symsize=0.3
if nbd gt 0 then oplot,[finalstr.param[bd].fparam[0]],[finalstr.param[bd].fparam[3]],$
  psym=6,color=2,thick=2,symsize=0.3
plot,[finalstr.param.fparam[0]],[finalstr.param.fparam[6]],psym=6,xtitle='Teff',ytitle='[alpha/Fe]',$
  thick=2,charthick=2,xrange=[6500,3500],yrange=[-1,1],symsize=0.3
if nbd gt 0 then oplot,[finalstr.param[bd].fparam[0]],[finalstr.param[bd].fparam[6]],$
  psym=6,color=2,thick=2,symsize=0.3
plot,[finalstr.param.fparam[0]],[finalstr.param.fparam[4]],psym=6,xtitle='Teff',ytitle='[C/Fe]',$
  thick=2,charthick=2,xrange=[6500,3500],yrange=[-1,1],symsize=0.3
if nbd gt 0 then oplot,[finalstr.param[bd].fparam[0]],[finalstr.param[bd].fparam[4]],$
  psym=6,color=2,thick=2,symsize=0.3
plot,[finalstr.param.fparam[0]],[finalstr.param.fparam[5]],psym=6,xtitle='Teff',ytitle='[N/Fe]',$
  thick=2,charthick=2,xrange=[6500,3500],yrange=[-1,1],symsize=0.3
if nbd gt 0 then oplot,[finalstr.param[bd].fparam[0]],[finalstr.param[bd].fparam[5]],$
  psym=6,color=2,thick=2,symsize=0.3
if tag_exist(finalstr.lib,'elem_symbol') then $
  for ielem=0,n_elements(finalstr.lib.elem_symbol)-1 do $
  plot,[finalstr.param.fparam[0]],[finalstr.param.elem[ielem]],psym=6,xtitle='Teff',ytitle=finalstr.lib.elem_symbol[ielem],$
    thick=2,charthick=2,xrange=[6500,3500],yrange=[-1,1],symsize=0.3
device,/close
ps2gif,path+'/'+field+'.eps',/eps,/delete

; table listing of all stars with flags, parameters, etc.
printf,html,'<br><IMG SRC='+field+'.gif>'
printf,html,'<br>Click on column labels to sort'
printf,html,'<TABLE BORDER=2 CLASS=sortable>'
if keyword_set(bestclass) then begin
  printf,html,'<TR><TD>Object<TD>S/N<TD>Best class<TD>Chi^2<TD>Teff<TD>log g<TD>vmicro<TD>[Fe/H]<TD>[C/Fe]<TD>[N/Fe]<TD>[alpha/Fe]' 
endif else begin
  printf,html,'<TR><TD>Object<TD>S/N<TD>Class<TD>Chi^2<TD>Teff<TD>log g<TD>vmicro<TD>[Fe/H]<TD>[C/Fe]<TD>[N/Fe]<TD>[alpha/Fe]' 
endelse
params=aspcap_params(npar=npar)
if n_elements(params) eq 8 then printf,html,'<TD>VSINI'
if n_elements(params) eq 9 then printf,html,'<TD>VSINI<TD>PARAM O'
if tag_exist(finalstr.lib,'elem_symbol') then $
  for i=0,n_elements(finalstr.lib.elem_symbol)-1 do printf,html,'<TD>'+finalstr.lib.elem_symbol[i]

for i=0,n_elements(finalstr.param)-1 do begin
  if keyword_set(bestclass) then plotdir='param/class_'+finalstr.param[i].class+'/plots/' else plotdir='plots/'
  if tag_exist(finalstr.param,'file') then $
  objid=finalstr.param[i].file else $
  objid=finalstr.param[i].apogee_id
  if tag_exist(finalstr.param,'starflag') then sflag=finalstr.param[i].starflag else sflag=0
  printf,html,'<TR><TD><a href='+plotdir+strtrim(objid,2)+'.jpg>'+objid,'</a><br>'+starflag(sflag)+'<br>'+aspcapflag(finalstr.param[i].aspcapflag,1)+'<br>'+aspcapflag(finalstr.param[i].aspcapflag,2)+'<br>'+field
  printf,html,'<TD>'+string(format='(f6.1)',finalstr.param[i].snr)
  if finalstr.param[i].aspcapflag gt 0 then $
    printf,html,'<TD><CENTER><FONT color=red>',finalstr.param[i].class,'</FONT></CENTER>'  else $
    printf,html,'<TD><CENTER>',finalstr.param[i].class,'</CENTER>' 
  printf,html,'<TD>',finalstr.param[i].param_chi2,format='(a,f8.1)'
  for ipar=0,n_elements(finalstr.param[i].fparam)-1 do begin
    printf,html,'<TD><TABLE>'
    if ipar eq 0 then formstr='(a,f8.0,a)' else formstr='(a,f8.2,a)'
    if ipar eq 0 then errstr='(a,f7.0)' else err='(a,f7.2)'
    ;if ipar eq 2 then par=10^finalstr.param[i].fparam[ipar] else par=finalstr.param[i].fparam[ipar]
    if params[ipar] eq 'LOG10VDOP' or params[ipar] eq 'LGVSINI' then par=10^finalstr.param[i].fparam[ipar] else par=finalstr.param[i].fparam[ipar]
    printf,html,'<TR><TD>',par,'<TD>+/-',format=formstr
    printf,html,'<TD>',sqrt(finalstr.param[i].fparam_cov[ipar,ipar]),format=errstr
    if ipar eq 2 then par=10^finalstr.param[i].param[ipar] else par=finalstr.param[i].param[ipar]
    printf,html,'<TR><TD>',par,'<TD>+/-',format=formstr
    printf,html,'<TD>',sqrt(finalstr.param[i].param_cov[ipar,ipar]),format=errstr
    printf,html,'</TABLE>'
  endfor
  if tag_exist(finalstr.lib,'elem_symbol') then $
    for ielem=0,n_elements(finalstr.lib.elem_symbol)-1 do begin
     printf,html,'<TD><TABLE>'
     printf,html,'<TR><TD><a href=elem/'+finalstr.param[i].file+'.'+finalstr.lib.elem_symbol[ielem]+'.jpg>'+string(format='(f8.2)',finalstr.param[i].felem[ielem])+'</A>'
     printf,html,'<TR><TD><a href=elem/'+finalstr.param[i].file+'.'+finalstr.lib.elem_symbol[ielem]+'.jpg>'+string(format='(f8.2)',finalstr.param[i].elem[ielem])+'</A>'
     printf,html,'</TABLE>'
    endfor
  if keyword_set(spectra) then printf,html,'<TD><A HREF='+plotdir+strtrim(objid,2)+'.jpg><IMG SRC='+plotdir+strtrim(objid,2)+'_1.jpg></A>'
  if keyword_set(elem) and tag_exist(finalstr.lib,'elem_symbol') then begin
    for iel=0,n_elements(finalstr.lib.elem_symbol)-1 do begin
      el =strtrim(finalstr.lib.elem_symbol[iel],2)
      printf,html,'<TD><A HREF=elem/'+strtrim(objid,2)+'.'+el+'.jpg><IMG SRC=elem/'+strtrim(objid,2)+'.'+el+'.jpg height=400></A>'
    endfor
  endif
endfor
printf,html,'</table>'

printf,html,'</body></html>'

free_lun,html

end
