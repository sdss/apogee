pro apmkallplan,mjdstart,mjdend,vers=vers

pipedir = getenv('APOGEEREDUCEPLAN_DIR')
plandir=pipedir+'/pro/'
dir1m=pipedir+'/data/1m/'

openw,all,'allplan.pro',/get_lun
; vers keyword can be used to override strict versioning! beware
if keyword_set(vers) then begin
  printf,all,'vers='''+vers+'''' 
  apsetver,vers=vers
endif else printf,all,'undefine,vers'
telescope=['apo25m','lco25m','apo1m']
inst=['apogee-n','apogee-s','apogee-n']
prefix=['ap','as','ap']
for itele=0,2 do begin
 printf,all,'apsetver,vers=vers,telescope="'+telescope[itele]+'"'
 for mjd=mjdstart,mjdend do begin
  cmjd=string(mjd,format='(i5.5)')
  if itele eq 2 then begin
    tmp=date_conv(float(mjd)+2400000.5,'F')
    name1m=strmid(tmp,2,2)+strmid(tmp,5,2)+strmid(tmp,8,2)
    if file_test(dir1m+name1m) then begin
      file_mkdir,'exposures/'+inst[itele]+'/'+cmjd+'/plan'
      printf,all,'mkplan1m,'+"'"+name1m+"'"
    endif
  endif else begin
   if file_test(plandir+telescope[itele]+'/'+telescope[itele]+'_'+cmjd+'.pro') then begin
    printf,all,'@'+telescope[itele]+'_'+cmjd
    file_mkdir,'exposures/'+inst[itele]+'/'+cmjd+'/plan'
    openw,plan,'exposures/'+inst[itele]+'/'+cmjd+'/plan/'+prefix[itele]+'MJD-'+cmjd+'.par',/get_lun
    printf,plan,'apred_vers  '+vers
    printf,plan,'telescope  '+telescope[itele]
    printf,plan,'mjd  '+cmjd
    free_lun,plan
   endif
  endelse
 endfor
endfor

free_lun,all

end
