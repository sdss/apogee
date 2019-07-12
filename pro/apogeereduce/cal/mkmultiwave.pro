;======================================================================
pro mkmultiwave,waveid,name=name,clobber=clobber,nowait=nowait,file=calfile

  if ~keyword_set(name) then name=string(waveid[0])
  dirs=getdir(apodir,caldir,spectrodir,vers)
  caldir=dirs.caldir
  file=dirs.prefix+string(format='("Wave-",i8.8)',name)
  ;if another process is alreadying make this file, wait!
  while file_test(caldir+'wave/'+file+'.lock') do begin
    if keyword_set(nowait) then begin
      print,' Wavecal file: ', file, ' already being made (.lock file exists)'
      return
    endif
    apwait,file,10
  endwhile
  ; does product already exist?
  if file_test(caldir+'/wave/'+file+'py.dat') and not keyword_set(clobber) then begin
    print,' Wavecal file: ', file+'py.dat', ' already made'
    return
  endif

  print,'making wave: ', waveid
  ; open .lock file
  openw,lock,/get_lun,caldir+'wave/'+file+'.lock'
  free_lun,lock

  ; process the frames and find lines
  for i=0,n_elements(waveid)-1,2 do makecal,wave=waveid[i],file=getenv('APOGEE_DIR')+'/data/cal/'+dirs.instrument+'-wave.par',/nofit

  ; new Python version!
  cmd=['apmultiwavecal','--name',name,'--vers',dirs.apred,'--plot','--hard','--inst',dirs.instrument] ;,'--verbose']
  for i=0,n_elements(waveid)-1 do cmd=[cmd,string(waveid[i])]
  spawn,cmd,/noshell
  if file_test(apogee_filename('Wave',num=waveid[0],chip='a')) then begin
    openw,1,caldir+'/wave/'+file+'py.dat'
    close,1
  endif else stop,'HALT:  failed to make wavecal',waveid
;
;  cmjd=getcmjd(waveid[0])
;  wavefile = dirs.expdir+cmjd+'/'+string(format='(i8.8)',waveid)
;  apmultiwavecal,wavefile,name=name,/save,/multi,clobber=clobber ;,/pl

  file_delete,caldir+'wave/'+file+'.lock'

end
