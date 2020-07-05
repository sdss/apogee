;======================================================================
pro mklsf,lsfid,waveid,darkid=darkid,flatid=flatid,psfid=psfid,fiberid=fiberid,clobber=clobber,full=full,newwave=newwave,$
          pl=pl,fibers=fibers,nowait=nowait

  if not keyword_set(newwave) then newwave=0

  dirs=getdir(apodir,caldir,spectrodir,vers)
  caldir=dirs.caldir
  file=apogee_filename('LSF',num=lsfid[0],/nochip)
  file=file_dirname(file)+'/'+file_basename(file,'.fits')

  ;if another process is alreadying make this file, wait!
  while file_test(file+'.lock') do begin
    if keyword_set(nowait) then begin
      print,' LSF file: ', file, ' already being made (.lock file exists)'
      return
    endif
    apwait,file,10
  endwhile

  ; does product already exist?
  if file_test(file+'.sav') and not keyword_set(clobber) then begin
    print,' LSF file: ',file+'.sav',' already made'
    return
  endif

  ; open .lock file
  openw,lock,/get_lun,file+'.lock'
  free_lun,lock

  cmjd=getcmjd(psfid)

  lsffile = apogee_filename('1D',num=lsfid[0],chip='c')

  mkpsf,psfid,darkid=darkid,flatid=flatid,fiberid=fiberid,/clobber
  w=approcess(lsfid,dark=darkid,flat=flatid,psf=psfid,flux=0,/doproc,/skywave,/clobber)
  cmd=['apskywavecal','dummy','--frameid',string(lsfid),'--waveid',string(waveid),'--apred',dirs.apred,'--telescope',dirs.telescope]
  spawn,cmd,/noshell

  lsffile = file_dirname(lsffile)+'/'+string(format='(i8.8)',lsfid)
  if size(waveid,/type) eq 7 then wavefile = caldir+'wave/'+waveid else $
    wavefile = caldir+'wave/'+string(format='(i8.8)',waveid)
  psffile = caldir+'/psf/'+string(format='(i8.8)',psfid)
  aplsf,lsffile,wavefile,psf=psffile,/gauss,pl=pl
  if keyword_set(full) then aplsf,lsffile,wavefile,psf=psffile,/clobber,pl=pl,fibers=fibers ;,porder=[1,1,1,1,1,0],/pl

;,porder=2

    file_delete,file+'.lock'
end
