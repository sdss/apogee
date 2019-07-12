;======================================================================
pro mklittrow,littrowid,cmjd=cmjd,darkid=darkid,flatid=flatid,sparseid=sparseid,fiberid=fiberid,clobber=clobber

  dirs=getdir()
  caldir=dirs.caldir

  file=dirs.prefix+string(format='("Littrow-b-",i8.8)',littrowid)
  ;if another process is alreadying making this file, wait!
  while file_test(caldir+'littrow/'+file+'.lock') do apwait,file,10
  ; does product already exist?
  if file_test(caldir+'/littrow/'+file+'.fits') and not keyword_set(clobber) then begin
    print,' littrow file: ',file+'.fits',' already made'
    return
  endif
  ; open .lock file
  openw,lock,/get_lun,caldir+'littrow/'+file+'.lock'
  free_lun,lock

  ; make empirical PSF with broader smoothing in columns so that Littrow ghost is not incorporated as much
  mkpsf,littrowid,darkid=darkid,flatid=flatid,sparseid=sparseid,fiberid=fiberid,average=200,/clobber
  ; process the frame with this PSF to get model that does not have Littrow ghost
  psfdir=apogee_filename('PSF',num=littrowid,chip='b')
  wavefile=0
  indir=apogee_filename('2D',num=littrowid,/dir,chip='b')
  ap2dproc,indir+'/'+string(format='(i8.8)',littrowid),$
             psfdir+'/'+string(format='(i8.8)',littrowid),4,wavefile=wavefile,/clobber

  ; read in the 2D file and the model, and use them to find the Littrow ghost
  a2=apread('2D',num=littrowid,chip='b')
  a2mod=apread('2Dmodel',num=littrowid,chip='b')
  a=a2.flux
  amask=a2.mask
  scat_remove,a,scat=1
  amod=a2mod.flux
  bad=where((amask and badmask()) gt 0)
  a[bad]=!values.f_nan
  l=where(median(a[1200:1500,*]-amod[1200:1500,*],20) gt 10,complement=nl)

  ; write out an integer mask
  litt=intarr(2048,2048)
  tmp=a[1200:1500,*]*0 & tmp[l]=1 & tmp[nl]=0
  litt[1250:1450,*]=tmp[50:250,*]
  file=apogee_filename('Littrow',num=littrowid,chip='b')

  MKHDR,head,litt   ;,/image
  leadstr = 'APMKLITTROW: '
  sxaddhist,leadstr+systime(0),head
  info = GET_LOGIN_INFO()
  sxaddhist,leadstr+info.user_name+' on '+info.machine_name,head
  sxaddhist,leadstr+'IDL '+!version.release+' '+!version.os+' '+!version.arch,head
  sxaddhist,leadstr+' APOGEE Reduction Pipeline Version: '+getvers(),head
  mwrfits,litt,file,head

  ; move PSFs to littrow directory since they are not a standard PSF!
  outdir=caldir+'littrow/'+cmjd+'/'
  file_mkdir,outdir
  files=file_search(caldir+'psf/*'+string(format='(i8.8)',littrowid)+'*.fits')
  file_move,files,outdir,/over
  files=file_search(file_dirname(apogee_filename('1D',num=littrowid,chip='b'))+'/*1D*'+string(format='(i8.8)',littrowid)+'*.fits')
  file_move,files,outdir,/over
  files=file_search(file_dirname(apogee_filename('2Dmodel',num=littrowid,chip='b'))+'/*2Dmodel*'+string(format='(i8.8)',littrowid)+'*.fits')
  file_move,files,outdir,/over

  file=dirs.prefix+string(format='("Littrow-b-",i8.8)',littrowid)
  file_delete,caldir+'/littrow/'+file+'.lock'
end
