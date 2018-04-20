;======================================================================
pro mkpersist,persistid,dark,flat,cmjd=cmjd,darkid=darkid,flatid=flatid,sparseid=sparseid,fiberid=fiberid,clobber=clobber,thresh=thresh

  if not keyword_set(thresh) then thresh=0.1

  dirs=getdir()
  caldir=dirs.caldir()

  file=dirs.prefix+string(format='("Persist-c-",i8.8)',persistid)
  ;if another process is alreadying making this file, wait!
  while file_test(caldir+'persist/'+file+'.lock') do apwait,file,10
  ; does product already exist?
  if file_test(caldir+'/persist/'+file+'.fits') and not keyword_set(clobber) then begin
    print,' persist file: ',file+'.fits',' already made'
    return
  endif
  ; open .lock file
  openw,lock,/get_lun,caldir+'persist/'+file+'.lock'
  free_lun,lock

  if keyword_set(cmjd) then begin
    d=approcess([dark,flat],cmjd=cmjd,darkid=darkid,flatid=flatid,psfid=psfid,nfs=1,/doap3dproc) 
  endif else begin
    d=approcess([dark,flat],darkid=darkid,flatid=flatid,psfid=psfid,nfs=1,/doap3dproc)
  endelse

  d=apread('2D',num=dark)
  f=apread('2D',num=flat)

  ; write out an integer mask
  chip=['a','b','c']
  for ichip=0,2 do begin
    persist=intarr(2048,2048)
    r=d[ichip].flux/f[ichip].flux
    bad=where(d[ichip].mask and badmask() or f[ichip].mask and badmask(),nbad)
    if nbad gt 0 then r[bad]=0.
    rz=zap(r,[10,10])
 ;   atv,rz,min=0,max=0.1,/linear
    print,median(rz)
    bad=where(rz gt thresh/4.,nbad)
    if nbad gt 0 then persist[bad]=4
    bad=where(rz gt thresh/2.,nbad)
    if nbad gt 0 then persist[bad]=2
    bad=where(rz gt thresh,nbad)
    if nbad gt 0 then persist[bad]=1
    file=caldir+'/persist/'+dirs.prefix+'Persist-'+chip[ichip]+'-'+string(format='(i8.8)',persistid)+'.fits'
    mwrfits,persist,file,/create
    mwrfits,rz,file
  endfor

  file=dirs.prefix+string(format='("Persist-c-",i8.8)',persistid)
  file_delete,caldir+'/persist/'+file+'.lock'
end
