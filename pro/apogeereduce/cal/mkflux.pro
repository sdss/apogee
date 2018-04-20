;======================================================================

pro mkflux,ims,cmjd=cmjd,darkid=darkid,flatid=flatid,psfid=psfid,waveid=waveid,littrowid=littrowid,persistid=persistid,$
    clobber=clobber,onedclobber=onedclobber,bbtemp=bbtemp,plate=plate,plugid=plugid,holtz=holtz,temp=temp

  dirs=getdir(apodir,caldir,spectrodir,vers)
  caldir=dirs.caldir

  file=dirs.prefix+string(format='("Flux-c-",i8.8)',ims[0])
  ;if another process is alreadying making this file, wait!
  while file_test(caldir+'flux/'+file+'.lock') do apwait,file,10
  ; does product already exist?
  if file_test(caldir+'/flux/'+file+'.fits') and not keyword_set(clobber) then begin
    print,' flux file: ',file+'.fits',' already made'
    if n_elements(temp) ne 0 then goto,response
    return
  endif
  ; open .lock file
  openw,lock,/get_lun,caldir+'flux/'+file+'.lock'
  free_lun,lock

  if not keyword_set(plate) then plate=0

  ; need to make sure extraction is done without flux calibration
  i1=ims[0]
  files=file_search(dirs.expdir+getcmjd(i1)+'/'+dirs.prefix+'1D-?-'+string(format='(i8.8)',i1)+'.fits')
  if files[0] ne '' then file_delete,files
  if keyword_set(cmjd) then $
    d=approcess(ims,cmjd=cmjd,darkid=darkid,flatid=flatid,psfid=psfid,littrowid=littrowid,persistid=persistid,/nocr,nfs=1,/doproc) $
  else $
    d=approcess(ims,darkid=darkid,flatid=flatid,psfid=psfid,littrowid=littrowid,persistid=persistid,/nocr,nfs=1,/doproc)
  cmjd=getcmjd(i1)
  inpfile = dirs.expdir+cmjd+'/'+string(format='(i8.8)',i1)
  fluxdir = caldir+'/flux/'
  APMKFLUXCAL,inpfile,outdir=fluxdir,/clobber

  ; clean up in case someone might want to reduce these files with flux calibration
  files=file_search(dirs.expdir+getcmjd(i1)+'/'+dirs.prefix+'1D-?-'+string(format='(i8.8)',i1)+'.fits')
  if files[0] ne '' then file_delete,files

if keyword_set(holtz) then begin
;  below is Holtz flux calibration method

  nframes=n_elements(ims)
  for ii=0,nframes-1 do begin
   i=ims[ii]
   ;if keyword_set(cmjd) then frame=apread(i,err,mask,head,cmjd=cmjd,/oned) $
   ;else frame=apread(i,err,mask,head,/oned) 
   frame=apread('1D',num=i)
   if ii eq 0 then begin
     head0=frame[0].hdr
     sz=size(frame[0].flux)
     flux=fltarr(sz[1],sz[2],3)
   endif
   for ichip=0,2 do flux[*,*,ichip]+=frame[ichip].flux
  endfor
  
  bad=-1
  if keyword_set(plate) then begin
    if not keyword_set(cmjd) then cmjd=getcmjd(ims[0])
    fiber=getfiber(plate,cmjd,plugid=plugid)
    bad=where(fiber.fiberid lt 0)
  endif

  sz=size(flux)
  chips=['a','b','c']
  resp=fltarr(2048,300,3)
  for ichip=0,2 do begin
    if keyword_set(bbtemp) then begin
      file=dirs.prefix+string(format='("Wave-",a,"-",i8.8)',chips[ichip],waveid)
      wavetab=mrdfits(caldir+'/wave/'+file+'.fits',1)
      refspec=fltarr(sz[1],sz[2])
      for ifiber=0,sz[1]-1 do refspec[ifiber,*]=planck(wavetab[ifiber,*],bbtemp)
    endif else begin
      refflux=reform(flux[*,150,ichip],sz[1],1)
      ;refspec=zap(refflux,[100,1])
      refspec=refflux/refflux
    endelse
    rows=intarr(sz[2])+1
    refimg=rows##refspec
    tmp=zap(flux[*,*,ichip],[100,1])
    if ichip eq 1 then norm=tmp[1024,150]
    resp(*,*,ichip)=refimg/tmp
    if (bad[0] ge 0) then for i=0,n_elements(bad)-1 do resp[*,bad[i]]=0.
  endfor
  ; normalize to center of green chip
  for ichip=0,2 do begin
    resp(*,*,ichip)*=norm
    file=dirs.prefix+string(format='("Flux-",a,"-",i8.8)',chips[ichip],i1)
    mwrfits,resp[*,*,ichip],caldir+'flux/'+file+'.fits',head0,/create
  endfor
endif

  file=dirs.prefix+string(format='("Flux-c-",i8.8)',ims[0])
  file_delete,caldir+'/flux/'+file+'.lock'

  response:
  ; extra block if we are calculating response function 
  if n_elements(temp) gt 0 then begin
    file=dirs.prefix+string(format='("Response-c-",i8.8)',ims[0])
    ;if another process is alreadying making this file, wait!
    while file_test(caldir+'flux/'+file+'.lock') do apwait,file,10
    ; does product already exist?
    if file_test(caldir+'/flux/'+file+'.fits') and not keyword_set(clobber) then begin
      print,' flux file: ',file+'.fits',' already made'
      return
    endif
    ; open .lock file
    openw,lock,/get_lun,caldir+'flux/'+file+'.lock'
    free_lun,lock
    chips=['a','b','c']
    wave=mrdfits(caldir+'/wave/'+dirs.prefix+'Wave-'+chips[1]+'-'+string(format='(i8.8)',waveid)+'.fits',2)
    flux=mrdfits(caldir+'/flux/'+dirs.prefix+'Flux-'+chips[1]+'-'+string(format='(i8.8)',ims[0])+'.fits',3)
    bbnorm = flux[1024] / PLANCK(wave[1024,150],temp)
    for i=0,2 do begin
      wave=mrdfits(caldir+'/wave/'+dirs.prefix+'Wave-'+chips[i]+'-'+string(format='(i8.8)',waveid)+'.fits',2)
      flux=mrdfits(caldir+'/flux/'+dirs.prefix+'Flux-'+chips[i]+'-'+string(format='(i8.8)',ims[0])+'.fits',3)
      bbflux = PLANCK( wave[*,150], temp) * bbnorm
      mkhdr,head,bbflux/flux
      leadstr = 'APMKFLAT: '
      sxaddhist,leadstr+systime(0),head
      info = GET_LOGIN_INFO()
      sxaddhist,leadstr+info.user_name+' on '+info.machine_name,head
      sxaddhist,leadstr+'IDL '+!version.release+' '+!version.os+' '+!version.arch,head
      sxaddhist,leadstr+' APOGEE Reduction Pipeline Version: '+getvers(),head
      file=dirs.prefix+string(format='("Response-",a,"-",i8.8)',chips[i],ims[0])
      mwrfits,bbflux/flux,caldir+'flux/'+file+'.fits',head
    endfor

    file=dirs.prefix+string(format='("Response-c-",i8.8)',ims[0])
    file_delete,caldir+'/flux/'+file+'.lock'
  endif
 
end
