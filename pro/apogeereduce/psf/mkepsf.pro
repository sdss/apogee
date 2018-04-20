;======================================================================
pro mkepsf,ims,cmjd=cmjd,darkid=darkid,flatid=flatid,sparseid=sparseid,clobber=clobber,dmax=dmax,sdmax=sdmax,darkims=darkims,outid=outid,average=average,maxread=maxread,filter=filter,thresh=thresh,scat=scat

  if not keyword_set(dmax) then dmax=7
  if not keyword_set(outid) then outid=ims[0]
  if not keyword_set(sparseid) then sparseid=0

  dirs=getdir()
  caldir=dirs.caldir

  file=dirs.prefix+string(format='("EPSF-c-",i8.8)',outid)
  ;if another process is alreadying make this file, wait!
  while file_test(caldir+'psf/'+file+'.lock') do apwait,file,10
  ; does product already exist?
  if not file_test(caldir+'/psf/'+file+'.fits') or keyword_set(clobber) then begin
    ; open .lock file
    openw,lock,/get_lun,caldir+'psf/'+file+'.lock'
    free_lun,lock

    if keyword_set(cmjd) then $
      d=approcess(ims,cmjd=cmjd,darkid=darkid,flatid=flatid,nfs=1,/nocr,/doap3dproc,maxread=maxread) $
    else $
      d=approcess(ims,darkid=darkid,flatid=flatid,nfs=1,/nocr,/doap3dproc,maxread=maxread)

    for n=0,n_elements(ims)-1 do begin
      ;if keyword_set(cmjd) then frame=apread(ims[n],err,mask,head,cmjd=cmjd,/domask) $
      ;else frame=apread(ims[n],err,mask,head,/domask)
      frame=apread('2D',num=ims[n])
      for ichip=0,2 do begin
        bad=where(frame[ichip].mask and badmask())
        frame[ichip].flux[bad]=!values.f_nan
        frame[ichip].err[bad]=!values.f_nan
      endfor

      if n eq 0 then red=frame.flux else red+=frame.flux
    endfor
    red/=n_elements(ims)

    if keyword_set(darkims) then begin
      d=approcess(darkims,darkid=darkid,/nocr,maxread=maxread,/doap3dproc)
      for n=0,n_elements(darkims)-1 do begin
        ;frame=apread(darkims[n],err,mask,head,/domask)
        frame=apread('2D',num=darkims[n],/domask)
        if n eq 0 then darksum=frame.flux else darksum+=frame.flux
      endfor
      darksum/=n_elements(darkims)
      red-=darksum
    endif
  
    tmp={id: outid,  ntrace: 0, bot: 0., mid: 0., top: 0.}
    tracelog=replicate(tmp,3)
  
    chips=['a','b','c']
    ntrace=intarr(3)
    for ichip=0,2 do begin
      ;apmkpsf_epsf,red[*,*,ichip],caldir,outid,ichip,dmax=dmax,sparseid=sparseid,scat=2,average=average

      if keyword_set(filter) then begin
        ; for sparse, we really don't want to have bad pixels, so try to replace them here, even with a
        ; broader filter than we want to use in EPSF construction
        tmp=red[*,*,ichip]
        filts=[10,50,100,200,300]
        for ifilt=0,n_elements(filts)-1 do begin
          bd = where(finite(tmp) eq 0,nbd)
          print,filts[ifilt],nbd
          if nbd gt 0 then begin
            zap=zap(tmp,[filts[ifilt],1])
            tmp[bd]=zap[bd]
          endif
        endfor
        red[*,*,ichip] = tmp
      endif

      apmkpsf_epsf,red[*,*,ichip],caldir,outid,ichip,dmax=dmax,sparseid=sparseid,average=average,thresh=thresh,scat=scat
    endfor
    sparse=dirs.prefix+string(format='("Sparse-",i8.8)',outid)
    mwrfits,red,caldir+'psf/'+sparse+'.fits',/create
  
    file_delete,caldir+'psf/'+file+'.lock'
  endif
  
end
