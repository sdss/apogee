;+
;
; MAKECAL
;
; This program contains routines to handle APOGEE calibration data
; It include routines:
;   makecal    will make one or ALL of the specified calibration product types
;              listed in the master calibration index file
;
; the actual construction of the calibration products is done with routines
;   called here, but located in apmkcal.pro
;
; Written by J.Holtzman Aug 2011
;-

;===========================================================================
;
;  makecal will make one or ALL of the specified calibration product types
;          listed in the master calibration index file
;
;   USAGE:  makecal,file=file,/dark,/flat,/wave,/lsf
;           OR
;           makecal,file=file,dark=darkid,flat=flatid,wave=waveid,lsf=lsfid
;
;   INPUT:  file=file: name of master calibration index file, if not specified
;               use default cal.par in calibration directory
;           /dark:  make all of the darks in the file
;           dark=darkid:  make the dark with name=darkid 
;           /flat:  make all of the flats in the file
;           flat=flatid:  make the flat with name=flatid 
;           /wave:  make all of the wavecals in the file
;           wave=waveid:  make the wavecal with name=waveid 
;           /lsf:  make all of the lsfs in the file
;           lsf=lsfid:  make the lsf with name=lsfid 
;
pro makecal,file=file,det=det,dark=dark,flat=flat,wave=wave,multiwave=multiwave,lsf=lsf,bpm=bpm,$
    psf=psf,flux=flux,sparse=sparse,fiber=fiber,$
    littrow=littrow,persist=persist,modelpersist=modelpersist,response=response,mjd=mjd,full=full,$
    newwave=newwave,nskip=nskip,average=average,clobber=clobber,vers=vers,telescope=telescope,nofit=nofit

  if keyword_set(vers) and keyword_set(telescope) then apsetver,vers=vers,telescope=telescope
  dirs=getdir(apo_dir,cal_dir,spectro_dir,apo_vers,lib_dir)
  
; get default file name if file not specified
  if keyword_set(file) then begin
    if  strpos(file,'/') lt 0 then file=file_dirname(dirs.calfile)+'/'+file 
  endif else file=dirs.calfile
  calfile=dirs.calfile

  if not keyword_set(full) then full=0
  if not keyword_set(newwave) then newwave=0
  if not keyword_set(nskip) then nskip=1

  ; read calibration master file into calibration structures
  readcal,file,darkstr,flatstr,sparsestr,fiberstr,badfiberstr,fixfiberstr,wavestr,lsfstr,bpmstr,fluxstr,detstr,littrowstr,persiststr,persistmodelstr,responsestr,multiwavestr
  ; make detector file as called for
  if keyword_set(det) then begin
    print,'makecal det: ', det
    if det gt 1 then begin
      i=where(detstr.name eq det)
      if i lt 0 then begin
        print,'No matching calibration line for ', det
        stop
      endif
      mkdet,detstr[i].name,detstr[i].linid
    endif 
  endif

  ; make darks as called for
  if keyword_set(dark) then begin
    print,'makecal dark: ', dark
    if dark gt 1 then begin
      i=where(darkstr.name eq dark)
      if i lt 0 then begin
        print,'No matching calibration line for ', dark
        stop
      endif
      ims=getnums(darkstr[i].frames)
      cmjd=getcmjd(ims[0],mjd=mjd)
      getcal,mjd,calfile,detid=detid
      makecal,det=detid
      mkdark,ims,clobber=clobber
    endif else begin
      if keyword_set(mjd) then  begin
        num=getnum(mjd) 
        red=where(darkstr.frames/10000L eq num)
      endif else red=indgen(n_elements(darkstr))
      if (red[0] ge 0) then begin
       for i=0,n_elements(red)-1 do begin
        ims=getnums(darkstr[red[i]].frames)
        cmjd=getcmjd(ims[0],mjd=mjd)
        getcal,mjd,calfile,detid=detid
        makecal,det=detid
        mkdark,ims,clobber=clobber
       endfor
      endif
    endelse
  endif

  ; make flats as called for
  if keyword_set(flat) then begin
    print,'makecal flat: ', flat
    if flat gt 1 then begin
      i=where(flatstr.name eq flat)
      if i lt 0 then begin
        print,'No matching calibration line for ', flat
        stop
      endif
      ims=getnums(flatstr[i].frames)
      cmjd=getcmjd(ims[0],mjd=mjd)
      getcal,mjd,calfile,darkid=darkid
      makecal,dark=darkid
      mkflat,ims,darkid=darkid,nrep=flatstr[i].nrep,dithered=flatstr[i].dithered,clobber=clobber
    endif else begin
      if keyword_set(mjd) then  begin
        num=getnum(mjd) 
        red=where(flatstr.frames/10000L eq num)
      endif else red=indgen(n_elements(darkstr))
      if (red[0] ge 0) then begin
       for i=0,n_elements(red)-1 do begin
        ims=getnums(flatstr[red[i]].frames)
        cmjd=getcmjd(ims[0],mjd=mjd)
        getcal,mjd,calfile,darkid=darkid
        makecal,dark=darkid
        mkflat,ims,darkid=darkid,nrep=flatstr[i].nrep,dithered=flatstr[i].dithered,clobber=clobber
       endfor
      endif
    endelse
  endif

  if keyword_set(bpm) then begin
    print,'makecal bpm: ', bpm
    if bpm gt 1 then begin
      i=where(bpmstr.name eq bpm)
      if i lt 0 then begin
        print,'No matching calibration line for ', bpm
        stop
      endif
      makecal,dark=bpmstr[i].darkid
      makecal,flat=bpmstr[i].flatid
      mkbpm,bpmstr[i].name,darkid=bpmstr[i].darkid,flatid=bpmstr[i].flatid,clobber=clobber
    endif else begin
      if keyword_set(mjd) then  begin
        num=getnum(mjd) 
        red=where(bpmstr.frames/10000L eq num)
      endif else red=indgen(n_elements(bpmstr))
      if (red[0] ge 0) then begin
       for i=0,n_elements(red)-1 do begin
        makecal,dark=bpmstr[i].darkid,clobber=clobber
        makecal,flat=bpmstr[i].flatid,clobber=clobber
        mkbpm,bpmstr[red[i]].name,darkid=bpmstr[i].darkid,flatid=bpmstr[i].flatid,clobber=clobber
       endfor
      endif
    endelse
  endif

  if keyword_set(sparse) then begin
    print,'makecal sparse: ', sparse
    if sparse gt 1 then begin
      i=where(sparsestr.name eq sparse)
      if i lt 0 then begin
        print,'No matching calibration line for ', sparse
        stop
      endif
      ims=getnums(sparsestr[i].frames)
      cmjd=getcmjd(ims[0],mjd=mjd)
      getcal,mjd,calfile,darkid=darkid,flatid=flatid,bpmid=bpmid
      makecal,dark=darkid
      makecal,flat=flatid
      makecal,bpm=bpmid
      darkims=getnums(sparsestr[i].darkframes)
      maxread=getnums(sparsestr[i].maxread)
      if n_elements(maxread) ne 3 then begin
        print,'sparse maxread does not have 3 elements! '
        stop
      endif
      mkepsf,ims,darkid=darkid,flatid=flatid,darkims=darkims,dmax=sparsestr[i].dmax,maxread=maxread,clobber=clobber,/filter,thresh=0.2,scat=2
    endif
  endif

  if keyword_set(fiber) then begin
    print,'makecal fiber: ', fiber
    if fiber gt 1 then begin
      cmjd=getcmjd(fiber,mjd=mjd)
      getcal,mjd,calfile,darkid=darkid,flatid=flatid,sparseid=sparseid
      mkpsf,fiber,darkid=darkid,flatid=flatid,sparseid=sparseid
    endif
  endif

  if keyword_set(psf) and ~keyword_set(flux) then begin
    print,'makecal psf: ', psf
    if psf gt 1 then begin
      cmjd=getcmjd(psf,mjd=mjd)
      getcal,mjd,calfile,darkid=darkid,flatid=flatid,sparseid=sparseid,fiberid=fiberid,littrowid=littrowid
      makecal,littrow=littrowid
      mkpsf,psf,darkid=darkid,flatid=flatid,sparseid=sparseid,fiberid=fiberid,littrowid=littrowid,clobber=clobber
    endif
  endif

  if keyword_set(littrow) then begin
    print,'makecal littrow: ', littrow
    if littrow gt 1 then begin
      cmjd=getcmjd(littrow,mjd=mjd)
      getcal,mjd,calfile,darkid=darkid,flatid=flatid,sparseid=sparseid,fiberid=fiberid
      makecal,flat=flatid
      mklittrow,littrow,cmjd=cmjd,darkid=darkid,flatid=flatid,sparseid=sparseid,fiberid=fiberid,clobber=clobber
    endif
  endif

  if keyword_set(persist) then begin
    print,'makecal persist: ', persist
    if persist gt 1 then begin
      i=where(persiststr.name eq persist)
      if i lt 0 then begin
        print,'No matching calibration line for ', persist
        stop
      endif
      cmjd=getcmjd(persist,mjd=mjd)
      getcal,mjd,calfile,darkid=darkid,flatid=flatid,sparseid=sparseid,fiberid=fiberid
      mkpersist,persist,persiststr[i].darkid,persiststr[i].flatid,thresh=persiststr[i].thresh,cmjd=cmjd,darkid=darkid,flatid=flatid,sparseid=sparseid,fiberid=fiberid,clobber=clobber
    endif
  endif

  if keyword_set(modelpersist) then begin
    print,'makecal modelpersist: ', modelpersist
    if modelpersist gt 1 then begin
      i=where(persistmodelstr.name eq modelpersist)
      if i lt 0 then begin
        print,'No matching calibration line for ', modelpersist
        stop
      endif
      cmjd=getcmjd(modelpersist,mjd=mjd)
      getcal,mjd,calfile,darkid=darkid,flatid=flatid,sparseid=sparseid,fiberid=fiberid
      mkpersistmodel,modelpersist
    endif
  endif

  if keyword_set(flux) then begin
    if not keyword_set(psf) then begin
      print,'for flux, need psf as well!'
      stop
    endif
    print,'makecal flux: ', flux
    if flux gt 1 then begin
      cmjd=getcmjd(flux,mjd=mjd)
      makecal,psf=psf
      getcal,mjd,calfile,darkid=darkid,flatid=flatid,littrowid=littrowid,waveid=waveid
      makecal,littrow=littrowid
      mkflux,flux,darkid=darkid,flatid=flatid,psfid=psf,littrowid=littrowid,waveid=waveid,clobber=clobber
    endif
  endif

  if keyword_set(response) then begin
    print,'makecal response: ', response
    if response gt 1 then begin
      i=where(responsestr.name eq response,nres)
      if nres eq 0 then begin
        print,'No matching calibration line for ', response
        stop
      endif else if nres gt 1 then i=i[0]
      file=apogee_filename('Response',chip='c',num=response)
      ; does product already exist?
      if file_test(file) and not keyword_set(clobber) then begin
        print,' flux file: ',file,' already made'
        return
      endif

      cmjd=getcmjd(response,mjd=mjd)
      getcal,mjd,calfile,darkid=darkid,flatid=flatid,littrowid=littrowid,waveid=waveid,fiberid=fiberid
      makecal,psf=responsestr[i].psf
      makecal,wave=waveid
      makecal,fiber=fiberid
      makecal,littrow=littrowid
      mkflux,response,darkid=darkid,flatid=flatid,psfid=responsestr[i].psf,littrowid=littrowid,waveid=waveid,temp=responsestr[i].temp,clobber=clobber
    endif
  endif

  if keyword_set(wave) then begin
    print,'makecal wave: ', wave
    if wave gt 1 then begin
      i=where(wavestr.name eq wave,nwave)
      if nwave le 0 then begin
        print,'No matching calibration line for ', wave
        return
        stop
      endif
      ims=getnums(wavestr[i[0]].frames)
      cmjd=getcmjd(ims[0],mjd=mjd)
      getcal,mjd,calfile,darkid=darkid,flatid=flatid,bpmid=bpmid,fiberid=fiberid
      makecal,bpm=bpmid
      makecal,fiber=fiberid
      mkwave,ims,name=wavestr[i[0]].name,darkid=darkid,flatid=flatid,psfid=wavestr[i[0]].psfid,fiberid=fiberid,clobber=clobber,nofit=nofit
    endif else begin
      if keyword_set(mjd) then  begin
        num=getnum(mjd) 
        red=where(wavestr.frames/10000L eq num)
      endif else red=indgen(n_elements(wavestr))
      if (red[0] ge 0) then begin
       for i=0,n_elements(red)-1,nskip do begin
        ims=getnums(wavestr[red[i]].frames)
        cmjd=getcmjd(ims[0],mjd=mjd)
        getcal,mjd,calfile,darkid=darkid,flatid=flatid,bpmid=bpmid,fiberid=fiberid
        makecal,bpm=bpmid
        makecal,fiber=fiberid
        mkwave,ims,name=wavestr[red[i]].name,darkid=darkid,flatid=flatid,psfid=wavestr[red[i]].psfid,fiberid=fiberid,clobber=clobber,/nowait,nofit=nofit
       endfor
      endif
    endelse
  endif
  if keyword_set(multiwave) then begin
    print,'makecal multiwave: ', multiwave
    if multiwave gt 1 then begin
      i=where(multiwavestr.name eq multiwave,nwave)
      if nwave le 0 then begin
        print,'No matching calibration line for ', multiwave
        stop
      endif
      ims=getnums(multiwavestr[i[0]].frames)
      mkmultiwave,ims,name=multiwavestr[i[0]].name,clobber=clobber,file=file
    endif else begin
      if keyword_set(mjd) then  begin
        num=getnum(mjd) 
        red=where(multiwavestr.frames/10000L eq num)
      endif else red=indgen(n_elements(multiwavestr))
      if (red[0] ge 0) then begin
       for i=0,n_elements(red)-1,nskip do begin
        ims=getnums(multiwavestr[red[i]].frames)
        mkmultiwave,ims,name=multiwavestr[red[i]].name,clobber=clobber,file=file,/nowait
       endfor
      endif
    endelse
  endif
  if keyword_set(lsf) then begin
    print,'makecal lsf: ', lsf
    if lsf gt 1 then begin
      i=where(lsfstr.name eq lsf,nlsf)
      if nlsf le 0 then begin
        print,'No matching calibration line for ', lsf
        stop
      endif
      ims=getnums(lsfstr[i[0]].frames)
      cmjd=getcmjd(ims[0],mjd=mjd)
      getcal,mjd,calfile,darkid=darkid,flatid=flatid,multiwaveid=waveid
      mklsf,ims,waveid,darkid=darkid,flatid=flatid,psfid=lsfstr[i[0]].psfid,full=full,newwave=newwave,clobber=clobber
    endif else begin
      if keyword_set(mjd) then  begin
        num=getnum(mjd) 
        red=where(lsfstr.frames/10000L eq num)
      endif else red=indgen(n_elements(lsfstr))
      if (red[0] ge 0) then begin
       for i=0,n_elements(red)-1,nskip do begin
        ims=getnums(lsfstr[red[i]].frames)
        cmjd=getcmjd(ims[0],mjd=mjd)
        getcal,mjd,calfile,darkid=darkid,flatid=flatid,waveid=waveid
        mklsf,ims,waveid,darkid=darkid,flatid=flatid,psfid=lsfstr[i].psfid,full=full,newwave=newwave,clobber=clobber
       endfor
      endif
    endelse
  endif
end
