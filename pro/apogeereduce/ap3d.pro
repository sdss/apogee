;+
;
; AP3D
;
; This program processes all of the APOGEE RAW datacubes for
; a single night.
;
; INPUTS:
;  planfiles  Input list of plate plan files
;  /verbose  Print a lot of information to the screen
;  /stp      Stop at the end of the prrogram
;
; OUTPUTS:
;  The RAW APOGEE 3D datacube files are processed and 2D 
;  images are output.
;
; USAGE:
;  IDL>ap3d,planfiles
;
; Written by D.Nidever  Feb. 2010
; Modifications J. Holtzman 2011+
;-

pro ap3d,planfiles,verbose=verbose,stp=stp,rogue=rogue,clobber=clobber,refonly=refonly

if n_elements(verbose) eq 0 then verbose=0  ; NOT verbose by default
if not keyword_set(clobber) then clobber=0  ; NOT clobber by default
t0 = systime(1)

nplanfiles = n_elements(planfiles)
; Not enough inputs
if nplanfiles eq 0 then begin
  print,'Syntax - ap3d,planfiles'
  return
endif

print,''
print,'RUNNING AP3D'
print,''
print,strtrim(nplanfiles,2),' PLAN files'

chiptag = ['a','b','c']

;--------------------------------------------
; Loop through the unique PLATE Observations
;--------------------------------------------
FOR i=0L,nplanfiles-1 do begin
  t1 = systime(1)

  planfile = planfiles[i]

  print,''
  print,'========================================================================='
  print,strtrim(i+1,2),'/',strtrim(nplanfiles,2),'  Processing Plan file ',planfile
  print,'========================================================================='

  ; Load the plan file
  ;--------------------
  print,'' & print,'Plan file information:'
  ;APLOADPLAN,planfile,planstr,/verbose,error=planerror,/expand,/newlog
  APLOADPLAN,planfile,planstr,/verbose,error=planerror,/newlog
  if n_elements(planerror) gt 0 then goto,BOMB
  dirs=getdir(apogee_dir,caldir,spectro_dir,apred_vers)
  prefix=dirs.prefix
  logfile=apogee_filename('Diag',plate=planstr.plateid,mjd=planstr.mjd)

  ; Try to make the required calibration files (if not already made)
  ; Then check if the calibration files exist
  ;--------------------------------------
   ; Det, BPM, Dark, Flat calib files
  apgundef,detfiles,bpmfiles,darkfiles,flatfiles,littrowfiles,persistfiles
  apgundef,persistmodelfiles,histfiles

  ; apDetector file : sets gain and readout noise
  if planstr.detid ne 0 then begin
    mkdet,planstr.detid
    detfiles = apogee_filename('Detector',num=planstr.detid,chip=chiptag)
    dettest = FILE_TEST(detfiles)
    if min(dettest) eq 0 then begin
      bd = where(dettest eq 0,nbd)
      if nbd gt 0 then stop,'halt: ', detfiles[bd],' NOT FOUND'
    endif
  endif

  ; apDark file  : dark frame
  if planstr.darkid ne 0 then begin
    makecal,dark=planstr.darkid
    darkfiles = apogee_filename('Dark',num=planstr.darkid,chip=chiptag)
    darktest = FILE_TEST(darkfiles)
    if min(darktest) eq 0 then begin
      bd = where(darktest eq 0,nbd)
      if nbd gt 0 then stop,'halt: ',darkfiles[bd],' NOT FOUND'
    endif
  endif

  ; apFlat file : flat field
  if planstr.flatid ne 0 then begin
    makecal,flat=planstr.flatid
    flatfiles = apogee_filename('Flat',num=planstr.flatid,chip=chiptag)
    flattest = FILE_TEST(flatfiles)
    if min(flattest) eq 0 then begin
      bd = where(flattest eq 0,nbd)
      if nbd gt 0 then stop,'halt: ',flatfiles[bd],' NOT FOUND'
    endif
  endif

  ; apBPM file : bad pixel mask
  if planstr.bpmid ne 0 then begin
    makecal,bpm=planstr.bpmid
    bpmfiles = apogee_filename('BPM',num=planstr.bpmid,chip=chiptag)
    bpmtest = FILE_TEST(bpmfiles)
    if min(bpmtest) eq 0 then begin
      bd = where(bpmtest eq 0,nbd)
      if nbd gt 0 then stop,'halt: ',bpmfiles[bd],' NOT FOUND'
    endif
  endif

  ; apLittrow file : littrow ghost pixel mask
  if tag_exist(planstr,'littrowid') eq 0 then add_tag,planstr,'littrowid',0,planstr
  if planstr.littrowid ne 0 then begin
    makecal,littrow=planstr.littrowid
    littrowfiles = apogee_filename('Littrow',num=planstr.littrowid,chip='b')
    littrowtest = FILE_TEST(littrowfiles)
    if min(littrowtest) eq 0 then begin
      bd = where(littrowtest eq 0,nbd)
      if nbd gt 0 then stop,'halt: ',littrowfiles[bd],' NOT FOUND'
    endif
  endif

  ; apPersist file : persistence pixel mask
  if tag_exist(planstr,'persistid') eq 0 then add_tag,planstr,'persistid',0,planstr
  if planstr.persistid ne 0 then begin
    makecal,persist=planstr.persistid
    persistfiles = apogee_filename('Persist',num=planstr.persistid,chip=chiptag)
    persisttest = FILE_TEST(persistfiles)
    if min(persisttest) eq 0 then begin
      bd = where(persisttest eq 0,nbd)
      if nbd gt 0 then stop,'halt: ',persistfiles[bd],' NOT FOUND'
    endif
  endif

  ; apPersistModel file : persistence model parameters
  if tag_exist(planstr,'persistmodelid') eq 0 then add_tag,planstr,'persistmodelid',0,planstr
  if planstr.persistmodelid ne 0 then begin
    ;makecal,modelpersist=persistmodelid
    ; The apPersistModel calibration files are created separately
    ;  corrections exist only for "b" (green) and "c" (blue) for now
    persistmodelfiles = apogee_filename('PersistModel',mjd=planstr.persistmodelid,chip=['b','c'])
    persistmodeltest = FILE_TEST(persistmodelfiles)
    if min(persistmodeltest) eq 0 then begin
      bd = where(persistmodeltest eq 0,nbd)
      if nbd gt 0 then stop,'halt: ',persistmodelfiles[bd],' NOT FOUND'
    endif
    makehist,planstr.mjd,dark=planstr.darkid
    histfiles = apogee_filename('Hist',mjd=planstr.mjd,chip=chiptag)
    histtest = FILE_TEST(histfiles)
    if min(histtest) eq 0 then begin
      bd = where(histtest eq 0,nbd)
      if nbd gt 0 then stop,'halt: ',histfiles[bd],' NOT FOUND'
    endif
  endif
  
  ; Are there enough files
  nframes = n_elements(planstr.apexp)
  if nframes lt 1 then begin
    print,'No frames to process'
    goto,BOMB
  endif

  if keyword_set(rogue) then begin
    ; rogue option uses alternate reduction routines as
    ;   an independent test
     objs=where(planstr.apexp.flavor eq 'object')
     ims=lonarr(n_elements(objs))
     single=-1 & singlename='' & mapname='header'
     if objs[0] ge 0 then begin
      for iframe=0,n_elements(objs)-1 do begin
       ;reads,strmid(planstr.apexp[objs[iframe]].name[0],6,8),im
       reads,planstr.apexp[objs[iframe]].name[0],im
       ims[iframe]=im
      endfor
      single=planstr.apexp[objs].single
      singlename=planstr.apexp[objs].singlename
      ;mapname=planstr.apexp[objs].mapname
      mapname=strmid(file_basename(planstr.plugmap),11)
      mapname=strmid(mapname,0,strlen(mapname)-4)
     endif 
     plateid=0 & darkid=0L & flatid=0L & psfid=0L & fluxid=0L & waveid=0L
     reads,planstr.plateid,plateid
     reads,file_basename(planstr.darkid),darkid
     reads,file_basename(planstr.flatid),flatid
     reads,file_basename(planstr.psfid),psfid
     reads,file_basename(planstr.fluxid),fluxid
     reads,file_basename(planstr.waveid),waveid
     apred_holtz,ims,plateid,darkid,flatid,psfid,fluxid,waveid,clobber=clobber,platetype=planstr.platetype,$
        starfiber=single,starnames=singlename,mapname=mapname,logfile=logfile

  endif else begin

  ; Process each frame
  ;-------------------
  For j=0L,nframes-1 do begin

    ; Make the filenames and check the files
    chipfiles = apogee_filename('R',chip=chiptag,num=planstr.apexp[j].name,mjd=planstr.mjd)
    info = APFILEINFO(chipfiles,/silent)
    framenum = info[0].fid8   ; the frame number
    ;okay = (info.exists AND info.rawfmt AND info.allchips AND (info.mjd5 eq planstr.mjd) AND $
    okay = (info.exists AND info.rawfmt AND info.allchips AND $
            ((info.naxis eq 3) OR (info.exten eq 1)))
    ;  must be 3D or have extensions
    if min(okay) lt 1 then begin
      bd = where(okay eq 0,nbd)
      print,'halt: There is a problem with files: ',strjoin((chipfiles)(bd),' ')
      stop
      goto,BOMB1
    endif

    print,''
    print,'------------------------------------------------------'
    print,strtrim(j+1,2),'/',strtrim(nframes,2),'  Processing files for Frame Number >>',strtrim(framenum,2),'<<'
    seq=strtrim(j+1,2)+'/'+strtrim(nframes,2)
    print,'------------------------------------------------------'

    ; Determine file TYPE
    ;----------------------
    ; dark - should be processed with 
    ; flat
    ; lamps
    ; object frame
    ;exptype = strtrim(info[0].exptype,2)
    exptype = planstr.apexp[j].flavor

    ;obstype = SXPAR(head,'OBSTYPE',count=nobs)
    if exptype eq '' or exptype eq '0' then begin
      error = 'NO OBSTYPE keyword found for '+base1
      print,error
      goto,BOMB1
    endif
    exptype = strlowcase(strtrim(exptype,2))

    ; This is a DARK frame
    ;----------------------
    if exptype eq 'dark' then begin
      stop,'halt: This is a DARK frame.  This should be processed with APMKSUPERDARK.PRO'
    end

    ; LOOP through the file types
    CASE exptype OF

      ;------------
      ; FLAT FRAME
      ;------------
      'psf': begin
        print,'This is a FLAT frame'

        ; Settings to use
        ;-----------------
        usedet = 1
        usebpm = 1
        usedark = 1
        useflat = 1
        uselittrow = 1
        usepersist = 1
        dopersistcorr = 0
 	nocr = 1
        crfix = 0
        criter = 0
        satfix = 1
        uptheramp = 0
        nfowler =  1
        rd3satfix = 0
        ;if nreads eq 3 then rd3satfix=1

      end ; lamp frame

      ;------------
      ; LAMP FRAME
      ;------------
      'lamp': begin
        print,'This is a LAMP frame'

        ; Settings to use
        ;-----------------
        usedet = 1
        usebpm = 1
        usedark = 1
        useflat = 1
        uselittrow = 1
        usepersist = 1
        dopersistcorr = 0
        nocr = 0
        crfix = 1
        satfix = 1
        uptheramp = 0
        nfowler = 1
        rd3satfix = 0
        ;if nreads eq 3 then rd3satfix=1

      end ; lamp frame
      'wave': begin
        print,'This is a WAVE frame'

        ; Settings to use
        ;-----------------
        usedet = 1
        usebpm = 1
        usedark = 1
        useflat = 1
        uselittrow = 1
        usepersist = 1
        dopersistcorr = 0
        nocr = 0
        crfix = 1
        criter = 0
        satfix = 1
        uptheramp = 0
        nfowler = 1
        rd3satfix = 0
        ;if nreads eq 3 then rd3satfix=1

      end ; lamp frame

      ;--------------
      ; OBJECT FRAME
      ;--------------
      'object': begin
        print,'This is an OBJECT frame'

        ; Settings to use
        ;-----------------
        usedet = 1
        usebpm = 1
        usedark = 1
        useflat = 1
        uselittrow = 1
        usepersist = 1
        dopersistcorr = 1
        nocr = 0
        if planstr.platetype eq 'single' then nocr=1
        crfix = 1
        criter = 0
        satfix = 1
        uptheramp = 1
        nfowler = 0
        rd3satfix = 0

      end ; object frame

      'flux': begin
        print,'This is an FLUX frame'

        ; Settings to use
        ;-----------------
        usedet = 1
        usebpm = 1
        usedark = 1
        useflat = 1
        uselittrow = 1
        usepersist = 1
        dopersistcorr = 0
        nocr = 1
        crfix = 0
        criter = 0
        satfix = 1
        uptheramp = 0
        nfowler = 1
        rd3satfix = 0

      end ; object frame

      else: begin
        print,exptype+' NOT SUPPORTED'
        goto,BOMB1
      endelse

    ENDCASE ; file type


    ;----------------------------------
    ; Looping through the three chips
    ;----------------------------------
    For k=0,2 do begin

      file = chipfiles[k]

      ; Check header
      head = headfits(file,errmsg=errmsg)
      if errmsg ne '' then begin
        error = 'There was an error loading the HEADER for '+file
        print,error
        goto,BOMB2
      endif

      ; Check that this is a data CUBE OR has extensions
      naxis = sxpar(head,'NAXIS')
      FITS_READ,file,dumim,dumhead,exten_no=1,message=read_message,/no_abort
      if naxis ne 3 and read_message ne '' then begin
        error = 'FILE must contain a 3D DATACUBE OR image extensions'
        print,error
        goto,BOMB2
      endif

      ; Chip specific calibration filenames
      apgundef,detcorr,bpmcorr,darkcorr,flatcorr,littrowcorr,persistcorr
      if usedet and planstr.detid ne 0 then detcorr = detfiles[k]
      if usebpm and planstr.bpmid ne 0 then bpmcorr = bpmfiles[k]
      if usedark and planstr.darkid ne 0 then darkcorr = darkfiles[k]
      if useflat and planstr.flatid ne 0 then flatcorr = flatfiles[k]
      if uselittrow and planstr.littrowid ne 0 and k eq 1 then littrowcorr = littrowfiles
      if usepersist and planstr.persistid ne 0 then persistcorr = persistfiles[k]
      if dopersistcorr and (planstr.persistmodelid ne 0) and (chiptag[k] eq 'c') then begin
      ;        and (chiptag[k] eq 'b' or (chiptag[k] eq 'c' and planstr.mjd lt 56860L)) then begin
        persistmodelcorr = persistmodelfiles[k-1]
        histcorr = histfiles[k]
      endif else apgundef,persistmodelcorr
      ; note this q3fix still fails for apo1m flat/PSF processing, which calls ap3dproc directly
      q3fix=0
      if tag_exist(planstr,'q3fix') then if k eq 2 and planstr.q3fix eq 1 then q3fix=1 
      if k eq 2 and planstr.mjd gt 56930L and planstr.mjd lt 57600L then q3fix=1 

      if tag_exist(planstr,'usereference') then usereference=planstr.usereference else usereference = 1
      if tag_exist(planstr,'maxread') then maxread=planstr.maxread else undefine,maxread

      print,''
      print,'-----------------------------------------'
      print,' Processing chip '+chiptag[k]+' - '+file_basename(file)
      print,'-----------------------------------------'
      print,''

      ; Output file
      outfile = apogee_filename('2D',chip=chiptag[k],num=framenum)
      ; Does the output directory exist?
      if file_test(file_dirname(outfile),/directory) eq 0 then FILE_MKDIR,file_dirname(outfile)

      ; PROCESS the file
      ;-------------------
      ;satfix = 0
      AP3DPROC,file,outfile,detcorr=detcorr,bpmcorr=bpmcorr,darkcorr=darkcorr,littrowcorr=littrowcorr,$
               persistcorr=persistcorr,flatcorr=flatcorr,nocr=nocr,crfix=crfix,criter=criter,satfix=satfix,$
               rd3satfix=rd3satfix,nfowler=nfowler,cleanuprawfile=1,$
               uptheramp=uptheramp,verbose=verbose,error=procerror,clobber=clobber,logfile=logfile,$
               fitsdir=getlocaldir(),q3fix=q3fix,maxread=maxread,persistmodelcorr=persistmodelcorr,histcorr=histcorr,$
               usereference=usereference,refonly=refonly,seq=seq

      BOMB2:

    Endfor ; chip loop

    BOMB1:

  Endfor ; frame loop

  endelse
  BOMB:

  writelog,logfile,'AP3D: '+file_basename(planfile)+string(format='(f8.2)',systime(1)-t1)
ENDFOR  ; plan file looop

print,'AP3D finished'
dt = systime(1)-t0
print,'dt = ',strtrim(string(dt,format='(F10.1)'),2),' sec'

;stop

if keyword_set(stp) then stop

end
