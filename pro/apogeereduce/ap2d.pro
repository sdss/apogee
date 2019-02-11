;+
;
; AP2D
;
; This program processes 2D APOGEE spectra.  It extracts the
; spectra.
;
; INPUTS:
;  planfiles  Input list of plate plan files
;  /verbose   Print a lot of information to the screen
;  /stp       Stop at the end of the prrogram
;
; OUTPUTS:
;  1D extracted spectra are output.  One file for each frame.
;
; USAGE:
;  IDL>ap2d
;
; Written by D.Nidever  Mar. 2010
; Modifications: J. Holtzman 2011+
;-

pro ap2d,planfiles,verbose=verbose,stp=stp,clobber=clobber,exttype=exttype,mapper_data=mapper_data

common savedepsf, savedepsffiles, epsfchip

savedepsffiles=[' ',' ',' ']
epsfchip=0

if n_elements(verbose) eq 0 then verbose=0  ; NOT verbose by default
; clobber will redo PSF, Flux and 1D frames (but not other fundamental calibration frames)
if not keyword_set(clobber) then clobber=0  ; NOT clobber by default,
if not keyword_set(exttype) then exttype=4

t0 = systime(1)

nplanfiles = n_elements(planfiles)
; Not enough inputs
if nplanfiles eq 0 then begin
  print,'Syntax - ap2d,planfiles'
  return
endif

print,''
print,'RUNNING AP2D'
print,''
print,strtrim(nplanfiles,2),' PLAN files'

chiptag = ['a','b','c']
apgundef,wavefile,responsefile

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
  APLOADPLAN,planfile,planstr,/verbose,error=planerror
  if n_elements(planerror) gt 0 then goto,BOMB
  logfile=apogee_filename('Diag',plate=planstr.plateid,mjd=planstr.mjd)

  ; don't extract dark frames
  if tag_exist(planstr,'platetype') then $
    if planstr.platetype eq 'dark' or planstr.platetype eq 'intflat' then goto,BOMB

  ; Try to make the required calibration files (if not already made)
  ; Then check if the calibration files exist
  ;--------------------------------------

  ; apPSF files 
  if planstr.sparseid ne 0 then makecal,sparse=planstr.sparseid
  if planstr.fiberid ne 0 then makecal,fiber=planstr.fiberid
  
  makecal,psf=planstr.psfid,clobber=clobber
  tracefiles = apogee_filename('PSF',num=planstr.psfid,chip=chiptag)
  tracefile = file_dirname(tracefiles[0])+'/'+string(format='(i8.8)',planstr.psfid)
  tracetest = FILE_TEST(tracefiles)  
  if min(tracetest) eq 0 then begin
    bd1 = where(tracetest eq 0,nbd1)
    if nbd1 gt 0 then stop,'halt: ',tracefiles[bd1],' NOT FOUND'
    for ichip=0,2 do begin
      p=mrdfits(tracefiles[ichip],1,/silent)
      if n_elements(p) ne 300 then begin
        print, 'halt: tracefile ', tracefiles[ichip],' does not have 300 traces'
      endif
    endfor
  endif

  ; apWave files : wavelength calibration
  waveid=planstr.waveid
  if tag_exist(planstr,'platetype') then if planstr.platetype eq 'cal' then waveid=0
  if waveid gt 0 then begin
    makecal,multiwave=waveid
  ; move 1dwavecal to new skycal below
  ;  wavefiles = apogee_filename('Wave',chip=chiptag,num=waveid)
  ;  wavefile = file_dirname(wavefiles[0])+'/'+string(format='(i8.8)',waveid)
  ;  wavetest = FILE_TEST(wavefiles)
  ;  if min(wavetest) eq 0 then begin
  ;    bd1 = where(wavetest eq 0,nbd1)
  ;    if nbd1 gt 0 then stop,'halt: ',wavefiles[bd1],' NOT FOUND'
  ;  endif
  endif

  ; apFlux files : since individual frames are usually made per plate,
  if planstr.fluxid ne 0 then begin
    makecal,flux=planstr.fluxid,psf=planstr.psfid,clobber=clobber
    fluxfiles = apogee_filename('Flux',chip=chiptag,num=planstr.fluxid)
    fluxfile = file_dirname(fluxfiles[0])+'/'+string(format='(i8.8)',planstr.fluxid)
    fluxtest = FILE_TEST(fluxfiles)  
    if min(fluxtest) eq 0 then begin
      bd1 = where(fluxtest eq 0,nbd1)
      if nbd1 gt 0 then stop,'halt: ',fluxfiles[bd1],' NOT FOUND'
    endif
  endif else fluxtest=0

  ; apResponse files 
  if tag_exist(planstr,'responseid') eq 0 then add_tag,planstr,'responseid',0,planstr
  if planstr.responseid ne 0 then begin
    makecal,response=planstr.responseid
    responsefiles = apogee_filename('Response',chip=chiptag,num=planstr.responseid)
    responsefile = file_dirname(responsefiles[0])+'/'+string(format='(i8.8)',planstr.responseid)
    responsetest = FILE_TEST(responsefiles)  
    if min(responsetest) eq 0 then begin
      bd1 = where(responsetest eq 0,nbd1)
      if nbd1 gt 0 then stop,'halt: ',responsefiles[bd1],' NOT FOUND'
    endif
  endif

  ; Load the Plug Plate Map file
  ;------------------------------
  if tag_exist(planstr,'platetype') then if planstr.platetype eq 'cal' or planstr.platetype eq 'single' then plugmap=0 else begin
    print,'' & print,'Plug Map file information:'
    plugfile = planstr.plugmap
    if tag_exist(planstr,'fixfiberid') then fixfiberid=planstr.fixfiberid
    if tag_exist(planstr,'badfiberid') then badfiberid=planstr.badfiberid
    plugmap=getplatedata(planstr.plateid,string(planstr.mjd,format='(i5.5)'),plugid=planstr.plugmap,fixfiberid=fixfiberid,badfiberid=badfiberid,mapper_data=mapper_data)
    if n_elements(plugerror) gt 0 then stop,'halt: error with plugmap: ',plugfile
    plugmap.mjd = planstr.mjd   ; enter MJD from the plan file
  endelse

  ; Are there enough files
  nframes = n_elements(planstr.apexp)
  if nframes lt 1 then begin
    print,'No frames to process'
    goto,BOMB
  endif

  ; Process each frame
  ;-------------------
  For j=0L,nframes-1 do begin

    ; Make the filenames and check the files

    rawfiles = apogee_filename('R',chip=chiptag,num=planstr.apexp[j].name)
    rawinfo = APFILEINFO(rawfiles,/silent)        ; this returns useful info even if the files don't exist
    framenum = rawinfo[0].fid8   ; the frame number
    files = apogee_filename('2D',chip=chiptag,num=framenum)
    inpfile = file_dirname(files[0])+'/'+framenum
    info = APFILEINFO(files,/silent)
    okay = (info.exists AND info.sp2dfmt AND info.allchips AND (info.mjd5 eq planstr.mjd) AND $
            ((info.naxis eq 3) OR (info.exten eq 1)))
    if min(okay) lt 1 then begin
      bd = where(okay eq 0,nbd)
      stop,'halt: There is a problem with files: ',strjoin((files)(bd),' ')
    endif

    print,''
    print,'-----------------------------------------'
    print,strtrim(j+1,2),'/',strtrim(nframes,2),'  Processing Frame Number >>',strtrim(framenum,2),'<<'
    print,'-----------------------------------------'

    ; Run AP2DPROC
    if tag_exist(planstr,'platetype') then if planstr.platetype eq 'cal' or planstr.platetype eq 'single' then skywave=0 else skywave=1
    if tag_exist(planstr,'platetype') then if planstr.platetype eq 'sky' then plugmap=0
    outdir=apogee_filename('1D',num=framenum,chip='a',/dir)
    if file_test(outdir,/directory) eq 0 then FILE_MKDIR,outdir
    if min(fluxtest) eq 0 or planstr.apexp[j].flavor eq 'flux' then $
      AP2DPROC,inpfile,tracefile,exttype,outdir=outdir,$
      wavefile=wavefile,skywave=skywave,plugmap=plugmap,clobber=clobber,/compress $
    else if waveid gt 0 then begin
      AP2DPROC,inpfile,tracefile,exttype,outdir=outdir,$
       fluxcalfile=fluxfile,responsefile=responsefile,$
       wavefile=wavefile,skywave=skywave,plugmap=plugmap,clobber=clobber,/compress 
    endif else $
      AP2DPROC,inpfile,tracefile,exttype,outdir=outdir,$
       fluxcalfile=fluxfile,responsefile=responsefile,$
       clobber=clobber,/compress 

    BOMB1:

  Endfor ; frame loop
  
  ; now add in wavelength calibration information, with shift from skylines
  if waveid gt 0 then begin
      if skywave then spawn,['apskywavecal',planfile],/noshell $
      else  spawn,['apskycal',planfile,'--nosky'],/noshell
  endif

  BOMB:

  ; compress 2D files
  nframes = n_elements(planstr.apexp)
  for j=0L,nframes-1 do begin
    files = apogee_filename('2D',num=planstr.apexp[j].name,chip=chiptag)
    modfiles = apogee_filename('2Dmodel',num=planstr.apexp[j].name,chip=chiptag)
    for jj=0,n_elements(files)-1 do begin
      if file_test(files[jj]) then begin
        file_delete,files[jj]+'.fz',/allow_nonexistent
 ;       SPAWN,['fpack','-D','-Y',files[jj]],/noshell
      endif
      if file_test(modfiles[jj]) then begin
        file_delete,modfiles[jj]+'.fz',/allow_nonexistent
        SPAWN,['fpack','-D','-Y',modfiles[jj]],/noshell
      endif
    endfor
  endfor


  writelog,logfile,'AP2D: '+file_basename(planfile)+string(format='(f8.2)',systime(1)-t1)

ENDFOR  ; plan file loop

apgundef,epsfchip

print,'AP2D finished'
dt = systime(1)-t0
print,'dt = ',strtrim(string(dt,format='(F10.1)'),2),' sec'

;stop

if keyword_set(stp) then stop

end
