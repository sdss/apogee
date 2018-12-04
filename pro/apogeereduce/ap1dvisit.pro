;+
;
; AP1DVISIT
;
; This program processes 1D APOGEE spectra.  It does dither
; combination, wavelength calibration and sky correction.
;
; INPUTS:
;  planfiles  Input list of plate plan files
;  /clobber   Don't use the apCframe files previously created (if it exists)
;  /verbose   Print a lot of information to the screen
;  /stp       Stop at the end of the prrogram
;
; OUTPUTS:
;  1D dither combined, wavelength calibrated, and sky corrected
;  spectra.  4Kx300x5 with planes:
;         1. flux (sky subtracted, absorption/flux corrected)
;         2. (wavelength, if needed)
;         3. flux error
;         4. sky estimate shifted to wavelength sampling of object
;         5. telluric/flux response shifted to wavelength sampling of object
;         6. flags
;  The names are apSpec-[abc]-PLATE4-MJD5.fits
;
; USAGE:
;  IDL>ap1dvisit,planfiles
;
; Written by D.Nidever  Mar. 2010
; Modifications J. Holtzman 2011+
;-

pro ap1dvisit,planfiles,clobber=clobber,verbose=verbose,stp=stp,newwave=newwave,test=test,mapper_data=mapper_data,halt=halt,dithonly=dithonly

common telluric,convolved_telluric

undefine,convolved_telluric
 
;setdisp,/silent
if n_elements(verbose) eq 0 then verbose=0  ; NOT verbose by default

t0 = systime(1)
t1 = systime(1)

nplanfiles = n_elements(planfiles)
; Not enough inputs
if nplanfiles eq 0 then begin
  print,'Syntax - ap1dvisit,planfiles,clobber=clobber,verbose=verbose,stp=stp'
  return
endif

print,''
print,'RUNNING AP1DVISIT'
print,''
print,strtrim(nplanfiles,2),' PLAN files'

chiptag = ['a','b','c']

;--------------------------------------------
; Loop through the unique PLATE Observations
;--------------------------------------------
FOR i=0L,nplanfiles-1 do begin

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

  ; Get APOGEE directories
  dirs=getdir(apogee_dir,cal_dir,spectro_dir,apred_vers=apred_vers,datadir=datadir)
  logfile=apogee_filename('Diag',plate=planstr.plateid,mjd=planstr.mjd)

  ; only process "normal" plates
  if tag_exist(planstr,'platetype') then $
    if planstr.platetype ne 'normal' and $
       planstr.platetype ne 'twilight' and $
       planstr.platetype ne 'sky' and $
       planstr.platetype ne 'single' and $
       planstr.platetype ne 'cal' $
       then goto,BOMB

  if tag_exist(planstr,'survey') then survey=planstr.survey else survey='apogee'

  ; Load the Plug Plate Map file
  ;------------------------------
  print,'' & print,'Plug Map file information:'
  if tag_exist(planstr,'fixfiberid') then fixfiberid=planstr.fixfiberid 
  if tag_exist(planstr,'badfiberid') then badfiberid=planstr.badfiberid 
  if planstr.platetype eq 'single' then begin
   plugfile=getenv('APOGEEREDUCE_DIR')+'/data/plPlugMapA-0001.par' 
   plugmap=getplatedata(planstr.plateid,string(planstr.mjd,format='(i5.5)'),obj1m=planstr.apexp[0].singlename,starfiber=planstr.apexp[0].single,fixfiberid=fixfiberid)
  endif else if planstr.platetype eq 'twilight' then begin
   plugmap=getplatedata(planstr.plateid,string(planstr.mjd,format='(i5.5)'),/twilight)
  endif else if planstr.platetype eq 'cal' then begin
   print,'no plugmap for cal frames'
  endif else begin
    plugfile = planstr.plugmap
    ;APLOADPLUGMAP,plugfile,plugmap,/verbose,error=plugerror,fixfiberid=fixfiberid
    plugmap=getplatedata(planstr.plateid,string(planstr.mjd,format='(i5.5)'),plugid=planstr.plugmap,fixfiberid=fixfiberid,badfiberid=badfiberid,mapper_data=mapper_data)
  endelse
  if keyword_set(stp) then stop

  if n_elements(plugerror) gt 0 then goto,BOMB

  if planstr.platetype ne 'cal' then begin
    plugmap.mjd = planstr.mjd   ; enter MJD from the plan file

    ; get objects
    if planstr.platetype eq 'single' then obj=where(plugmap.fiberdata.objtype ne 'SKY' and plugmap.fiberdata.spectrographid eq 2) else $
    obj=where(plugmap.fiberdata.objtype ne 'SKY' and plugmap.fiberdata.spectrographid eq 2 and plugmap.fiberdata.mag[1] gt 7.5)

  endif

  ; Check if the calibration files exist
  ;--------------------------------------
  makecal,lsf=planstr.lsfid,/full
  wavefiles = apogee_filename('Wave',chip=chiptag,num=planstr.waveid)
  wavetest = FILE_TEST(wavefiles)
  lsffiles = apogee_filename('LSF',chip=chiptag,num=planstr.lsfid)
  lsftest = FILE_TEST(lsffiles)
  if min(wavetest) eq 0 or min(lsftest) eq 0 then begin
    bd1 = where(wavetest eq 0,nbd1)
    if nbd1 gt 0 then print,wavefiles[bd1],' NOT FOUND'
    bd2 = where(lsftest eq 0,nbd2)
    if nbd2 gt 0 then print,lsffiles[bd2],' NOT FOUND'
    goto,BOMB
  endif

  ; Do the output directories exist?
  plate_dir=apogee_filename('Plate',mjd=planstr.mjd,plate=planstr.plateid,chip='a',/dir)
  if file_test(plate_dir,/directory) eq 0 then FILE_MKDIR,plate_dir
  s=strsplit(plate_dir,'/',/extract)
  if dirs.telescope ne 'apo1m' and planstr.platetype ne 'cal' then begin
   if ~file_test(spectro_dir+'/plates/'+s[-2]) then file_link,'../'+s[-5]+'/'+s[-4]+'/'+s[-3]+'/'+s[-2],spectro_dir+'/plates/'+s[-2]
   cloc=strtrim(string(format='(i)',plugmap.locationid),2)
   file_mkdir,spectro_dir+'/location_id/'+dirs.telescope
   if ~file_test(spectro_dir+'/location_id/'+dirs.telescope+'/'+cloc) then file_link,'../../'+s[-5]+'/'+s[-4]+'/'+s[-3],spectro_dir+'/location_id/'+dirs.telescope+'/'+cloc
  endif

  ; Are there enough files
  nframes = n_elements(planstr.apexp)
  ;if nframes lt 2 then begin
  ;  print,'Need 2 OBSERVATIONS/DITHERS to Process'
  ;  goto,BOMB
  ;endif

  ; Start the plots directory
  plots_dir = plate_dir+'/plots/'
  if file_test(plots_dir,/directory) eq 0 then FILE_MKDIR,plots_dir

  ; do we need to iterate for telluric?
  niter=1
  if tag_exist(planstr,'telliter') then if planstr.telliter gt 0 then niter=2

  undefine,visitstr

  for iter=1,niter  do begin

  undefine,alltellstar
  undefine,allframes

  ; do we already have apPlate file?
  file=apogee_filename('Plate',chip='c',mjd=planstr.mjd,plate=planstr.plateid)
  if file_test(file) and not keyword_set(clobber) then begin
    print,'File already exists: ', file
    goto,dorv
  endif

  ; Process each frame
  ;-------------------
  shiftstr = REPLICATE({index:-1L,framenum:'',$
               shift:999999.0,shifterr:999999.0,pixshift:0.,sn:-1.},nframes)

  ; assume no dithering until we see that a dither has been commanded from the
  ;   header cards
  nodither = 1
  ntellerror=0

  For j=0L,nframes-1 do begin

    t1=systime(1)

    ; for ASDAF plate, fix up the plugmap structure
    if planstr.platetype eq 'asdaf' then begin
      all=where(plugmap.fiberdata.spectrographid eq 2)
      plugmap.fiberdata[all].objtype = 'SKY'
      star=where(plugmap.fiberdata.spectrographid eq 2 and $
                 plugmap.fiberdata.fiberid eq planstr.apexp[j].single) 
      plugmap.fiberdata[star].objtype = 'STAR'
      plugmap.fiberdata[star].tmass_style= planstr.apexp[j].singlename
      plugmap.fiberdata[star].mag=fltarr(5)+99.
      if tag_exist(planstr,'hmag') then plugmap.fiberdata[star].mag[1]=planstr.hmag else plugmap.fiberdata[star].mag[1]=5

    endif

    ; Make the filenames and check the files
    rawfiles = apogee_filename('R',chip=chiptag,num=planstr.apexp[j].name)
    rawinfo = APFILEINFO(rawfiles,/silent)        ; this returns useful info even if the files don't exist
    framenum = rawinfo[0].fid8   ; the frame number
    files = apogee_filename('1D',chip=chiptag,num=framenum)
    info = APFILEINFO(files,/silent)
    okay = (info.exists AND info.sp1dfmt AND info.allchips AND (info.mjd5 eq planstr.mjd) AND $
            ((info.naxis eq 3) OR (info.exten eq 1)))
    if min(okay) lt 1 then begin
      bd = where(okay eq 0,nbd)
      stop,'halt: There is a problem with files: ',strjoin((files)(bd),' ')
    endif

    print,''
    print,'-----------------------------------------'
    print,strtrim(j+1,2),'/',strtrim(nframes,2),'  Processing Frame Number >>',strtrim(framenum,2),'<<'
    print,'-----------------------------------------'

    
    ;------------------------------------------
    ; Correcting and Calibrating the ap1D files
    ;------------------------------------------
    cfiles = apogee_filename('Cframe',chip=chiptag,num=framenum,plate=planstr.plateid,mjd=planstr.mjd)
    if keyword_set(clobber) or iter eq 2 or $
         not file_test(cfiles[0]) or $
         not file_test(cfiles[1]) or $
         not file_test(cfiles[2]) then begin

      writelog,logfile,' 1d processing '+file_basename(files[0])+string(format='(f8.2)',systime(1)-t1)
      ; Load the 1D files
      ;--------------------
      APLOADFRAME,files,frame0,/exthead  ; loading frame 1

      ; Fix INF and NAN
      for k=0,2 do begin
        bdnan = where(finite(frame0.(k).flux) eq 0 or finite(frame0.(k).err) eq 0,nbdnan)
        if nbdnan gt 0 then begin
          frame0.(k).flux[bdnan] = 0.0
          frame0.(k).err[bdnan] = baderr()
          frame0.(k).mask[bdnan] = 1   ; bad
        endif

        ; Fix ERR=0
        bdzero = where(frame0.(k).err le 0,nbdzero)
        if nbdzero gt 0 then begin
          frame0.(k).flux[bdzero] = 0.0
          frame0.(k).err[bdzero] = baderr()
          frame0.(k).mask[bdzero] = 1
        endif
      endfor

      ; Add Wavelength and LSF information to the frame structure
      ;---------------------------------------------------------
      ; Loop through the chips
      for k=0,2 do begin
        chstr = frame0.(k)
        ; Get the LSF calibration data
        FITS_READ,lsffiles[k],lsfcoef,lhead

        if keyword_set(newwave) then begin
          remove_tags, chstr,'WCOEF',newstr
          chstr=newstr
        endif
        ; Add to the chip structure
        ; Wavelength calibration data already added by ap2dproc with ap1dwavecal
        ;if tag_exist(frame0.(0),'WCOEF') and not keyword_set(newwave) then begin
        if tag_exist(chstr,'WCOEF') and not keyword_set(newwave) then begin
          chstr = CREATE_STRUCT(temporary(chstr),'LSFFILE',lsffiles[k],'LSFCOEF',$
                                lsfcoef,'WAVE_DIR',plate_dir,'WAVEFILE',wavefiles[k])
        ; Need wavelength information
        endif else begin
          FITS_READ,wavefiles[k],wcoef,whead,exten=1
          chstr = CREATE_STRUCT(temporary(chstr),'LSFFILE',lsffiles[k],'LSFCOEF',$
                                lsfcoef,'WAVEFILE',wavefiles[k],'WCOEF',wcoef,'WAVE_DIR',plate_dir)
        endelse

        ; Now add this to the final FRAME structure
        if k eq 0 then begin
          frame = CREATE_STRUCT('chip'+chiptag[k],chstr)
        endif else begin
          frame = CREATE_STRUCT(frame,'chip'+chiptag[k],chstr)
        endelse
      endfor
      apgundef,frame0   ; free up memory
      if keyword_set(stp) then stop

      ;----------------------------------
      ; STEP 1:  Measure dither Shift
      ;----------------------------------
      print,'STEP 1: Measuring the DITHER SHIFT with APDITHERSHIFT'
      ; Not first frame, measure shift relative to 1st frame
      dither_commanded = sxpar(frame.(0).header,'DITHPIX')
      print,'dither_commanded: ',dither_commanded
      if j gt 0 then print,'ref_dither_commanded: ',ref_dither_commanded
      print,'nodither: ', nodither
      if j gt 0 then $
        if dither_commanded ne 0 and abs(dither_commanded-ref_dither_commanded) gt 0.002 then nodither=0
      if (j gt 0) and not nodither then begin
        ashift=[0.0,0.0] & ashifterr=0.0
        ;APDITHERSHIFT,ref_frame,frame,ashift,ashifterr
        if tag_exist(planstr,'platetype') then begin
          if planstr.platetype eq 'sky' or planstr.platetype eq 'cal' then begin
            plot=1
            pfile=plate_dir+'/plots/dithershift-'+framenum 
          endif else begin
            ;pfile=0 & plot=0
            plot=1
            pfile=plate_dir+'/plots/dithershift-'+framenum 
          endelse
        endif
           
        if planstr.platetype eq 'single' then nofit=1 else nofit=0
        APDITHERSHIFT,ref_frame,frame,shift,shifterr,/xcorr,pfile=pfile,plot=plot,plugmap=plugmap,nofit=nofit,mjd=planstr.mjd
        if keyword_set(stp) then stop
        print,'Measured dither shift: ',ashift,shift

      ; First frame, reference frame
      endif else begin
        ; note reference frame wants to include sky and telluric!
        ref_frame = frame
        shift = [0.0,0.0] & ashift=[0.0,0.0]
        shifterr = 0.0 & ashifterr=0.0
        if dither_commanded ne 0 then ref_dither_commanded = dither_commanded
        print,'Shift = 0.0'
      endelse
      apaddpar,frame,'APDITHERSHIFT: Measuring the dither shift',/history
      if shift[0] eq 0.0 then apaddpar,frame,'APDITHERSHIFT: This is the REFERENCE FRAME',/history
      apaddpar,frame,'DITHSH',shift[0],' Measured dither shift (pixels)'
      apaddpar,frame,'DITHSLOP',shift[1],' Measured dither shift slope (pixels/fiber)'
      apaddpar,frame,'EDITHSH',shifterr,' Dither shift error (pixels)'
      ;apaddpar,frame,'ADITHSH',ashift,' Measured dither shift (pixels)'
      ;apaddpar,frame,'AEDITHSH',ashifterr,' Dither shift error (pixels)'
      
      writelog,logfile,'  dithershift '+string(format='(f8.2)',systime(1)-t1)+string(format='(f8.2)',systime(1)-t0)

      if tag_exist(planstr,'platetype') then $
        if planstr.platetype ne 'normal' and planstr.platetype ne 'single' and planstr.platetype ne 'twilight' then goto,BOMB1

      ;----------------------------------
      ; STEP 2:  Wavelength Calibrate
      ;----------------------------------
      ; THIS IS NOW DONE AS PART OF AP2DPROC, CAN RERUN AS A TEST OR NEEDED IF NEWWAVE IS SPECIFIED...
      print,'STEP 2: Wavelength Calibrating with AP1DWAVECAL'
      plotfile = plate_dir+'/plots/pixshift_chip-'+framenum 
      if keyword_set(dithonly) then AP1DWAVECAL_REFIT,frame,frame_wave,plugmap=plugmap,/verbose,/plot,pfile=plotfile
      plotfile = plate_dir+'/plots/pixshift-'+framenum 
      if planstr.platetype eq 'twilight' then $
      AP1DWAVECAL,frame,frame_wave,/verbose,/plot,pfile=plotfile else $
      AP1DWAVECAL,frame,frame_wave,plugmap=plugmap,/verbose,/plot,pfile=plotfile

      apgundef,frame  ; free up memory
      writelog,logfile,'  wavecal '+string(format='(f8.2)',systime(1)-t1)+string(format='(f8.2)',systime(1)-t0)

      if keyword_set(dithonly) then goto, BOMB1
      ;stop

      ;----------------------------------
      ; STEP 3:  Airglow Subtraction
      ;----------------------------------
      print,'STEP 3: Airglow Subtraction with APSKYSUB'
      APSKYSUB,frame_wave,plugmap,frame_skysub,subopt=1,error=skyerror 
      if n_elements(skyerror) gt 0 and planstr.platetype ne 'twilight' then begin
        stop,'halt: APSKYSUB Error: ',skyerror
        apgundef,frame_wave,frame_skysub,skyerror
      endif
      apgundef,frame_wave  ; free up memory
      writelog,logfile,'  airglow '+string(format='(f8.2)',systime(1)-t1)+string(format='(f8.2)',systime(1)-t0)
      if keyword_set(stp) then stop

      if tag_exist(planstr,'platetype') then $
        if planstr.platetype ne 'normal' and planstr.platetype ne 'single' and planstr.platetype ne 'twilight' then goto,BOMB1

      ;----------------------------------
      ; STEP 4:  Telluric Correction
      ;----------------------------------
      print,'STEP 4: Telluric Correction with APTELLURIC'
      if planstr.platetype eq 'single' then begin
        starfit=2  & single=1 & pltelstarfit=1
      endif else if planstr.platetype eq 'twilight' then begin
        starfit=0  
      endif else begin
        starfit=1 & single=0 & pltelstarfit=0 & visitstr=0
      endelse
      if tag_exist(planstr,'pltelstarfit') then $
         pltelstarfit=planstr.pltelstarfit
      if tag_exist(planstr,'usetelstarfit') then $
         usetelstarfit=1 else usetelstarfit=0
      if tag_exist(planstr,'maxtellstars') then $
         maxtellstars=planstr.maxtellstars else maxtellstars=0
      if tag_exist(planstr,'tellzones') then $
         tellzones=planstr.tellzones else tellzones=0
      APTELLURIC,frame_skysub,plugmap,frame_telluric,tellstar,starfit=starfit,$
        single=single,pltelstarfit=pltelstarfit,usetelstarfit=usetelstarfit,$
        maxtellstars=maxtellstars,tellzones=tellzones,specfitopt=1,$
        plots_dir=plots_dir,error=telerror,/save,/preconv,visitstr=visitstr,test=test
      tellstar.im=planstr.apexp[j].name
      ADD_TAG,frame_telluric,'TELLSTAR',tellstar,frame_telluric
      if n_elements(alltellstar) eq 0 then alltellstar=tellstar else alltellstar=[alltellstar,tellstar]
      if n_elements(telerror) gt 0 and planstr.platetype ne 'single' then begin
        print,'not halted: APTELLURIC Error: ',telerror
        ntellerror+=1
        apgundef,frame_skysub,frame_telluric,telerror
        goto, BOMB1
      endif
      apgundef,frame_skysub  ; free up memory
      writelog,logfile,'  telluric '+string(format='(f8.2)',systime(1)-t1)+string(format='(f8.2)',systime(1)-t0)
      ;-----------------------
      ; Output apCframe files
      ;-----------------------
      print,'Writing output apCframe files'
      outfiles = apogee_filename('Cframe',chip=chiptag,num=framenum,plate=planstr.plateid,mjd=planstr.mjd)
      APVISIT_OUTCFRAME,frame_telluric,plugmap,outfiles,/silent

    endif  ; correcting and calibrating ap1D files


    ;---------------------------------------------
    ; Using the apCframe files previously created
    ;---------------------------------------------

    ; Make the filenames and check the files
    ; Cframe files
    cfiles = apogee_filename('Cframe',chip=chiptag,num=framenum,plate=planstr.plateid,mjd=planstr.mjd)
    cinfo = APFILEINFO(cfiles,/silent)
    okay = (cinfo.exists AND cinfo.allchips AND (cinfo.mjd5 eq planstr.mjd) AND $
            ((cinfo.naxis eq 3) OR (cinfo.exten eq 1)))
    if min(okay) lt 1 then begin
      bd = where(okay eq 0,nbd)
      stop,'halt: There is a problem with files: ',strjoin((cfiles)(bd),' ')
    endif

    print,'Using apCframe files: '+cfiles

    ; Load the apCframe file
    APLOADCFRAME,cfiles,frame_telluric,/exthead

    ; Get the dither shift information from the header
    if j eq 0 then begin
      ref_dither_commanded = sxpar(frame_telluric.(0).header,'DITHPIX')
    endif else begin
      dither_commanded = sxpar(frame_telluric.(0).header,'DITHPIX')
      if dither_commanded ne 0 and abs(dither_commanded-ref_dither_commanded) gt 0.002 then nodither=0
    endelse
    shift = sxpar(frame_telluric.(0).header,'DITHSH')
    shifterr = sxpar(frame_telluric.(0).header,'EDITHSH')
    pixshift = sxpar(frame_telluric.(0).header,'MEDWSH')

    ; Add to the ALLFRAMES structure
    ;--------------------------------
    if n_elements(allframes) eq 0 then allframes=frame_telluric else allframes=[allframes,frame_telluric]

    ; Update SHIFTSTR
    shiftstr[j].index = j
    shiftstr[j].framenum = framenum
    shiftstr[j].shift = shift
    shiftstr[j].shifterr = shifterr
    shiftstr[j].pixshift = pixshift
    ; get S/N of brightest non-saturated object, just for sorting by S/N
    if planstr.platetype eq 'single' then obj=where(plugmap.fiberdata.objtype ne 'SKY' and plugmap.fiberdata.spectrographid eq 2) else $
    obj=where(plugmap.fiberdata.objtype ne 'SKY' and plugmap.fiberdata.spectrographid eq 2 and plugmap.fiberdata.mag[1] gt 7.5 and plugmap.fiberdata.fiberid ne 195)
    hmag=plugmap.fiberdata[obj].mag[1]
    isort=sort(hmag)
    ibright=obj[isort[0]]
    fbright=300-plugmap.fiberdata[ibright].fiberid
    shiftstr[j].sn = median(frame_telluric.(1).flux[*,fbright]/frame_telluric.(1).err[*,fbright])
    if keyword_set(stp) then stop

    ; instead of using S/N of brightest object, which is subject to any issue with that object,
    ;  use median frame zeropoint instead
    if planstr.platetype ne 'single' then begin
      fiberloc=300-plugmap.fiberdata[obj].fiberid
      zero=median(hmag+2.5*alog10(median(frame_telluric.(1).flux[*,fiberloc],dim=1)))
      shiftstr[j].sn=zero
    endif
    apgundef,frame_telluric  ; free up memory

    BOMB1:


  ENDFOR  ; frame loop

  if keyword_set(dithonly) then return


  ; write summary telluric file
  if planstr.platetype eq 'single' then begin

    tellstarfile=$
      planstr.plate_dir+'/apTellstar-'+strtrim(planstr.mjd,2)+'-'+strtrim(planstr.name,2)+'.fits'  
    mwrfits,alltellstar,tellstarfile,/create
  endif else if planstr.platetype eq 'normal'  then begin
    tellstarfile=apogee_filename('Tellstar',plate=planstr.plateid,mjd=planstr.mjd)
    mwrfits,alltellstar,tellstarfile,/create
  endif
  t1=systime(1)

  if tag_exist(planstr,'platetype') then $
    if planstr.platetype ne 'normal' and planstr.platetype ne 'single' and planstr.platetype ne 'twilight' then goto,BOMB

  ; Remove frames that had problems from SHIFTSTR
  bdframe = where(shiftstr.index eq -1,nbdframe)
  if nodither eq 1 then minframes=1 else minframes=2
  if (nframes-nbdframe) lt minframes then begin
    print,'halt: Error: dont have two good frames to proceed'
    stop
    goto, BOMB
  endif
  if nbdframe gt 0 then REMOVE,bdframe,shiftstr
  ; stop with telluric errors
  if ntellerror gt 0 then begin
    print, ntellerror,' frames had APTELLURIC errors'
    if keyword_set(halt) then stop,'halt: ', ntellerror,' frames had APTELLURIC errors'
  endif

  ;----------------------------------
  ; STEP 5:  Dither Combining
  ;----------------------------------
  print,'STEP 5: Combining DITHER FRAMES with APDITHERCOMB'
;  APDITHERCOMB,allframes,shiftstr,pairstr,plugmap,combframe
  APDITHERCOMB,allframes,shiftstr,pairstr,plugmap,combframe,/median,/newerr,npad=50,nodither=nodither,/verbose

  writelog,logfile,' dithercomb '+file_basename(planfile)+string(format='(f8.2)',systime(1)-t1)+string(format='(f8.2)',systime(1)-t0)
  if n_elements(pairstr) eq 0 and nodither eq 0 then begin
    print,'halt: Error: no dither pairs'
    stop
    goto, BOMB
  endif

  ;----------------------------------
  ; STEP 6:  Flux Calibration
  ;----------------------------------
  print,'STEP 6: Flux Calibration with AP1DFLUXING'
  if planstr.platetype ne 'single' then begin
    fiberloc=300-plugmap.fiberdata[obj].fiberid
    zero=median(hmag+2.5*alog10(median(combframe.(1).flux[*,fiberloc],dim=1)))
  endif else zero=0.
  AP1DFLUXING,combframe,plugmap,finalframe

  ;----------------------------------------------------
  ; Output apPlate frames and individual visit spectra, and load apVisit headers with individual star info
  ;----------------------------------------------------
  print,'Writing output apPlate and apVisit files'
  if planstr.platetype eq 'single' then single=1 else single=0
  undefine,mjdfrac
  if tag_exist(planstr,'mjdfrac') then if planstr.mjdfrac eq 1 then $
    mjdfrac=sxpar(finalframe.(0).header,'JD-MID')-2400000.5 
  APVISIT_OUTPUT,finalframe,plugmap,shiftstr,pairstr,$
    /silent,single=single,iter=iter,mjdfrac=mjdfrac,survey=survey
  writelog,logfile,' output '+file_basename(planfile)+string(format='(f8.2)',systime(1)-t1)+string(format='(f8.2)',systime(1)-t0)

  ;---------------
  ; Radial velocity measurements for this visit
  ;--------------
  dorv:
  if tag_exist(planstr,'platetype') then $
    if planstr.platetype ne 'normal' and planstr.platetype ne 'single' then goto,BOMB
  print,'Radial velocity measurements'
  locid=plugmap.locationid
  visitstrfile = apogee_filename('VisitSum',plate=planstr.plateid,mjd=planstr.mjd,reduction=plugmap.fiberdata[obj].tmass_style)
  if tag_exist(planstr,'mjdfrac') then if planstr.mjdfrac eq 1 then begin
    cmjd=strtrim(string(planstr.mjd),2)
    s=strsplit(visitstrfile,cmjd,/extract,/regex)
    visitstrfile=s[0]+string(format='(f8.2)',mjdfrac)+s[1]
  endif
  if file_test(visitstrfile) and not keyword_set(clobber) and niter ne 2 then begin
    print,'File already exists: ', visitstrfile
    return
  endif
  outdir = file_dirname(visitstrfile)
  if file_test(outdir,/directory) eq 0 then FILE_MKDIR,outdir

  objind = where(plugmap.fiberdata.spectrographid eq 2 and $
                 plugmap.fiberdata.holetype eq 'OBJECT' and $
                 plugmap.fiberdata.objtype ne 'SKY',nobjind)
  objdata = plugmap.fiberdata[objind]
  obj = plugmap.fiberdata[objind].tmass_style
  targid = plugmap.fiberdata[objind].object
  targ1 = plugmap.fiberdata[objind].target1
  targ2 = plugmap.fiberdata[objind].target2
  targ3 = plugmap.fiberdata[objind].target3

  plate = plugmap.plateid
  mjd = plugmap.mjd
  platemjd5 = strtrim(plate,2)+'-'+strtrim(mjd,2)
  if keyword_set(stp) then stop

  apgundef,allvisitstr
  ; Plots directory
  for istar=0,n_elements(objind)-1 do begin
    visitfile=apogee_filename('Visit',plate=planstr.plateid,mjd=planstr.mjd,fiber=objdata[istar].fiberid,reduction=obj)
    if tag_exist(planstr,'mjdfrac') then if planstr.mjdfrac eq 1 then begin
      cmjd=strtrim(mjd,2)
      s=strsplit(visitfile,cmjd,/extract,/regex)
      visitfile=s[0]+cmjd+s[1]+string(format='(f8.2)',mjdfrac)+s[2]
    endif

    visitstr = {apogee_id:'',target_id:'',file:'',fiberid:0,plate:'0',mjd:0L,telescope:'',$
           survey:'',field:'',programname:'',$
           location_id:0,ra:0.0d0,dec:0.0d0,glon:0.0d0,glat:0.0d0,$
           j:0.0,j_err:0.0,h:0.0,h_err:0.0,k:0.0,k_err:0.0,$
           ak_targ:-99., ak_targ_method: 'NONE', ak_wise: -99., sfd_ebv: -99.,$
           ra_targ: 0.0d0, dec_targ: 0.0d0, $
           apogee_target1:0L,apogee_target2:0L,apogee_target3:0L,$
           targflags:'',snr: 0.0, starflag:0L,starflags: '',$
           dateobs:'',jd:0.0d0,BC:0.0,vtype:0,$
           VREL:999999.0,vrelerr:999999.0,VHELIO:999999.0,Vlsr:999999.0,Vgsr:999999.0,$
           chisq:0.0,rv_teff:0.0,rv_feh:99.9,rv_logg:99.9,rv_alpha: 99.9,rv_carb: 99.9,synthfile:''}

    visitstr.apogee_id=obj[istar]
    visitstr.target_id=targid[istar]
    visitstr.fiberid=objdata[istar].fiberid
    visitstr.plate=planstr.plateid
    visitstr.mjd=planstr.mjd
    visitstr.location_id=locid
    visitstr.telescope=dirs.telescope
    visitstr.apogee_target1=targ1[istar]
    visitstr.apogee_target2=targ2[istar]
    visitstr.apogee_target3=targ3[istar]
    visitstr.targflags=targflag(targ1[istar],targ2[istar],targ3[istar],survey=survey)
    visitstr.survey=survey
    visitstr.field=plugmap.field
    visitstr.programname=plugmap.programname
    visitstr.ak_targ=plugmap.fiberdata[objind[istar]].ak_targ
    visitstr.ak_targ_method=plugmap.fiberdata[objind[istar]].ak_targ_method
    visitstr.ak_wise=plugmap.fiberdata[objind[istar]].ak_wise
    visitstr.sfd_ebv=plugmap.fiberdata[objind[istar]].sfd_ebv
    aprv,visitfile,visitstr,/save,dir_plots=plots_dir
    MWRFITS,visitstr,visitfile,/silent
    PUSH,allvisitstr,visitstr
  endfor
  writelog,logfile,' aprv '+file_basename(planfile)+string(format='(f8.2)',systime(1)-t1)+string(format='(f8.2)',systime(1)-t0)

  ; Save all RV info for all stars to apVisitSum file
  ;---------------------------
  ; HDU0 - header only
  MKHDR,head0,0
  sxaddpar,head0,'PLATEID',planstr.plateid
  sxaddpar,head0,'MJD',planstr.mjd
  sxaddpar,head0,'EXPTIME',sxpar(finalframe.(0).header,'EXPTIME'),'Total visit exptime per dither pos'
  sxaddpar,head0,'JD-MID',sxpar(finalframe.(0).header,'JD-MID'),' JD at midpoint of visit'
  sxaddpar,head0,'UT-MID',sxpar(finalframe.(0).header,'UT-MID'),' Date at midpoint of visit'
  ncombine = sxpar(finalframe.(0).header,'NCOMBINE',count=num_ncombine)
  if num_ncombine eq 0 then ncombine=1
  sxaddpar,head0,'NCOMBINE',ncombine
  sxaddpar,head0,'ZEROPT',zero
  for j=0,ncombine-1 do sxaddpar,header,'FRAME'+strtrim(j+1,2),sxpar(finalframe.(0).header,'FRAME'+strtrim(j+1,2)),'Constituent frame'
  sxaddpar,head0,'NPAIRS',sxpar(finalframe.(0).header,'NPAIRS'),' Number of dither pairs combined'

  leadstr = 'AP1DVISIT: '
  sxaddpar,head0,'V_APRED',getvers()
  sxaddhist,leadstr+systime(0),head0
  info = GET_LOGIN_INFO()
  sxaddhist,leadstr+info.user_name+' on '+info.machine_name,head0
  sxaddhist,leadstr+'IDL '+!version.release+' '+!version.os+' '+!version.arch,head0
  sxaddhist,leadstr+dirs.prefix+'Visit information for '+strtrim(n_elements(allvisitstr),2)+' Spectra',head0
  FITS_WRITE,visitstrfile,0,head0
  ; HDU1 - structure
  MWRFITS,allvisitstr,visitstrfile,/silent

 endfor
 BOMB:
ENDFOR

print,'AP1DVISIT finished'
writelog,logfile,'AP1DVISIT '+file_basename(planfile)+string(format='(f8.2)',systime(1)-t1)+string(format='(f8.2)',systime(1)-t0)
dt = systime(1)-t0
print,'dt = ',strtrim(string(dt,format='(F10.1)'),2),' sec'

;stop

if keyword_set(stp) then stop

end
