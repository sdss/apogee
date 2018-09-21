;+
;
; APQL_WRAPPER_MANUAL
;
; This is the wrapper to run the APOGEE quicklook program manually,
; without the actor.
;
; INPUTS:
;   inlist: a csv file containing a list of exposures to be processed.  The format is: mjd, exp, plateId, plugfile
;   outfile: if this is specified, then the quicklook output (except for large arrays) will be written to this file.
;             The format will be a list of data in CSV format for the Quicklook, Quicklook60, Quicklook60_imbinzoom,
;             Quicklook60_repspec and Quicklook_prediction tables.
;   no_dbinsert: if this is not set, then the Apogee Quicklook DB will be updated  
; OUTPUTS:
;   A DB insert if no_dbinsert != 1.
;   An output file, if specified by the outfile input parameter.  See INPUTS for decription.
; USAGE:
;   e.g.,
;   apql_wrapper_manual,'inlist.txt', outfile='output.txt', /no_dbinsert
;-
;
pro apql_wrapper_manual, inlist, no_dbinsert=no_dbinsert, outfile=outfile, data_dir=data_dir, spectro_dir=spectro_dir


   RESOLVE_ROUTINE,'apwavecal_chip'

   ; define the IDL environment variables used by idlsql to point to the desired database
   ; this used to be called APQL_DEFSYS
   SDSS_DB_PARAMS

   npix = 2048L
   nchips = 3L
   nbin = 8
   npixbin = npix/nbin

   ;read in data from inlist
   if FILE_TEST(inlist) eq 0 then begin
      msg = "Input file not found: "+inlist
      print,msg
      return
   endif

   data = read_csv(inlist)
   mjds = strtrim(data.field1,2)
   frameids = strtrim(data.field2,2)
   plates = strtrim(data.field3,2)
   plugfiles = strtrim(data.field4,2)
   nframes = n_elements(frameids)

   ; Get APOGEE directories if the environment variables exists
   if not keyword_set(data_dir) then data_dir = APGETDIR('APQLDATA_DIR',/exists,error=direrr)
   if not keyword_set(spectro_dir) then spectro_dir = APGETDIR('APQLSPECTRO_DIR',/exists,error=direrr)

   ; we need data_dir and spectro_dir 
   if strlen(data_dir) eq 0 then begin
      msg = 'Missing data_dir -> aborted'
      print,msg
      return
   endif
   if strlen(spectro_dir) eq 0 then begin
      msg = 'Missing spectro_dir -> aborted'
      print,msg
      return
   endif

   chiptag = ['a','b','c']
   linelist_dir = spectro_dir+'/lib/linelists/'
   datadir = data_dir
   psfdir = spectro_dir+'/cal/psf/'
   bpmdir = spectro_dir+'/cal/bpm/'

   ; Get Psf file
   tinfo = APFILEINFO(FILE_SEARCH(psfdir+'apPSF-a-*.fits'),/silent)
   gdpsf = where(tinfo.exists eq 1 and tinfo.allchips eq 1 and tinfo.suffix ne '',ngdpsf)
   if ngdpsf eq 0 then begin
      msg='No good psf file found in '+psfdir
      print,msg
      return
   endif else begin
      ; sort through the files and use the latest one
      tinfogd = tinfo[gdpsf]       ; the good ones
      si = reverse(sort(tinfogd.mtime))
      psffiles = psfdir+'apPSF-'+chiptag+'-'+tinfogd[si[0]].suffix+'.fits'
   endelse

   ; Get BPM file
   binfo = APFILEINFO(FILE_SEARCH(bpmdir+'apBPM-a-*.fits'),/silent)
   gdbpm = where(binfo.exists eq 1 and binfo.allchips eq 1 and binfo.suffix ne '',ngdbpm)
   if ngdbpm eq 0 then begin
      msg='No good BPM file found in '+bpmdir
      print,msg
      ;return
   endif else begin
      ; sort through the files and use the latest one
      binfogd = binfo[gdbpm]       ; the good ones
      si = reverse(sort(binfogd.mtime))
      bpmfiles = bpmdir+'apBPM-'+chiptag+'-'+binfogd[si[0]].suffix+'.fits'
      bpmid = bpmdir+binfogd[si[0]].suffix
   endelse



   ; Initialize some arrays/structures
   ; Making sure the heap memory is released before deleting the variables
   ; to prevent a ever increasing memory usage while the program runs
   ;----------------------------------------------------------------------
   if size(str,/tname) eq 'STRUCT' then HEAP_FREE, str, /verbose
   if size(allstr,/tname) eq 'STRUCT' then HEAP_FREE, allstr, /verbose
   if size(prevstr,/tname) eq 'STRUCT' then HEAP_FREE, prevstr, /verbose
   if size(all_predictstr,/tname) eq 'STRUCT' then HEAP_FREE, all_predictstr, /verbose
   apgundef,str,allstr,prevstr,all_predictstr


   ; Create !APQL system variable
   ;------------------------------
   DEFSYSV,'!apql',exists=apql_exists
   if not apql_exists then begin
     apqlstr = {plugmap:{filename:'',plate:'',mjd:'',datastr:PTR_NEW()},$
                firstread:{frameid:'',readnum:0L,plate:'',data:PTR_NEW(),header:strarr(1000)},$
                psf:{frameid:'',mjd:'',data:REPLICATE({chipnum:0L,filename:'',tracestr:PTR_NEW(),psfim:PTR_NEW()},nchips)},$
                airglow:{filename:'',data:PTR_NEW()},current_filename:'', timePerRead:0.0}
     DEFSYSV,'!apql',apqlstr
   endif

   ; Load the psf files
   ;----------------------
   ;  This only needs to be loaded ONCE per night
   if !apql.psf.frameid eq '' then begin

     psfbase = file_basename(psffiles,'.fits')
     dum = strsplit(psfbase[0],'-',/extract)
     psfframeid = dum[2]

     psfhead = headfits(psffiles[0])
     !apql.psf.frameid = psfframeid
     !apql.psf.mjd = strtrim(sxpar(psfhead,'MJD'),2)

     for i=0,2 do begin
       tracestr = MRDFITS(psffiles[i],1,/silent)
       psfim = MRDFITS(psffiles[i],2,/silent)
       !apql.psf.data[i].chipnum = i+1
       !apql.psf.data[i].filename = psffiles[i]
       if PTR_VALID(!apql.psf.data[i].tracestr) then PTR_FREE,!apql.psf.data[i].tracestr
       !apql.psf.data[i].tracestr = PTR_NEW(tracestr)
       if PTR_VALID(!apql.psf.data[i].psfim) then PTR_FREE,!apql.psf.data[i].psfim
       !apql.psf.data[i].psfim = PTR_NEW(psfim)
     end

   endif

   ; Load the AIRGLOW linelist
   ;---------------------------
   if (!apql.airglow.filename eq '') then begin
     airstr = IMPORTASCII(linelist_dir+'airglow.txt',/header,/silent)
     nairstr = n_elements(airstr)
     !apql.airglow.filename = linelist_dir+'airglow.txt'
     if PTR_VALID(!apql.airglow.data) then PTR_FREE,!apql.airglow.data
     !apql.airglow.data = PTR_NEW(airstr)
   end

 
  ;for each exposure/frame
  for i=0,nframes-1 do begin
      plate = plates[i]
      mjd = mjds[i]
      frameid = frameids[i]
      plugfile = plugfiles[i]

      ; clear the pointers to free heap memory
      if size(str,/tname) eq 'STRUCT' then HEAP_FREE, str, /verbose
      apgundef,str

      if !apql.plugmap.plate ne plate or !apql.plugmap.mjd ne mjd then begin

        ; we got a new plugmap file
        APLOADPLUGMAP,plugfile,plugmap
        !apql.plugmap.filename = plugfile
        !apql.plugmap.plate = plugmap.plateid
        !apql.plugmap.mjd = plugmap.mjd
        ; clean up the heap for new data
        if PTR_VALID(!apql.plugmap.datastr) then PTR_FREE,!apql.plugmap.datastr
        !apql.plugmap.datastr = PTR_NEW(plugmap)

        ; Sort the plugmap data by fiberid
        fiberid = (*!apql.plugmap.datastr).fiberdata.fiberid
        si = sort(fiberid)
        (*!apql.plugmap.datastr).fiberdata = (*!apql.plugmap.datastr).fiberdata[si]

        ; with a new plate load the fitskeywords_errortype from the database
        get_sql_col,'SELECT pk, name,codename FROM apogeeqldb.fitskeywords_errortype', kw_pk, kw_name, kw_code, /string
        FITSKW_ERR = {pk:kw_pk, name:kw_name, code:kw_code}

        ; get the latest hmag_standard, snr_standard_goal and version from the database
        get_sql_col,'SELECT hmag_standard, snr_standard_goal, version from apogeeqldb.apogee_snr_goals '+ $
            'order by version DESC limit 1',hmag,snrGoal,version

        if n_elements(hmag) eq 0 then begin
          hmag=12.0
          snrGoal=30.0
          version=0
        endif else begin
          hmag=hmag[0]
          snrGoal=snrGoal[0]
          version=long(version[0])
        endelse
        SNR_GOALS = {hmag:hmag, snr_goals:snrGoal, version:version}

        ; get the latest list of required fits kweywords
        get_sql_col, 'SELECT MAX(version) FROM apogeeqldb.required_fitskeywords', latest_version,/string
        get_sql_col, 'SELECT pk,name,version,datatype,lowval,highval FROM apogeeqldb.required_fitskeywords WHERE version='+$
        latest_version[0], pk,name,version,dtype,lowval,highval,/string
        REQUIRED_FITSKW = {pk:pk, name:name, version:version, dtype:dtype, lowval:lowval, highval:highval}

       endif
   
      rawfiles = file_search(datadir+mjd+'/apRaw-'+frameid+'*.fits')
      get_sql_col,'select e.pk from exposure e where e.exposure_no ='+frameid,exp_pk,/long
      if n_elements(exp_pk) eq 0 then exp_pk='None' else exp_pk=exp_pk[0]

      ;for each read
      for j=0,n_elements(rawfiles)-1 do begin

        filename = rawfiles[j] ;strtrim(vals[0],2)
        basename = file_basename(filename,'.fits')
        !apql.current_filename = basename
        dum = strsplit(basename,'-',/extract)
        readnum = long(dum[2])
              
        ; Get the header information
        head = headfits(filename)
        plate = sxpar(head,'PLATEID')
        exptype = strtrim(strupcase(sxpar(head,'EXPTYPE')),2)  ; object, calibration, flat
        if exptype eq 'OBJECT' and !apql.plugmap.plate ne plate then begin
          msg = string(plate,format='("Wrong plugmap for this plate (",I0,")")',/print)
          print,msg
          msg = string(!apql.plugmap.plate,format='("Saved plugmap is for plate ",I0)',/print)
          print,msg
        endif

        ; Check if this is an Any-Star-Down-Any-Fiber (ASDAF) exposure
        ;  ASDAF if it's object and ra/dec coordinates do match those in plugmap
        if exptype eq 'OBJECT' and PTR_VALID(!apql.plugmap.datastr) then begin
          rahd = sxpar(head,'RA',count=nra)
          dechd = sxpar(head,'DEC',count=ndec)
          rapl = (*!apql.plugmap.datastr).racen
          decpl = (*!apql.plugmap.datastr).deccen
          ; we have valid coordinates
          if nra gt 0 and ndec gt 0 and min(valid_num([rahd,dechd,rapl,decpl])) eq 1 then begin
            dist = sphdist(rahd,dechd,rapl,decpl,/deg)
            if dist gt 2 then exptype='ASDAF'
            if exptype eq 'ASDAF' then print,'This is an ASDAF exposure'
          endif
        endif

        ; Load the first read if necessary
        ;----------------------------------
        if !apql.firstread.frameid ne frameid or readnum eq 1 then begin

          exp_filenames=filename
          ; Load the first read
          if file_test(filename) eq 0 then begin
            msg = string(filename,format='("First read file ",A," NOT FOUND")',/print)
            print,msg
            break
          endif
          ; It takes ~0.05 sec to load these images
          FITS_READ,filename,im0,head0,message=message0,/no_abort
          if message0 ne '' then begin
            msg = string(filename,format='("can not read file",A)',/print)
            print,msg
            break
          endif
          im0 = long(im0)
          ; Stuff it into the system variable structure
          nhead0 = n_elements(head0)
          !apql.firstread.header[0:nhead0-1] = head0
          if PTR_VALID(!apql.firstread.data) then PTR_FREE,!apql.firstread.data
          !apql.firstread.data = PTR_NEW(im0, /NO_COPY)
          !apql.firstread.plate = sxpar(head0,'PLATEID')
          !apql.firstread.frameid = frameid
          timePerRead = sxpar(head0,'INTOFF',count=count)
          if count gt 0 then begin
            !apql.timePerRead = timePerRead[0]/1000.0  ; INTOFF is in msec
          endif else begin
            !apql.timePerRead = 0.0
          endelse
          if size(allstr,/tname) eq 'STRUCT' then HEAP_FREE, allstr, /verbose
          apgundef,allstr  ; initialize allstr

          if size(all_predictstr,/tname) eq 'STRUCT' then HEAP_FREE, all_predictstr, /verbose
          apgundef,all_predictstr  ; initialize all_predictstr

          ; Initialize PREVSTR if this is a new plate, or new EXPTYPE, or new ASDAF
          nprevexp_plate = 0 & prev_exptype='' & new_asdaf=0
          nprevstr = n_elements(prevstr)

          if nprevstr gt 0 then begin
            dum = where(prevstr.plate eq plate,nprevexp_plate)
            prev_exptype = first_el(prevstr.exptype,/last)

            ; Check if this is the same ASDAF as before
            if exptype eq 'ASDAF' then begin
              rahd1 = sxpar(prevstr[nprevstr-1].header,'RA',count=nra1)
              dechd1 = sxpar(prevstr[nprevstr-1].header,'DEC',count=ndec1)
              rahd2 = sxpar(head,'RA',count=nra2)
              dechd2 = sxpar(head,'DEC',count=ndec2)
              ; we have good coordinates
              if min([nra1,ndec1,nra2,ndec2]) gt 0 and min(valid_num([rahd1,dechd1,rahd2,dechd2])) eq 1 then begin
                dist = sphdist(rahd1,dechd1,rahd2,dechd2,/deg)
                if dist gt 0.1 then new_asdaf=1
              endif
            endif ; asdaf
        endif ; nprevstr>0
        
        if nprevexp_plate eq 0 or (prev_exptype ne exptype) or (new_asdaf eq 1) then begin
          if size(prevstr,/tname) eq 'STRUCT' then HEAP_FREE,prevstr, /verbose
            apgundef,prevstr
          endif
        ; NREAD>1
        endif else if long(readnum) gt 1 then begin

          ; Start processing each UTR frame with APQL
          exp_filenames=[exp_filenames,filename]
          APQL,filename,str,allstr,prevstr,predict_str=predict_str, SNR_GOALS=snr_goals, $
            FITSKW_ERR=fitskw_err, REQUIRED_FITSKW=required_fitskw

          ; check if this was a FLAT and if missingFibers are available
          if strpos(strupcase(str.exptype),'FLAT') ge 0 then begin
            reply = string(long(str.frameid),long(str.readnum),str.nmissingfibers, $
                     format='("missingFibers=",A,",",I0,",",I0)',/print)
            if str.nmissingFibers gt 0 and PTR_VALID(str.missingFibers) then begin
              for k=0,str.nmissingFibers-1 do begin
                reply = string(reply,(*str.missingFibers)[k],format='(A,",",I0)',/print)
              endfor
            endif
            print, reply
          endif

          ; Add this to the ALLSTR structure
          ; free up pointers
          tempstr = {filename:'',frameid:'',readnum:'',exptime:0.0,header:strarr(1000),$
                      exptype:'',plate:'',date:'',jd:0.0d0,mjd5:'',snr_goals_version:0,hmag_standard:0.0,$
                      snr_standard_goal:0.0,snr_standard:0.0,delta_snr2_standard:0.0,$
                      logsnr_hmag_coef:fltarr(2),logsnr_hmag_coef_goal:fltarr(2),$
                      snr2_time_coef:fltarr(2),skyvar_meddev_perc:0.0,skyvar_stddev_perc:0.0,$
                      skyvar_contflux:0.0,skyvar_contflux_rate:0.0,skyvar_avglineflux:0.0,$
                      skyvar_avglineflux_rate:0.0,sky_status:0,dither_prevexp_measured:0.0,$
                      dither_prevexp_header:0.0,dither_relative:0.0,dither_status:0, $
                      fitsheader_status:0,exp_finished_status:0,visit_finished_status:0,$
                      medsky:str.medsky,medskylinestr:PTR_NEW(),expected_total_readnum:0.0,$
                      expected_total_readnum_recent:0.0,snr2_time_recent_coef:fltarr(2),wavefit_rms:0.0, $
                      wavelength_status:-1,arraydisplay:{data:bytarr(npixbin,npixbin*nchips),bscale:0.0,bzero:0.0,zscale:[0.0,0.0]},$
                      arraydisplay_sub:REPLICATE({data:PTR_NEW(),bscale:0.0,bzero:0.0,zscale:[0.0,0.0],yrange:[0L,0L]},3),$
                      representative_spectra:REPLICATE({data:PTR_NEW(),bscale:0.0,bzero:0.0,yrange:[0.0,0.0],$
                      objtype:'',fiberid:0,mag:fltarr(3),medsnr:0.0},9)}
          STRUCT_ASSIGN,str,tempstr
          if PTR_VALID(str.medskylinestr) then tempstr.medskylinestr = PTR_NEW(*str.medskylinestr)
          PUSH,allstr,tempstr
          apgundef,tempstr

          ; Erase information in STR that we don't need for DBINSERT
          if PTR_VALID(str.rawimage) then PTR_FREE,str.rawimage
          if PTR_VALID(str.cds_image) then PTR_FREE,str.cds_image
          if PTR_VALID(str.frame) then PTR_FREE,str.frame
          if PTR_VALID(str.medskylinestr) then PTR_FREE,str.medskylinestr
          if PTR_VALID(str.skyfiberstr) then PTR_FREE,str.skyfiberstr
          if PTR_VALID(str.skylinestr) then PTR_FREE,str.skylinestr

          ; Estimate the timePerRead from the first two UTR read that comes in
          ; in case fits keyword is missing
          if long(readnum) eq 2 and !apql.timePerRead eq 0.0 then begin
            date_obs=sxpar(!apql.firstread.header,'DATE-OBS',count=count)
            jd0=date2jd(date_obs[0])
            date_obs=sxpar(str.header,'DATE-OBS',count=count)
            jd1=date2jd(date_obs[0])
            !apql.timePerRead = (jd1-jd0)*(24.0d*60.0d*60.0d)  ; in seconds
          endif

          ; Update the database with IDL-SQL
          ; apogeedb.quicklook and apogeedb.quicklook60 tables

          print,'Starting DataBase INSERT'
          ; for now save the structure to an IDL save file and pass the filename
          ; to be consistent with apql_wrapper.pro
          savefile = filepath('apql_db_'+str.frameid+'.sav',/TMP)

          if n_elements(predict_str) gt 0 then begin
            ; PREDICT_STR has the prediction information every 6th read
            ; add filename to predict_str and add to all_predictstr
            ext_predict_str = REPLICATE(CREATE_STRUCT(predict_str[0],'filename',''),n_elements(predict_str))
            STRUCT_ASSIGN, predict_str, ext_predict_str  
            ext_predict_str.filename = file_basename(str.filename)
            PUSH,all_predictstr,ext_predict_str

            save, str, predict_str, fitskw_err, file=savefile 
          endif else begin
            save, str, fitskw_err, file=savefile
          endelse
 
          fitskw_version = required_fitskw.version[0]

          if NOT keyword_set(no_dbinsert) then begin 
            apql_dbinsert, savefile=savefile, fitskw_version=fitskw_version, exp_pk=exp_pk
            ;print,'DB INSERT'
          endif

          ; Free up the STR structure
          HEAP_FREE,str, /verbose
          apgundef,str

        endif ; readnum>1 

      endfor
  
      ; write allstr data to an output file
      if keyword_set(outfile) then begin
        apql_writeoutputfile, outfile, allstr, all_predictstr, exp_pk=exp_pk, fitskw_version=fitskw_version
      endif
      ;Update PREVSTR
      if n_elements(allstr) gt 0 then $
        PUSH,prevstr,first_el(allstr,/last)
  endfor ;exposure loop 
end
