;+
;
; APQL_WRAPPER
;
; This is the wrapper to the APOGEE quicklook program.
; It sets up the socket communication with the apogeeql actor and waits for 
; commands/messages to pass to the quicklook program.
; It also returns messages from the QuickLook software back to the actor.
;
; The APOGEE QUICKLOOK ACTOR spawns an IDL session and calls this 
; routine which opens the socket and sends a message back to the actor.
;
; This program is normally started by the APOGEE quicklook actor as:
; 
; idl -quiet -e apql_wrapper -args localhost 10038 data_dir=/data/apogee spectro_dir=/data/apogee/spectro 
;
; INPUTS:
;  apqlhost: name of the host to which the APOGE QuickLook ACTOR socket is connected
;  apqlport: port to which the APOGE QuickLook ACTOR socket is connected on the remote machine
;
; OUTPUTS:
;  Passes messages/info back and forth through the sockets
;  the three arrays.
;
; USAGE:
;
;-
;
;pro apql_wrapper, apqlhost, apqlport, data_dir=data_dir, spectro_dir=spectro_dir, no_dbinsert=no_dbinsert
pro apql_wrapper, apqlhost, apqlport, obs

   args=COMMAND_LINE_ARGS()
   apqrhost=args[0]
   apqrport=args[1]
   obs=args[2]

   if obs eq 'APO' then apsetver,vers='quicklook',telescope='apo25m'
   if obs eq 'LCO' then apsetver,vers='quicklook',telescope='lco25m'

   RESOLVE_ROUTINE,'apwavecal_chip'

   ; define the IDL environment variables used by idlsql to point to the desired database
   ; this used to be called APQL_DEFSYS
   SDSS_DB_PARAMS,obs=obs

   npix = 2048L
   nchips = 3L
   ditherPos = 0.0
   namedDitherPos = '?'

   ; to run the IDL_BRIDGE on a machine with no display attached, we need to
   ; fake a display with Xvfb
   ; SPAWN,'Xvfb :1 -screen 0 1600x1200x24 &'
   setenv,'DISPLAY=:1'

   ; Get APOGEE directories if the environment variables exists
   ;data_dir = APGETDIR('APQLDATA_DIR',/exists,error=direrr)
   ;spectro_dir = APGETDIR('APQLSPECTRO_DIR',/exists,error=direrr)
   ;archive_dir = APGETDIR('APQLARCHIVE_DIR',/exists,error=direrr)

   ; look for command line arguments first
   ;args=command_line_args(count=count)
   ;if count gt 0 then apqlhost=args[0]
   ;if count gt 1 then apqlport=fix(args[1])
   ;; overwrite the environment variables if they were provided on the command line
   ;; we can't guarantee the order in which the parameters were provided
   ;if count gt 2 then begin
   ;   if strpos(args[2],'data_dir') ge 0 then begin
   ;      data_dir = strmid(args[2],9)
;
;      endif else if strpos(args[2],'spectro_dir') ge 0 then begin
;         spectro_dir = strmid(args[2],12)
;      endif else if strpos(args[2],'archive_dir') ge 0 then begin
;         spectro_dir = strmid(args[2],12)
;      endif
;   endif
;   if count gt 3 then begin
;      if strpos(args[3],'data_dir') ge 0 then begin
;         data_dir = strmid(args[3],9)
;      endif else if strpos(args[3],'spectro_dir') ge 0 then begin
;         spectro_dir = strmid(args[3],12)
;      endif else if strpos(args[3],'archive_dir') ge 0 then begin
;         spectro_dir = strmid(args[3],12)
;      endif
;   endif
;   if count gt 4 then begin
;      if strpos(args[4],'data_dir') ge 0 then begin
;         data_dir = strmid(args[4],9)
;      endif else if strpos(args[4],'spectro_dir') ge 0 then begin
;         spectro_dir = strmid(args[4],12)
;      endif else if strpos(args[4],'archive_dir') ge 0 then begin
;         spectro_dir = strmid(args[4],12)
;      endif
;   endif
;
   ; if no arguments provided, use the default hosts and ports
   if n_elements(apqlhost) eq 0 then apqlhost='127.0.0.1' ; hubhost='10.25.1.1'
   if n_elements(apqlport) eq 0 then apqlport=10038

   ; we need data_dir and spectro_dir
;   if strlen(data_dir) eq 0 then begin
;      msg = 'Missing data_dir -> aborted'
;      print,msg
;      return
;   endif
;   if strlen(spectro_dir) eq 0 then begin
;      msg = 'Missing spectro_dir -> aborted'
;      print,msg
;      return
;   endif
;   if strlen(archive_dir) eq 0 then begin
;      msg = 'Missing archive_dir -> aborted'
;      print,msg
;      return
;   endif

   ; Initialize the message log file
   jd = systime(/julian)
   caldat,jd,month,day,year,hour,minute,second
   fm2='(I02)' & fm4='(I04)'
   messagelogfile = '/data-ql/logs/toapql_messages.'+string(month,fo=fm2)+string(day,fo=fm2)+$
       string(year,fo=fm4)+string(hour,fo=fm2)+string(minute,fo=fm2)+string(second,fo=fm2)+'.log'

   chiptag = ['a','b','c']
   ;linelist_dir = spectro_dir+'/lib/linelists/'
   ;datadir = data_dir+'/raw/'
   ;psfdir = spectro_dir+'/cal/psf/'
   ;bpmdir = spectro_dir+'/cal/bpm/'
   dirs=getdir()
   linelist_dir = dirs.libdir+'/skylines/'
   datadir = apogee_filename('Raw',num=0,read=1,/dir)
   psfdir = apogee_filename('PSF',num=0,chip='a',/dir)
   bpmdir = apogee_filename('BPM',num=0,chip='a',/dir)

   ; Get Psf file
   tinfo = APFILEINFO(FILE_SEARCH(psfdir+'*PSF-a-*.fits'),/silent)
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
   apgundef,str,allstr,prevstr


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

   ; connect to the apogeeql socket
   ;---------------------------
   print,'Connecting from apql_wrapper to apogeeql through: ',apqlhost,apqlport
   SOCKET, apqlactor_lun, apqlhost, apqlport, /get_lun, error=connect_error
   if connect_error ne 0 then begin
      print,'Error connecting to apogeeql ('+strtrim(apqlhost,2)+', '+strtrim(apqlport,2)+')'
      print, !ERROR_STATE.MSG
      return
   endif

   printf,apqlactor_lun, 'Hello from apql_wrapper'
   imagedir = ''

   ; create a new thread for inserting info in the database
   ; we assume for now that only one thread will be needed
   ;----------------------------------------------------
   dbBridge = obj_new('IDL_IDLBridge', OUTPUT='', CALLBACK='apql_dbinsert_callback')
   ; execute the startup file if we had one
   dbBridge->EXECUTE, '@' + PREF_GET('IDL_STARTUP')

   ;=================
   ; LARGE LOOP
   ;=================
   flag = 0
   num_cmds = 0L
   WHILE (flag eq 0) do begin

      ; clear the pointers to free heap memory
      if size(str,/tname) eq 'STRUCT' then HEAP_FREE, str, /verbose
      apgundef,str

     ; Get information from apogeeql actor through the socket
     ;------------------------------------------------------
     result = FILE_POLL_INPUT(apqlactor_lun)
     astring=''
     ; a True results from FILE_POLL_INPUT could indicate the presence of data but also an error. 
     ; skip processing to EOD if we encounter an error reading
     ON_IOERROR, EOD
     READF,apqlactor_lun, astring

     GOTO, NO_IOERROR

     EOD:
        print,'Error reading apqlactor_lun -> testing the socket connection'
        ; send a message to the actor to make sure it's still alive
        ON_IOERROR, NULL
        MESSAGE, /RESET
        printf, apqlactor_lun, 'PONG'

        if strpos(!ERROR_STATE.SYS_MSG,'Broken') ge 0 then begin
            ; the socket was closed at the other end -> quit the program
            print,'socket is close -> EXITING apql_wrapper'
            OBJ_DESTROY, dbBridge
            return
        endif

     NO_IOERROR:
     ; data was succesfully read from apqlactor_lun (cancel ON_IOERROR)
     ON_IOERROR, NULL
     print,strtrim(num_cmds+1,2),' ql Message = ',astring
     ; Write message to message log file
     WRITELINE,messagelogfile,astring,/append


     ; checks if a shutdown request was sent
     keywords = strsplit(astring,/extract)
     ; print,'Wrapper received: ',keywords
     pos = strpos(keywords[0],'=')
     if pos ge 0 then thiscmd=strupcase(strmid(keywords[0],0,pos)) else thiscmd=strupcase(keywords[0])


     CASE thiscmd OF
        'STARTING': begin
           ; first message from the apogeeql actor
           printf, apqlactor_lun, 'Listening ...'
           snrAxisRange = [0.0,35.0]
           printf, apqlactor_lun, snrAxisRange,format='("snrAxisRange=",2(F0.2,","))'
           printf, apqlactor_lun, 'STARTED'
           ; reset the missingFibers info for STUI - will get updated when a FLAT is taken
           printf, apqlactor_lun, "missingFibers=0,0,nan"
           end
        'PING': begin
           ; aliveness test
           printf, apqlactor_lun, 'PONG'
           end
        'QUIT': begin
           printf, apqlactor_lun, 'Quitting the IDL QuickLook handler'
           free_lun, apqlactor_lun
           OBJ_DESTROY, dbBridge
           return
           end
        'DITHERPOSITION': begin
           ; a new dither position was detected
           ; expecting:   ditherPosition=13.9977,A
           values=strsplit(strmid(astring,15),',',/extract)
           ditherPos = float(values[0])
           namedDitherPos = values[1]
           print,'0 - namedDitherPos=',namedDitherPos
           end
        'PLUGMAPINFO': begin
           ; a new plugmapfile was detected
           ; expecting the plateId, fscan_mjd, fscan_id, filename
           ; plugMapInfo=4934,55704,1,/data-ql/plugmaps/plPlugMapA-4934-55704-01.par

           ; the following command could be used to extract the information directly from the database
           ; instead of writing and reading a temporary file but the yanny_read only deals with text files.
           ; get_sql_col,'select file from pl_plugmap_m as p1, plugging as p2, plate as p3, '+ $
           ;  'cartridge as c1 where c1.number='+strtrim(string(cartridge),2)+$
           ;  ' and p2.cartridge_pk=c1.pk and p2.plate_pk=p3.pk and p1.plugging_pk=p2.pk and '+ $
           ;  'p3.plate_id='+strtrim(string(plate),2)+' and p1.pointing_name='+"'"+pointing+"'",plugmap,/string

           values=strsplit(strmid(astring,12),',',/extract)
           plate = values[0]
           fscan_mjd = values[1]
           fscan_id = values[2]
           filename = values[3]

           if !apql.plugmap.plate ne plate or !apql.plugmap.mjd ne fscan_mjd then begin
              plugfile = file_search(filename,count=nplugfile)
              if nplugfile eq 0 then begin
                 msg = 'NO plugmap file in '+filename
                 print,msg
                 printf, apqlactor_lun, 'Error: '+msg
                 ; return
              endif else begin
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
              endelse

              ; with a new plate load the fitskeywords_errotype from the database
              get_sql_col,'SELECT pk, name,codename FROM apogeeqldb.fitskeywords_errortype', kw_pk, kw_name, kw_code, /string
              FITSKW_ERR = {pk:kw_pk, name:kw_name, code:kw_code}

              ; get the latest hmag_standard, snr_standard_goal and version from the database
              get_sql_col,'SELECT hmag_standard, snr_standard_goal, version from apogeeqldb.apogee_snr_goals '+ $
                 'order by version DESC limit 1',hmag,snrGoal,version

              if n_elements(hmag) eq 0 then begin
                 ;print,'Error finding rows in apogee_snr_goals'
                 ;print,'Using default of hmag=12.0, snr_goal=30.0, version=0 '
                 hmag=12.0
                 snrGoal=30.0
                 version=0
              endif else begin
                 hmag=hmag[0]
                 snrGoal=snrGoal[0]
                 version=long(version[0])
              endelse
              SNR_GOALS = {hmag:hmag, snr_goals:snrGoal, version:version}
              ; send these values to the actor to send to STUI
              printf, apqlactor_lun, snrGoal,hmag,format='("snrGoal=",F0.2,",",F0.2)'

              ; get the latest list of required fits kweywords
              get_sql_col, 'SELECT MAX(version) FROM apogeeqldb.required_fitskeywords', latest_version,/string
              get_sql_col, 'SELECT pk,name,version,datatype,lowval,highval FROM apogeeqldb.required_fitskeywords WHERE version='+$
                 latest_version[0], pk,name,version,dtype,lowval,highval,/string
              REQUIRED_FITSKW = {pk:pk, name:name, version:version, dtype:dtype, lowval:lowval, highval:highval}

              ; reset the missingFibers info for STUI - will get updated when a FLAT is taken
              printf, apqlactor_lun, "missingFibers=0,0,nan"

              ; get the list of previous OBJECT exposures with this plate if any
              ; if we are here, we assume the loaded plate is for APOGEE (tested in apogeeql)
              plate_prev_exp = apql_get_prevexp(plate, snrgoal=snrGoal, count=count,/silent)

              ;   tempstr = {plateId:'',expNum:0, expName:'', exptime:0.0, numReads:0, snrGoal:float(snrgoal), $
              ;   ditherPos:0.0,snr:0.0,netExpTime:0.0,netSnr:0.0, expType:'Object', namedDitherPos:''}

              if count gt 0 then begin
                  ; send the list of previous exposures for this plate back to the actor
                  for i=0,count-1 do begin
                      reply = string(plate_prev_exp[i].plateid, plate_prev_exp[i].expNum,  plate_prev_exp[i].expName, $
                          plate_prev_exp[i].exptime, plate_prev_exp[i].numReads, plate_prev_exp[i].snrgoal, $
                          plate_prev_exp[i].ditherPos, plate_prev_exp[i].snr, plate_prev_exp[i].netExpTime, $
                          plate_prev_exp[i].netSnr, plate_prev_exp[i].exptype, plate_prev_exp[i].namedDitherPos, $
                          format='("exposureData=",I0,",",I0,",",A,",",F0.2,",",A,",",5(F0.2,","),A,",",A1)',/print)
                      printf, apqlactor_lun, reply
                  endfor
              endif

           endif
           end

        'UTR': begin
           if strpos(astring,'UTR=DONE') ge 0 then begin

              ; Update PREVSTR
              if n_elements(allstr) gt 0 then $
                PUSH,prevstr,first_el(allstr,/last)


              ; the UTR exposure just completed (or got aborted)
              ; the quickreduction is now handled in its own thread from the python actor
              if n_elements(allstr) gt 0 then begin
                 ; send the new exposureData to the actor
                 ;   use the prevstr structure
                 netExpTime = total(prevstr.exptime)
                 netSnrH12 = sqrt( total(prevstr.snr_standard^2) )  ; add in quadrature

                 p = (n_elements(allstr) -1)>0
                 if allstr[p].exptype ne 'OBJECT' then begin
                     ; if not an object simply reset the STUI table
                     numExpForThisPlate=1   
                     reply = string(allstr[p].plate, numExpForThisPlate, strtrim((allstr[p].frameid),2), allstr[p].exptime, $
                         allstr[p].readnum, allstr[p].snr_standard_goal, allstr[p].dither_prevexp_measured, $
                         allstr[p].snr_standard, netExpTime, netSnrH12, allstr[p].exptype, namedDitherPos, $
                         format='("exposureData=",I0,",",I0,",",A,",",F0.2,",",A,",",5(F0.2,","),A,",",A1)',/print)
                     printf, apqlactor_lun, reply
                 endif else begin
                     ; if we have an object exposure, add to the table loaded when the plate was loaded
                     if n_elements(plate_prev_exp) eq 0 then begin
                         ; this can only happen if we're here before we got a plugmapInfo
                         plate_prev_exp = apql_get_prevexp(allstr[p].plate, count=count,/silent)
                     endif
                     count = n_elements(plate_prev_exp)
                     if count eq 1 then begin
                         ; check if we have a valid entry or an empty structure
                         if long(plate_prev_exp[0].plateId) ne long(allstr[p].plate) then begin
                             ; either an empty list or a new plate
                             pos = 0
                         endif else begin
                             ; one valid previous exposure
                             plate_prev_exp = [plate_prev_exp, plate_prev_exp[0]]
                             pos = n_elements(plate_prev_exp)-1
                         endelse
                     endif else begin
                         if long(plate_prev_exp[0].plateId) ne long(allstr[p].plate) then begin
                             ; a new plate -> start with an empty structure
                             plate_prev_exp = plate_prev_exp[0]
                             pos = 0
                         endif else begin
                             ; add a new row to the list
                             plate_prev_exp = [plate_prev_exp, plate_prev_exp[0]]
                             pos = n_elements(plate_prev_exp)-1
                         endelse
                     endelse
                     help,plate_prev_exp,pos
                     plate_prev_exp[pos].plateId = allstr[p].plate
                     plate_prev_exp[pos].expNum  = pos+1
                     plate_prev_exp[pos].expName = strtrim((allstr[p].frameid),2)
                     plate_prev_exp[pos].exptime = allstr[p].exptime
                     plate_prev_exp[pos].numReads = long(allstr[p].readnum)
                     plate_prev_exp[pos].ditherPos = float(allstr[p].dither_prevexp_measured)
                     plate_prev_exp[pos].snr = float(allstr[p].snr_standard)
                     plate_prev_exp[pos].netExpTime = total(plate_prev_exp.exptime) 
                     plate_prev_exp[pos].netSnr = sqrt(total(plate_prev_exp.snr^2.0)) 
                     ; STUI wants expType capitalized
                     exptype = strtrim(allstr[p].exptype,2)
                     plate_prev_exp[pos].expType = strupcase(strmid(exptype,0,1))+strlowcase(strmid(exptype,1))
                     plate_prev_exp[pos].namedDitherPos = namedDitherPos 

                     ; because of the way STUI displays these we need to print everyone of them
                     for i=0,pos do begin
                         reply = string(plate_prev_exp[i].plateid, plate_prev_exp[i].expNum,  plate_prev_exp[i].expName, $
                              plate_prev_exp[i].exptime, plate_prev_exp[i].numReads, plate_prev_exp[i].snrgoal, $
                              plate_prev_exp[i].ditherPos, plate_prev_exp[i].snr, plate_prev_exp[i].netExpTime, $
                              plate_prev_exp[i].netSnr, plate_prev_exp[i].exptype, plate_prev_exp[i].namedDitherPos, $
                              format='("exposureData=",I0,",",I0,",",A,",",F0.2,",",A,",",5(F0.2,","),A,",",A1)',/print)
                         printf, apqlactor_lun, reply
                     endfor
                 endelse
              endif


           ; Normal read
           endif else if strpos(astring,'UTR=') ge 0 then begin

              ; we got a new filename for the UTR exposure (full path)
              ; expecting /data-ql/MJD5/apRaw-DDDDXXXX-RRR.fits,2374,2,50
              ; where MJD5 is the MJD of the day
              ; where DDDD is the day number starting at 0 from Jan.1 2011
              ; where XXXX is the exposure number of that day
              ; RRR is the UTR read number for this exposure
              ; message is UTR=absolute_filename, exposure_pk, readnum, nreadsCommanded
              vals = strsplit(strmid(astring,4),',',/extract)
              filename = strtrim(vals[0],2)
              exposure_pk = strtrim(vals[1],2)
              readnum = long(vals[2])
              numReadsCommanded = long(vals[3])
              ;filename = strmid(astring,4)
              ; get the filename without the full path
              basename = file_basename(filename,'.fits')
              mjd5 = file_basename(file_dirname(filename))
              !apql.current_filename = basename
              dum = strsplit(basename,'-',/extract)
              frameid = strtrim(dum[1],2)
              readnum = long(dum[2])

              ; Get the header information
              head = headfits(filename)
              plate = sxpar(head,'PLATEID')
              exptype = strtrim(strupcase(sxpar(head,'EXPTYPE')),2)  ; object, calibration, flat
              if exptype eq 'OBJECT' and !apql.plugmap.plate ne plate then begin
                 msg = string(plate,format='("Wrong plugmap for this plate (",I0,")")',/print)
                 printf, apqlactor_lun, 'Error: '+msg
                 print,msg
                 msg = string(!apql.plugmap.plate,format='("Saved plugmap is for plate ",I0)',/print)
                 printf, apqlactor_lun, 'Error: '+msg
                 print,msg
                 ; return
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

              ; print,'EXPTYPE=',exptype

              ; Load the first read if necessary
              ;----------------------------------
              if !apql.firstread.frameid ne frameid or readnum eq 1 then begin

                   exp_filenames=filename
                   ; Load the first read
                   if file_test(filename) eq 0 then begin
                     msg = string(filename,format='("First read file ",A," NOT FOUND")',/print)
                     printf, apqlactor_lun, 'Error: '+msg
                     print,msg
                     break
                     ; return
                   endif
                   ; It takes ~0.05 sec to load these images
                   FITS_READ,filename,im0,head0,message=message0,/no_abort
                   if message0 ne '' then begin
                     msg = string(filename,format='("can not read file",A)',/print)
                     printf, apqlactor_lun, 'Error: '+msg
                     print,msg
                     break
                     ; return
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

                   ; load list of required fits keyword errors if necessary
                   if n_elements(fitskw_err) eq 0 then begin
                       get_sql_col,'SELECT pk, name,codename FROM apogeeqldb.fitskeywords_errortype', $
                           kw_pk, kw_name, kw_code, /string
                       FITSKW_ERR = {pk:kw_pk, name:kw_name, code:kw_code}
                   endif

                   ; load list of required fits keywords if necessary
                   if n_elements(required_fitskw) eq 0 then begin
                     ; get the latest list of required fits kweywords
                     get_sql_col, 'SELECT MAX(version) FROM apogeeqldb.required_fitskeywords', latest_version,/string
                     get_sql_col, 'SELECT pk,name,version,datatype,lowval,highval FROM apogeeqldb.required_fitskeywords WHERE version='+$
                       latest_version[0], pk,name,version,dtype,lowval,highval,/string
                     REQUIRED_FITSKW = {pk:pk, name:name, version:version, dtype:dtype, lowval:lowval, highval:highval}
                   endif

              ; NREAD>1
              endif else if long(readnum) gt 1 then begin

                   ; Start processing each UTR frame with APQL
                   exp_filenames=[exp_filenames,filename]
                   APQL,filename,str,allstr,obs=obs,prevstr,predict_str=predict_str, SNR_GOALS=snr_goals, $
                      FITSKW_ERR=fitskw_err, REQUIRED_FITSKW=required_fitskw

                   ; check if this was a FLAT and if missingFibers are available
                   if strpos(strupcase(str.exptype),'FLAT') ge 0 then begin
                       reply = string(long(str.frameid),long(str.readnum),str.nmissingfibers, $
                                     format='("missingFibers=",A,",",I0,",",I0)',/print)
                       if str.nmissingFibers gt 0 and PTR_VALID(str.missingFibers) then begin
                           for i=0,str.nmissingFibers-1 do begin
                               reply = string(reply,(*str.missingFibers)[i],format='(A,",",I0)',/print)
                           endfor
                       endif
                      printf, apqlactor_lun, reply
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
                             fits_header_status:0,exp_finished_status:0,visit_finished_status:0,$
                             medsky:str.medsky,medskylinestr:PTR_NEW()}
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
                   status = dbBridge->Status(error=errmsg) 
                   ;if status eq 0 then begin
                   if fix(status) ne 1 then begin
                       ; print,'Starting DataBase INSERT'
                       ; for now save the structure to an IDL save file and pass the filename
                       ; to the routine (can't pass pointers to threads or shared memory)
                       savefile = filepath('apql_db_'+str.frameid+'.sav',/TMP)
                       if n_elements(predict_str) gt 0 then begin
                           ; PREDICT_STR has the prediction information every 6th read
                           save, str, predict_str, fitskw_err, file=savefile 
                       endif else begin
                         save, str, fitskw_err, file=savefile
                       endelse
                       dbBridge->SetVar, 'savefile', savefile
                       dbBridge->SetVar, 'exp_pk', exposure_pk
                       dbBridge->SetVar, 'obs', obs
                       dbBridge->SetVar, 'fitskw_version', required_fitskw.version[0]
                       if NOT keyword_set(no_dbinsert) then $
                           dbBridge->Execute,'apql_dbinsert, obs=obs, savefile=savefile, fitskw_version=fitskw_version, exp_pk=exp_pk', /NOWAIT
                   endif else begin
                      print,'******---------*******'
                      print,'Error: dbBridge status=',status,'  errmsg='+errmsg
                      print,'******---------*******'
                   endelse

                   ; Send information for STUI back to the actor
                   if n_elements(str) gt 0 then begin
                      ; form the string to return to the actor
                      ; change the bits from 1=valid to 0=valid and 1=bad
                      ; The status values in STR are 0-good, 1-bad
                      validBits = ishft((str.wavelength_status eq 0),3)
                      validBits += ishft((str.sky_status eq 0),2)
                      validBits += ishft((str.dither_status eq 0),1)
                      validBits += (str.fitsheader_status eq 0)
                      ; the coeffs are per seconds and not per readnum
                      reply = string(long(str.frameid), long(str.readnum), str.snr_standard, $
                                     str.snr2_time_coef[0], str.snr2_time_coef[1] * !apql.timePerRead, $
                                     str.snr2_time_recent_coef[0], str.snr2_time_recent_coef[1] * !apql.timePerRead, validBits,  $
                                     str.dither_prevexp_measured, str.dither_prevexp_header,  $
                                     avg(str.waverange_diff[*,0]), $
                                     str.expected_total_readnum * !apql.timePerRead, $
                                     str.expected_total_readnum, $
                                     numReadsCommanded, str.delta_snr2_standard,str.expType, namedDitherPos, $
                                     format='("utrData=",A,",",I0,",",F0.2,",",4(F0.4,","),"0x",z02,",",'+$
                                     '4(F0.2,","),F0.2,",",I0,",",F0.2,",",A,",",A1)',/print)
                      printf, apqlactor_lun, reply

                      ; send the URL (will have to figure out where they're coming from 
                      printf,apqlactor_lun,"rootURL=www.apo.nmsu.edu/apogeeql"
                      printf,apqlactor_lun,"exposureURL=exposures"
                      printf,apqlactor_lun,"fitsURL=fits"
                      printf,apqlactor_lun,"dithersURL=dither"
                      printf,apqlactor_lun,"skyURL=sky"
                      printf,apqlactor_lun,"waveURL=wave"
                   endif

                   ; send the predicted exposures (only calculated every 6 UTR)
                   count=n_elements(predict_str)
                   if count gt 0 then begin
                       ; form the predictedExposure message
                       format='("predictedExposure=",I0,",",I0,",",A,",",F0.2,",",I0,",",F0.2,",Object,",I0,",",A1)'
                       for i=0,count-1 do begin
                           if predict_str[i].ditherexposure eq 1 then begin
                               namedPos='B' 
                           endif else if predict_str[i].ditherexposure eq 2 then begin
                               namedPos='A' 
                           endif else begin
                               namedPos='?' 
                           endelse
                           reply=string(plate,i+1,predict_str[i].frameid,predict_str[i].exptime,$
                               predict_str[i].nreads, predict_str[i].snr_standard, $
                               predict_str[i].ditherexposure,namedPos,$
                               format=format)
                           printf, apqlactor_lun, reply
                       endfor
                   endif

                   ; Free up the STR structure
                   HEAP_FREE,str, /verbose
                   apgundef,str

               endif ; readnum>1 
           endif ; new UTR filename
           end

        else: begin
           ; echo the input string back to the sender
           nchar = strlen(astring)
           printf, apqlactor_lun, 'Received '+strtrim(string(nchar),2)+' characters'
           printf, apqlactor_lun, astring
        end

     ENDCASE


     ; Increment counter
     num_cmds++

   ENDWHILE  ; large loop

END
