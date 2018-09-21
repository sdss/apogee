;+
;
; APQR_WRAPPER
;
; This is the wrapper to the APOGEE quickreduce program.
; It sets up the socket communication with the apogeeql actor and waits for 
; commands/messages to pass to the apquickred program.
; It also returns messages from the apquickred software back to the actor.
;
; The APOGEE QUICKLOOK ACTOR spawns an IDL session and calls this 
; routine which opens the socket and sends a message back to the actor.
;
; This program is normally started by the APOGEE quicklook actor as:
; 
; idl -quiet -e apqr_wrapper -args localhost 10039
;
; INPUTS:
;  apqrhost: name of the host to which the APOGE QuickLook ACTOR socket is connected
;  apqrport: port to which the APOGE QuickLook ACTOR socket is connected on the remote machine
;
; OUTPUTS:
;  Passes messages/info back and forth through the sockets
;  the three arrays.
;
; USAGE:
;
;-
;
;===========================================================================================
pro qrBridge_cb, status, error, qrBridge, apqractor_lun
   ; return error messages back to the actor
   print,'apqractor_lun=',apqractor_lun

   if status eq 2 then begin
      print,  'APQUICKRED Completed normally'
      printf, apqractor_lun, 'APQUICKRED Completed normally'
   endif else if status eq 3 then begin
      print, 'APQUICKRED Error: ',error
      printf, apqractor_lun, 'APQUICKRED Error: ',error
      ; do a RETALL to get back to the main in case of an error
      qrBridge->EXECUTE, 'RETALL'
   endif else if status eq 4 then begin
      print, 'APQUICKRED Aborted: ',error
      printf, apqractor_lun, 'APQUICKRED Aborted: ',error
      ; do a RESET_SESSION to start fresh (will read the idl_startup)
      qrBridge->EXECUTE, '.RESET_SESSION'
   endif

end

;===========================================================================================
;
pro apqr_wrapper, apqrhost, apqrport, data_dir=data_dir, spectro_dir=spectro_dir


   RESOLVE_ROUTINE,'apwavecal_chip'

   MAX_NUM_BRIDGES = 5
   npix = 2048L
   nchips = 3L

   ; define the IDL environment variables used by idlsql to point to the desired database
   ; this used to be called APQL_DEFSYS
   SDSS_DB_PARAMS

   ; to run the IDL_BRIDGE on a machine with no display attached, we need to
   ; fake a display with Xvfb
   ; SPAWN,'Xvfb :1 -screen 0 800x640x8 &'
   setenv,'DISPLAY=:1'

   ; Get APOGEE directories if the environment variables exists
   data_dir = APGETDIR('APQLDATA_DIR',/exists,error=direrr)
   spectro_dir = APGETDIR('APQLSPECTRO_DIR',/exists,error=direrr)
   archive_dir = APGETDIR('APQLARCHIVE_DIR',/exists,error=direrr)
   quickred_dir = APGETDIR('APQLQUICKRED_DIR',/exists,error=direrr)

   ; look for command line arguments first
   args=command_line_args(count=count)
   if count gt 0 then apqrhost=args[0]
   if count gt 1 then apqrport=fix(args[1])
   ; overwrite the environment variables if they were provided on the command line
   ; we can't guarantee the order in which the parameters were provided
   if count gt 2 then begin
      if strpos(args[2],'data_dir') ge 0 then begin
         data_dir = strmid(args[2],9)

      endif else if strpos(args[2],'spectro_dir') ge 0 then begin
         spectro_dir = strmid(args[2],12)
      endif else if strpos(args[2],'archive_dir') ge 0 then begin
         spectro_dir = strmid(args[2],12)
      endif
   endif
   if count gt 3 then begin
      if strpos(args[3],'data_dir') ge 0 then begin
         data_dir = strmid(args[3],9)
      endif else if strpos(args[3],'spectro_dir') ge 0 then begin
         spectro_dir = strmid(args[3],12)
      endif else if strpos(args[3],'archive_dir') ge 0 then begin
         spectro_dir = strmid(args[3],12)
      endif
   endif
   if count gt 4 then begin
      if strpos(args[4],'data_dir') ge 0 then begin
         data_dir = strmid(args[4],9)
      endif else if strpos(args[4],'spectro_dir') ge 0 then begin
         spectro_dir = strmid(args[4],12)
      endif else if strpos(args[4],'archive_dir') ge 0 then begin
         spectro_dir = strmid(args[4],12)
      endif
   endif

   ; if no arguments provided, use the default hosts and ports
   if n_elements(apqrhost) eq 0 then apqrhost='127.0.0.1' ; hubhost='10.25.1.1'
   if n_elements(apqrport) eq 0 then apqrport=10039

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
   if strlen(archive_dir) eq 0 then begin
      msg = 'Missing archive_dir -> aborted'
      print,msg
      return
   endif

   ;datadir = data_dir+'/raw/'
   psfdir = spectro_dir+'/cal/psf/'
   bpmdir = spectro_dir+'/cal/bpm/'

   ; Initialize the message log file
   jd = systime(/julian)
   caldat,jd,month,day,year,hour,minute,second
   fm2='(I02)' & fm4='(I04)'
   messagelogfile = '/data-ql/logs/toapqr_messages.'+string(month,fo=fm2)+string(day,fo=fm2)+$
       string(year,fo=fm4)+string(hour,fo=fm2)+string(minute,fo=fm2)+string(second,fo=fm2)+'.log'


   ; Get Psf file
   chiptag = ['a','b','c']
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
      ;psffiles = psfdir+'apPSF-'+chiptag+'-'+tinfogd[si[0]].suffix+'.fits'
      psfid = psfdir+tinfogd[si[0]].suffix
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

   ; connect to the apogeeql socket
   ;---------------------------
   print,'Connecting from apqr_wrapper to apogeeql through: ',apqrhost,apqrport
   SOCKET, apqractor_lun, apqrhost, apqrport, /get_lun, error=connect_error
   if connect_error ne 0 then begin
      print,'Error connecting to apogeeqr ('+strtrim(apqrhost,2)+', '+strtrim(apqrport,2)+')'
      print, !ERROR_STATE.MSG
      return
   endif

   printf,apqractor_lun, 'Hello from apqr_wrapper'

   ; create a single new thread (to start with) for the quick reduction program
   ;----------------------------------------------------
   qrBridge = obj_new('IDL_IDLBridge', OUTPUT='', CALLBACK='qrBridge_cb', USERDATA=apqractor_lun)
   ; execute the startup file if we had one
   qrBridge->EXECUTE, '@' + PREF_GET('IDL_STARTUP')

   ;=================
   ; LARGE LOOP
   ;=================
   flag = 0
   num_cmds = 0
   WHILE (flag eq 0) do begin

     ; Get information from apogeeql actor through the socket
     ;------------------------------------------------------
     result = FILE_POLL_INPUT(apqractor_lun)
     astring=''

     ; a True results from FILE_POLL_INPUT could indicate the presence of data but also an error. 
     ; skip processing to EOD if we encounter an error reading
     ON_IOERROR, EOD
     READF,apqractor_lun, astring

     GOTO, NO_IOERROR

     EOD:
        print,'Error reading apqractor_lun -> testing the socket connection'
        ; send a message to the actor to make sure it's still alive
        ON_IOERROR, NULL
        MESSAGE, /RESET
        printf, apqractor_lun, 'PONG'

        if strpos(!ERROR_STATE.SYS_MSG,'Broken') ge 0 then begin
            ; the socket was closed at the other end -> quit the program
            print,'socket to python actor is close -> EXITING apqr_wrapper'
            OBJ_DESTROY, qrBridge
            return
        endif

     NO_IOERROR:
     ; data was succesfully read from apqractor_lun (cancel the ON_IOERROR)
     ON_IOERROR, NULL
     print,strtrim(num_cmds+1,2),' qr Message = ',astring
     WRITELINE,messagelogfile,astring,/append


     ; checks if a shutdown request was sent
     keywords = strsplit(astring,/extract)
     ; print,'Wrapper received: ',keywords
     pos = strpos(keywords[0],'=')
     if pos ge 0 then thiscmd=strupcase(strmid(keywords[0],0,pos)) else thiscmd=strupcase(keywords[0])

     CASE thiscmd OF
        'STARTING': begin
           ; first message from the apogeeql actor
           printf, apqractor_lun, 'Listening ...'
           printf, apqractor_lun, 'STARTED'
           end
        'PING': begin
           ; aliveness test
           printf, apqractor_lun, 'PONG'
           end
        'QUIT': begin
           printf, apqractor_lun, 'Quitting the IDL Quickred handler'
           free_lun, apqractor_lun
           OBJ_DESTROY, qrBridge
           return
           end

        'PLUGMAPINFO': begin
           ; a new plugmapfile was detected
           ; expecting the plateId, fscan_mjd, fscan_id, filename
           ; plugMapInfo=4027,55625,2,/tmp/plPlugMapM-4027-55625-02.par
           values=strsplit(strmid(astring,12),',',/extract)
           print,'apqr_wrapper  PLUGMAPINFO  values=',values
           plateId = values[0]
           plateScanMJD = values[1]
           plateScanId = values[2]
           filename = values[3]
           plugfile = file_search(filename,count=nplugfile)
        end

        'UTR': begin
           if strpos(astring,'UTR=DONE') ge 0 then begin
              ; the UTR exposure just completed (or got aborted)
              ; expecting   UTR=DONE,frameid,mjd5,exp_pk
              ; as in:      UTR=DONE,01280012,55690,4724
              res=strsplit(astring,',',/extract)
              frameid=res[1]
              mjd5   = res[2]
              exp_pk = long(res[3])
              rawdir = data_dir+strtrim(mjd5,2)+'/'
              ;outdir = archive_dir+strtrim(mjd5,2)+'/'
              bundledir = archive_dir+strtrim(mjd5,2)+'/'
              quickreddir = quickred_dir+strtrim(mjd5,2)+'/'

              ; You can run another execute command even after there was
              ; an error.  You can do the same after there was an abort as well.
              ; You CANNOT run another execute command while
              ; the previous command is still running.

              ; look for the first available idl session
              pos=-1
              avail_pos=-1
              print,'Number of running qrBridge processes: ',n_elements(qrBridge)
              res=OBJ_VALID(qrBridge)
              p=where(res eq 1,count)
              print,'Number of valid qrBridge objects: ',count,'  res=',res
              if count ne n_elements(qrBridge) then begin
                  ; we have some dead qrBridge
                  p=where(res eq 0,complement=goodObj, badcount)
                  if badcount gt 0 then begin
                      OBJ_DESTROY,qrBridge[res[p]]
                      if n_elements(goodObj) gt 0 then begin
                          qrBridge = qrBridge[res[goodObj]]
                      endif else begin
                          apgundef,qrBridge
                      endelse
                  endif
              endif

              for proc=0,n_elements(qrBridge)-1 do begin
                  if qrBridge[proc]->Status() ne 1 then begin
                      if avail_pos[0] eq -1 then avail_pos=proc else avail_pos=[avail_pos,proc]
                      qrBridge[proc]->Execute,'retall'
                  endif
              endfor

              ; if more than 2 bridges available, kill the others
              if n_elements(avail_pos) gt 2 then begin
                  for i=2,n_elements(avail_pos)-1 do begin
                      print, 'Destroying qrBridge object ',avail_pos[i]
                      OBJ_DESTROY, qrBridge[avail_pos[i]]
                  endfor
                  res=OBJ_VALID(qrBridge)
                  p=where(res eq 1,count)
                  if count gt 0 then begin
                      qrBridge = qrBridge[res[p]]
                      ; find the first available qrbridge
                  endif else begin
                      ; all the bridges are dead -> destroy them all and start a new one
                      for i=0,n_elements(qrBridge)-1 do OBJ_DESTROY, qrBridge[i]
                      myobj = obj_new('IDL_IDLBridge', OUTPUT='', CALLBACK='qrBridge_cb', USERDATA=apqractor_lun)
                      qrBridge = myobj
                  endelse
              endif

              ; find the firts available bridge
              pos=-1
              res = OBJ_VALID(qrBridge)
              valid = where(res eq 1, count)
              for p=0,count-1 do begin
                  if qrBridge[valid[p]]->Status() ne 1 then begin
                      pos=valid[p]
                      break
                  endif
              endfor

              ; create a new bridge if all busy (up to a maximum)
              if pos eq -1 then begin
                  if n_elements(qrBridge) lt MAX_NUM_BRIDGES then begin
                      ; all processes are still busy - create a new one (within limits)
                      print,'---------> Starting a new apquickred Bridge'
                      myobj = obj_new('IDL_IDLBridge', OUTPUT='', CALLBACK='qrBridge_cb', USERDATA=apqractor_lun)
                      ; execute the startup file if we had one
                      myobj->EXECUTE, '@' + PREF_GET('IDL_STARTUP')
                      qrBridge = [qrBridge, myobj]
                      pos = n_elements(qrBridge)-1
                  endif else begin
                      print,''
                      print,'****************'
                      print,' ERROR: reached maximum number of apquickred processes (',MAX_NUM_BRIDGES,')'
                      print,'        skipping apquickred for ',frameid
                      print,'****************'
                  endelse
              endif

              if pos ge 0 then begin
                 ; qrBridge is not busy -> assume we can send it commands
                 ; print,'Starting Quick Reduction 1'
                 ; status=3 indicates an error so try to retall
                 if not (OBJ_VALID(qrBridge[pos]))[0] then begin
                     ; the object is invalid -> destroy it a start a new bridge
                     OBJ_DESTROY, qrBridge[pos]
                     print,'---------> Starting a new apquickred Bridge'
                     myobj = obj_new('IDL_IDLBridge', OUTPUT='', CALLBACK='qrBridge_cb', USERDATA=apqractor_lun)
                     ; execute the startup file if we had one
                     myobj->EXECUTE, '@' + PREF_GET('IDL_STARTUP')
                     qrBridge[pos] = myobj
                 endif
                 if qrBridge[pos]->Status() eq 3 then qrBridge[pos]->Execute,'retall'
                 qrBridge[pos]->SetVar, 'frameid', frameid
                 qrBridge[pos]->SetVar, 'rawdir', rawdir
                 qrBridge[pos]->SetVar, 'bundledir', bundledir
                 qrBridge[pos]->SetVar, 'quickreddir', quickreddir
                 qrBridge[pos]->SetVar, 'exp_pk', exp_pk
                 qrBridge[pos]->SetVar, 'bpmid', bpmid
                 qrBridge[pos]->SetVar, 'psfid', psfid
                 qrBridge[pos]->SetVar, 'mjd5', mjd5
                 ; pass the whole structure
                 if n_elements(plugfile) gt 0 then $
                   qrBridge[pos]->SetVar, 'plugfile', plugfile
                 qrBridge[pos]->Execute, $
                    'apquickred,frameid,plugfile=plugfile,rawdir=rawdir,bundledir=bundledir,'+$
                    'quickreddir=quickreddir,bpmid=bpmid,psfid=psfid,exp_pk=exp_pk,mjd5=mjd5',/NOWAIT
                 break
              endif
           endif
           end

        else: begin
           ; echo the input string back to the sender
           nchar = strlen(astring)
           printf, apqractor_lun, 'Received '+strtrim(string(nchar),2)+' characters'
           printf, apqractor_lun, astring
        end

     ENDCASE


     ; Increment counter
     num_cmds++

   ENDWHILE  ; large loop

END
