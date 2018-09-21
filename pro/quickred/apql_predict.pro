;+
;
; APQL_PREDICT
;
; This is part of the APOGEE quicklook software and predicts
; what exposures to take for the rest of this plate visit.
; This should only be run every ~60 sec.
; This uses a "simulation" based on current seeing/weather
; conditions and the rules for when to spot an exposure to
; make predictions.
;
; INPUTS:
;  str           Structure with information on the current red
;  allstr        Same as STR but for previous reads of this exposure
;  prevstr       Structure of information for all previous exposures
;                  for this plate.  This is essentially the last
;                  element of ALLSTR for each exposure.
;  /summary      Give summary of exposures.
;  /verbose      Print a lot to the screen.
;  /silent       Don't print anything to the screen.
;
; OUTPUTS:
;  predict_str   Summary structure of the predicted exposures from
;                  the simulation.
;  allstatus_str The array of STATUS_STR structures from APQL_STATUSCHECK
;                  for all the simulated reads.
;  =error        The error message if one occurred.
;
; USAGE:
;  IDL>apql_predict,str,allstr,prevstr,predict_str,allstatus_str
;
; By D.Nidever  April 2011
;-
pro apql_predict,str,allstr,prevstr,predict_str,allstatus_str,verbose=verbose,$
                 silent=silent,error=error,summary=summary

   ; Error Handling
   ;------------------
   ; Establish error handler. When errors occur, the index of the
   ; error is returned in the variable Error_status:  
   ;CATCH, Error_status 

   ;This statement begins the error handler:  
   ;if (Error_status ne 0) then begin 
   ;   error = !ERROR_STATE.MSG  
   ;   if not keyword_set(silent) then print,error
   ;   CATCH, /CANCEL 
   ;   return
   ;endif

   apgundef,predict_str,allstatus_str

   ; Not enough inputs
   if n_elements(str) eq 0 then begin
     error = 'Not enough inputs'
     if not keyword_set(silent) then $
       print,'Syntax - apql_predict,str,allstr,prevstr,predict_str,allstatus_str,'+$
           'verbose=verbose,summary=summary,silent=silent,error=error'
     return
   endif

   ;timestep = 10.6  ; seconds per read
   timestep = !apql.timePerRead
   if timestep eq 0.0 then timestep=10.6

   ; Template structure for SIM_STR and others
   tempstr = {frameid:'',readnum:'',exptime:0.0,snr_standard_goal:0.0,snr_standard:0.0,$
              exp_finished_status:0,visit_finished_status:0, exptype:''}

   ; Do we have any signal
   if str.snr_standard eq 0.0 or str.logsnr_hmag_coef[1] eq 0.0 or FINITE(str.snr_standard) eq 0 then begin
     error = 'APQL_PREDICT: No signal detected.  Cannot do prediction'
     if not keyword_set(silent) then print,error
     return
   endif

   ; Initialize SIM_ALLSTR
   nallstr = n_elements(allstr)
   if nallstr gt 0 then begin
     sim_allstr = REPLICATE(tempstr,nallstr)
     for i=0,nallstr-1 do begin
       sim_allstr[i].frameid = allstr[i].frameid
       sim_allstr[i].readnum = allstr[i].readnum
       sim_allstr[i].exptime = allstr[i].exptime
       sim_allstr[i].snr_standard_goal = allstr[i].snr_standard_goal
       sim_allstr[i].snr_standard = allstr[i].snr_standard
       sim_allstr[i].exptype = allstr[i].exptype
     end
   endif

   ; Initialize SIM_PREVSTR
   nprevstr = n_elements(prevstr)
   if nprevstr gt 0 then begin
     sim_prevstr = REPLICATE(tempstr,nprevstr)
     for i=0,nprevstr-1 do begin
       sim_prevstr[i].frameid = prevstr[i].frameid
       sim_prevstr[i].readnum = prevstr[i].readnum
       sim_prevstr[i].exptime = prevstr[i].exptime
       sim_prevstr[i].snr_standard_goal = prevstr[i].snr_standard_goal
       sim_prevstr[i].snr_standard = prevstr[i].snr_standard
       sim_prevstr[i].exptype = prevstr[i].exptype
     end
   end


   ; Make "fake" structure that only have the information in it that
   ; we need to check the status
   snr2_time_coef = str.snr2_time_coef   ; SNR^2 vs. readnum
   ;if snr2_time_coef[0] lt 0.0 then snr2_time_coef[0] = 0   ; need non-negative y-offset

   framenum0 = long(str.frameid)   ; starting frame number


   ; Exposure loop
   ;----------------
   exposure_count = 0
   exp_endflag = 0
   WHILE (exp_endflag eq 0) do begin

     ; Get frame number
     framenum = framenum0 + exposure_count

     if keyword_set(verbose) then begin
       print,'EXPOSURE = ',strtrim(framenum,2)
       print,'================'
       print,'Frameid  Readnum  S/N       Exptime Visittime  Dither Tstatus_exp Tstatus_visit Fstatus_exp Fstatus_visit'
       print,'-----------------------------------------------------------------------------------------------------------'
     endif

     ; Read loop
     ;----------
     if exposure_count eq 0 then readnum = long(str.readnum) $ ; current exposure
       else readnum = 2      ; start with 2nd read
     read_endflag = 0
     While (read_endflag eq 0) do begin

       ; Initialize SIM_STR
       sim_str = tempstr
       sim_str.frameid = string(framenum,format='(I08)')
       sim_str.readnum = string(readnum,format='(I03)')
       sim_str.snr_standard_goal = str.snr_standard_goal
       sim_str.exptype = str.exptype
       ; current exposure
       if exposure_count eq 0 then begin
         sim_str.exptime = str.exptime + (readnum-long(str.readnum))*timestep
         snr2_starting_meas = str.snr_standard^2                                ; actual measured starting SNR2
         snr2_starting_poly = poly( long(str.readnum)*timestep,snr2_time_coef)  ; polyfit starting SNR2
         snr2_current_poly = poly( readnum*timestep, snr2_time_coef)            ; polyfit current SNR2
         snr_standard = sqrt( snr2_starting_meas + (snr2_current_poly - snr2_starting_poly) )
         sim_str.snr_standard = snr_standard
       ; future exposures
       endif else begin
         sim_str.exptime = readnum*timestep
         snr_standard = sqrt( poly(readnum*timestep,snr2_time_coef) )
         if finite(snr_standard) eq 0 then snr_standard=0
         sim_str.snr_standard = snr_standard > 0   ; must be >=0
       endelse


       ; What we need to update
       ; str.exptime
       ; str.snr_standard
       ; Nprevstr to get which exp of the dither it is
       ; prevstr.snr_standard for S/N of previous exposure

       ; Check the status
       APQL_STATUSCHECK,sim_str,sim_allstr,sim_prevstr,status_str

       ; Add to SIM_ALLSTR
       PUSH,sim_allstr,sim_str

       ; Add to ALLSTATUS_STR
       PUSH,allstatus_str,status_str

       ; Finished
       ; to prevent an infinite loop if there's an error in APQL_STATUSCHECK
       ; changed to ge 0 (since the status is initialized to -1)
       ;if status_str.exp_finished_status ge 0 then read_endflag=1
       if status_str.exp_finished_status gt 0 then read_endflag=1

       ; Printing info
       if keyword_set(verbose) then begin
         print,sim_str.frameid,' ',sim_str.readnum,sim_str.snr_standard,$
               status_str.exposure_time, status_str.visit_time, status_str.ditherexposure,$
               status_str.tstatus_exp, status_str.tstatus_visit, status_str.exp_finished_status, status_str.visit_finished_status
       endif

       ;stop

       readnum++  ; increment read number

     Endwhile

     ; Add to SIM_PREVSTR
     nsim_allstr = n_elements(sim_allstr)
     PUSH,sim_prevstr,sim_allstr[nsim_allstr-1]

     ; Finished
     if status_str.visit_finished_status gt 0 then exp_endflag=1


     if keyword_set(verbose) then begin
       print,'-----------------------------------------------------------------------------------------------------------'
     endif


     exposure_count++  ; increment exposures counter

     ;stop

   ENDWHILE


   ; Print past exposures
   if n_elements(prevstr) gt 0 and (keyword_set(verbose) or keyword_set(summary)) then begin
     print,'PREVIOUS EXPOSURES'
     print,'------------------------------------------------------------------'
     print,'NUM  Frameid Nreads Exptime    S/N   Dither EXPStatus VISITStatus'
     print,'------------------------------------------------------------------'
     for i=0,n_elements(prevstr)-1 do begin
       format = '(I2,A10,I5,F10.2,F8.2,I6,I8,I8)'

       dithpix = sxpar(prevstr[i].header,'DITHPIX')
       if i eq 0 then begin
         ditherexposure = 1
         dithpix0 = dithpix
       endif else begin
         delta_dither = dithpix-dithpix0
         if abs(delta_dither) gt 0.3 then ditherexposure=2 else ditherexposure=1
       endelse

       print,i+1,prevstr[i].frameid,long(prevstr[i].readnum),prevstr[i].exptime,prevstr[i].snr_standard,$
             ditherexposure,prevstr[i].exp_finished_status,prevstr[i].visit_finished_status,format=format
     end
     print,'------------------------------------------------------------------'
   endif

   ; Make a "summary" structure of the predictions
   ;------------------------------------------------
   ;  this includes CURRENT and FUTURE exposures, not past ones
   ui = uniq(allstatus_str.frameid,sort(allstatus_str.frameid))
   uframeid = allstatus_str[ui].frameid
   nexposures = n_elements(uframeid)
   predict_str = REPLICATE({frameid:'',nreads:0L,exptime:0.0,snr_standard:0.0,tmin_exp:0.0,tmax_exp:0.0,tmin_visit:0.0,$
                            tmax_visit:0.0,ditherexposure:0,exp_finished_status:0,visit_finished_status:0},nexposures)
   if keyword_set(verbose) or keyword_set(summary) then begin
     print,'CURRENT AND PREDICTED EXPOSURES'
     print,'------------------------------------------------------------------'
     print,'NUM  Frameid Nreads Exptime    S/N   Dither EXPStatus VISITStatus'
     print,'------------------------------------------------------------------'
   end
   for i=0,nexposures-1 do begin
     ind = where(allstatus_str.frameid eq uframeid[i],nind)
     lastind = max(ind)
     lastread = allstatus_str[lastind]
     predict_str[i].frameid = uframeid[i]
     predict_str[i].nreads = long(lastread.readnum)
     predict_str[i].exptime = lastread.exposure_time
     predict_str[i].snr_standard = lastread.snr_standard
     predict_str[i].tmin_exp = lastread.tmin_exp
     predict_str[i].tmax_exp = lastread.tmax_exp
     predict_str[i].tmin_visit = lastread.tmin_visit
     predict_str[i].tmax_visit = lastread.tmax_visit
     predict_str[i].ditherexposure = lastread.ditherexposure
     predict_str[i].exp_finished_status = lastread.exp_finished_status
     predict_str[i].visit_finished_status = lastread.visit_finished_status

     if keyword_set(verbose) or keyword_set(summary) then begin
       format = '(I2,A10,I5,F10.2,F8.2,I6,I8,I8)'
       print,i+1,predict_str[i].frameid,predict_str[i].nreads,predict_str[i].exptime,predict_str[i].snr_standard,$
             predict_str[i].ditherexposure,predict_str[i].exp_finished_status,predict_str[i].visit_finished_status,format=format
     endif
   endfor

   if keyword_set(verbose) or keyword_set(summary) then $
     print,'------------------------------------------------------------------'

   visit_exptime = total(predict_str.exptime)
   ;visit_snr = sqrt( total(predict_str.snr_standard^2) )
   visit_snr = sqrt( total(sim_prevstr.snr_standard^2) ) ; includes previous exposures
   if keyword_set(verbose) then begin
     print,'Visit exptime = ',stringize(visit_exptime,ndec=2,/nocomma),' sec'
     print,'Visit S/N     = ',stringize(visit_snr,ndec=2)
   endif

   ;stop

end
