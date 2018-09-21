;+
;
; APQL_STATUSCHECK
;
; Check the completion status of the current APOGEE exposure using
; a set of "rules".
;
; INPUTS:
;  str          Structure with information on the current red
;  allstr       Same as STR but for previous reads of this exposure
;  prevstr      Structure of information for all previous exposures
;                 for this plate.  This is essentially the last
;                 element of ALLSTR for each exposure.
;  /silent  Don't print anything to the screen.
;
; OUTPUTS:
;  status_str   A structure giving information on the completion status
;  =error   The error message if one occurred.
;
; USAGE:
;  IDL>apql_statuscheck,str,allstr,prevstr,status_str
;
; By D.Nidever  April 2011
;-
pro apql_statuscheck,str,allstr,prevstr,status_str,silent=silent,error=error


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

; Not enough inputs
if n_elements(str) eq 0 then begin
  error = 'Not enough inputs'
  if not keyword_set(silent) then $
    print,'Syntax - apql_statuscheck,str,allstr,prevstr,status_str,silent=silent,error=error'
  return
endif


; Quicklook Exposure Time Estimation Algorithm
; 
; This capability of the quicklook is really it's main purpose and since all of the quicklook and STUI software
; is now coming together we need to decide on a working (but not necessarily final) strategy for the quicklook.
; 
; This is of course difficult because we are trying to optimize several things:
;
; 1. Don't want to be readnoise dominated, so we want to get enough reads to beat that down. That effectively sets
;     a lower exptime limit.
; 2. Don't want the airglow lines to saturate since we will use them for wavelength zeropoint calibration and
;     probably for LSF determination.
; 3. Want each exposure to have a corresponding dither pair. Don't want "orphan" exposures if we can help it.
; 4. Want each exposure of a dither pair to have similar S/N
; 5. Want the total S/N in the plate "visit" (all exposures combined) to reach a certain level. It's still not
;     clear yet what that is.
; 6. We are limited to a minimum of 50 min. and a maximum of 60 min. for each plate "visit" on MARVELS co-observing
;     plates.
; 7. Latency considerations
;
; I was thinking about this last night and came up with some ideas on potentially how to do this. These are not
; really well fleshed out and just something to get the discussion started.
;
; We need a set of rules (decision tree) that quicklook can use to determine what to do for any exposure. Then
; at any given point during the visit it can set out a "plan" for what to do for the rest of that plate visit.
; It can do this by assuming that the (S/N)2/time rate will be the same (and constant) as what it is currently
; (averaged of last ~1 min.) and then do a little "simulation" of running the time forward, checking the rules,
; and tabulating the "decisions" (of whether to continue exposing or terminate and start a new exposure). Then
; it can summarize the "results" in a table as the prediction (at the current point) of what to do for the rest
; of the visit. Maybe only run this simulation every minute or so because with variable seeing/cloud conditions 
; this might make the prediction dance around too much, but maybe the averaging might keep it steady. Hopefully
; it's clear what I'm saying.
;
; There should probably be fairly strict time restrictions that are part of the "rules":
;
; T1. TMIN_EXP - the minimum exposure time for an individual exposure set by the need to not be readnoise dominated.
;                 ~5 min?? 
; T2. TMAX_EXP - the maximum exposure time for an individual exposure set by the need to not saturate the airglow
;                 lines since we need them for wavelength zeropoint calibration and likely LSF calibration. ~30 min.
; T3. TMIN_VISIT - the minimum time for a plate visit. For MARVELS shared plates this is 50 min. Not sure what it
;                 is for APOGEE-only plates.
; T4. TMAX_VISIT - the maximum time for a plate visit. For MARVELS shared plates this is 60 min. Not sure what
;                 it is for APOGEE-only plates.
;
; Here's a first pass at the rules:
;
; R1. If this is the first exposure of a dither pair then go as long as needed to get the the goal S/N at the
;       fiducial H magnitude, bounded by the time restrictions. Make sure that there is enough time at the end
;       of this exposure (using TMAX_VISIT) to be able to take a second exposure (in the dither pair) of equal
;       length. This ensures that we get the dithers in pairs. This doesn't take the S/N consideration into account
;       (for the second exposure of the dither pair), but it's a start.
; R2. If this is the second exposure of a dither pair, then go as long as needed to match the S/N of the first
;       exposure of this dither pair, bounded by the time restrictions.
; R3. If we have passed TMIN_VISIT _and_ (we have reached TMAX_VISIT _or_ there is not enough time to take another
;       dither pair (maybe using the average length of the previous exposures plus corrections for any seeing/cloud
;       condition changes) ), then call this visit done.
;
; This is just a first draft, but I think it's headed in the right direction.



; Initialize STATUS_STR
status_str = {frameid:str.frameid,readnum:str.readnum,exposure_time:0.0,snr_standard:str.snr_standard,visit_time:0.0,$
              tmin_exp:0L,tmax_exp:0L,tmin_visit:0L,tmax_visit:0L,$
              ditherexposure:-1, snr_finished:-1,tstatus_exp:-1,tstatus_visit:-1,secondexp_limit:-1,$
              exp_finished_status:-1, visit_finished_status:-1}


;===================
; OBJECT EXPOSURE
;===================
If str.exptype eq 'OBJECT' then begin



  ;====================
  ; TIME CONSTRAINTS
  ;====================

  ;--------------------------------------------------------
  ; T1. TMIN_EXP - the minimum exposure time for an individual exposure set by the need to not be readnoise dominated.
  ;                 ~5 min??
  TMIN_EXP = 300

  ;--------------------------------------------------------
  ; T2. TMAX_EXP - the maximum exposure time for an individual exposure set by the need to not saturate the airglow
  ;                 lines since we need them for wavelength zeropoint calibration and likely LSF calibration. ~30 min.
  TMAX_EXP = 1800

  ;--------------------------------------------------------
  ; T3. TMIN_VISIT - the minimum time for a plate visit. For MARVELS shared plates this is 50 min. Not sure what it   
  ;                 is for APOGEE-only plates.
  TMIN_VISIT = 3000

  ;--------------------------------------------------------
  ; T4. TMAX_VISIT - the maximum time for a plate visit. For MARVELS shared plates this is 60 min. Not sure what
  ;                 it is for APOGEE-only plates.
  TMAX_VISIT = 3600

  ; But them into a structure
  status_str.tmin_exp = tmin_exp
  status_str.tmax_exp = tmax_exp
  status_str.tmin_visit = tmin_visit
  status_str.tmax_visit = tmax_visit


  ;============
  ; RULES
  ;============


  ; ditherexposure = 1  first exp of dither pair
  ; ditherexposure = 2  second exp of dither pair

  ; Figure out which exposure of the pair this is
  ;------------------------------------------------
  nprevstr = n_elements(prevstr)

  ; First exposure of plate
  if n_elements(prevstr) eq 0 then begin
    ditherexposure = 1

  ; Not first exposure of plate
  endif else begin

    ; use dither_prevexp_header to determine this

    ; can't we just use odd-first, even-second???
    ; just need to check for "bad" exposures
    if odd(nprevstr+1) eq 1 then ditherexposure=1 else ditherexposure=2

  endelse ; not first exposure of plate
  status_str.ditherexposure = ditherexposure


  ; EXP_FINISHED_STATUS bit values:
  ; 1-reached S/N goal
  ; 2-hit maximum exposure time limit
  ; 4-hit maximum visit time limit
  ; 8-hit time constraint for second exposure  (only for 1st exposures
  ;                                             of a dither pair)
  exp_finished_status = 0

  ; VISIT_FINISHED_STATUS bit values:
  ; 1-hit maximum visit limit
  ; 2-not enough time for another dither pair
  visit_finished_status = 0


  ; Check time restrictions
  ;-------------------------
  ; tstatus: 0-short, 1-in "good" window, 2-too long
  ;  exposure time
  exposure_time = str.exptime
  if exposure_time lt tmin_exp then tstatus_exp=0
  if exposure_time ge tmin_exp and exposure_time le tmax_exp then tstatus_exp=1 
  if exposure_time gt tmax_exp then tstatus_exp=2 
  status_str.tstatus_exp = tstatus_exp
  status_str.exposure_time = exposure_time
  ;  visit time
  if nprevstr gt 0 then visit_time = TOTAL(prevstr.exptime)+str.exptime else visit_time=exposure_time
  if visit_time lt tmin_visit then tstatus_visit=0
  if visit_time ge tmin_visit and visit_time le tmax_visit then tstatus_visit=1 
  if visit_time gt tmax_visit then tstatus_visit=2 
  status_str.tstatus_visit = tstatus_visit
  status_str.visit_time = visit_time

  ; Time left in this visit
  visit_time_left = tmax_visit - visit_time


  ;=================================
  ; First exposure of a dither pair
  ;=================================

  ;-------------------------------------------
  ; R1. If this is the first exposure of a dither pair then go as long as needed to get the goal S/N at the
  ;       fiducial H magnitude, bounded by the time restrictions. Make sure that there is enough time at the end
  ;       of this exposure (using TMAX_VISIT) to be able to take a second exposure (in the dither pair) of equal
  ;       length. This ensures that we get the dithers in pairs. This doesn't take the S/N consideration into account
  ;       (for the second exposure of the dither pair), but it's a start.
  IF (ditherexposure eq 1) THEN BEGIN


    ; Check S/N for fiducial magnitude
    ;-------------------------------------
    if str.snr_standard ge str.snr_standard_goal then snr_finished=1 else snr_finished=0
    status_str.snr_finished = snr_finished

    ; Check time for a second exposure
    ;-----------------------------------
    ;beg_visit_time = visit_time - exposure_time  ; vist_time at beg of this exposure
    if (exposure_time ge visit_time_left) then secondexp_limit=1 else secondexp_limit=0
    status_str.secondexp_limit = secondexp_limit
  
    ; Finished with this EXPOSURE
    ;-----------------------------
    ; We reached our S/N goal and we are long enough
    if snr_finished eq 1 and tstatus_exp eq 1 then exp_finished_status += 1

    ; We're hitting up against the MAXIMUM time limits
    if tstatus_exp eq 2 then exp_finished_status += 2
    if tstatus_visit eq 2 then exp_finished_status += 4

    ; We're hitting the second exposure limit
    if secondexp_limit eq 1 then exp_finished_status += 8
    status_str.exp_finished_status = exp_finished_status

    ; Finished with this VISIT
    ;-------------------------
    ; We reached the maximum time limit
    if tstatus_visit eq 2 then visit_finished_status = 1
    status_str.visit_finished_status = visit_finished_status

    ;stop


  ;==================================
  ; Second exposure of a dither pair
  ;==================================

  ;-------------------------------------------
  ; R2. If this is the second exposure of a dither pair, then go as long as needed to match the S/N of the first
  ;       exposure of this dither pair, bounded by the time restrictions.
  ENDIF ELSE BEGIN


    ; Check S/N compared to first exposure of the dither pair
    ;---------------------------------------------------------
    snr_standard_prev = prevstr[nprevstr-1].snr_standard
    if str.snr_standard ge snr_standard_prev then snr_finished=1 else snr_finished=0
    status_str.snr_finished = snr_finished


    ; Finished with this EXPOSURE
    ;-----------------------------
    ; We reached our S/N goal and we are long enough
    if snr_finished eq 1 and tstatus_exp eq 1 then exp_finished_status += 1

    ; We're hitting up against the MAXIMUM time limits
    if tstatus_exp eq 2 then exp_finished_status += 2
    if tstatus_visit eq 2 then exp_finished_status += 4
    status_str.exp_finished_status = exp_finished_status


    ; Finished with this VISIT
    ;-------------------------
    ; We reached the maximum time limit
    if tstatus_visit eq 2 then visit_finished_status = 1

    ;-------------------------------------------
    ; R3. If we have passed TMIN_VISIT _and_ (we have reached TMAX_VISIT _or_ there is not enough time to take another  
    ;       dither pair (maybe using the average length of the previous exposures plus corrections for any seeing/cloud
    ;       condition changes) ), then call this visit done.
    if tstatus_visit eq 1 and ( visit_time_left lt 2*tmin_exp) then visit_finished_status += 2

    status_str.visit_finished_status = visit_finished_status

    ;stop


  ENDELSE ; second exposure of dither pair

  ; Update STR with the vital information
  str.exp_finished_status = status_str.exp_finished_status
  str.visit_finished_status = status_str.visit_finished_status


; NON-OBJECT Exposures
;---------------------
Endif else begin

  ; Just use the S/N criterion
  ;---------------------------
  if str.snr_standard ge str.snr_standard_goal then snr_finished=1 else snr_finished=0
  status_str.snr_finished = snr_finished
  str.exp_finished_status = snr_finished
  ; for non-objects we only want one exposure, so visit=exposure
  str.visit_finished_status = snr_finished

Endelse


;stop

end
