;+
;
; APQL_GET_PREVEXP
;
; This routine is called by the quicklook program to get the list of 
; OBJECT exposures already taken with the specifed plateId.
;
; The program returns a structure with the information needed by   
; the exposureData message sent back to the actor then back to STUI.
; 
;-
;
FUNCTION apql_get_prevexp, plateid, snrgoal=snrgoal, count=count, silent=silent

   t0 = systime(1)  

   count=0
   if n_elements(snrgoal) eq 0 then snrgoal=30.0

   tempstr = {plateId:'',expNum:0, expName:'', exptime:0.0, numReads:0, snrGoal:float(snrgoal), $
       ditherPos:0.0,snr:0.0,netExpTime:0.0,netSnr:0.0, expType:'Object', namedDitherPos:''}

   if n_elements(plateid) eq 0 or n_elements(plateid) gt 1 then begin
      if not keyword_set(silent) then $
          print,'Usage:  result = apql_get_prevexp(plateId, snrgoal=snrgoal, count=count)'
      return, tempstr
   endif

   ; make sure the plateId is for APOGEE-2
   valid=0
   get_sql_col,"SELECT plt.pk FROM ((platedb.plate AS plt "+$
     "INNER JOIN platedb.plate_to_survey AS p2s ON (plt.pk=p2s.plate_pk)) "+$
     "INNER JOIN platedb.survey AS sv ON (p2s.survey_pk=sv.pk)) "+$
     "WHERE sv.label LIKE 'APOGEE-2' AND plt.plate_id="+strtrim(string(plateid),2), valid,/long
   if valid eq 0 then begin
       message,'not a valid APOGEE-2 plate',/cont
       return, tempstr
   endif

   exp_no='-1'
   get_sql_col, "SELECT plt.plate_id, exp.exposure_no, exp.exposure_time, ql.readnum, "+$
        "qr.dither_pixpos, qr.snr_standard "+$
        "FROM ((((((platedb.plate AS plt INNER JOIN platedb.plate_pointing AS pltg ON (plt.pk=pltg.plate_pk))  "+$
        "INNER JOIN platedb.observation AS obs ON (obs.plate_pointing_pk=pltg.pk))  "+$
        "INNER JOIN platedb.exposure AS exp ON (exp.observation_pk=obs.pk)) "+$
        "INNER JOIN platedb.exposure_flavor as expflv ON (exp.exposure_flavor_pk=expflv.pk)) "+$
        "INNER JOIN apogeeqldb.quickred as qr ON (qr.exposure_pk=exp.pk)) "+$
        "INNER JOIN apogeeqldb.quicklook as ql ON (qr.last_quicklook_pk=ql.pk)) "+$
        "WHERE plt.plate_id="+strtrim(string(plateid),2)+" AND expflv.label='Object' ORDER BY exp.exposure_no", $
        plt, exp_no, exp_time, readNum, ditherpos, snr, /string

    if exp_no[0] eq '-1' then begin
        ; nothing was returned from the query
        count=0
        return, tempstr
    endif

    count = n_elements(exp_no)
    outstr=replicate(tempstr,count)
    netExpTime=0.0
    ditherA = 13.0 ; dither position A corresponds to this pixel value

    for i=0,count-1 do begin
        outstr[i].plateId = plt[i]
        outstr[i].expNum  = i+1
        outstr[i].expName = strtrim(string(long(exp_no[i]),format='(i08)'),2)
        outstr[i].exptime = float(exp_time[i])
        outstr[i].numReads = long(readNum[i])
        outstr[i].ditherPos = float(ditherPos[i])
        outstr[i].snr = float(snr[i])
        netExpTime += exp_time[i]
        outstr[i].netExpTime = netExpTime 
        netSnr = sqrt(total(snr[0:i]^2.0))
        outstr[i].netSnr = netSnr 
        if abs(ditherPos[i]-ditherA) lt 0.1 then $
            outstr[i].namedDitherPos='A' $
        else $
            outstr[i].namedDitherPos='B'
    endfor

    if not keyword_set(silent) then print,'apql_get_prevexp completed'
    dt = systime(1)-t0
    if not keyword_set(silent) then print,'dt = ',strtrim(dt,2),' sec.'
    return, outstr

END
