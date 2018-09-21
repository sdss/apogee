;+
;
; UPDATE_QR_LASTQL
;
; This routine is meant to update the value in the quickred table for the
; last_quicklook_pk column.  A bug was found in the code where it was using the
; last quicklook pk entry but when the query was executed, many more quicklook
; might have happen before for a different exposure.
; 
; The code was modified to get the last pk for this exposure but we need to fix
; the data already in the database.
;
;-
;
pro update_qr_lastql

   ; get he list of all the QUICKRED entries
   get_sql_col, "SELECT pk, exposure_pk FROM apogeeqldb.quickred ORDER BY pk",qr_pk, exp_pk, /string

   count=n_elements(qr_pk)
   print,''
   print,'Number of rows to update: ',count
   print,''

   for i=0L, count-1L do begin
      print,strtrim(string(i+1L),2)+' / '+strtrim(string(count),2)+'  qr_pk=',qr_pk[i],'  exp_pk=',exp_pk[i]
      get_sql_col, "SELECT pk FROM apogeeqldb.quicklook WHERE exposure_pk="+strtrim(exp_pk[i],2)+$
          " ORDER BY pk DESC LIMIT 1",last_qlpk,/string
      if n_elements(last_qlpk) ne 1 then begin
          print,'    n_elements(last_qlpk)=',n_elements(last_qlpk)
          continue
      endif else begin
          print,"    update last_qlpk=",last_qlpk[0]
          exec_sql,'UPDATE apogeeqldb.quickred SET last_quicklook_pk='+last_qlpk[0]+' WHERE pk='+qr_pk[i]
      endelse
   endfor


END
