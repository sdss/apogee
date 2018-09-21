pro apql_dbinsert_callback, status, error, dbBridge, userdata

; This informs the actor that the apql database insert completed

; Inform the actor that the dbinsert finished
; send back the quicklook_pk value

; If new predictions were added infor the actor
if n_elements(readnum) gt 0 then if (long(readnum) mod 6 eq 0) then begin
  reply = 'NEW PREDICTIONS'
  printf, apqlactor_lun,reply
endif

end
