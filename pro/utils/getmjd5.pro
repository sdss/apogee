function getmjd5,head,error=error,silent=silent,sdss=sdss

;+
;
; GETMJD5
;
; This function returns the MJD5 5-digit Modified
; Julian Date given an observation header.  The
; header must contain DATE-OBS.
;
; INPUTS:
;  head     The observation header.  Must be a string array
;  /sdss    Use the SDSS version which adds 0.3 days to MJD
;  /silent  Don't print anything to the screen.
;  /stp     Stop at the end of the program.
;
; OUTPUTS:
;  MJD5, the 5-digit Modified Julian date (number).
;  The trailing decimal is clipped.
;  =error  The error message if one occured.
;
; USAGE:
;  IDL>mjd5 = getmjd5(head)
;
; By D.Nidever  Feb 2010
;-

apgundef,error

; Not enough inputs
if n_elements(head) eq 0 then begin
  error = 'Syntax - mjd5 = getmjd5(head,sdss=sdss,error=error,silent=silent)'
  print,error
  return,-1
endif

; Head must be string array
nhead = size(head,/n_elements)
type = size(head,/type)
if nhead lt 2 or type ne 7 then begin
  error = 'HEAD must be a string array'
  if not keyword_set(silent) then print,error
  return,-1
endif

; Get DATE-OBS
dateobs = SXPAR(head,'DATE-OBS',count=count)
if count eq 0 then begin
  error = 'NO DATE-OBS keyword in header'
  if not keyword_set(silent) then print,error
  return,-1
endif
if size(dateobs,/type) ne 7 then begin
  error = 'DATE-OBS is NOT a string'
  if not keyword_set(silent) then print,error
  return,-1
endif


; Parse DATE-OBS
;----------------
;  DATE-OBS= '2007-08-17T02:05:36.7' 
arr = strsplit(dateobs,'T',/extract)
if n_elements(arr) ne 2 then begin
  error = 'DATE-OBS = '+dateobs+' NOT in the correct format. YYYY-MM-DDTHH:MM:SS.S'
  if not keyword_set(silent) then print,error
  return,-1
endif
datearr = strsplit(arr[0],'-',/extract)
timearr = strsplit(arr[1],':',/extract)
if n_elements(datearr) lt 3 or n_elements(timearr) lt 3 then begin
  error = 'DATE-OBS = '+dateobs+' NOT in the correct format. YYYY-MM-DDTHH:MM:SS.S'
  if not keyword_set(silent) then print,error
  return,-1
endif

; Check each number
year = datearr[0]
if VALID_NUM(year,/integer) eq 0 then begin
  error = 'YEAR = '+year+' is NOT a valid integer'
  if not keyword_set(silent) then print,error
  return,-1
endif
year_num = long(year)

month = datearr[1]
if VALID_NUM(month,/integer) eq 0 then begin
  error = 'MONTH = '+month+' is NOT a valid integer'
  if not keyword_set(silent) then print,error
  return,-1
endif
month_num = long(month)
if month_num lt 1 or month_num gt 12 then begin
  error = 'MONTH = '+month+' MUST be 1-12'
  if not keyword_set(silent) then print,error
  return,-1
endif

day = datearr[2]
if VALID_NUM(day,/integer) eq 0 then begin
  error = 'DAY = '+day+' is NOT a valid integer'
  if not keyword_set(silent) then print,error
  return,-1
endif
day_num = long(day)
if day_num lt 1 or day_num gt 31 then begin
  error = 'DAY = '+day+' Must be 1-31'
  if not keyword_set(silent) then print,error
  return,-1
endif

hour = timearr[0]
if VALID_NUM(hour,/integer) eq 0 then begin
  error = 'hour = '+hour+' is NOT a valid integer'
  if not keyword_set(silent) then print,error
  return,-1
endif
hour_num = long(hour)
if hour_num lt 0 or hour_num gt 24 then begin
  error = 'HOUR = '+hour+' Must be 1-24'
  if not keyword_set(silent) then print,error
  return,-1
endif

minute = timearr[1]
if VALID_NUM(minute,/integer) eq 0 then begin
  error = 'MIN = '+minute+' is NOT a valid integer'
  if not keyword_set(silent) then print,error
  return,-1
endif
minute_num = long(minute)
if minute_num lt 0 or minute_num gt 60 then begin
  error = 'MIN = '+minute+' Must be 1-60'
  if not keyword_set(silent) then print,error
  return,-1
endif

sec = timearr[2]
if VALID_NUM(sec) eq 0 then begin
  error = 'SEC = '+sec+' is NOT a valid number'
  if not keyword_set(silent) then print,error
  return,-1
endif
sec_num = float(sec)
if sec_num lt 0 or sec_num gt 60 then begin
  error = 'SEC = '+sec+' Must be 1-60'
  if not keyword_set(silent) then print,error
  return,-1
endif


; Getting Julian Date
;---------------------
jd = JULDAY(month_num, day_num, year_num, hour_num, minute_num, sec_num)

; Modified Julian Date
;----------------------
; The modified Julian date is MJD = JD - 2400000.5
; The Julian day starts at NOON, while MJD starts at midnight
MJD = JD - 2400000.5
if keyword_set(sdss) then MJD+=0.3    ; use SDSS version
mjd5 = long(mjd)  ; clip the decimals


if keyword_set(stp) then stop

return,mjd5

end
