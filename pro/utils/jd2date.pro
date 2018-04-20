function jd2date,jd

;+
;
; JD2DATE
;
; Convert a Julain date to a date string
;
; INPUTS:
;  jd      The Julian date
;
; OUTPUTS:
;  output  The date string in YYYY-MM-DDTHH:MM:SS.S format.
;
; USAGE:
;  IDL>date = jd2date(394830.00)
;
; By D.Nidever  Dec 2010
;-

;date will be in this format
;  2004-08-08T06:55:20

; Not enough inputs
if n_elements(jd) eq 0 then begin
  print,'Syntax - date = jd2date(jd)'
  return,-1
endif

; Convert JD to date
CALDAT,jd,month,day,year,hour,min,sec

date = string(year,format='(I4)')+'-'+string(month,format='(I02)')+'-'+string(day,format='(I02)')+'T'
date = date+ string(hour,format='(I02)')+':'+string(min,format='(I02)')+':'+strtrim(string(sec,format='(F04.1)'),2)

return,date
end
