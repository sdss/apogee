;+
; NAME: 
;  addslash
; PURPOSE: 
;  Append a trailing / to string (if needed).
; DESCRIPTION:
; CATEGORY:
;  Utility
; CALLING SEQUENCE:
;  addslash,name
; INPUTS:
;  name - string to modify
; OPTIONAL INPUT PARAMETERS:
; KEYWORD INPUT PARAMETERS:
; OUTPUTS:
;  name - string with trailing slash
; KEYWORD OUTPUT PARAMETERS:
; COMMON BLOCKS:
; SIDE EFFECTS:
; RESTRICTIONS:
; PROCEDURE:
; MODIFICATION HISTORY:
;  95/06/08, Written by Marc W. Buie, Lowell Observatory
;  99/02/08, MWB, added clause for Mac style directories.
;-
function addslash,in_name
;   if badpar(in_name,7,0,CALLER='ADDSLASH: (name) ') then return,''
   name = in_name
   if strlen(name) eq 0 then name = './'
   if !version.os_family eq 'Windows' then begin
      if (strmid(name,strlen(name)-1,1) ne '/') and $
         (strmid(name,strlen(name)-1,1) ne '\') then name = name + '\'
   endif else if !version.os_family eq 'MacOS' then begin
      if strmid(name,strlen(name)-1,1) ne ':' then name = name + ':'
   endif else begin
      if strmid(name,strlen(name)-1,1) ne '/' then name = name + '/'
   endelse
   return,name
end
