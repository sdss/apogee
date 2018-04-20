pro apaddpar,frame,name,value,comment,history=history,stp=stp

;+
;
; APADDPAR
;
; This adds keywords to the chip headers in an APOGEE FRAME
; structure.  It behaves almost the same as SXADDPAR
; and SXADDHIST
;
; For some reason SXADDPAR can't operate on the header
; in the structure, so we need to copy it to a string array,
; run SXADDPAR, and then copy it back.  The header string
; array in the frame structure must be large enough to
; handle the new header
;
; If an array of frame structure are input, then only
; the first element (frame) will be updated.
;
; INPUTS:
;  frame     The APOGEE frame structure
;  name      The name of the parameter to add.  If /history is
;              set then this is history string.
;  value     The value for the parameter.
;  comment   Comment field.  If /history is set then setting
;              /comment will cause "COMMENT " to be written
;              to the header instead of "HISTORY ". 
;  /history  This should be added as a history comment line
;  /stp      Stop at the end of the program.
;
; OUTPUTS:
;  The chip headers in the frame structure are updated with
;  the appropriate keyword.
;
; USAGE:
;  IDL>apaddpar,frame,'GAIN',1.0
;   or
;  IDL>adaddpar,frame,'I DID THIS',/history
;
; By D.Nidever  July 2010
;-

; Not enough inputs
if n_elements(frame) eq 0 then begin
  print,'Syntax - apaddpar,frame,name,value,history=history,stp=stp'
  return
endif

; Check that this is a structre
type = size(frame,/type)
if type ne 8 then begin
  print,'FRAME MUST be a structure'
  return
endif

; Check that there are three tags
tags = tag_names(frame)
if n_elements(tags) ne 3 then begin
  print,'FRAME must have three chip sub-structure'
  return
endif


; Loop through the Chips
;------------------------
FOR i=0,2 do begin

  ; Check that header tag exists
  tags2 = tag_names(frame.(i))
  indhd = where(tags2 eq 'HEADER',nindhd)
  if nindhd eq 0 then begin
    print,'HEADER array NOT FOUND'
    return
  endif

  ; Get the original header
  head = frame[0].(i).header
  nhead0 = n_elements(head)

  ; Remove trailing blank lines
  indend = where(stregex(head,'^END',/boolean) eq 1,nindend)
  if indend[0] eq -1 then indend=nhead0-1
  head = head[0:indend[0]]


  ; REGULAR keyword parameter to add
  ;---------------------------------
  if not keyword_set(history) then begin
    if keyword_set(comment) then SXADDPAR,head,name,value,comment else $
    SXADDPAR,head,name,value

  ; History parameter
  endif else begin
    SXADDHIST,name,head,comment=comment
  endelse


  ; Check that the array in the frame structure is large enough
  nhead = n_elements(head)
  if nhead gt nhead0 then begin
    print,'New header is longer than structure header array'
    return
  endif

  ; Stuff it back in
  frame[0].(i).header = ''   ; wipe it clean
  frame[0].(i).header[0:nhead-1] = head

ENDFOR

if keyword_set(stp) then stop

end
