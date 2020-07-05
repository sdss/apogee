function struct_rename_tags, struct, oldtags, newtags, addtag=addtag, addval=addval
;Change the names of existing structure fields
;
;Inputs:
;
;Input/Output:
; struct (structure) structure to be modified.
;
; J. Holtzman, June 2013

if n_params() lt 3 then begin
  print, 'syntax: struct_rename_tags, struct, tags, newtags'
  return,-1
endif

if n_elements(oldtags) ne n_elements(newtags) then begin
  print,'Number of new tags must match number of old tags'
  stop
endif

;Get list of structure tags.
tags = tag_names(struct)
ntags = n_elements(tags)

; check to see if oldtags exist in structure
for i=0,n_elements(oldtags)-1 do  begin
  oldtags[i] = strupcase(strtrim(oldtags[i],2))
  newtags[i] = strupcase(strtrim(newtags[i],2))
  j=where(tags eq oldtags[i],nj)
  if nj eq 0 then print, oldtags[i]+' does not appear in existing structure'
endfor

;Check that input is a structure.
if size(struct, /tname) ne 'STRUCT' then begin
  message, 'first argument is not a structure'
endif

; create the base structure with the new names
for itag=0,ntags-1 do begin
  j=where(oldtags eq strupcase(strtrim(tags[itag],2)),nmatch)
  if itag eq 0 then $
   if nmatch eq 0 then new = create_struct(tags[itag], struct[0].(itag)) else new=create_struct(newtags[j], struct[0].(itag)) $
  else $
   if nmatch eq 0 then new = create_struct(new, tags[itag], struct[0].(itag)) else new=create_struct(new, newtags[j], struct[0].(itag))
endfor

; add new tags if requested
if n_elements(addtag) eq 1 then new = create_struct(new, addtag, addval)
if n_elements(addtag) gt 1 then for i=0,n_elements(addtag)-1 do new = create_struct(new, addtag[i], addval[i])

; populate all structure elements
out=replicate(new,n_elements(struct))
for i=0,n_elements(struct)-1 do for itag=0,ntags-1 do out[i].(itag)=struct[i].(itag)

return,out

end

