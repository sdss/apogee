function apsetpar,planstr,tag,default

; function to get parameter from structure, or set default

if tag_exist(planstr,tag) then begin
  tnames=tag_names(planstr)
  tindex=where(strcmp(tnames,strupcase(tag)) eq 1)
  val=planstr.(tindex)
endif
if n_elements(val) eq 0 then val=default

return,val
end

