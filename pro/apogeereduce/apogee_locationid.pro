function apogee_locationid,field,survey

common com_apogeereduce_field, plans

if(n_tags(plans) eq 0) then $
   plans= yanny_readone(getenv('PLATELIST_DIR')+'/platePlans.par')

if n_elements(survey) gt 0 then begin
  if survey eq 'manga-apogee2' then $
    f=where(plans.comments eq field and plans.survey eq survey,nf) $
  else $
    f=where(plans.name eq field and plans.survey eq survey,nf)
  print,nf,' fields found matching ',field,' and ',survey
endif else begin
  f=where(plans.name eq field and (strpos(plans.survey,'apogee') ge 0),nf)
  if nf eq 0 then f=where(plans.name eq 'APG_'+field and (strpos(plans.survey,'apogee') ge 0),nf)
  if nf eq 0 then f=where(plans.name eq 'APGS_'+field and (strpos(plans.survey,'apogee') ge 0),nf)
  if nf eq 0 then f=where(plans.comments eq field and (strpos(plans.survey,'apogee') ge 0),nf)
  print,nf,' fields found matching ',field
  if nf eq 0 then stop
  if nf gt 0 then survey=plans[f[0]].survey
endelse

if nf gt 0 then return,plans[f[0]].locationid 

; if we didn't find it, maybe we were given a locationid as a character string, so return the int (old style)
; this will crash if it's really a character string without a match!
return,long(field)

end

