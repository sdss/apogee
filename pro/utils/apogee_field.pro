function apogee_field,loc,lplate,survey,programname,addloc=addloc

common com_apogeereduce_field, plans

dirs=getdir()

if(n_tags(plans) eq 0) then $
   plans= yanny_readone(getenv('PLATELIST_DIR')+'/platePlans.par')


; 1m case
if dirs.telescope eq 'apo1m' then begin
  survey='apo1m'
  programname=lplate
  if loc eq 1 and size(lplate,/type) eq 7 then return,lplate
  if loc eq 0 and size(lplate,/type) eq 7 then return,lplate
endif

if loc eq 0 and lplate eq 0 then return,''

; special error case from bad plugmap file
if loc eq 4386 and lplate eq 5915 then loc=4368

f=where(plans.plateid gt 4800 and (loc eq 0 or plans.locationid eq loc) and (lplate eq 0 or plans.plateid eq lplate),nf)

if nf gt 0 then begin

   for i=0,nf-1 do begin
     cnames= plans[f[i]].name
     survey = plans[f[i]].survey
     programname = plans[f[i]].programname
     if survey eq 'manga-apogee2' then cnames = plans[f[i]].comments

     ;; take only first pointing
     words= strsplit(cnames, /extr)
     cname= words[0]

     if survey eq 'manga-apogee2' and keyword_set(addloc) then cname = cname+'_loc'+strtrim(string(format='(i4)',loc),2)

     ;; remove APG_ if necessary
     if(strmatch(cname, 'APG_*')) then $
        cname= strmid(cname, 4)

     if(strmatch(cname, 'APGS_*')) then $
        cname= strmid(cname, 5)

     case cname of
        'BULGE_04+00': cname='004+00'
        'HALO_M13': cname='M13'
        'DISK_27+00': cname='027+00'
        'DISK_30+08': cname='030+08'
        else:
     endcase
     if (i gt 0) then begin
       print, 'More than one match for plate: '+strtrim(string(lplate),2)+' loc: '+strtrim(string(loc),2) +$
              ' in platePlans.par!'
       ; location 5042 has both apogee2 and apogee2s in platePlans
       if loc eq 5042 then survey='apogee2'
       if cname ne cname0 then print,'Field names do not match!',plans[f].name, cname0, cname
     endif
     cname0=cname
   endfor
   return, cname 

endif else begin
   return,'Unknown'
endelse

end

