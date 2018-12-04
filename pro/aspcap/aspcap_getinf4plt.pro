pro aspcap_getinf4plt,file,objtype,objid,vrad,vhel,vraderr,llongi,llati,Jmag,Hmag,Kmag

pepe=readfits(file,ext=0,head,/silent)
objtype=strcompress(strmid(aspcap_extrackeyw(head,'OBJTYP'),1,8),/remove_all)
objid=strcompress(strmid(aspcap_extrackeyw(head,'OBJID'),1,strlen(aspcap_extrackeyw(head,'OBJID'))),/remove_all)
if strmid(objid,0,4) eq '2M2M' then objid=strmid(objid,2,strlen(objid))
vrad=double(aspcap_extrackeyw(head,'VRAD'))
vhel=double(aspcap_extrackeyw(head,'VHELIO'))
vraderr=double(aspcap_extrackeyw(head,'VERR'))
llongi=double(aspcap_extrackeyw(head,'GLON'))
llati=double(aspcap_extrackeyw(head,'GLAT'))
;Hmag=double(aspcap_extrackeyw(head,'H'))
;Jmag=double(aspcap_extrackeyw(head,'J'))
;Kmag=double(aspcap_extrackeyw(head,'K'))
Hmag=sxpar(head,'H')
Jmag=sxpar(head,'J')
Kmag=sxpar(head,'K')
end
