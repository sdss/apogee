function getplatedata,plate,cmjd,plugid=plugid,asdaf=asdaf,mapa=mapa,obj1m=obj1m,fixfiberid=fixfiberid,noobject=noobject,stop=stop,skip=skip,twilight=twilight,badfiberid=badfiberid,apogees=apogees,mapper_data=mapper_data,starfiber=starfiber

; getplatedata loads up a structure with plate information and information about the 300 APOGEE fibers
;  This is obtained from a plPlugMapA file or from a 
;  plPlugMapM+plateHolesSorted combination
; returned structure includes:
;    fiberid, ra, dec, eta, zeta, hmag, objtype, obj (name)
;  for each of 300 APOGEE (spectrographid=2) files

dirs=getdir(apodir,datadir=datadir)
if n_elements(mapper_data) eq 0 then mapper_data=dirs.mapperdir
if size(plate,/type) eq 7 then begin
  cplate=plate 
  platenum=0L
endif else begin
  cplate=strtrim(string(format='(i6.4)',plate),2)
  platenum=long(plate)
endelse
mjd=0L
reads,cmjd,mjd

; create the output fiber structure
tmp={fiberid: 0, ra: 0.d0, dec: 0.d0, eta: 0.d0, zeta: 0.d0, hmag: 0., objtype: 'none', holetype: 'OBJECT', object: '', tmass_style: '', target1: 0L, target2: 0L, target3: 0L, spectrographid: 2, mag: fltarr(5), ak_targ: -99., ak_targ_method : 'none', ak_wise: -99., sfd_ebv: -99.}
guide=replicate(tmp,16)
loc=0L
if keyword_set(obj1m) then begin
  if keyword_set(fixfiberid) then begin
    if fixfiberid eq 1 then begin
      fiberid=[218,220,222,223,226,227,228,229,230,231]
      if ~keyword_set(starfiber) then starfiber=229
    endif
  endif else begin
    fiberid=[218,219,221,223,226,228,230]
    if ~keyword_set(starfiber) then starfiber=223
  endelse
  ;fiber=replicate(tmp,9)
  fiber=replicate(tmp,n_elements(fiberid))
  platedata={plate: platenum, mjd: mjd, plateid: cplate, locationid: 1L, field: ' ', programname: '', $
             cmjd: cmjd, ha: [-99.,-99.,-99.], fiberdata: fiber, guidedata: guide}
  platedata.field = cplate
  ;fiber.fiberid=[218,219,220,221,223,226,228,230,231]
  ;fiber.objtype=['SKY','SKY','SKY','SKY','STAR','SKY','SKY','SKY','SKY']
  ;if n_elements(fixfiberid) gt 0 then begin
  ;  if fixfiberid eq 1 then begin
  ;    fiber.fiberid=[218,220,222,223,226,227,228,229,230,231]
  ;    fiber.objtype=['SKY','SKY','SKY','SKY','SKY','SKY','SKY','SKY','SKY','SKY']
  ;    if ~keyword_set(starfiber) then starfiber=229
  ;  endif
  ;endif else begin
  ;  fiber.fiberid=[218,219,221,223,226,228,230]
  ;  fiber.objtype=['SKY','SKY','SKY','STAR','SKY','SKY','SKY']
  ;  fiber.objtype=['SKY','SKY','SKY','SKY','SKY','SKY','SKY']
  ;  if ~keyword_set(starfiber) then starfiber=223
  ;endelse
  fiber.fiberid=fiberid
  fiber.objtype=replicate('SKY',n_elements(fiberid))
  j=where(fiber.fiberid eq starfiber)
  fiber[j].objtype = 'STAR'
  j=where(fiber.objtype eq 'SKY')
  fiber[j].target2 = 2L^4
  obj=mrdfits(getenv('APOGEEREDUCEPLAN_DIR')+'/data/1m/'+cplate+'.fits',1,status=status)

  if status eq 0 then begin
   j=where(strtrim(obj.name,2) eq strtrim(obj1m,2),nj)
   if nj gt 0 then begin
    ifiber=where(fiber.fiberid eq starfiber)
    fiber[ifiber].object = obj[j].name
    fiber[ifiber].tmass_style = strtrim(obj[j].name,2)
    fiber[ifiber].hmag = obj[j].h
    fiber[ifiber].mag[1] = obj[j].h
    fiber[ifiber].ra = obj[j].ra
    fiber[ifiber].dec = obj[j].dec
    fiber[ifiber].target2 = 2L^22
   endif
  endif else stop,'halt: no file found with object information!'
  platedata.fiberdata = fiber  
  return,platedata
endif
if keyword_set(twilight) then begin
  fiber=replicate(tmp,300)
  platedata={plate: platenum, mjd: mjd, plateid: cplate, locationid: 1L, field: ' ', programname: '', cmjd: cmjd, ha: [-99.,-99.,-99.], fiberdata: fiber, guidedata: guide}
  platedata.field = cplate
  ;j=where(indgen(300) mod 10 eq 0)
  ;fiber[j].target2 = 2L^4
  ;fiber[j].objtype = 'SKY'
  fiber.hmag=10.
  fiber.mag=[10.,10.,10.,10.,10.]
  fiber.objtype='STAR'
  fiber.fiberid = indgen(300)+1
  platedata.fiberdata = fiber  
  return,platedata
endif
fiber=replicate(tmp,300)
reads,cplate,platenum
platedata={plate: platenum, mjd: mjd, plateid: cplate, locationid: 0L, field: ' ', programname: '', cmjd: cmjd, ha: [-99.,-99.,-99.], fiberdata: fiber, guidedata: guide}
platedata.field = apogee_field(loc,platenum,survey,programname)
platedata.locationid = loc
platedata.programname = programname

; do we want to use a plPlugMapA file with the matching already done?
havematch=0
if keyword_set(mapa) then root = 'plPlugMapA' else root='plPlugMapM'
if keyword_set(plugid) then begin
  tmp=file_basename(plugid,'.par')
  if strpos(tmp,'plPlug') ge 0 then tplugid=strmid(tmp,11) else tplugid=tmp
  plugfile=root+'-'+tplugid+'.par'
endif else begin
  tplugid=root+'-'+cplate
  plugfile=root+'-'+cplate+'.par'
endelse
if keyword_set(mapa) then plugdir=datadir+cmjd+'/' else begin
  tmp=strsplit(tplugid,'-',/extract)
  plugmjd=tmp[1]
  plugdir=mapper_data+'/'+plugmjd+'/'
endelse

; Does the plugfile exist? If so, load it
if file_test(plugdir+'/'+plugfile) then aploadplugmap,plugdir+'/'+plugfile,plugmap,fixfiberid=fixfiberid else $
   if keyword_set(skip) then return,0 else stop,'halt: cannot find plugmap file '+plugdir+'/'+plugfile

platedata.locationid = plugmap.locationid
platedata.ha[0] = plugmap.ha[0]
platedata.ha[1] = plugmap.ha_observable_min[0]
platedata.ha[2] = plugmap.ha_observable_max[0]
if not keyword_set(mapa) then begin
  ; get the plateHolesSorted file for thie plate and read it
  ;reads,cplate,plate
  platestr=string(format='(i6.6)',platenum)
  platedir=getenv('PLATELIST_DIR')+'/plates/'+strmid(platestr,0,4)+'XX/'+platestr
  holefile='plateHolesSorted-'+platestr+'.par'
  print,'yanny_read,'+platedir+'/'+holefile
  yanny_read,platedir+'/'+holefile,pdata,/anon,hdr=hdr
  p=*pdata
  yanny_free,pdata
  ; use locationid from plateHoles files as there are a few cases where plugmapM is
  ;  wrong
  j=where(strpos(strupcase(hdr),'LOCATIONID') ge 0)
  tmp=strsplit(hdr[j],/ext)
  loc=0L
  reads,tmp[1],loc
  platedata.locationid = loc
endif

; load guide stars
for i=0,15 do begin
  m=where(plugmap.fiberdata.holetype eq 'GUIDE' and $
          plugmap.fiberdata.fiberid eq i,nm)
  guide[i].fiberid = plugmap.fiberdata[m].fiberid
  guide[i].ra = plugmap.fiberdata[m].ra 
  guide[i].dec = plugmap.fiberdata[m].dec 
  guide[i].eta = plugmap.fiberdata[m].eta
  guide[i].zeta = plugmap.fiberdata[m].zeta
  guide[i].spectrographid = plugmap.fiberdata[m].spectrographid
endfor
platedata.guidedata = guide

; find matching plugged entry for each spectrum and load up the output information from correct source(s)
;if dirs.instrument eq  'apogee-s' and mjd lt 57810 then plugmap.fiberdata.fiberid += 1
for i=0,299 do begin
  fiber[i].spectrographid = -1
  m=where(plugmap.fiberdata.holetype eq 'OBJECT' and $
          plugmap.fiberdata.spectrographid eq 2 and $
          plugmap.fiberdata.fiberid eq 300-i,nm)
  if keyword_set(badfiberid) then begin
    j=where(badfiberid eq 300-i,nbad)
    if nbad gt 0 then begin
      print,'fiber index ',i,' declared as bad'
      nm=0
    endif 
  endif
  if nm gt 1 then begin
      print,'halt: more than one match for fiber id !! MARVELS??'
      print,plugmap.fiberdata[m].fiberid,plugmap.fiberdata[m].primtarget,plugmap.fiberdata[m].sectarget
      stop
  endif
  if nm eq 1 then begin
    fiber[i].fiberid = plugmap.fiberdata[m].fiberid
    fiber[i].ra = plugmap.fiberdata[m].ra 
    fiber[i].dec = plugmap.fiberdata[m].dec 
    fiber[i].eta = plugmap.fiberdata[m].eta
    fiber[i].zeta = plugmap.fiberdata[m].zeta
    fiber[i].target1 = plugmap.fiberdata[m].primTarget
    fiber[i].target2 = plugmap.fiberdata[m].secTarget
    fiber[i].spectrographid = plugmap.fiberdata[m].spectrographid

    if keyword_set(asdaf) then begin
      ; special for asdaf object plates
      if 300-i eq asdaf then begin
        fiber[i].objtype = 'STAR'
        fiber[i].hmag = 0.
      endif else begin
        fiber[i].objtype = 'SKY'
        fiber[i].hmag = -99.999
      endelse
    endif else begin
      fiber[i].objtype = plugmap.fiberdata[m].objtype
      ; fix up objtype
      fiber[i].objtype = 'STAR'
      fiber[i].holetype = plugmap.fiberdata[m].holetype
      if keyword_set(mapa) then begin
        ; HMAG's are correct from plPlugMapA files
        fiber[i].hmag = plugmap.fiberdata[m].mag[1]
        fiber[i].object = strtrim(plugmap.fiberdata[m].tmass_style,2)
        fiber[i].tmass_style = strtrim(plugmap.fiberdata[m].tmass_style,2)
        if is_bit_set(fiber[i].sectarget,9) eq 1 then fiber[i].objtype='HOT_STD'
        if is_bit_set(fiber[i].sectarget,4) eq 1 then fiber[i].objtype='SKY'
      endif else begin
        ; get matching stars from coordinate match
        match = where(abs(p.target_ra-fiber[i].ra) lt 0.00002 and $
                      abs(p.target_dec-fiber[i].dec) lt 0.00002,nmatch)
        if nmatch gt 0 then begin
          if tag_exist(p,'apogee2_target1') and platenum gt 7500 then begin
            fiber[i].target1 = p[match].apogee2_target1
            fiber[i].target2 = p[match].apogee2_target2
            fiber[i].target3 = p[match].apogee2_target3
            apogee2=1
          endif else begin
            fiber[i].target1 = p[match].apogee_target1
            fiber[i].target2 = p[match].apogee_target2
            apogee2=0
          endelse
          if is_bit_set(fiber[i].target2,9) eq 1 then fiber[i].objtype='HOT_STD'
          if is_bit_set(fiber[i].target2,4) eq 1 then begin
            object='SKY' 
            hmag=99.99
            fiber[i].mag = [hmag,hmag,hmag,hmag,hmag]
            fiber[i].objtype='SKY'
          endif else begin
            tmp=strtrim(p[match].targetids,2)
            len=strlen(tmp)
            object=strmid(tmp,len-16)
            if strpos(tmp,'A') eq 0 then $
             object='AP'+object else object='2M'+object
            hmag = p[match].tmass_h
            fiber[i].mag = [p[match].tmass_j,p[match].tmass_h,p[match].tmass_k,0.,0.]
            ; adopt PM un-adjusted  coordinates
            fiber[i].ra-=p[match].pmra/1000./3600./cos(fiber[i].dec*!pi/180.)*(p[match].epoch-2000.)
            fiber[i].dec-=p[match].pmdec/1000./3600.*(p[match].epoch-2000.)
          endelse
          fiber[i].hmag = hmag
          fiber[i].object = object
          fiber[i].tmass_style = object
        endif else stop,'no match found in plateHoles!',fiber[i].ra,fiber[i].dec, i
      endelse
    endelse
  endif else begin
    fiber[i].fiberid = -1
    print,'no match for fiber index',i
  endelse
endfor

; load apogeeObject file to get proper name and coordinates
; get apogeeObject catalog info for this field
if apogee2 then apogeeobject='apogee2Object' else apogeeobject='apogeeObject'
if ~keyword_set(noobject) then begin
 targetdir=getenv('APOGEE_TARGET')
 objectfile=targetdir+'/'+apogeeobject+'/'+apogeeobject+'_'+strtrim(apogee_field(platedata.locationid,platenum),2)+'.fits'
 if not file_test(objectfile) then begin
  print,'cant find apogeeObject file: ', objectfile
 endif else begin
  print,'reading apogeeObject file: ',objectfile
  objects=mrdfits(objectfile,1)
  spherematch,objects.ra,objects.dec,fiber.ra,fiber.dec,10./3600.,match1,match2,dist,maxmatch=1
  for i=0,299 do begin
   if fiber[i].objtype eq 'STAR' or fiber[i].objtype eq 'HOT_STD' then begin
     j=where(match2 eq i,nj)
     if nj gt 0 then begin
      if strtrim(fiber[i].object,2) ne strtrim(objects[match1[j]].apogee_id) then begin
        print,'apogeeObject differs from plateHoles: '
        print,fiber[i].object+' '+objects[match1[j]].apogee_id
        print,fiber[i].ra,' ',objects[match1[j]].ra
        print,fiber[i].dec,' ',objects[match1[j]].dec
      endif
      fiber[i].tmass_style=objects[match1[j]].apogee_id
      fiber[i].ra=objects[match1[j]].ra
      fiber[i].dec=objects[match1[j]].dec
      if finite(objects[match1[j]].ak_targ) then fiber[i].ak_targ=objects[match1[j]].ak_targ
      fiber[i].ak_targ_method=objects[match1[j]].ak_targ_method
      if finite(objects[match1[j]].ak_wise) then fiber[i].ak_wise=objects[match1[j]].ak_wise
      if finite(objects[match1[j]].sfd_ebv) then fiber[i].sfd_ebv=objects[match1[j]].sfd_ebv
     endif else print,'not halted, but no match in object file for ',fiber[i].object
   endif 
  endfor
 endelse
endif

platedata.fiberdata = fiber
if keyword_set(stop) then stop
return,platedata
end
