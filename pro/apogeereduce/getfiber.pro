function getfiber,plate,cmjd,plugid=plugid,single=single,mapa=mapa,write=write

; fiberid loads up and returns structure with information about 300 APOGEE
;  fibers. This is obtained from a plPlugMapA file or from a 
;  plPlugMapM+plateHolesSorted combination
; returned structure includes:
;    fiberid, ra, dec, eta, zeta, hmag, objtype, obj (name)
;  for each of 300 APOGEE (spectrographid=2) files
; with /write, writes a plPlugMapA file from plPlugMapM+plateHoles

getdir,apodir,datadir=datadir

; do we want to use a plPlugMapA file with the matching already done?
havematch=0
cplate=strtrim(string(format='(i6.4)',plate),2)

if keyword_set(plugid) and keyword_set(mapa) then $
  plugfile='plPlugMapA-'+plugid+'.par' else $
  plugfile='plPlugMapA-'+cplate+'.par'

; Does the plPlugMapA file exist? If so, load it
if file_test(datadir+cmjd+'/'+plugfile) then begin
  print,'usingplPlugMapA'
  aploadplugmap,datadir+cmjd+'/'+plugfile,plugmap
  havematch=1  
endif else begin
  ; otherwise, get the plateHolesSorted file for thie plate and read it
  spawn,'cp '+datadir+'/plates/*/'+string(format='(i6.6)',plate)+'/plateHolesSort* '+datadir+cmjd
  holefile=string(format='("plateHolesSorted-",i6.6,".par")',plate)
  print,'yanny_read,'+datadir+cmjd+'/'+holefile
  yanny_read,datadir+cmjd+'/'+holefile,pdata
  p=*pdata
  yanny_free,pdata
  ;p=yanny_readone(datadir+cmjd+'/'+holefile)

  ; get the plPlugMapM file 
  if keyword_set(plugid) then $
    plugfile='plPlugMapM-'+plugid+'.par' else $
    plugfile='plPlugMapM-'+cplate+'.par'
  aploadplugmap,datadir+cmjd+'/'+plugfile,plugmap

  ; with write keyword, construct a plPlugMapA file using info from plateHoles
  if keyword_set(write) then begin
    if keyword_set(plugid) then $
      outfile='plPlugMapN-'+plugid+'.par' else $
      outfile='plPlugMapN-'+cplate+'.par'
    yanny_read,datadir+cmjd+'/'+plugfile,plugdata,hdr=hdr,enums=enums,structs=structs,stnames=stnames
    ptmp=*plugdata
    new=create_struct(ptmp[0],'tmass_style','-')
    out=replicate(new,n_elements(ptmp))
    tags=tag_names(ptmp)
    ntags=n_elements(tags)
    skymask = 16L
    hotmask = 512L
    extmask = 1024L
    starmask = skymask or hotmask

    for i=0,n_elements(ptmp)-1 do begin
      for itag=0,ntags-1 do out[i].(itag)=ptmp[i].(itag)
      out[i].tmass_style='-'
      if ptmp[i].holetype eq 'OBJECT' and $
         ptmp[i].spectrographid eq 2 then begin
        match = where(abs(p.target_ra-ptmp[i].ra) lt 0.00001 and $
                      abs(p.target_dec-ptmp[i].dec) lt 0.00001)
        if match ge 0 then begin
          if (out[i].sectarget and skymask) gt 0 then out[i].objtype = 'SKY' else begin
            name = strtrim(p[match].targetids,2)
            out[i].tmass_style = strtrim('2M'+strmid(p[match].targetids,7),2)
            out[i].mag[0] = p[match].tmass_j
            out[i].mag[1] = p[match].tmass_h
            out[i].mag[2] = p[match].tmass_k
            if (out[i].sectarget and hotmask) gt 0 then out[i].objtype = 'HOT_STD'
            if (out[i].primtarget and extmask) gt 0 then out[i].objtype = 'EXTOBJ'
            if (out[i].sectarget and starmask) eq 0 and (out[i].primtarget and extmask) eq 0 then out[i].objtype = 'STAR'
          endelse
        endif
      endif 
    endfor
    *plugdata=out
    ns=n_elements(structs)
    s=[structs[0:ns-2],' char tmass_style[30];',structs[ns-1]]

    yanny_write,outfile,plugdata,hdr=hdr,enums=enums,structs=s,stnames=stnames
    yanny_free,plugdata
  endif
endelse

; create the output fiber structure
tmp={fiberid: 0, ra: 0.d0, dec: 0.d0, eta: 0.d0, zeta: 0.d0, hmag: 0., objtype: 'none', object: '', primtarg: 0L, sectarg: 0L}
fiber=replicate(tmp,300)

; find matching plugged entry for each spectrum and load up the output information from correct source(s)
for i=0,299 do begin
  m=where(plugmap.fiberdata.holetype eq 'OBJECT' and $
          plugmap.fiberdata.spectrographid eq 2 and $
          plugmap.fiberdata.fiberid eq 300-i)
  if m[0] ge 0 then begin
    if n_elements(m) gt 1 then begin
      print,'more than one match for fiber id !! MARVELS??'
      print,plugmap.fiberdata[m].fiberid,plugmap.fiberdata[m].primtarget,plugmap.fiberdata[m].sectarget
      stop
    endif
    fiber[i].fiberid = plugmap.fiberdata[m[0]].fiberid 
    fiber[i].ra = plugmap.fiberdata[m[0]].ra 
    fiber[i].dec = plugmap.fiberdata[m[0]].dec 
    fiber[i].eta = plugmap.fiberdata[m[0]].eta
    fiber[i].zeta = plugmap.fiberdata[m[0]].zeta
    fiber[i].primtarg = plugmap.fiberdata[m[0]].primTarget
    fiber[i].sectarg = plugmap.fiberdata[m[0]].secTarget

    if keyword_set(single) then begin
      ; special for single object plates
      if 300-i eq single then begin
        fiber[i].objtype = 'STAR'
        fiber[i].hmag = 0.
      endif else begin
        fiber[i].objtype = 'SKY'
        fiber[i].hmag = -99.999
      endelse
    endif else begin
      fiber[i].objtype = plugmap.fiberdata[m[0]].objtype
      if havematch then begin
        ; HMAG's are correct from plPlugMapA files
        fiber[i].hmag = plugmap.fiberdata[m[0]].mag[1]
        fiber[i].object = strtrim(plugmap.fiberdata[m[0]].tmass_style,2)
      endif else begin
        ; get matching stars from coordinate match
        match = where(abs(p.target_ra-fiber[i].ra) lt 0.00001 and $
                    abs(p.target_dec-fiber[i].dec) lt 0.00001)
        if match ge 0 then begin
          fiber[i].hmag = p[match].tmass_h
          fiber[i].object = strtrim(p[match].targetids,2)
        endif
      endelse
    endelse
  endif else begin
    fiber[i].fiberid = -1
    print,'no match for fiber ',i
  endelse
endfor

return,fiber
end
