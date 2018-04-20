pro aploadplugmap,plugfile,plugmap,verbose=verbose,silent=silent,$
                  stp=stp,error=error,fixfiberid=fixfiberid

;+
;
; APLOADPLUGMAP
;
; This program loads an APOGEE plugmap file.
;
; INPUTS:
;  plugfile  The absolute path of the plugmap file
;
; OUTPUTS:
;  plugmap   The plugmap structure with all the
;              relevant information
;  /verbose  Print a lot of information to the screen
;  /silent   Don't print anything to the screen
;  /stp      Stop at the end of the program
;  =error    The error message if one occurred
;
; USAGE;
;  IDL>aploadplugmap,plugfile,plugmap
;
; By D.Nidever  May 2010
;-

apgundef,error,plugmap

; Not enough inputs
nplugfile = n_elements(plugfile)
if nplugfile eq 0 then begin
  if not keyword_set(silent) then begin
    print,'Syntax - aploadplugmap,plugfile,plugmap,verbose=verbose,silent=silent,'
    print,'                       stp=stp,error=error'
  endif
  error = 'Not enough inputs'
  return
endif

; More than one file
if nplugfile gt 1 then begin
  error = 'More than one plugmap file. Only ONE allowed'
  if not keyword_set(silent) then print,error
  return
endif


; Check that the plug file exists
if file_test(plugfile) eq 0 then begin
  error = plugfile+' NOT FOUND'
  if not keyword_set(silent) then print,error
  return
endif

;----------
; Strip path from plug file name, and change to that directory
thisplug = file_basename(plugfile)
thispath = file_dirname(plugfile)
CD, thispath, current=origdir

; There is a READPLUGMAP.PRO in idlspec2d, but most of the stuff
; it does is unnecessary for APOGEE.

;----------
; Find the PLUGMAPOBJ structure
;  the loading is pretty slow
YANNY_READ, thisplug, pdata, hdr=hdr, errcode=errcode ;,/anonymous
if errcode ne 0 then begin
  error = 'ERROR LOADING '+plugfile
  if not keyword_set(silent) then print,error
  return
endif
if n_elements(pdata) eq 0 then begin
  error = 'No data in the plugmap file'
  if not keyword_set(silent) then print,error
  return
endif
for i=0, N_elements(pdata)-1 do begin
   if (tag_names(*pdata[i], /structure_name) EQ 'PLUGMAPOBJ') then $
    plugdata = *pdata[i]
endfor
YANNY_FREE, pdata

if (N_elements(plugdata) EQ 0) then begin
  stop
  error = 'ABORT: No PLUGMAPOBJ structure in plugmap file ' + plugfile
  if not keyword_set(silent) then print,error
  CD, origdir
  return
endif
CD,origdir

; Add ETA/ZETA to plugmap structure
;-----------------------------------
ADD_TAG,plugdata,'ZETA',0.0d0,plugdata
ADD_TAG,plugdata,'ETA',0.0d0,plugdata


; Make the PLUGMAP structure
;----------------------------
;plugmap = {hdr:hdr,fiberdata:plugdata,reddeningMed:0.0,tileid:'',raCen:0.0d0,decCen:0.0d0,plateid:'',$
;           temp:0.0,hamin:0.0,hamax:0.0,mjd:''}
plugmap = {hdr:hdr,fiberdata:plugdata,mjd:''}
tags = TAG_NAMES(plugmap)
apgundef,hdr,plugdata


; This explains the plPlugMap format in more detail.
; http://www.sdss.org/dr7/dm/flatFiles/plPlugMap.html

; Get the header keywords
gdline = where(strtrim(strmid(plugmap.hdr,0,1),2) ne '#' and $
               strtrim(plugmap.hdr,2) ne '',ngdline)
keywords = strarr(ngdline)
keytype = lonarr(ngdline)
for i=0,ngdline-1 do begin
  line = strtrim(plugmap.hdr[gdline[i]],2)
  dum = strsplit(line,' ',/extract)
  keywords[i] = idl_validname(dum[0],/convert_all)
endfor


;----------
; Find keywords from the header
;keywords = ['reddeningMed','tileId','raCen','decCen','plateId','locationId','temp','haMin','haMax']
;keytype = [4,7,5,5,7,7,4,4,4]
nkeywords = n_elements(keywords)
for j=0,nkeywords-1 do begin
 if j gt 0 then dup = where(keywords[0:j-1] eq keywords[j]) else dup=-1
 if dup lt 0 then begin

  ; NOTE, the keywords are CASE SENSITIVE
  value = YANNY_PAR(plugmap.hdr,keywords[j],count=count)

  ; NO match, check without case sensitivity
  if count eq 0 then begin
    ; must be at the beginning and have a space afterwards
    ind = where(stregex(strtrim(plugmap.hdr,2),'^'+keywords[j]+' ',/boolean,/fold_case) eq 1,nind)
    if nind gt 0 then begin
      line = strtrim(plugmap.hdr[ind[0]],2)
      len = strlen(keywords[j])
      value = strtrim(strmid(line,len),2)
      count = 1
    end
  endif

  ; No match
  if count eq 0 then begin
    ;error = 'KEYWORD '+keywords[j]+' NOT FOUND IN '+plugfile
    ;if not keyword_set(silent) then print,error
    ;return
    if not keyword_set(silent) then print,'KEYWORD '+keywords[j]+' NOT FOUND IN '+plugfile
    value='missing'
  endif

  ; Get type
  len = strlen(strtrim(value,2))
  type = 7
  if min(valid_num(value)) eq 1 and max(len) lt 7 then type=4
  if min(valid_num(value)) eq 1 and max(len) ge 7 then type=5
  if min(valid_num(value,/integer)) eq 1 then type=3   ; use min in case value is an array

  ;type = keytype[j]
  plugmap = CREATE_STRUCT(plugmap,keywords[j],fix(value,type=type))
  tags = tag_names(plugmap)
  ;tagind = where(tags eq strupcase(keywords[j]))
  ;type = size(plugmap.(tagind),/type)
  ;plugmap.(tagind) = fix(value,type=type)
  ; Print out the keyword/value pair
  if keyword_set(verbose) then print,keywords[j],' = ',fix(value,type=type)
 endif
end

; Get ETA/ZETA using plate center
;---------------------------------
ROTSPHCEN,plugmap.fiberdata.ra,plugmap.fiberdata.dec,plugmap.raCen,plugmap.decCen,zeta,eta,/gnomic
plugmap.fiberdata.zeta = zeta
plugmap.fiberdata.eta = eta

; Fix bit 6 erroneously set in plugmap files for tellurics 
ind = where(plugmap.fiberdata.spectrographid eq 2 and $
                     plugmap.fiberdata.holetype eq 'OBJECT' and $
                     plugmap.fiberdata.objtype eq 'HOT_STD',nstar)
m=0L & m='FFFFFFDF'X
if nstar gt 0 then plugmap.fiberdata[ind].sectarget = long(plugmap.fiberdata[ind].sectarget and m)

; remove duplicate plateId entries if they exist
if n_elements(plugmap.plateid) gt 1 then begin
  for i=1,n_elements(plugmap.plateid)-1 do if plugmap.plateid[0] ne plugmap.plateid[i] then $
     stop,'halt: multiple plate ids in plugmap do not match!'
  plateid=plugmap.plateid[0]
  remove_tags,plugmap,'plateid',plugtmp
  plugmap=create_struct(plugtmp, 'plateid', plateid)
endif

; custom errors in mapping?
if n_elements(fixfiberid) gt 0 then begin
 if fixfiberid eq 1 then begin
  star=where(plugmap.fiberdata.spectrographid eq 2)
  for istar=0,n_elements(star)-1 do begin
    fiberid=plugmap.fiberdata[star[istar]].fiberid
    if fiberid ge 0 then begin
      subid=(fiberid - 1) mod 30
      bundleid=(fiberid-subid)/30
      plugmap.fiberdata[star[istar]].fiberid = (9-bundleid)*30 + subid +1
    endif
print,istar,fiberid,subid,bundleid,plugmap.fiberdata[star[istar]].fiberid
  endfor
 endif
 if fixfiberid eq 2 then begin
   ; MTP#2 rotated
   star=where(plugmap.fiberdata.spectrographid eq 2 and plugmap.fiberdata.holetype eq 'OBJECT')
   fiberid=plugmap.fiberdata[star].fiberid
   j=where(fiberid ge 31 and fiberid le 36)
   plugmap.fiberdata[star[j]].fiberid =  fiberid[j] + 23
   j=where(fiberid ge 37 and fiberid le 44)
   plugmap.fiberdata[star[j]].fiberid =  fiberid[j] + 8
   j=where(fiberid ge 45 and fiberid le 52)
   plugmap.fiberdata[star[j]].fiberid =  fiberid[j] - 8
   j=where(fiberid ge 54 and fiberid le 59)
   plugmap.fiberdata[star[j]].fiberid = fiberid[j] - 23

   ; missing fibers from unpopulated 2 of MTP
   j=where(plugmap.fiberdata[star].fiberid eq 53 or plugmap.fiberdata[star].fiberid eq 60)
   plugmap.fiberdata[star[j]].fiberid = -1
 endif
endif

if keyword_set(stp) then stop

end
