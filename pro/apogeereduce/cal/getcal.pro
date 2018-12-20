;+
;
; GETCAL
;
; This program contains routines to handle APOGEE calibration data
; It include routines:
;   getcal   : given a cal structure and MJD, returns info from appropriate
;              MJD for requested date
;   readcalstr returns the name of calibration product in input structure to
;              be used for specified date
;   getnums    parses a string with file names into array of file numbers
;
; the actual construction of the calibration products is done with routines
;   called here, but located in apmkcal.pro
;
; Written by J.Holtzman Aug 2011
;-
;======================================================================
;
; readcalstr returns the name of calibration product in input structure to
;            be used for specified date
;   USAGE:
;            name=readcalstr(str,mjd)
;   INPUT
;            str: name of a calibration structure with minimum tags
;                 mjd1, mjd2, name
;            mjd: desired MJD
;   OUTPUT:
;            name: name of calibration product valid for input MJD
;

function readcalstr,str,mjd

  n=0
  for i=0,n_elements(str)-1 do begin
    if mjd ge str[i].mjd1 and mjd le str[i].mjd2 then begin
      ret=i & n+=1
    endif
  endfor
  if n eq 0 then begin
;    print,'No cal product found for mjd ', mjd
;    stop
    return,0L
  endif else if n gt 1 then begin
    print,'Multiple cal products found for mjd ', mjd, ' will use last: ',str[ret].name
    stop
  endif
  return,str[ret].name
end

;======================================================================

pro getcal,mjd,file,darkid=darkid,flatid=flatid,sparseid=sparseid,bpmid=bpmid,waveid=waveid,multiwaveid=multiwaveid,lsfid=lsfid,fluxid=fluxid,detid=detid,fiberid=fiberid,badfiberid=badfiberid,fixfiberid=fixfiberid,littrowid=littrowid,persistid=persistid,persistmodelid=persistmodelid,responseid=responseid

  ; get the calibration files for desired date (mjd) from master calibration index (file)
  readcal,file,darkstr,flatstr,sparsestr,fiberstr,badfiberstr,fixfiberstr,wavestr,lsfstr,bpmstr,fluxstr,detstr,littrowstr,persiststr,persistmodelstr,responsestr,multiwavestr
  darkid=readcalstr(darkstr,mjd)
  flatid=readcalstr(flatstr,mjd)
  sparseid=readcalstr(sparsestr,mjd)
  fiberid=readcalstr(fiberstr,mjd)
  badfiberid=getnums(readcalstr(badfiberstr,mjd))
  fixfiberid=getnums(readcalstr(fixfiberstr,mjd))
  bpmid=readcalstr(bpmstr,mjd)
  waveid=readcalstr(wavestr,mjd)
  multiwaveid=readcalstr(multiwavestr,mjd)
  lsfid=readcalstr(lsfstr,mjd)
  fluxid=readcalstr(fluxstr,mjd)
  detid=readcalstr(detstr,mjd)
  littrowid=readcalstr(littrowstr,mjd)
  persistid=readcalstr(persiststr,mjd)
  persistmodelid=readcalstr(persistmodelstr,mjd)
  responseid=readcalstr(responsestr,mjd)
  print,'dark: ', darkid
  print,'flat: ', flatid
  print,'bpm: ', bpmid
  print,'sparse: ', sparseid
  print,'fiber: ', fiberid
  print,'badfiber: ', badfiberid
  print,'fixfiber: ', fixfiberid
  print,'wave: ', waveid
  print,'multiwave: ', multiwaveid
  print,'lsf: ', lsfid
  print,'flux: ', fluxid
  print,'det: ', detid
  print,'littrow: ', littrowid
  print,'persist: ', persistid
  print,'persistmodel: ', persistmodelid
  print,'response: ', responseid
end

