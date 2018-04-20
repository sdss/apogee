;+
;
; MAKECAL
;
; This program contains routines to handle APOGEE calibration data
; It include routines:
;   readcal  : reads a master list file specifying calibration frames
;              to be used on different dates, and the information needed
;              to create these frames, returns info in structures for
;              each cal product
;   getcal   : given a cal structure and MJD, returns info from appropriate
;              MJD for requested date
;   readcalstr returns the name of calibration product in input structure to
;              be used for specified date
;   getnums    parses a string with file names into array of file numbers
;   makecal    will make one or ALL of the specified calibration product types
;              listed in the master calibration index file
;
; the actual construction of the calibration products is done with routines
;   called here, but located in apmkcal.pro
;
; Written by J.Holtzman Aug 2011
;-

;======================================================================

; READCAL
; 
; INPUTS:
;   file     Input name of master calibration list file (e.g., cal.par)
;
; OUTPUTS:
;   darkstr,rnstr,flatstr,gainstr,wavestr,lsfstr:
;            structures containing information to allow construction
;            of each cal product, and giving range of MJD over which
;            product is to be used
; structures are as follows:
; darkstr = { mjd1: mjd1, mjd2: mjd2, name: name, frames: frames }
; flatstr = { mjd1: mjd1, mjd2: mjd2, name: name, frames: frames, nrep: nrep, dithered: dithered}
; wavestr = { mjd1: mjd1, mjd2: mjd2, name: name, frames: frames, psfid: psfid}
; lsfstr = { mjd1: mjd1, mjd2: mjd2, name: name, frames: frames , psfid: psfid}
; detstr = { mjd1: mjd1, mjd2: mjd2, name: name, linid: linid}
; bpmstr = { mjd1: mjd1, mjd2: mjd2, name: name, darkid: darkid, flatid: flatid}
; rnstr = { frame1: frame1, frame2: frame2}
; gainstr = { frame1: frame1, frame2: frame2}
; littrowstr = { mjd1: mjd1, mjd2: mjd2, name: name, psfid: psfid}
; persiststr = { mjd1: mjd1, mjd2: mjd2, name: name, darkid: darkid, flatid: flatid, thresh: thresh}
; responsestr = { mjd1: mjd1, mjd2: mjd2, name: name, fluxid: fluxid, psfid: psfid, temp: temp}
; sparsestr = { mjd1: mjd1, mjd2: mjd2, name: name, frames: frames, darkframes: darkframes, dmax: dmax, maxread: maxread}
;         

pro readcal,file,darkstr,flatstr,sparsestr,fiberstr,badfiberstr,fixfiberstr,wavestr,lsfstr,bpmstr,fluxstr,detstr,littrowstr,persiststr,persistmodelstr,responsestr

; read all the lines in the calibration listfile
openr, lun, /get_lun, file
aline = ''  
i=0
while ~EOF(lun) do begin
  readf, lun, aline  
  if i eq 0 then line=aline else line=[line,aline]
  i+=1
endwhile
free_lun,lun

; establish generic variable types
mjd1=0L & mjd2=0L & name=0L & frames='' & darkframes='' & nframes=0 & nrep=0 & dithered=0 & maxread=''
darkid=0L & flatid=0L & psfid=0L & fluxid=0L & linid=0L & fixfiberid=0L

; extract the darks information and load dark structure
darks=where(strpos(line,'dark') eq 0)
for i=0,n_elements(darks)-1 do begin
  if darks[i] ge 0 then begin
    fields=strsplit(line(darks[i]),/extract)
    reads,fields[1],mjd1
    reads,fields[2],mjd2
    reads,fields[3],name
    reads,fields[4],frames
    str = { mjd1: mjd1, mjd2: mjd2, name: name, frames: frames }
    if i eq 0 then darkstr=str else darkstr = struct_append(darkstr,str)
  endif
endfor

; extract the flats information
flats=where(strpos(line,'flat') eq 0)
for i=0,n_elements(flats)-1 do begin
  if flats[i] ge 0 then begin
    fields=strsplit(line(flats[i]),/extract)
    reads,fields[1],mjd1
    reads,fields[2],mjd2
    reads,fields[3],name
    reads,fields[4],frames
    reads,fields[5],nrep
    reads,fields[6],dithered
    str = { mjd1: mjd1, mjd2: mjd2, name: name, frames: frames, nrep: nrep, dithered: dithered}
    if i eq 0 then flatstr=str else flatstr = struct_append(flatstr,str)
  endif
endfor

; extract the sparsepack PSF information
sparse=where(strpos(line,'sparse') eq 0)
for i=0,n_elements(sparse)-1 do begin
  if sparse[i] ge 0 then begin
    fields=strsplit(line(sparse[i]),/extract)
    reads,fields[1],mjd1
    reads,fields[2],mjd2
    reads,fields[3],name
    reads,fields[4],frames
    reads,fields[5],darkframes
    reads,fields[6],dmax
    reads,fields[7],maxread
    str = { mjd1: mjd1, mjd2: mjd2, name: name, frames: frames, darkframes: darkframes, dmax: dmax, maxread: maxread}
    if i eq 0 then sparsestr=str else sparsestr = struct_append(sparsestr,str)
  endif
endfor

; extract the fiber location information
fiber=where(strpos(line,'fiber') eq 0)
for i=0,n_elements(fiber)-1 do begin
  if fiber[i] ge 0 then begin
    fields=strsplit(line(fiber[i]),/extract)
    reads,fields[1],mjd1
    reads,fields[2],mjd2
    reads,fields[3],name
    str = { mjd1: mjd1, mjd2: mjd2, name: name }
    if i eq 0 then fiberstr=str else fiberstr = struct_append(fiberstr,str)
  endif
endfor

; extract the badfiber location information
badfiber=where(strpos(line,'badfiber') eq 0)
for i=0,n_elements(badfiber)-1 do begin
  if badfiber[i] ge 0 then begin
    fields=strsplit(line(badfiber[i]),/extract)
    reads,fields[1],mjd1
    reads,fields[2],mjd2
    reads,fields[3],frames
    str = { mjd1: mjd1, mjd2: mjd2, name: frames}
    if i eq 0 then badfiberstr=str else badfiberstr = struct_append(badfiberstr,str)
  endif
endfor

; extract the fixfiber location information
fixfiber=where(strpos(line,'fixfiber') eq 0)
for i=0,n_elements(fixfiber)-1 do begin
  if fixfiber[i] ge 0 then begin
    fields=strsplit(line(fixfiber[i]),/extract)
    reads,fields[1],mjd1
    reads,fields[2],mjd2
    reads,fields[3],name
    str = { mjd1: mjd1, mjd2: mjd2, name: name}
    if i eq 0 then fixfiberstr=str else fixfiberstr = struct_append(fixfiberstr,str)
  endif
endfor

; extract the wavecal information
waves=where(strpos(line,'wave') eq 0)
for i=0,n_elements(waves)-1 do begin
  if waves[i] ge 0 then begin
    fields=strsplit(line(waves[i]),/extract)
    reads,fields[1],mjd1
    reads,fields[2],mjd2
    reads,fields[3],name
    reads,fields[4],frames
    reads,fields[5],psfid
    str = { mjd1: mjd1, mjd2: mjd2, name: name, frames: frames, psfid: psfid}
    if i eq 0 then wavestr=str else wavestr = struct_append(wavestr,str)
  endif
endfor

; extract the lsf information
lsfs=where(strpos(line,'lsf') eq 0)
for i=0,n_elements(lsfs)-1 do begin
  if lsfs[i] ge 0 then begin
    fields=strsplit(line(lsfs[i]),/extract)
    reads,fields[1],mjd1
    reads,fields[2],mjd2
    reads,fields[3],name
    reads,fields[4],frames
    reads,fields[5],psfid
    str = { mjd1: mjd1, mjd2: mjd2, name: name, frames: frames , psfid: psfid}
    if i eq 0 then lsfstr=str else lsfstr = struct_append(lsfstr,str)
  endif
endfor

; extract the detector information
dets=where(strpos(line,'det') eq 0)
for i=0,n_elements(dets)-1 do begin
  if dets[i] ge 0 then begin
    fields=strsplit(line(dets[i]),/extract)
    reads,fields[1],mjd1
    reads,fields[2],mjd2
    reads,fields[3],name
    reads,fields[4],linid
    str = { mjd1: mjd1, mjd2: mjd2, name: name, linid: linid}
    if i eq 0 then detstr=str else detstr = struct_append(detstr,str)
  endif
endfor

; extract the BPM information
bpms=where(strpos(line,'bpm') eq 0)
for i=0,n_elements(bpms)-1 do begin
  if bpms[i] ge 0 then begin
    fields=strsplit(line(bpms[i]),/extract)
    reads,fields[1],mjd1
    reads,fields[2],mjd2
    reads,fields[3],name
    reads,fields[4],darkid
    reads,fields[5],flatid
    str = { mjd1: mjd1, mjd2: mjd2, name: name, darkid: darkid, flatid: flatid}
    if i eq 0 then bpmstr=str else bpmstr = struct_append(bpmstr,str)
  endif
endfor

; extract the LITTROW information
littrows=where(strpos(line,'littrow') eq 0)
for i=0,n_elements(littrows)-1 do begin
  if littrows[i] ge 0 then begin
    fields=strsplit(line(littrows[i]),/extract)
    reads,fields[1],mjd1
    reads,fields[2],mjd2
    reads,fields[3],name
    reads,fields[4],psfid
    str = { mjd1: mjd1, mjd2: mjd2, name: name, psfid: psfid}
    if i eq 0 then littrowstr=str else littrowstr = struct_append(littrowstr,str)
  endif
endfor

; extract the PERSISTENCE information
persists=where(strpos(line,'persist ') eq 0)
for i=0,n_elements(persists)-1 do begin
  if persists[i] ge 0 then begin
    fields=strsplit(line(persists[i]),/extract)
    reads,fields[1],mjd1
    reads,fields[2],mjd2
    reads,fields[3],name
    reads,fields[4],darkid
    reads,fields[5],flatid
    reads,fields[6],thresh
    str = { mjd1: mjd1, mjd2: mjd2, name: name, darkid: darkid, flatid: flatid, thresh: thresh}
    if i eq 0 then persiststr=str else persiststr = struct_append(persiststr,str)
  endif
endfor

; extract the PERSISTENCE model information
persists=where(strpos(line,'persistmodel') eq 0)
for i=0,n_elements(persists)-1 do begin
  if persists[i] ge 0 then begin
    fields=strsplit(line(persists[i]),/extract)
    reads,fields[1],mjd1
    reads,fields[2],mjd2
    reads,fields[3],name
    str = { mjd1: mjd1, mjd2: mjd2, name: name}
    if i eq 0 then persistmodelstr=str else persistmodelstr = struct_append(persistmodelstr,str)
  endif
endfor

; extract the RESPONSE information
responses=where(strpos(line,'response') eq 0)
for i=0,n_elements(responses)-1 do begin
  if responses[i] ge 0 then begin
    fields=strsplit(line(responses[i]),/extract)
    reads,fields[1],mjd1
    reads,fields[2],mjd2
    reads,fields[3],name
    reads,fields[4],fluxid
    reads,fields[5],psfid
    reads,fields[6],temp
    str = { mjd1: mjd1, mjd2: mjd2, name: fluxid, psf: psfid, temp: temp}
    if i eq 0 then responsestr=str else responsestr = struct_append(responsestr,str)
  endif
endfor

; extract the readnoise information
frame1=0L & frame2=0L
rns=where(strpos(line,'rn') eq 0)
for i=0,n_elements(rns)-1 do begin
  if rns[i] ge 0 then begin
    fields=strsplit(line(rns[i]),/extract)
    reads,fields[1],frame1
    reads,fields[2],frame2
    str = { frame1: frame1, frame2: frame2}
    if i eq 0 then rnstr=str else rnstr = struct_append(rnstr,str)
  endif
endfor
; extract the gain information
gains=where(strpos(line,'gain') eq 0)
for i=0,n_elements(gains)-1 do begin
  if gains[i] ge 0 then begin
    fields=strsplit(line(gains[i]),/extract)
    reads,fields[1],frame1
    reads,fields[2],frame2
    str = { frame1: frame1, frame2: frame2}
    if i eq 0 then gainstr=str else gainstr = struct_append(gainstr,str)
  endif
endfor

end

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

;==========================================================================
; getnums
;       parses a string with file names into array of file numbers
;         aaa-bbb      parses to list of numbers bewteen aaa and bbb (e.g., for consecutive frames)
;         aaa,bbb,ccc  parses to list of aaa bbb ccc (e.g. for non-consecutive frames)
;         aaa          parses to aaa (e.g., for single frames)
;
function getnums,list

  undefine,out
  if strpos(list,',') gt 0 then files=strsplit(list,',',/extract) else files=list
  for i=0,n_elements(files)-1 do begin
    if strpos(files[i],'-') gt 0 then begin
      elements=strsplit(files[i],'-',/extract)
      startval=0L & endval=0L
      reads,elements[0],startval
      reads,elements[1],endval
      n=endval-startval+1
      ims=startval+indgen(n)
    endif else begin
      ims=0L
      reads,files[i],ims
    endelse
    if n_elements(nums) eq 0 then nums=ims else nums=[nums,ims]
  endfor
  return,nums

  if strpos(files,'-') gt 0 then begin
    elements=strsplit(files,'-',/extract)
    startval=0L & endval=0L
    reads,elements[0],startval
    reads,elements[1],endval
    n=endval-startval+1
    return,startval+indgen(n)
  endif else if strpos(files,',') gt 0 then begin
    elements=strsplit(files,',',/extract)
    nums=lonarr(n_elements(elements))
    nnn=0L
    for i=0,n_elements(elements)-1 do begin
      reads,elements[i],nnn
      nums[i]=nnn
    endfor
    return,nums
  endif else begin
    nnn=0L
    reads,files,nnn
    return,nnn
  endelse
  
end

;======================================================================

pro getcal,mjd,file,darkid=darkid,flatid=flatid,sparseid=sparseid,bpmid=bpmid,waveid=waveid,lsfid=lsfid,fluxid=fluxid,detid=detid,fiberid=fiberid,badfiberid=badfiberid,fixfiberid=fixfiberid,littrowid=littrowid,persistid=persistid,persistmodelid=persistmodelid,responseid=responseid

  ; get the calibration files for desired date (mjd) from master calibration index (file)
  readcal,file,darkstr,flatstr,sparsestr,fiberstr,badfiberstr,fixfiberstr,wavestr,lsfstr,bpmstr,fluxstr,detstr,littrowstr,persiststr,persistmodelstr,responsestr
  darkid=readcalstr(darkstr,mjd)
  flatid=readcalstr(flatstr,mjd)
  sparseid=readcalstr(sparsestr,mjd)
  fiberid=readcalstr(fiberstr,mjd)
  badfiberid=getnums(readcalstr(badfiberstr,mjd))
  fixfiberid=getnums(readcalstr(fixfiberstr,mjd))
  bpmid=readcalstr(bpmstr,mjd)
  waveid=readcalstr(wavestr,mjd)
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
  print,'lsf: ', lsfid
  print,'flux: ', fluxid
  print,'det: ', detid
  print,'littrow: ', littrowid
  print,'persist: ', persistid
  print,'persistmodel: ', persistmodelid
  print,'response: ', responseid
end

;===========================================================================
;
;  makecal will make one or ALL of the specified calibration product types
;          listed in the master calibration index file
;
;   USAGE:  makecal,file=file,/dark,/flat,/wave,/lsf
;           OR
;           makecal,file=file,dark=darkid,flat=flatid,wave=waveid,lsf=lsfid
;
;   INPUT:  file=file: name of master calibration index file, if not specified
;               use default cal.par in calibration directory
;           /dark:  make all of the darks in the file
;           dark=darkid:  make the dark with name=darkid 
;           /flat:  make all of the flats in the file
;           flat=flatid:  make the flat with name=flatid 
;           /wave:  make all of the wavecals in the file
;           wave=waveid:  make the wavecal with name=waveid 
;           /lsf:  make all of the lsfs in the file
;           lsf=lsfid:  make the lsf with name=lsfid 
;
pro makecal,file=file,det=det,dark=dark,flat=flat,wave=wave,lsf=lsf,bpm=bpm,$
    psf=psf,flux=flux,sparse=sparse,fiber=fiber,$
    littrow=littrow,persist=persist,modelpersist=modelpersist,response=response,mjd=mjd,full=full,$
    newwave=newwave,nskip=nskip,average=average,clobber=clobber,vers=vers,telescope=telescope

  if keyword_set(vers) and keyword_set(telescope) then apsetver,vers=vers,telescope=telescope
  dirs=getdir(apo_dir,cal_dir,spectro_dir,apo_vers,lib_dir)

  ; get default file name if file not specified
  if keyword_set(file) then file=cal_dir+'/'+file $
  else file=dirs.calfile
  calfile=dirs.calfile
 
  if not keyword_set(full) then full=0
  if not keyword_set(newwave) then newwave=0
  if not keyword_set(nskip) then nskip=1

  ; read calibration master file into calibration structures
  readcal,file,darkstr,flatstr,sparsestr,fiberstr,badfiberstr,fixfiberstr,wavestr,lsfstr,bpmstr,fluxstr,detstr,littrowstr,persiststr,persistmodelstr,responsestr
  ; make detector file as called for
  if keyword_set(det) then begin
    print,'makecal det: ', det
    if det gt 1 then begin
      i=where(detstr.name eq det)
      if i lt 0 then begin
        print,'No matching calibration line for ', det
        stop
      endif
      mkdet,detstr[i].name,detstr[i].linid
    endif 
  endif

  ; make darks as called for
  if keyword_set(dark) then begin
    print,'makecal dark: ', dark
    if dark gt 1 then begin
      i=where(darkstr.name eq dark)
      if i lt 0 then begin
        print,'No matching calibration line for ', dark
        stop
      endif
      ims=getnums(darkstr[i].frames)
      cmjd=getcmjd(ims[0],mjd=mjd)
      getcal,mjd,calfile,detid=detid
      makecal,det=detid
      mkdark,ims,clobber=clobber
    endif else begin
      if keyword_set(mjd) then  begin
        num=getnum(mjd) 
        red=where(darkstr.frames/10000L eq num)
      endif else red=indgen(n_elements(darkstr))
      if (red[0] ge 0) then begin
       for i=0,n_elements(red)-1 do begin
        ims=getnums(darkstr[red[i]].frames)
        cmjd=getcmjd(ims[0],mjd=mjd)
        getcal,mjd,calfile,detid=detid
        makecal,det=detid
        mkdark,ims,clobber=clobber
       endfor
      endif
    endelse
  endif

  ; make flats as called for
  if keyword_set(flat) then begin
    print,'makecal flat: ', flat
    if flat gt 1 then begin
      i=where(flatstr.name eq flat)
      if i lt 0 then begin
        print,'No matching calibration line for ', flat
        stop
      endif
      ims=getnums(flatstr[i].frames)
      cmjd=getcmjd(ims[0],mjd=mjd)
      getcal,mjd,calfile,darkid=darkid
      makecal,dark=darkid
      mkflat,ims,darkid=darkid,nrep=flatstr[i].nrep,dithered=flatstr[i].dithered,clobber=clobber
    endif else begin
      if keyword_set(mjd) then  begin
        num=getnum(mjd) 
        red=where(flatstr.frames/10000L eq num)
      endif else red=indgen(n_elements(darkstr))
      if (red[0] ge 0) then begin
       for i=0,n_elements(red)-1 do begin
        ims=getnums(flatstr[red[i]].frames)
        cmjd=getcmjd(ims[0],mjd=mjd)
        getcal,mjd,calfile,darkid=darkid
        makecal,dark=darkid
        mkflat,ims,darkid=darkid,nrep=flatstr[i].nrep,dithered=flatstr[i].dithered,clobber=clobber
       endfor
      endif
    endelse
  endif

  if keyword_set(bpm) then begin
    print,'makecal bpm: ', bpm
    if bpm gt 1 then begin
      i=where(bpmstr.name eq bpm)
      if i lt 0 then begin
        print,'No matching calibration line for ', bpm
        stop
      endif
      makecal,dark=bpmstr[i].darkid
      makecal,flat=bpmstr[i].flatid
      mkbpm,bpmstr[i].name,darkid=bpmstr[i].darkid,flatid=bpmstr[i].flatid,clobber=clobber
    endif else begin
      if keyword_set(mjd) then  begin
        num=getnum(mjd) 
        red=where(bpmstr.frames/10000L eq num)
      endif else red=indgen(n_elements(bpmstr))
      if (red[0] ge 0) then begin
       for i=0,n_elements(red)-1 do begin
        makecal,dark=bpmstr[i].darkid,clobber=clobber
        makecal,flat=bpmstr[i].flatid,clobber=clobber
        mkbpm,bpmstr[red[i]].name,darkid=bpmstr[i].darkid,flatid=bpmstr[i].flatid,clobber=clobber
       endfor
      endif
    endelse
  endif

  if keyword_set(sparse) then begin
    print,'makecal sparse: ', sparse
    if sparse gt 1 then begin
      i=where(sparsestr.name eq sparse)
      if i lt 0 then begin
        print,'No matching calibration line for ', sparse
        stop
      endif
      ims=getnums(sparsestr[i].frames)
      cmjd=getcmjd(ims[0],mjd=mjd)
      getcal,mjd,calfile,darkid=darkid,flatid=flatid,bpmid=bpmid
      makecal,dark=darkid
      makecal,flat=flatid
      makecal,bpm=bpmid
      darkims=getnums(sparsestr[i].darkframes)
      maxread=getnums(sparsestr[i].maxread)
      if n_elements(maxread) ne 3 then begin
        print,'sparse maxread does not have 3 elements! '
        stop
      endif
      mkepsf,ims,darkid=darkid,flatid=flatid,darkims=darkims,dmax=sparsestr[i].dmax,maxread=maxread,clobber=clobber,/filter,thresh=0.2,scat=2
    endif
  endif

  if keyword_set(fiber) then begin
    print,'makecal fiber: ', fiber
    if fiber gt 1 then begin
      cmjd=getcmjd(fiber,mjd=mjd)
      getcal,mjd,calfile,darkid=darkid,flatid=flatid,sparseid=sparseid
      mkpsf,fiber,darkid=darkid,flatid=flatid,sparseid=sparseid
    endif
  endif

  if keyword_set(psf) and ~keyword_set(flux) then begin
    print,'makecal psf: ', psf
    if psf gt 1 then begin
      cmjd=getcmjd(psf,mjd=mjd)
      getcal,mjd,calfile,darkid=darkid,flatid=flatid,sparseid=sparseid,fiberid=fiberid,littrowid=littrowid
      makecal,littrow=littrowid
      mkpsf,psf,darkid=darkid,flatid=flatid,sparseid=sparseid,fiberid=fiberid,littrowid=littrowid,clobber=clobber
    endif
  endif

  if keyword_set(littrow) then begin
    print,'makecal littrow: ', littrow
    if littrow gt 1 then begin
      cmjd=getcmjd(littrow,mjd=mjd)
      getcal,mjd,calfile,darkid=darkid,flatid=flatid,sparseid=sparseid,fiberid=fiberid
      makecal,flat=flatid
      mklittrow,littrow,cmjd=cmjd,darkid=darkid,flatid=flatid,sparseid=sparseid,fiberid=fiberid,clobber=clobber
    endif
  endif

  if keyword_set(persist) then begin
    print,'makecal persist: ', persist
    if persist gt 1 then begin
      i=where(persiststr.name eq persist)
      if i lt 0 then begin
        print,'No matching calibration line for ', persist
        stop
      endif
      cmjd=getcmjd(persist,mjd=mjd)
      getcal,mjd,calfile,darkid=darkid,flatid=flatid,sparseid=sparseid,fiberid=fiberid
      mkpersist,persist,persiststr[i].darkid,persiststr[i].flatid,thresh=persiststr[i].thresh,cmjd=cmjd,darkid=darkid,flatid=flatid,sparseid=sparseid,fiberid=fiberid,clobber=clobber
    endif
  endif

  if keyword_set(modelpersist) then begin
    print,'makecal modelpersist: ', modelpersist
    if modelpersist gt 1 then begin
      i=where(persistmodelstr.name eq modelpersist)
      if i lt 0 then begin
        print,'No matching calibration line for ', modelpersist
        stop
      endif
      cmjd=getcmjd(modelpersist,mjd=mjd)
      getcal,mjd,calfile,darkid=darkid,flatid=flatid,sparseid=sparseid,fiberid=fiberid
      mkpersistmodel,modelpersist
    endif
  endif

  if keyword_set(flux) then begin
    if not keyword_set(psf) then begin
      print,'for flux, need psf as well!'
      stop
    endif
    print,'makecal flux: ', flux
    if flux gt 1 then begin
      cmjd=getcmjd(flux,mjd=mjd)
      makecal,psf=psf
      getcal,mjd,calfile,darkid=darkid,flatid=flatid,littrowid=littrowid,waveid=waveid
      makecal,littrow=littrowid
      mkflux,flux,darkid=darkid,flatid=flatid,psfid=psf,littrowid=littrowid,waveid=waveid,clobber=clobber
    endif
  endif

  if keyword_set(response) then begin
    print,'makecal response: ', response
    if response gt 1 then begin
      i=where(responsestr.name eq response,nres)
      if nres eq 0 then begin
        print,'No matching calibration line for ', response
        stop
      endif else if nres gt 1 then i=i[0]
      file=apogee_filename('Response',chip='c',num=response)
      ; does product already exist?
      if file_test(file) and not keyword_set(clobber) then begin
        print,' flux file: ',file,' already made'
        return
      endif

      cmjd=getcmjd(response,mjd=mjd)
      getcal,mjd,calfile,darkid=darkid,flatid=flatid,littrowid=littrowid,waveid=waveid,fiberid=fiberid
      makecal,psf=responsestr[i].psf
      makecal,wave=waveid
      makecal,fiber=fiberid
      makecal,littrow=littrowid
      mkflux,response,darkid=darkid,flatid=flatid,psfid=responsestr[i].psf,littrowid=littrowid,waveid=waveid,temp=responsestr[i].temp,clobber=clobber
    endif
  endif


  if keyword_set(wave) then begin
    print,'makecal wave: ', wave
    if wave gt 1 then begin
      i=where(wavestr.name eq wave,nwave)
      if nwave le 0 then begin
        print,'No matching calibration line for ', wave
        stop
      endif
      ims=getnums(wavestr[i[0]].frames)
      cmjd=getcmjd(ims[0],mjd=mjd)
      getcal,mjd,calfile,darkid=darkid,flatid=flatid
      mkwave,ims,darkid=darkid,flatid=flatid,psfid=wavestr[i[0]].psfid
    endif else begin
      if keyword_set(mjd) then  begin
        num=getnum(mjd) 
        red=where(wavestr.frames/10000L eq num)
      endif else red=indgen(n_elements(wavestr))
      if (red[0] ge 0) then begin
       for i=0,n_elements(red)-1,nskip do begin
        ims=getnums(wavestr[red[i]].frames)
        cmjd=getcmjd(ims[0],mjd=mjd)
        getcal,mjd,calfile,darkid=darkid,flatid=flatid
        mkwave,ims,darkid=darkid,flatid=flatid,psfid=wavestr[red[i]].psfid
       endfor
      endif
    endelse
  endif
  if keyword_set(lsf) then begin
    print,'makecal lsf: ', lsf
    if lsf gt 1 then begin
      i=where(lsfstr.name eq lsf,nlsf)
      if nlsf le 0 then begin
        print,'No matching calibration line for ', lsf
        stop
      endif
      ims=getnums(lsfstr[i[0]].frames)
      cmjd=getcmjd(ims[0],mjd=mjd)
      getcal,mjd,calfile,darkid=darkid,flatid=flatid,waveid=waveid
      mklsf,ims,waveid,darkid=darkid,flatid=flatid,psfid=lsfstr[i[0]].psfid,full=full,newwave=newwave
    endif else begin
      if keyword_set(mjd) then  begin
        num=getnum(mjd) 
        red=where(lsfstr.frames/10000L eq num)
      endif else red=indgen(n_elements(lsfstr))
      if (red[0] ge 0) then begin
       for i=0,n_elements(red)-1,nskip do begin
        ims=getnums(lsfstr[red[i]].frames)
        cmjd=getcmjd(ims[0],mjd=mjd)
        getcal,mjd,calfile,darkid=darkid,flatid=flatid,waveid=waveid
        mklsf,ims,waveid,darkid=darkid,flatid=flatid,psfid=lsfstr[i].psfid,full=full,newwave=newwave
       endfor
      endif
    endelse
  endif
end
