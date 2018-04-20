;+
;
; APFILEINFO
;
; This function returns a structure that contains lots of information
; on APOGEE files (similar to the IDL FILE_INFO.PRO function).  This
; information can be used to figure out if a file is okay to use for
; a given task/program.
;
; For now this is mainly meant to check RAW files, but hopefully in
; the future it will be more multi-purpose.
;
; INPUTS:
;  files    An array of file names
;  /silent  Don't print anything to the screen.
;  /stp     Stop at the end of the program.
;
; OUTPUTS:
;  info     A structure (one element for each element of "files") that
;             contains information on what's in each file.
;
; USAGE:
;  IDL>info = apfileinfo(files)
;
; By D.Nidever  Feb. 2010
;-

function apfileinfo,files,silent=silent,stp=stp

; Not enough inputs
if n_elements(files) eq 0 then begin
  print,'Syntax - info = apfileinfo(files,silent=silent,stp=stp)'
  return,{exists:0}
endif

; Files is not a string
if size(files,/type) ne 7 then begin
  print,'FILES must be a string (array)'
  return,{exists:0}
endif

nfiles = n_elements(files)

chiptag = ['a','b','c']

; Check the input files
filestr = REPLICATE({file:'',exists:0,filesize:0LL,atime:0LL,ctime:0LL,mtime:0LL,$
                     fpack:0,exten:0,dir:'',base:'',ext:'',rawfmt:0,darkfmt:0,$ 
                     flatfmt:0,linfmt:0,sp2dfmt:0,sp1dfmt:0,chip:'',fid8:'',fmjd5:'',$
                     suffix:'',allchips:0,mjd5:0L,naxis:0,size:lonarr(4),$
                     sz2048:0,nreads:0L,dateobs:'',jd:0.0d0,plateid:0L,cartid:0L,$
                     object:'',obstype:'',exptype:'',imagetyp:'',exptime:0.0,$
                     ra:0.0d,dec:0.0d0,dithpix:0.0,collpist:0.0,colpitch:0.0,collyaw:0.0,$
                     tcammid:0.0,tdetbase:0.0,lampqrtz:-1,lampune:-1,lampthar:-1},nfiles)

tags = tag_names(filestr)
filestr.file = files

for i=0,nfiles-1 do begin

  ; Break the filename into dir/filename parts
  fil = files[i]
  filfull = FILE_EXPAND_PATH(fil)
  fil_dir = FILE_DIRNAME(filfull)
  fil_base = FILE_BASENAME(filfull)
  filestr[i].dir = fil_dir
  filestr[i].base = fil_base
  if strpos(fil_base,'.') ne -1 then $
    ext = first_el(strsplit(fil_base,'.',/extract),/last) else ext=''
  filestr[i].ext = ext

  ; Get info 
  info = FILE_INFO(files[i])

  ; if file doesn't exist, check for compressed version
  if ( not info.exists) then begin
    info = FILE_INFO(files[i]+'.fz')
    if  info.exists then filestr[i].fpack = 1
    files[i]=files[i]+'.fz'
  endif
  if strpos(files[i],'.fz') ge 0 then filestr[i].fpack=1

  ; Stuff in the times
  filestr[i].atime = info.atime
  filestr[i].ctime = info.ctime
  filestr[i].mtime = info.mtime

  ;; Does the filename have the correct format
  ;;  the format needs to be apR-[abc]-XXXXXXXX.fits where the Xs are numbers
  ;filestr[i].apfile = 1   ; by default until we find something wrong
  ;if strmid(fil_base,0,4) ne 'apR-' then filestr[i].apfile = 0
  ;if stregex('abc',strmid(fil_base,4,1),/boolean) eq 0 then filestr[i].apfile = 0
  ;if strmid(fil_base,5,1) ne '-' then filestr[i].apfile = 0
  ;if valid_num(strmid(fil_base,6,8)) eq 0 then filestr[i].apfile = 0
  ;if ext ne 'fits' then filestr[i].apfile = 0
  ;if filestr[i].apfile eq 0 then begin
  ;  if not keyword_set(silent) then print,files[i],' is NOT in the APOGEE format. apR-[abc]-XXXXXXXX.fits'
  ;  goto,BOMB
  ;endif

  ; Check standard APOGEE file formats
  ; RAW files: apR-[abc]-XXXXXXXX.fits  (where the Xs are 8 integer digits)
  filestr[i].rawfmt = stregex(files[i],'a[ps]R-[abc]-[0-9]{8}.fits',/boolean)
  if ext eq 'apz' then $
    filestr[i].rawfmt = stregex(files[i],'a[ps]R-[abc]-[0-9]{8}.apz',/boolean)
  ; DARK files: apDark-[abc]-XXXXX.fits (where the Xs are 5 integer digits)
  ;filestr[i].darkfmt = stregex(files[i],'apDark-[abc]-[0-9]{5}.fits',/boolean)
  filestr[i].darkfmt = stregex(files[i],'a[ps]Dark-[abc]-[0-9]{8}.fits',/boolean)
  ; FLAT files: apFlat-[abc]-XXXXX.fits (where the Xs are 5 integer digits)
  ;filestr[i].flatfmt = stregex(files[i],'apFlat-[abc]-[0-9]{5}.fits',/boolean)
  filestr[i].flatfmt = stregex(files[i],'a[ps]Flat-[abc]-[0-9]{8}.fits',/boolean)
  ; LIN files: apLin-[abc]-XXXXX.fits (where the Xs are 5 integer digits)
  filestr[i].linfmt = stregex(files[i],'a[ps]Lin-[abc]-[0-9]{5}.fits',/boolean)
  ; SP2D files: ap2D-[abc]-XXXXXXXX.fits  (where the Xs are 8 integer digits)
  filestr[i].sp2dfmt = stregex(files[i],'a[ps]2D-[abc]-[0-9]{8}.fits',/boolean)
  ; SP1D files: ap1D-[abc]-XXXXXXXX.fits  (where the Xs are 8 integer digits)
  filestr[i].sp1dfmt = stregex(files[i],'a[ps]1D-[abc]-[0-9]{8}.fits',/boolean)

  ; cal, bpm, trace, wave, lsf, flux, cframe, plate, visit, star

  ; Chip
  if filestr[i].rawfmt eq 1 then filestr[i].chip = strmid(fil_base,4,1)
  if filestr[i].darkfmt eq 1 then filestr[i].chip = strmid(fil_base,7,1)
  if filestr[i].flatfmt eq 1 then filestr[i].chip = strmid(fil_base,7,1)
  if filestr[i].linfmt eq 1 then filestr[i].chip = strmid(fil_base,6,1)
  if filestr[i].sp2dfmt eq 1 then filestr[i].chip = strmid(fil_base,5,1)
  if filestr[i].sp1dfmt eq 1 then filestr[i].chip = strmid(fil_base,5,1)

  ; ID8, only for RAW, SP2D, SP1D files
  if filestr[i].rawfmt eq 1 then filestr[i].fid8 = strmid(fil_base,6,8)
  if filestr[i].sp2dfmt eq 1 then filestr[i].fid8 = strmid(fil_base,7,8)
  if filestr[i].sp1dfmt eq 1 then filestr[i].fid8 = strmid(fil_base,7,8)

  ; MJD5, only for Dark, Flat and Lin files
  if filestr[i].darkfmt eq 1 then filestr[i].fmjd5 = strmid(fil_base,9,5)
  if filestr[i].flatfmt eq 1 then filestr[i].fmjd5 = strmid(fil_base,9,5)
  if filestr[i].linfmt eq 1 then filestr[i].fmjd5 = strmid(fil_base,8,5)

  ; Get suffix for this file, if in the right format
  ;   PREFIX-?-SUFFIX.EXT 
  base = file_basename(fil_base,'.'+filestr[i].ext)
  dum = strsplit(base,'-',/extract)
  if n_elements(dum) ge 3 then begin                      ; PREFIX-?-SUFFIX.EXT 
    prefix = dum[0]+'-'+dum[1]+'-'
    prelen = strlen(prefix)
    suffix = strmid(base,prelen)  ; get everything after the prefix
    filestr[i].suffix = suffix
  endif
  ;if n_elements(dum) eq 3 then filestr[i].suffix = dum[2]  ; PREFIX-?-SUFFIX.EXT 
  if n_elements(dum) eq 2 then filestr[i].suffix = dum[1]  ; PREFIX-SUFFIX.EXT, i.e. apFlux files

  ; Do all three chips exist?
  prefix=strmid(filestr[i].base,0,2)
  chiptest = 0
  suffix='.fits'
  if filestr[i].fpack eq 1 then suffix=suffix+'.fz'
  if filestr[i].rawfmt eq 1 or filestr[i].darkfmt eq 1 or filestr[i].flatfmt eq 1 or $
    filestr[i].linfmt eq 1 or filestr[i].sp2dfmt eq 1 or filestr[i].sp1dfmt eq 1 then begin
    for j=0,2 do begin
      if filestr[i].rawfmt eq 1 then chfil=filestr[i].dir+'/'+prefix+'R-'+chiptag[j]+'-'+filestr[i].fid8+suffix
      if ext eq 'apz' then $
        if filestr[i].rawfmt eq 1 then chfil=filestr[i].dir+'/'+prefix+'R-'+chiptag[j]+'-'+filestr[i].fid8+'.apz'
      ;if filestr[i].darkfmt eq 1 then chfil=filestr[i].dir+'/apDark-'+chiptag[j]+'-'+filestr[i].fmjd5+suffix
      ;if filestr[i].flatfmt eq 1 then chfil=filestr[i].dir+'/apFlat-'+chiptag[j]+'-'+filestr[i].fmjd5+suffix
      if filestr[i].darkfmt eq 1 then chfil=filestr[i].dir+'/'+prefix+'Dark-'+chiptag[j]+'-'+filestr[i].suffix+suffix
      if filestr[i].flatfmt eq 1 then chfil=filestr[i].dir+'/'+prefix+'Flat-'+chiptag[j]+'-'+filestr[i].suffix+suffix
      if filestr[i].linfmt eq 1 then chfil=filestr[i].dir+'/'+prefix+'Lin-'+chiptag[j]+'-'+filestr[i].fmjd5+suffix
      if filestr[i].sp2dfmt eq 1 then chfil=filestr[i].dir+'/'+prefix+'2D-'+chiptag[j]+'-'+filestr[i].fid8+suffix
      if filestr[i].sp1dfmt eq 1 then chfil=filestr[i].dir+'/'+prefix+'1D-'+chiptag[j]+'-'+filestr[i].fid8+suffix
      test = FILE_TEST(chfil)
      chiptest += test
    end
  ; Non-standard filename
  endif else begin
    ; If the file has this type of filename PREFIX-?-SUFFIX.EXT 
    ; then we can check for allchips
    dum = strsplit(fil_base,'-',/extract)
    if n_elements(dum) ge 3 and filestr[i].suffix ne '' then begin
      pre = dum[0]
      ;suf = dum[2]
      suf = filestr[i].suffix+'.fits'
      for j=0,2 do begin
        chfil=filestr[i].dir+'/'+pre+'-'+chiptag[j]+'-'+suf
        test = FILE_TEST(chfil)
        chiptest += test
      end
    endif
  endelse

  if chiptest eq 3 then filestr[i].allchips=1 else filestr[i].allchips=0
  if filestr[i].allchips eq 0 then begin
    if not keyword_set(silent) then print,'ONLY ',strtrim(chiptest,2),' chip files found for ',$
        filestr[i].file,'.  NEED 3 CHIP files'
    ;goto,BOMB
  endif

  ; Does the file exist
  filestr[i].exists = info.exists
  if info.exists eq 0 then begin
    if not keyword_set(silent) then print,files[i],' NOT FOUND'
    goto,BOMB
  endif

  ; File size
  filestr[i].filesize = info.size
  if info.size eq 0 then begin
    if not keyword_set(silent) then print,files[i],' IS EMPTY'
    goto,BOMB
  endif

  ; Does it have extensions
  dum = headfits(files[i],exten=1,errmsg=read_messsage,/silent)
  if n_elements(read_message) eq 0 then filestr[i].exten=1
  ;FITS_READ,files[i],exten_no=1,message=read_message,/no_abort
  ;if read_message eq '' then filestr[i].exten=1

  ; NAXIS and SIZE
  head = headfits(files[i],errmsg=errmsg)
  if ext eq 'apz' then head = headfits(files[i],exten=1,errmsg=errmsg)  ; get exten=1 header for compressed files
  if errmsg ne '' then begin
    print,'Error loading header for ',files[i]
    if not keyword_set(silent) then print,errmsg
    ;goto,BOMB
  end
  naxis = sxpar(head,'NAXIS')
  filestr[i].naxis = naxis
  for j=0,naxis-1 do filestr[i].size[j] = sxpar(head,'NAXIS'+strtrim(j+1,2))

  ; Get MJD5
  mjd5 = GETMJD5(head,error=mjd5err,silent=silent,/sdss)
  filestr[i].mjd5 = mjd5
  if mjd5 eq -1 then begin
    if not keyword_set(silent) then print,'No DATE/TIME information for ',files[i]
  endif
  
  ; Get JD
  dateobs = sxpar(head,'DATE-OBS',count=ndateobs)
  if ndateobs gt 0 then begin
    jd = DATE2JD(dateobs)
    filestr[i].jd = jd
  endif

  ; Is the array size correct?
  naxis1 = sxpar(head,'NAXIS1')
  naxis2 = sxpar(head,'NAXIS2')
  naxis3 = sxpar(head,'NAXIS3')
  filestr[i].sz2048 = 1  ; okay bey default until we find otherwise
  if (naxis1 ne 2048 or naxis2 ne 2048) then begin
    filestr[i].sz2048 = 0
    if not keyword_set(silent) then print,files[i],' is NOT 2048x2048'
  endif

  ; Nreads
  if naxis eq 3 then filestr[i].nreads = naxis3

  ; Try NFRAMES or NFOWLER
  if naxis eq 2 then begin
    nreads = sxpar(head,'NFOWLER',count=num_nfowler)
    if num_nfowler eq 0 then nreads = sxpar(head,'NFRAMES',count=num_nframes)
    filestr[i].nreads = nreads
  endif

  ; DATEOBS, OBSTYPE, EXPTIME, DITHER, CDITHER, PLATE
  ;dateobs = sxpar(head,'DATE-OBS',count=ndateobs)
  ;if ndateobs gt 0 then filestr[i].dateobs = strtrim(dateobs,2)
  ;object = sxpar(head,'OBJECT',count=nobject)
  ;if nobject gt 0 then filestr[i].object = strtrim(object,2)
  ;obstype = sxpar(head,'OBSTYPE',count=nobstype)
  ;if nobstype gt 0 then filestr[i].obstype = strtrim(obstype,2)
  ;imagetyp = sxpar(head,'IMAGETYP',count=nimagetyp)
  ;if nimagetyp gt 0 then filestr[i].imagetyp = strtrim(imagetyp,2)
  ;filestr[i].exptime = sxpar(head,'EXPTIME')
  ;filestr[i].dither = sxpar(head,'DITHER')
  ;filestr[i].cdither = sxpar(head,'CDITHER')
  ;plate = sxpar(head,'PLATE',count=nplate)
  ;if nplate gt 0 then filestr[i].plate = strtrim(plate,2)

  ; Other header keyword values to get
  hdkeywords = ['DATE-OBS','PLATEID','CARTID','OBJECT','OBSTYPE','EXPTYPE','IMAGETYP','EXPTIME','RA','DEC',$
                'DITHPIX','COLLPIST','COLPITCH','COLLYAW','TCAMMID','TDETBASE','LAMPQRTZ','LAMPUNE','LAMPTHAR']
  strtags = ['DATEOBS','PLATEID','CARTID','OBJECT','OBSTYPE','EXPTYPE','IMAGETYP','EXPTIME','RA','DEC',$
             'DITHPIX','COLLPIST','COLPITCH','COLLYAW','TCAMMID','TDETBASE','LAMPQRTZ','LAMPUNE','LAMPTHAR']
  nhdkeywords = n_elements(hdkeywords)

  ; Loop through header keywords
  for j=0,nhdkeywords-1 do begin

    value = sxpar(head,hdkeywords[j],count=count)
    strind = where(tags eq strupcase(strtags[j]),nstrind)
    type = size(filestr[i].(strind[0]),/type)
    if count gt 0 then begin
      if type eq 7 then filestr[i].(strind[0]) = strtrim(fix(value,type=type),2) else $
        filestr[i].(strind[0]) = fix(value,type=type)
    end

  end


  BOMB:

end

if keyword_set(stp) then stop

return,filestr

end
