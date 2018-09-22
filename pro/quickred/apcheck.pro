;+
; apcheck.pro
;
; script to be run daily before transfer to SAS to:
;      1. check to see if each data-ics frame has a corresponding annotated frame;
;         log missing frames in missingannotated file
;      2. insure that all annotated files have a corresponding compressed bundled file in
;         the archive directory, and if not, attempt the bundling; log missing
;         frames in missingbundled file
;      3. copy the plugmap files to the archive directory, if desired
;           (NOT NEEDED SINCE ACTOR IS DOING THIS NOW)
;      4. make the HTML logs using image header information and user comments
;           from the apcomments/MJD5.comment file
;      5. do a checksum of all files in the archive directory to make .md5sum file
;
;  USAGE:
;      apcheck,[mjd=mjd|/today],/clobber,/nomd5,/all,/copyplug,/do1m,/loop
;
;     apcheck,/today,/all
;       should be run, e.g. by crontab, just before the MJD rollover
;     apcheck,/today,/nomd5
;       can be run periodically through the night to update HTML logs on the fly. In
;       this case, only want to check for bundled files for annotated files that
;       have been around for a while (30 minutes), because actor may be in the process
;       of bundling more recent files
;
;   INPUT:
;      mjd=mjd : works on specified MJD rather than MJD at runtime
;      /today  : works on MJD at runtime
;      nomd5   : suppress md5sum calculations
;      /all    : checks for bundled files for ALL annotated files for the MJD,
;                not just for ones more than 30 minutes old
;      /clobber : force rebundling even if compressed bundled file exists
;      /nolog : suppress log file creation
;      /copyplug : copies plugmaps (if available) for all plates listed in
;                 NAME card of any FITS header from the night
;      /do1m   : works on 1m directories
;      /loop   : goes into infinite loop checking for files
;   OUTPUT:
;      creates compressed bundled files in archive directory if they are missing
;      creates MJD5.md5sum file in archive directory
;
;-

pro apcheck,mjd=mjd,today=today,clobber=clobber,nomd5=nomd5,all=all,copyplug=copyplug,nolog=nolog,do1m=do1m,loop=loop,force=force,file=file,obs=obs,apred=apred

if n_elements(apred) eq 0 then apred='quickred'
if n_elements(obs) eq 0 then obs='apo'
if obs eq 'apo' then apsetver,vers=apred,telescope='apo25m'
if obs eq 'lco' then apsetver,vers=apred,telescope='lco25m'

; use today's date with /today, or specified date with mjd=
if keyword_set(today) then begin
  get_juldate,jd
  mjd=floor(jd-2400000.5+0.3)
endif else if not keyword_set(mjd) then begin
  print,'You must specify an MJD with MJD= or /today'
  return
endif
cmjd=string(format='(i5.5)',mjd)

; input directory for annotated frames
icsdir='/data-ics/'+string(format='(i4.4)',mjd-55562)+'/'
if not file_test(icsdir,/dir) and not keyword_set(force) then return

; input directory for annotated frames
rawdir=apogee_filename('Raw',num=(mjd-55562)*10000,read=0,/dir)
;if keyword_set(do1m) then rawdir='/data-ql/data/'+cmjd+'/1m/' else $
;rawdir='/data-ql/data/'+cmjd+'/'
if not file_test(rawdir,/dir) then file_mkdir,rawdir

; output directory to go to archive
outdir=apogee_filename('R',num=(mjd-55562)*10000,chip='a',/dir)
;if keyword_set(do1m) then outdir='/data/apogee/archive1m/'+cmjd+'/' else $
;outdir='/data/apogee/archive/'+cmjd+'/'
if not file_test(outdir,/dir) then file_mkdir,outdir

condition=1
while (condition) do begin
if not keyword_set(loop) then condition=0

; check for complete list of annotated frames
print,'icsdir: ', icsdir
files=file_search(icsdir+'*.fits') 
if files[0] eq '' and not keyword_set(force) then return
openw,1,rawdir+'/missingannotated'
openw,2,outdir+'/'+cmjd+'.missingannotated'
for i=0,n_elements(files)-1 do begin
;print,files[i]
  if not file_test(rawdir+file_basename(files[i])) then begin
    printf,1,files[i]
    printf,2,files[i]
  endif
endfor
close,1
close,2

; find all files with the first read if /all is specified
; find only files older than 30 minutes otherwise
print,'rawdir: ', rawdir
if keyword_set(all) then $
  files=file_search(rawdir+'*-001.fits') $
else if keyword_set(file) then $
  files=file_search(rawdir+'*Raw-'+file+'-001.fits') $
else $
  spawn,'find '+rawdir+' -name \*-001.fits -mmin +30',files

platenames=strarr(n_elements(files))
openw,1,rawdir+'/missingbundled'
openw,2,outdir+'/'+cmjd+'.missingbundled'
for i=0,n_elements(files)-1 do begin
 if files[i] ne '' then begin
  ; parse out filenumber
  name=file_basename(files[i])
  words=strsplit(name,'-',/extract)
  filenumber=words[1]
  ; if .apz file exists, assume all is already done OK (unless /clobber specified)
  ; otherwise, run quickred
  outfiles=apogee_filename('R',num=filenumber,chip=['a','b','c'])
  if file_test(outfiles[0]) eq 0 or $
     file_test(outfiles[1]) eq 0 or $
     file_test(outfiles[2]) eq 0 or $
     keyword_set(clobber) then begin
    printf,1,file_basename(files[i])
    printf,2,file_basename(files[i])
    if not file_test(outdir,/dir) then file_mkdir,outdir
    APQBUNDLE,filenumber,indir=rawdir,outdir=outdir,error=bundle_error,/clobber
    ;cubefiles = outdir+'apR-'+['a','b','c']+'-'+filenumber+'.fits'
    cubefiles = outdir+file_basename(outfiles,'.apz')+'.fits'
    apzip,cubefiles
    FILE_DELETE,cubefiles,/verbose,/allow
  endif else print,' done ',outfiles

  ; get header for PLATE id of this file
  if file_test(outfiles[2]) then begin
  head=headfits(outfiles[2],exten=1)
  ;platelist[i]=sxpar(head,'PLATEID')
  platenames[i]=sxpar(head,'NAME')
  endif
 endif
endfor
close,1
close,2

spawn,'df /data-ics /data-ql /data',dfdata
openw,1,rawdir+'/df'
;openw,2,outdir+'/'+cmjd+'.df'
printf,1,dfdata
;printf,2,dfdata
close,1
;close,2

; with /copyplug, copy over the plates used on this night
if keyword_set(copyplug) then begin
  plates=platenames[uniq(platenames,sort(platenames))]
  for iplate=0,n_elements(plates)-1 do begin
    pdate=strsplit(file_basename(plates[iplate]),'-',/extract)
    if n_elements(pdate) gt 1 then begin
      date=pdate[1]
      print,'/scan/'+date+'/'+'plPlugMapM-'+plates[iplate]+'.par'
      if file_test('/scan/'+date+'/'+'plPlugMapM-'+plates[iplate]+'.par') then begin
        file_copy,'/scan/'+date+'/'+'plPlugMapM-'+plates[iplate]+'.par',outdir,/overwrite
        file_chmod,outdir+'plPlugMapM-'+plates[iplate]+'.par',/g_write
      endif
      if file_test ('/data-ql/plugmaps/plPlugMapA-'+plates[iplate]+'.par') then begin
       file_copy,'/data-ql/plugmaps/plPlugMapA-'+plates[iplate]+'.par',outdir,/overwrite
       file_chmod,outdir+'plPlugMapA-'+plates[iplate]+'.par',/g_write
      endif
    endif
  endfor
endif

; move quickred files
reddir='/data/apogee/quickred/'+cmjd+'/'
if keyword_set(all) then spawn,'mv '+outdir+'/ap?D*.fits* '+reddir

; make HTML log using apglist
if not keyword_set(nolog) then begin
  files=file_search(outdir+'*R-a*.apz')
  if file_test('/home/apogee/apcomments/'+cmjd+'.comments') then $
  apglist,files,['DATE-OBS','NREAD','EXPTYPE','LAMPQRTZ','LAMPTHAR','LAMPUNE','PLATEID','SECZ','SEEING','OBSCMNT','COMMENT','COLLPIST','COLPITCH','DITHPIX','TCAMMID','TLSDETB'],outname=outdir+cmjd+'.log',mjd=mjd,comfile='/home/apogee/apcomments/'+cmjd+'.comments' $
  else $
  apglist,files,['DATE-OBS','NREAD','EXPTYPE','LAMPQRTZ','LAMPTHAR','LAMPUNE','PLATEID','SECZ','SEEING','OBSCMNT','COMMENT','COLLPIST','COLPITCH','DITHPIX','TCAMMID','TLSDETB'],outname=outdir+cmjd+'.log',mjd=mjd
  spawn,'cp '+outdir+cmjd+'.log.html /home/holtz/apogee_logs/'
  spawn,'cp /home/apogee/apcomments/'+cmjd+'.comments '+outdir
endif

if condition eq 1 then wait,loop
endwhile

; remove any files from failed bundling
files=file_search(outdir+'apzip.*')
if n_elements(files) gt 0 and files[0] ne '' then file_delete,files

; make a checksum of all of the files
if not keyword_set(nomd5) then begin
  file_delete,outdir+cmjd+'.md5sum',/allow_nonexistent
  ;spawn,'md5sum '+outdir+'*',result
  spawn,'(cd '+outdir+'; md5sum *)',result
  print,'checksum to ',outdir+cmjd+'.md5sum'
  openw,1,outdir+cmjd+'.md5sum'
  printf,1,result
  close,1
endif
end
