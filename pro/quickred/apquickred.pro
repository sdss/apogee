;+
;
; APQUICKRED
;
; APOGEE Quick Reduction for on the mountain
;
; INPUTS:
;  frameid       The frame number
;  =plugmap      The plugmag structure for this plate.
;  =rawdir       The base directory for the raw file.  The default is
;                  APODATA_DIR+'raw/'
;  =bundledir    The final output directory for the compressed bundled
;                  files.  The default is APQLARCHIVE_DIR+MJD5
;  =quickreddir  The final output directory for the quick reduced
;                  ap2D and ap1D files.  The default is APQLQUICKRED_DIR+MJD5
;  =bpmid        The directory+ID8 concatenated for the bad pixel mask
;                  to use for these exposures.
;  =psfid        The directory+ID8 concatenated for the PSF
;                  calibration file to use for these exposures
;  /no_compress  Whether to compress the output data.  The default is
;                  to compress, so to not compress use /no_compress
;  /no_dbinsert  Do not update the database.  The default is to update
;                  the database.
;  =snr_goals    A structure giving the S/N goals, normally
;                  obtained from the apogeeql database.
;  =exp_pk       The exposure primary key from the platedb databse.
;  =outfile      If present, will write the quickred outputs to the
;                named output file (excluding long arrays) 
;  /stp          Stop at the end of the program.
;
; OUTPUTS:
;  The bundled files are created and compressed.  The
;  datacubes are collapsed and extracted.
;  =dbstr        The database structure.
;
; USAGE:
;  IDL>apquickred,118
;
; By D.Nidever  Dec 2010
;
;
;-
;
pro apquickred,frameid,plugmap=plugmap,rawdir=rawdir,bundledir=bundledir,quickreddir=quickreddir,$
               bpmid=bpmid,psfid=psfid,no_compress=no_compress,snr_goals=snr_goals,stp=stp,$
               no_dbinsert=no_dbinsert,dbstr=dbstr,exp_pk=exp_pk,plugfile=plugfile,$
               mjd5=mjd5,outfile=outfile

print,'in APQUICKRED'

t0 = systime(1)

; Error Handling
;------------------
; Establish error handler. When errors occur, the index of the
; error is returned in the variable Error_status:  
;CATCH, Error_status 
;
;This statement begins the error handler:  
;if (Error_status ne 0) then begin 
;   error = !ERROR_STATE.MSG  
;   if not keyword_set(silent) then print,error
;   CATCH, /CANCEL 
;   return
;endif

apgundef,dbstr

; setup SDSS database parameters
SDSS_DB_PARAMS

; Processing steps:
; 1) Bundle
; 2) Collapse datacube
; 3) Extract
; 4) Diagnostics
; 5) Database update
; 6) Compress


; Not enough inputs
nframeid = n_elements(frameid)
if nframeid ne 1 then begin
  print,'Syntax - apquickred,frameid,plugmap=plugmap,rawdir=rawdir,bundledir=bundledir,quickreddir=quickreddir,'
  print,'               bpmid=bpmid,psfid=psfid,no_compress=no_compress,snr_goals=snr_goals,stp=stp,'
  print,'               no_dbinsert=no_dbinsert,exp_pk=exp_pk,dbstr=dbstr,plugfile=plugfile'
  return
endif

; Default S/N goal parameters
if keyword_set(snr_goals) then begin
   hmag_standard = snr_goals.hmag
   snr_standard_goal = snr_goals.snr_goals
   hmag_standard_version= snr_goals.version
endif else begin
   hmag_standard = 12.0       ; Hmag used for getting S/N
   snr_standard_goal = 30.0   ; S/N goal at standard Hmag
   hmag_standard_version = 0
endelse

print,''
print,'==============================='
print,'Running APOGEE QUICK REDUCTION'
print,'==============================='
print,systime(0)


; Get APOGEE directories
;data_dir = APGETDIR('APQLDATA_DIR',/exists,error=direrr)
;if n_elements(direrr) gt 0 then return
;spectro_dir = APGETDIR('APQLSPECTRO_DIR',/exists,error=direrr)
;if n_elements(direrr) gt 0 then return
;archive_dir = APGETDIR('APQLARCHIVE_DIR',/exists,error=direrr)
;if n_elements(direrr) gt 0 then return
;quickred_dir = APGETDIR('APQLQUICKRED_DIR',/exists,error=direrr)
;if n_elements(direrr) gt 0 then return
;
;raw_dir = data_dir
;bundle_dir = archive_dir
;bpmdir = spectro_dir+'cal/bpm/'
;psfdir = spectro_dir+'cal/psf/'
;quickred_dir=quickred_dir+'/'+string(mjd5,format='(I0)')+'/'

; Using input directories
;if n_elements(rawdir) gt 0 then raw_dir=addslash(rawdir)
;if n_elements(bundledir) gt 0 then bundle_dir=addslash(bundledir)
;if n_elements(quickreddir) gt 0 then quickred_dir=addslash(quickreddir)

; Creating quickred_dir directory
;if file_test(quickred_dir,/directory) eq 0 then begin
;  print,'Directory ',quickred_dir,' does NOT exist.  Creating it.'
;  FILE_MKDIR,quickred_dir
;endif


chiptag = ['a','b','c']
npix = 2048L
nchips = 3

; Frameid
;frameid = strtrim(string(long(frameid),format='(I04)'),2)
;framenumber=strmid(frameid,0,4)

; Get all of the files for this frame
raw_dir=apogee_filename('Raw',num=frameid,read=1,/dir)
rawfiles = file_search(raw_dir+'apRaw-'+strtrim(frameid,2)+'-???.fits',count=nrawfiles)
if nrawfiles eq 0 then begin
  ;print,'NO RAW files for ',raw_dir+'apRaw-'+strtrim(frameid,2)+'-???.fits'
  print,'NO RAW files for ',rawfiles
  return
endif
base = file_basename(rawfiles,'.fits')
num = strmid(base,15,3)
nreads = nrawfiles


;-----------------------
; Get BPM file to use
;-----------------------

; Input BPM file
if n_elements(bpmid) gt 0 then begin
  ;bpmfiles = file_dirname(bpmid)+'/apBPM-'+chiptag+'-'+file_basename(bpmid)+'.fits'
  bpmfiles = apogee_filename('BPM',num=bpmid,chip=chiptag)
  bpm_test = file_test(bpmfiles)
  bd_bpmtest = where(bpm_test eq 0,nbd_bpmtest)
  if nbd_bpmtest gt 0 then begin
    print,bpmfiles[bd_bpmtest],' NOT FOUND'
    apgundef,bpmfiles,bpmcorr
  endif else begin
    print,'Using BPM = ',bpmid
    bpmcorr = file_dirname(bpmid)+'/apBPM-'+chiptag+'-'+file_basename(bpmid)+'.fits'
  endelse
endif

; Find BPM file to use
if n_elements(bpmcorr) eq 0 then begin

  bpmdir = apogee_filename('BPM',num=0,chip=chiptag,/dir)
  bpmfiles = file_search(bpmdir+'*BPM-a-*.fits',count=nbpmfiles)
  if nbpmfiles gt 0 then begin
    bpminfo = file_info(bpmfiles)
    info = file_info(rawfiles[0])
    bestbpm = first_el(minloc( abs(bpminfo.mtime-info.mtime) ))
    bpmframeids = strmid(file_basename(bpmfiles,'.fits'),8)
    ;bpmframeids = strmid(file_basename(bpmfiles,'.fits'),8,8)
    ;bpmframenums = long(bpmframeids)
    ;bestbpm = first_el(minloc( abs(bpmframenums-long(iframeid)) ))
    bpmcorr = bpmdir+'/apBPM-'+chiptag+'-'+bpmframeids[bestbpm]+'.fits'
    print,'Using previously created BPM files ',bpmdir+bpmframeids[bestbpm]
  endif else begin
    print,'NO BPM calibration file found.'
  endelse

endif

;-----------------------
; Get PSF file to use
;-----------------------

; Input PSF file
if n_elements(psfid) gt 0 then begin
  psffiles = file_dirname(psfid)+'/apPSF-'+chiptag+'-'+string(long(file_basename(psfid)),format='(I08)')+'.fits'
  psf_test = file_test(psffiles)
  bd_psftest = where(psf_test eq 0,nbd_psftest)
  if nbd_psftest gt 0 then begin
    print,psffiles[bd_psftest],' NOT FOUND'
    apgundef,psffiles,psfcorr
  endif else begin
    print,'Using PSF = ',psfid
    psfcorr = psfid
  endelse
endif

; Find PSF file to use
if n_elements(psfcorr) eq 0 then begin

  psfdir = apogee_filename('PSF',num=0,chip='a',/dir)
  psffiles = file_search(psfdir+'apPSF-a-*.fits',count=npsffiles)
  if npsffiles gt 0 then begin
    psfframeids = strmid(file_basename(psffiles,'.fits'),8,8)
    psfframenums = long(psfframeids)
    bestpsf = first_el(minloc( abs(psfframenums-long(frameid)) ))
    psfcorr = psfdir+psfframeids[bestpsf]
    print,'Using previously created PSF files ',psfcorr
  endif else begin
    print,'NO PSF calibration file found. CANNOT extract the spectra.'
  endelse

endif

print,''
print,'-----------------------------------------------------------'
print,'Processing Frame = ',frameid,'  Nreads=',strtrim(nreads,2)
print,'-----------------------------------------------------------'

; Maybe make sure that the file has completely copied over before
; reading it in

; Make sure the final directory exists
bundle_dir=apogee_filename('R',num=frameid,chip='a',/dir)
if file_test(bundle_dir,/directory) eq 0 then begin
  print,'Directory ',bundle_dir,' does NOT exist.  Creating it.'
  FILE_MKDIR,bundle_dir
endif

;stop

;----------------
; Step 1 - Bundle
;-----------------
print,''
print,'Step 1 - Bundling'
print,'-----------------'
print,''
APQBUNDLE,frameid,indir=raw_dir,outdir=bundle_dir,error=bundle_error,/clobber
if n_elements(bundle_error) gt 0 then begin
  print,'APQBUNDLE Error'
  return
endif

;stop

;-----------------------------
; Step 2 - Collapse datacube
;-----------------------------
print,''
print,'Step 2 - Collapsing the datacube'
print,'--------------------------------'
print,''
;cubefiles = bundle_dir+'apR-'+chiptag+'-'+frameid+'.fits'
;imfiles = quickred_dir+'ap2D-'+chiptag+'-'+frameid+'.fits'
cubefiles=apogee_filename('R',num=frameid,chip=chiptag)
; datamodel uses .apz, but here we want .fits, since compressing is done later
cubefiles=file_dirname(cubefiles)+'/'+file_basename(cubefiles,'.apz')+'.fits'
imfiles=apogee_filename('2D',num=frameid,chip=chiptag)
apgundef,output2d
; Loop through the chips
for i=0,2 do begin

  apgundef,ibpmcorr,idarkcorr,ap3d_error,output2d_chip

  ; Bundled file exists
  if file_test(cubefiles[i]) eq 1 then begin

    if n_elements(bpmcorr) gt 0 then ibpmcorr=bpmcorr[i]
    AP3DQUICK,cubefiles[i],imfiles[i],nfowler=10,bpmfile=ibpmcorr,error=ap3d_error,/clobber,/outlong,$
                 output=output2d_chip

    ; Add to output2d
    PUSH,output2d,output2d_chip   ; this is an array of the flux
    print,''

  ; Bundled file NOT found
  endif else begin
    print,cubefiles[i],' NOT FOUND'
  endelse

endfor  ; chip loop

;stop

;---------------------------
; Step 3 - Extract Spectra
;---------------------------
print,''
print,'Step 3 - Extracting the spectra'
print,'-------------------------------'

; Does it make sense to extract this exposure type
head = headfits(cubefiles[0])
exptype = strtrim(strupcase(sxpar(head,'EXPTYPE')),2)
;exptype = 'OBJECT'
doextract = 0
;if exptype eq 'OBJECT' or exptype eq 'FLAT' or exptype eq 'SKY' or exptype eq 'CALIB' or $
;   exptype eq 'LOCALFLAT' or exptype eq 'SUPERFLAT' then doextract=1
if exptype eq 'OBJECT' or exptype eq 'FLAT' or exptype eq 'SKY' or exptype eq 'CALIB' or $
   exptype eq 'SUPERFLAT' or exptype eq 'QUARTZFLAT' or $
   exptype eq 'ARCLAMP' or exptype eq 'DOMEFLAT' or exptype eq 'BLACKBODY' then doextract=1

; Extract
apgundef,output1d
if n_elements(psfcorr) gt 0 and keyword_set(doextract) then begin
  extract_type = 3 ; 2 ;1 ;3  ;1
  quickred_dir=apogee_filename('1D',num=frameid,chip='a',/dir)
  specfiles=apogee_filename('1D',num=frameid,chip=chiptag)
  AP2DPROC,quickred_dir+string(format='(i8.8)',frameid),psfcorr,extract_type,outdir=quickred_dir,fixbadpix=0,/clobber,/outlong,$
           output=output1d
  ; output1d is a structure
  ;specfiles = quickred_dir+'ap1D-'+chiptag+'-'+frameid+'.fits'

; Don't extract
endif else begin
  if n_elements(psfcorr) eq 0 then print,'NO PSF FILE. CANNOT extract the spectra'
  if not keyword_set(doextract) then print,'Cannot do extraction for EXPTYPE='+exptype
endelse

;---------------------------------------------
; Step 4 - Calculate diagnostic information
;---------------------------------------------
; Quicklook-like diagnostics

print,''
print,'Step 4 - Calculate diagnostic information'
print,'--------------------------------------------'
print,''

; Initialize the DBSTR structure
info = apfileinfo(cubefiles[0],/silent)
nbin = 8
npixbin = npix/nbin
yloarr = [500,1000,1500]
yhiarr = yloarr+99
nsub = n_elements(yloarr)
dbstr = {frameid:frameid,header:head,exptype:exptype,exptime:info.exptime,plateid:info.plateid,$
       date:info.dateobs,mjd:info.mjd5,dithpix:sxpar(head,'DITHPIX'),arraydisplay_nbin:nbin,$
       arraydisplay:{data:bytarr(npixbin,npixbin*nchips),bscale:0.0,bzero:0.0,zscale:[0.0,0.0]},$
       arraydisplay_sub:REPLICATE({data:bytarr(npixbin*nchips,npixbin),bscale:0.0,bzero:0.0,$
       zscale:[0.0,0.0],yrange:[0L,0L]},nsub)}
dbstr.arraydisplay_sub.yrange[0] = yloarr
dbstr.arraydisplay_sub.yrange[1] = yhiarr
; Add 1D spectra to DBSTR
if n_elements(output1d) gt 0 then begin
  ; nfibers = n_elements(output1d[0].flux[0,*])
  nfibers = n_elements(output1d.(0).flux[0,*])
  dbstr = CREATE_STRUCT(dbstr,{frame:output1d,qr_spectrum:REPLICATE({spectrum:uintarr(npix,nchips),$
                         fiberid:0, bzero:0.0, bscale:0.0, medsnr:0.0},nfibers)})
endif

; Get S/N only for OBJECT spectra with a PLUGMAP
if n_elements(output1d) gt 0 and exptype eq 'OBJECT' and (n_elements(plugmap) gt 0 or n_elements(plugfile) gt 0) then begin

   dbstr = CREATE_STRUCT(dbstr,{snr:fltarr(nfibers),hmag:fltarr(nfibers),objtype:strarr(nfibers),$
                             logsnr_hmag_coef:fltarr(2),hmag_standard:hmag_standard,$
                             hmag_standard_version:hmag_standard_version,snr_standard:-1.0})
   if n_elements(plugmap) eq 0 then begin
      ; we only got the filename and need to read it in
      APLOADPLUGMAP,plugfile,plugmap
   endif
   APQUICKRED_SNRMAG,dbstr,plugmap=plugmap

print,'SNR: ',dbstr.snr_standard
openu,1,'snr.dat'
printf,1,dbstr.snr_standard
close,1


endif


;-------------------------------
; Step 5 - Update the database
;-------------------------------

  print,''
  print,'Step 5 - Database prep and insert'
  print,'--------------------------'
  print,''

; Add the binned images and spectra
APQUICKRED_DBPREP,output2d,dbstr

; Update the quickred tables in the apogeeql schema
if not keyword_set(no_dbinsert) then begin
    ;savefile = filepath('apqr_db_'+dbstr.frameid+'.sav',/TMP)
    ;print,'apqr_dbinsert savefile  -----------------> ',savefile
    ;save, dbstr, file=savefile 
    ;print,' *************  n_elements(dbstr)=',n_elements(dbstr)
    ;help,dbstr,/st
    APQUICKRED_DBINSERT,instruct=dbstr, exp_pk=exp_pk
endif 

; If an output file has been specified, write the quickred data to it.
if keyword_set(outfile) then begin
    print,'Writing quickred data to ',outfile
    APQUICKRED_WRITEOUTPUTFILE, outfile, dbstr, exp_pk=exp_pk
endif



;--------------------
; Step 6 - Compress
;--------------------
if NOT keyword_set(no_compress) then begin
  print,''
  print,'Step 6 - Compressing files'
  print,'--------------------------'
  print,''

  print,'Compressing the bundled files'
  APZIP,cubefiles

  ; Compressing the 2D image files with fpack
  print,'Compressing the 2D images with FPACK'
  FILE_DELETE,imfiles+'.fz',/allow,/quiet  ; delete compressed files if they already exist
  for i=0,2 do begin
    print,'Compressing ',imfiles[i]
    SPAWN,['fpack','-D',imfiles[i]],out,errout,/noshell   ; delete original after compressing
  end
  ; Deleting the original 2D image files


  ; Compressing the 1D spectral files with fpack
  if n_elements(specfiles) gt 0 then begin
     print,'Compressing the 1D spectral files with FPACK'
     FILE_DELETE,specfiles+'.fz',/allow,/quiet  ; delete compressed files if they already exist
     for i=0,2 do begin
       print,'Compressing ',specfiles[i]
       SPAWN,['fpack','-D',specfiles[i]],out,errout,/noshell   ; delete original after compressing
     endfor
  endif

  ; Remove the bundled files
  print,'Deleting Bundled files'
  FILE_DELETE,cubefiles,/verbose,/allow

endif


dt = systime(1)-t0
print,'dt = ',strtrim(dt,2),' sec'

;stop

if keyword_set(stp) then stop

end
