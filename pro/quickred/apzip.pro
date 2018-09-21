pro apzip,input,delete=delete,silent=silent,error=error

;+
;
; APZIP
;
; This program compresses the raw APOGEE files
; using various techniques.
;
; If the output compressed file already exists then
; it is automatically overwritten!
;
; This program is specificially designed to compress
; ONLY raw APOGEE data.  It must be in this EXACT format:
;  HDU0: header but NO data
;  HDU1: header, read1 image  as UNSIGNED INTEGERS (BITPIX=16 or UINT)
;  HDU2: header, read2 image  as UNSIGNED INTEGERS (BITPIX=16 or UINT)
;  and so on for all the reads.
;
; INPUTS:
;  input      A list of input raw bundled APOGEE fits files.
;  /delete    Delete the original file after successfully compressing.
;  /silent    Don't print anything to the screen
;
; OUTPUTS:
;  The files are compressed and have filenames with
;  extensions of ".apz".
;  =error     The error message, if one occurred..
;
; USAGE:
;  IDL>apzip,'apR-a-00000085.fits'
;
; By D.Nidever  August 2010
;-

t0 = systime(1)

; Error Handling
;------------------
; Establish error handler. When errors occur, the index of the  
; error is returned in the variable Error_status:  
;CATCH, Error_status 

;This statement begins the error handler:  
;if (Error_status ne 0) then begin 
;   error = !ERROR_STATE.MSG  
;   if not keyword_set(silent) then print,error
;   CATCH, /CANCEL 
;   return
;endif


; Not enough inputs
ninput = n_elements(input)
if ninput eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax - apzip,input,delete=delete,silent=silent,error=error'
  return
endif

; Get the inputs
LOADINPUT,input,files,count=nfiles

; More than one file input
if nfiles gt 1 then begin
  if not keyword_set(silent) then print,strtrim(nfiles,2),' Files input'
  for i=0,nfiles-1 do begin
    if not keyword_set(silent) then print,strtrim(i+1,2),'/',strtrim(nfiles,2)
    if not keyword_set(silent) then print,''
    apzip,files[i],delete=delete,silent=silent,error=error
    if not keyword_set(silent) then print,''
  end
  return
endif

; Does the file exist
if file_test(files) eq 0 then begin
  error = files+' NOT FOUND'
  if not keyword_set(silent) then print,error
  return
endif

; Check that "fpack" is available
spawn,'fpack -H',out,errout
if errout[0] ne '' then begin
  error = 'FPACK not found'
  if not keyword_set(silent) then print,error
  return
endif

; Check that the extension is ".fits"
dir = file_dirname(files)+'/'
fil = file_basename(files)
dum = strsplit(fil,'.',/extract)
ext = dum[n_elements(dum)-1]
if ext ne 'fits' then begin
  error = 'Extension must be .fits'
  if not keyword_set(silent) then print,error
  return
endif

; Test that we can read the file
FITS_READ,files,im,head,exten=0,message=message,/no_abort
if message ne '' then begin
  error = 'ERROR reading '+files+' '+message
  if not keyword_set(silent) then print,error
  return
endif

; verify the validity of the checksum and datasum from the header
res=FITS_TEST_CHECKSUM(head, im, errmsg=errmsg)
if res eq -1 then begin
  ; the checksum doesn't match -> send a warning and keep going
  print,'APZIP ERROR:  BAD checksum for file (ext=0) '+files
  print,'    '+errmsg
  return
endif

; Temporary directory
;  use /tmp/ if possible otherwise the directory that the file is in
;tempdir = '/tmp/'
;if FILE_TEST(tempdir,/directory) eq 0 then
tempdir = FILE_DIRNAME(files)


infoin = FILE_INFO(files)
if not keyword_set(silent) then $
  print,'Compressing >>',files,'<< (',strtrim(string(infoin.size/1e6,format='(F10.2)'),2),' MB)'


; Check format and get number of reads
;--------------------------------------
; Check primary header, no data allowed
head0 = headfits(files,exten=0,errmsg=errmsg)
if sxpar(head0,'NAXIS') ne 0 then begin
  error = 'Primary HDU has data in it.  This is not allowed!'
  if not keyword_set(silent) then print,error
  return
endif

flag = 0
nreads = 0
while flag eq 0 do begin
  nreads++

  ; Getting header
  head = headfits(files,exten=nreads,errmsg=message)

  ; Checking data format
  if message eq '' then begin
    bitpix = sxpar(head,'BITPIX')
    if bitpix ne 16 then begin
      error = 'Error: BITPIX='+strtrim(bitpix,2)+'  READS must be UNSIGNED INTEGERS (BITPIX=16)'
      if not keyword_set(silent) then print,error
      return
    endif
  endif

  if message ne '' then flag=1  ; there was an error, stop
endwhile
nreads--  ; last one was too much

if not keyword_set(silent) then $
  print,'Nreads = ',strtrim(nreads,2)


; There is data to compress, Nreads>0
;-------------------------------------
if nreads ge 1 then begin

  ; Load first read
  FITS_READ,files,im,head,exten=1,/no_abort

  ; verify the validity of the checksum and datasum from the header
  res=FITS_TEST_CHECKSUM(head, im, errmsg=errmsg)
  if res eq -1 then begin
    ; the checksum doesn't match -> send a warning and keep going
    print,'APZIP ERROR:  BAD checksum for file (ext=1) '+files
    print,'    '+errmsg
    return
  endif

  sz1 = size(im)
  npix = sz1[1]
  im1 = long(im)  ; The first image


  ; Step I: Make dCounts temporary file
  ;------------------------------------
  if not keyword_set(silent) then $
    print,'Step I: Making dCounts temporary file'

  dcounts_tempfile = MKTEMP('apzip',outdir=tempdir)
  dcounts_tempfile = dcounts_tempfile[0]

  ; Initialize the dCounts temporary file
  head0 = headfits(files,exten=0)      ; get exten=0 header
  sxaddpar,head0,'SIMPLE','T',''
  sxdelpar,head0,'CHECKSUM'
  sxdelpar,head0,'DATASUM'
  FITS_ADD_CHECKSUM, head0, /no_timestamp
  MWRFITS,0,dcounts_tempfile,head0,/create,/no_comment  ; exten=0 is blank


  ; Loop through the reads
  ;  start with 2nd read
  lastim = im1
  tot_dcounts = im1*0.0
  for i=2,nreads do begin

    ; Load the next READ
    FITS_READ,files,im,head,exten=i,/no_abort
    ; verify the validity of the checksum and datasum from the header
    res=FITS_TEST_CHECKSUM(head, im, errmsg=errmsg)
    if res eq -1 then begin
      ; the checksum doesn't match -> send a warning and keep going
      print,'APZIP ERROR:  BAD checksum for file (ext='+strtrim(i,2)+') '+files
      print,'    '+errmsg
      return
    endif
    im = long(im)
    sz = size(im)

    ; Check that the image dimension is correct
    if sz[1] ne sz1[1] or sz[2] ne sz1[2] then begin
      error = 'Images dimensions of READ1 (in exten=1) and READ'+strtrim(i,2)+' (in exten='+strtrim(i,2)+') do NOT MATCH'
      if not keyword_set(silent) then print,error
      FILE_DELETE,dcounts_tempfile,/allow  ; delete temporary file
      return
    endif

    ; Make dCounts
    dcounts = im - lastim

    ; Fix the header
    sxaddpar,head,'BITPIX',32  ; needs to be LONG
    sxaddpar,head,'BZERO',0
    sxdelpar,head,'SIMPLE'      ; delete SIMPLE if present, only allowed in PDU
    sxdelpar,head,'CHECKSUM'
    sxdelpar,head,'DATASUM'
    FITS_ADD_CHECKSUM, head, dcounts, /no_timestamp

    ; Write to the temporary dCounts file
    MWRFITS,dcounts,dcounts_tempfile,head,/silent

    ; Save last read
    lastim = im

    tot_dcounts += dcounts  ; add to the sum of all dCounts

  endfor

  ; Calculate average dCounts
  avg_dcounts = ROUND( tot_dcounts/(nreads-1) )  ; must be an integer


  ; Initialize the final (pre-compressed) file
  ;--------------------------------------------
  outfile_precmp = MKTEMP('apzip',outdir=tempdir)
  outfile_precmp = outfile_precmp[0]

  ; Put Average dCounts in HDU0 with the original header
  head0 = headfits(files,exten=0)
  sxaddpar,head0,'BITPIX',32   ; needs to be LONG
  SXADDPAR,head0,'NAXIS',size(avg_dcounts,/n_dim),'Dimensionality',after='BITPIX'
  SXADDPAR,head0,'NAXIS1',n_elements(avg_dcounts[*,0]),after='NAXIS'
  SXADDPAR,head0,'NAXIS2',n_elements(avg_dcounts[0,*]),after='NAXIS1'
  sxaddpar,head0,'BZERO',0,after='NAXIS2'
  sxaddpar,head0,'BSCALE',1,after='BZERO'
  sxdelpar,head0,'CHECKSUM'
  sxdelpar,head0,'DATASUM'
  FITS_ADD_CHECKSUM, head0, avg_dcounts, /no_timestamp
  MWRFITS,avg_dcounts,outfile_precmp,head0,/create,/no_comment

  ; Put first read in exten=1
  read0 = uint(im1)
  head1 = headfits(files,exten=1)
  sxaddpar,head1,'BITPIX',16  ; leave as UINT
  sxaddpar,head1,'BZERO',32768
  sxdelpar,head1,'SIMPLE'      ; delete SIMPLE if present, only allowed in PDU
  sxdelpar,head1,'CHECKSUM'
  sxdelpar,head1,'DATASUM'
  FITS_ADD_CHECKSUM, head1, read0, /no_timestamp
  MWRFITS,read0,outfile_precmp,head1,/silent

  ; Step II: Load in dCounts and subtract AVG dCounts
  ;---------------------------------------------------
  if not keyword_set(silent) then $
    print,'Step II: Subtracting average dCounts'
  for i=1,nreads-1 do begin

    ; Load dCounts image (use mrdfits to keep header intact for checksum)
    dcounts = MRDFITS(dcounts_tempfile,i,head,/silent)

    res=FITS_TEST_CHECKSUM(head, dcounts, errmsg=errmsg)
    if res eq -1 then begin
      ; the checksum doesn't match -> send a warning and keep going
      print,'APZIP ERROR:  BAD checksum for file (ext='+strtrim(i,2)+') '+dcounts_tempfile
      print,'    '+errmsg
      return
    endif

    ; Subtract the average dcounts
    resid = dcounts - avg_dcounts

    ; Get the header for 2nd read of this pair
    ;  read=2 for first dcounts
    head = headfits(files,exten=i+1)

    ; Difference images minus Mean count rate
    sxaddpar,head,'BITPIX',32  ; needs to be LONG
    sxaddpar,head,'BZERO',0
    sxdelpar,head,'SIMPLE'      ; delete SIMPLE if present, only allowed in PDU
    sxdelpar,head,'CHECKSUM'
    sxdelpar,head,'DATASUM'
    FITS_ADD_CHECKSUM, head, resid, /no_timestamp
    MWRFITS,resid,outfile_precmp,head,/silent

  endfor

  ; Delete temporary file
  FILE_DELETE,dcounts_tempfile,/allow   ; delete temporary file

; No data to compress, Nreads=0
;-------------------------------
Endif else begin

  ; Making pre-compressed temporary filename
  outfile_precmp = MKTEMP('apzip',outdir=tempdir)
  outfile_precmp = outfile_precmp[0]
  FILE_COPY,files,outfile_precmp,/over,/allow  ; just copy file

Endelse


; Step III: Compress the file with fpack
;--------------------------------------
if not keyword_set(silent) then $
  print,'Step III: Compressing with fpack'
FILE_DELETE,outfile_precmp+'.fz',/allow
SPAWN,['fpack','-C',outfile_precmp],/noshell,out,errout   ; -C suppresses checksum update
outbd = where(stregex(out,'error',/fold_case,/boolean) eq 1,noutbd)
if errout[0] ne '' or noutbd gt 0 then begin
  if errout[0] ne '' then error = 'fpack error '+errout
  if nbdout gt 0 then error=error+out
  if not keyword_set(silent) then print,error
  FILE_DELETE,outfile_precmp,/allow   ; delete temporary file
  return
endif


; Make final output filename
dir = FILE_DIRNAME(files)+'/'
base = FILE_BASENAME(files,'.fits')
finalfile = dir+base+'.apz'

; Rename the compressed file
if file_test(finalfile) eq 1 and not keyword_set(silent) then print,'Overwriting ',finalfile
FILE_MOVE,outfile_precmp+'.fz',finalfile,/overwrite

; Delete temporary files
FILE_DELETE,outfile_precmp,/allow   ; delete temporary file

; Delete original file
if keyword_set(delete) then begin
  if not keyword_set(silent) then $
    print,'Deleting Original file ',files
  FILE_DELETE,files
endif

; Final compression
infoout = file_info(finalfile)
if not keyword_set(silent) then begin
  print,'Input file size = ',strtrim(infoin.size,2),' bytes'
  print,'Output file size = ',strtrim(infoout.size,2),' bytes'
  print,'Compression ratio = ',strtrim(float(infoin.size)/infoout.size,2)
endif

; Time elapsed
dt = systime(1)-t0
if not keyword_set(silent) then $
  print,'dt = ',strtrim(dt,2),' sec'

if keyword_set(stp) then stop

end
