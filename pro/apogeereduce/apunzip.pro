pro apunzip,input,clobber=clobber,delete=delete,silent=silent,error=error,no_checksum=no_checksum,fitsdir=fitsdir,nohalt=nohalt

;+
;
; APUNZIP
;
; This program uncompresses the raw APOGEE files
; that were compressed with APZIP
;
; This program is specificially designed to compress
; ONLY raw APOGEE data.  It assumes that the data is
; in this format:
;  HDU0: header but NO data
;  HDU1: header, read1 image  as UNSIGNED INTEGERS (BITPIX=16 or UINT)
;  HDU2: header, read2 image  as UNSIGNED INTEGERS (BITPIX=16 or UINT)
;  and so on for all the reads.
; The uncompression process returns the data to this exact format.
;
; INPUTS:
;  input     A list of input compressed raw bundled APOGEE fits files
;              with endings of .apz.
;  /clobber  If output file exists then overwrite it.
;  /delete   Delete compressed file after successfully uncompressing
;  /silent   Don't print anything to the screen.
;  /no_checksum If specified, will skip the checksum validation
;
; OUTPUTS:
;  The files are uncompressed and have filenames with
;  extensions of ".apz".
;  =error    The error message, if one occurred.
;
; USAGE:
;  IDL>apunzip,'apR-a-00000085.apz'
;
; By D.Nidever  August 2010
;
; Modified:
;   SBeland  Aug 2011 - Added the checksum
;-

t0 = systime(1)

; Error Handling
;------------------
; Establish error handler. When errors occur, the index of the  
; error is returned in the variable Error_status:  
CATCH, Error_status 

;This statement begins the error handler:  
if (Error_status ne 0) then begin 
   error = !ERROR_STATE.MSG  
   if not keyword_set(silent) then print,error
   CATCH, /CANCEL 
   return
endif

; Not enough inputs
ninput = n_elements(input)
if ninput eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax - apunzip,input,clobber=clobber,delete=delete,silent=silent,error=error,no_checksum=no_checksum'
  return
endif

; Get the inputs
LOADINPUT,input,files,count=nfiles


; More than one file input
if nfiles gt 1 then begin
  for i=0,nfiles-1 do apunzip,files[i],clobber=clobber,delete=delete,silent=silent,error=error,no_checksum=no_checksum
  return
endif

if file_test(files) eq 0 then begin
  error = files+' NOT FOUND'
  if not keyword_set(silent) then print,error
  return
endif

; Check that "funpack" is available
spawn,['funpack','-H'],out,errout,/noshell
if errout[0] ne '' then begin
  error = 'FUNPACK not found'
  if not keyword_set(silent) then print,error
  return
endif

; Check that the extension is ".apz"
dir = file_dirname(files)+'/'
fil = file_basename(files)
dum = strsplit(fil,'.',/extract)
ext = dum[n_elements(dum)-1]
if ext ne 'apz' then begin
  error = 'Extension must be .apz'
  if not keyword_set(silent) then print,error
  return
endif
base = file_basename(files,'.apz')

; Temporary directory
;  use /tmp/ if possible otherwise the directory that the file is in
;tempdir = '/tmp/'
;if FILE_TEST(tempdir,/directory) eq 0 then
if keyword_set(fitsdir) then tempdir = fitsdir else $
  tempdir = FILE_DIRNAME(files)

; Getting file info
info = FILE_INFO(files)
if not keyword_set(silent) then $
  print,'Uncompressing >>',files,'<< (',strtrim(string(info.size/1e6,format='(F10.2)'),2),' MB)'

; Final output filename
if keyword_set(fitsdir) then finalfile=fitsdir+base+'.fits' else $
  finalfile = dir+base+'.fits'

; if another process is working already on this file, wait until done,
;    then return
if file_test(finalfile+'.lock') then begin
  while file_test(files+'.lock') do apwait,files+'.lock',10
  return
endif

; Does the file exist?
; open .lock file
openw,lock,/get_lun,finalfile+'.lock'
free_lun,lock

if file_test(finalfile) eq 1 and keyword_set(clobber) then begin
  if not keyword_set(silent) then $
    print,'Overwriting ',finalfile
  FILE_DELETE,finalfile,/allow
endif
if file_test(finalfile) eq 1 and ~keyword_set(clobber) then begin
  if not keyword_set(silent) then $
    print,finalfile,' exists already.  Writing compressed file to ',finalfile+'.1'
  finalfile = finalfile+'.1'
endif


; uncompress the input file to a temporary file
; get a unique filename (and delete the created empty file)
outfile_uncmp = MKTEMP('apzip',outdir=tempdir)
outfile_uncmp = outfile_uncmp[0]
FILE_DELETE,outfile_uncmp,/allow


; Step I: Uncompress the file with funpack
;-------------------------------------------
if not keyword_set(silent) then $
  print,'Step I: Uncompress with funpack'
spawn,['funpack','-O',outfile_uncmp,'-C',files],out,errout,/noshell  ; -C suppresses checksum update
if n_elements(errout) gt 1 or errout[0] ne '' then begin
  error = 'halt:    fpack error '+errout
  if not keyword_set(silent) then if keyword_set(nohalt) then print,error else stop,error
  return
endif

; Get number of reads
flag = 0
nreads = 0
fits_open,outfile_uncmp,fcb
nreads=fcb.nextend
fits_close,fcb

if not keyword_set(silent) then $
  print,'        Nreads = ',strtrim(nreads,2)

; There is data to uncompress, Nreads>0
;---------------------------------------
if nreads ge 1 then begin


  ; Step II: Reconstructing the reads
  ;-----------------------------------
  if not keyword_set(silent) then $
    print,'Step II: Reconstructing the original reads'

  ; Load average dCounts image
  FITS_READ,outfile_uncmp,avg_dcounts,head0,exten=0,/no_abort,/noscale
  junk=sxpar(head0,'CHECKSUM',count=dcount)
  sz0 = size(avg_dcounts)
  if dcount gt 0 and not keyword_set(no_checksum) then begin
     print,'checking checksum 0'
     res = FITS_TEST_CHECKSUM(head0,avg_dcounts,errmsg=errmsg)
     if res eq -1 then begin
         error = '        Checksum failed on '+files+' (ext=0)'
         if not keyword_set(silent) then print,error
         return
     endif
  endif

  ; Load read=1 (first one)
  FITS_READ,outfile_uncmp,read1,head1,exten=1,/no_abort,/noscale
  sz1 = size(read1)
  junk=sxpar(head1,'CHECKSUM',count=dcount)
  if dcount gt 0 and not keyword_set(no_checksum) then begin
     print,'checking checksum 1'
     res = FITS_TEST_CHECKSUM(head1,read1,errmsg=errmsg)
     if res eq -1 then begin
         error = '        Checksum failed on '+files+' (ext=1)'
         if not keyword_set(silent) then print,error
         return
     endif
  endif

  ; Check that image dimensions of AVG_DCOUNTS and READ1 match
  if sz0[1] ne sz1[1] or sz0[2] ne sz1[2] then begin
    error = '         Images dimensions of AVERAGE DCOUNTS (in exten=0) and READ1 (in exten=1) do NOT MATCH'
    if not keyword_set(silent) then print,error
    FILE_DELETE,outfile_uncmp,/allow  ; delete temporary file
    return
  endif

  ; Write primary HDU
  sxaddpar,head0,'SIMPLE','T',''
  sxaddpar,head0,'BITPIX',16,''
  sxaddpar,head0,'NAXIS',0
  sxdelpar,head0,'NAXIS1'
  sxdelpar,head0,'NAXIS2'
  sxdelpar,head0,'PCOUNT'
  sxdelpar,head0,'GCOUNT'
  sxdelpar,head0,'CHECKSUM'
  sxdelpar,head0,'DATASUM'
  sxdelpar,head0,'BZERO'
  sxdelpar,head0,'BSCALE'
  fits_add_checksum, head0, /no_timestamp
  MWRFITS,0,finalfile,head0,/silent,/create,/no_comment

  ; Write first read
  sxaddpar,head1,'XTENSION','IMAGE','',before='SIMPLE'
  bitpix = sxpar(head1,'BITPIX',count=bxcount)
  if bxcount gt 0 then sxaddpar,head1,'BITPIX',bitpix,''
  sxaddpar,head1,'NAXIS',2
  sxaddpar,head1,'NAXIS1',n_elements(read1[*,0]),'', after='NAXIS'
  sxaddpar,head1,'NAXIS2',n_elements(read1[0,*]),'', after='NAXIS1'
  sxaddpar,head1,'PCOUNT',0,'', after='NAXIS2'
  sxaddpar,head1,'GCOUNT',1,'', after='PCOUNT'
  sxdelpar,head1,'SIMPLE'      ; delete SIMPLE if present, only allowed in PDU
  sxdelpar,head1,'CHECKSUM'
  sxdelpar,head1,'DATASUM'
  FITS_ADD_CHECKSUM, head1, read1, /no_timestamp
  MWRFITS,read1,finalfile,head1,/silent, /no_comment       ; write first read

  ; Loop through extensions and add them together
  lastim = read1
  For i=2,nreads do begin

    ; Read in "residual" image
    FITS_READ,outfile_uncmp,residim,head,exten=i,/no_abort,/noscale
    junk=sxpar(head,'CHECKSUM',count=dcount)
    sz = size(residim)
    if dcount gt 0 and not keyword_set(no_checksum) then begin
       print,'checking checksum',i
       res = FITS_TEST_CHECKSUM(head,residim,errmsg=errmsg)
       if res eq -1 then begin
           error = '         Checksum failed on '+files+' (ext='+strtrim(i,2)+')'
           if not keyword_set(silent) then print,error
           return
       endif
    endif

    ; Check that the image dimension is correct
    if sz[1] ne sz1[1] or sz[2] ne sz1[2] then begin
      error = '         Images dimensions of READ1 (in exten=1) and RESID'+strtrim(i-1,2)+' (in exten='+strtrim(i,2)+') do NOT MATCH'
      if not keyword_set(silent) then print,error
      FILE_DELETE,outfile_uncmp,/allow  ; delete temporary file
      return
    endif

    ; Re-construct the original counts
    ;----------------------------------
    ;  This is how the dcounts/resid were created:
    ;    dcounts[i] = read[i+1]-read[i]
    ;    resid[i] = dcounts[i] - avg_dcounts
    ;  So, adding avg_dcounts to resid gives back dcounts
    ;  and you just keep adding dCounts to the last read to
    ;  reconstruct all of the reads.
    origim = long(lastim) + residim + avg_dcounts
    origim = uint(origim)            ; must be unsigned integer

    ; Fix header
    sxaddpar,head,'BITPIX',16,''     ; unsigned integer
    sxaddpar,head,'XTENSION','IMAGE','',before='SIMPLE'
    sxaddpar,head,'NAXIS',2
    sxaddpar,head,'NAXIS1',n_elements(origim[*,0]),'', after='NAXIS'
    sxaddpar,head,'NAXIS2',n_elements(origim[0,*]),'', after='NAXIS1'
    sxaddpar,head,'PCOUNT',0,'', after='NAXIS2'
    sxaddpar,head,'GCOUNT',1,'', after='PCOUNT'
    sxdelpar,head,'SIMPLE'      ; delete SIMPLE if present, only allowed in PDU
    sxaddpar,head,'BZERO',32768,''
    sxaddpar,head,'BSCALE',1,''
    sxdelpar,head,'CHECKSUM'
    sxdelpar,head,'DATASUM'
    FITS_ADD_CHECKSUM, head, origim, /no_timestamp

    ; Now write the original read
    MWRFITS,origim,finalfile,head,/silent,/no_comment

    ; Save last read
    lastim = origim

  End


; No data to uncompress, Nreads=0
;---------------------------------
Endif else begin

  ; Just copy the file
  FILE_COPY,outfile_uncmp,finalfile,/over,/allow

Endelse


; Delete temporary file
FILE_DELETE,outfile_uncmp,/allow


; Delete original file
if keyword_set(delete) then begin
  if not keyword_set(silent) then $
    print,'Deleting Original file ',files
  FILE_DELETE,files
endif

; remove lock file
file_delete,finalfile+'.lock'

; Time elapsed
dt = systime(1)-t0
if not keyword_set(silent) then $
  print,'dt = ',strtrim(dt,2),' sec'


if keyword_set(stp) then stop

end
