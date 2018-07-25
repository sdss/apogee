function apg_checksum,input,silent=silent,fitsdir=fitsdir

;+
;
; APG_CHECKSUM
;
; This program uncompresses the APOGEE bundled and compressed apz files
; and verifies the checksum (and datasum) for every extension.
;
; This program assumes the files were created with apzip with the 
; fpack flag -C which doesn't update the value of the checksum and
; datasum to reflect the new compressed data.
;
; INPUTS:
;  input     A list of input compressed raw bundled APOGEE fits files
;
; OUTPUTS:
;  The program returns the resulting value of the FITS_TEST_CHECKSUM
;  for every input file (one per file).
;
; USAGE:
;  IDL> results = APG_CHECKSUM('apR-a-00000085.apz')
;
; By S.Beland   August 2011
;-

t0 = systime(1)
if not keyword_set(silent) then silent=0

; Not enough inputs
ninput = n_elements(input)
if ninput eq 0 then begin
  message,/con, 'Not enough inputs'
  message,/con,'Syntax - apunzip,input,clobber=clobber,delete=delete,silent=silent,error=error'
  return,-2
endif

; Get the inputs
LOADINPUT,input,files,count=nfiles


; More than one file input
results = intarr(nfiles)
if nfiles gt 1 then begin
  for i=0,nfiles-1 do begin
      results[i] = apg_checksum(files[i],silent=silent)
      if not silent then print,files[i]+'  ->  ',strtrim(results[i],2)
  endfor
  return,results
endif

; Does the file exist
if file_test(files) eq 0 then begin
  message,/con,files+' NOT FOUND'
  return,-2
endif

; Check that "funpack" is available
spawn,['funpack','-H'],out,errout,/noshell
if errout[0] ne '' then begin
  message,/con,'FUNPACK not found'
  return,-2
endif

; Check that the extension is ".apz"
dir = file_dirname(files)+'/'
fil = file_basename(files)
dum = strsplit(fil,'.',/extract)
ext = dum[n_elements(dum)-1]
if ext ne 'apz' then begin
  message,/con,'Extension must be .apz'
  return,-2
endif

; Temporary directory
if keyword_set(fitsdir) then tempdir = fitsdir else $
  tempdir = FILE_DIRNAME(files)

; Final output filename
tempfile = MKTEMP('apzip',outdir=tempdir)
FILE_DELETE,tempfile,/allow

; Step I: Uncompress the file with funpack
;-------------------------------------------
spawn,['funpack','-O',tempfile,'-C',files],out,errout,/noshell  ; -C suppresses checksum update
if errout[0] ne '' then begin
  message,/con,'fpack error '+errout
  return,-2
endif

; Get number of extensions
fits_open,tempfile,fcb
nreads=fcb.nextend
fits_close,fcb

; Read every extensions and calculate the checksum and datasum
;-------------------------------------------------------------
res = 1
for ext=0,nreads do begin
    fits_read,tempfile,d0,h0,exten_no=ext,/no_pdu,/noscale
    res = min([res,FITS_TEST_CHECKSUM(h0,d0,errmsg=errmsg)])
endfor

; cleanup after ourself
FILE_DELETE,tempfile,/allow

; Time elapsed
dt = systime(1)-t0
if not silent then print,'dt = ',strtrim(dt,2),' sec'

return, res

end
