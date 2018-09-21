pro apqbundle,frameid0,indir=indir,outdir=outdir,compress=compress,clobber=clobber,$
             outfiles=outfiles,error=error,stp=stp

;+
;
; APQBUNDLE
;
; This bundles the raw APOGEE frames and optionally compresses them
;
; The data will be output in this format:
;  HDU0: header but NO data
;  HDU1: header, read1 image  as UNSIGNED INTEGERS (BITPIX=16 or UINT)
;  HDU2: header, read2 image  as UNSIGNED INTEGERS (BITPIX=16 or UINT)
;  and so on for all the reads.
; 
; If /compress is set then APZIP will be used to compress the raw
; APOGEE data.  Use APUNZIP to uncompress the data.
;
; INPUTS:
;  frameid    The 8 character frame ID.
;  =indir     The input directory with the apRaw files.  By default
;               indir=/net/stream/apogee/data/raw/
;  =outdir    The output directory for the apR-[abc]-frameid.fits files
;               By default outdir=indir.
;  /compress  Compress the output files.
;  /clobber   Overwrite files if they already exist
;  /stp       Stop at the end of the program.
;
; OUTPUTS:
;  The raw bundles APOGEE files are names outdir/apR-[abc]-frameid.fits
;  If /compress is set then compressed files will also be created and
;  have names of outdir/apR-[abc]-frameid.apz.
;  =outfiles   The output files names.  Normally this is the bundles
;                files, but if /compress is set then this is the
;                compressed filenames.
;
; USAGE:
;  IDL>apqbundle,'00000117'
;
; By D.Nidever  Nov. 2010
;-

t0 = systime(1)

; Error Handling
;------------------
; Establish error handler. When errors occur, the index of the  
; error is returned in the variable Error_status:  
; CATCH, Error_status 

;This statement begins the error handler:  
; if (Error_status ne 0) then begin 
;    error = !ERROR_STATE.MSG  
;    if not keyword_set(silent) then print,error
;    CATCH, /CANCEL 
;    return
; endif


; Not enough inputs
if n_elements(frameid0) eq 0 then begin
  print,'Syntax - apqbundle,frameid,indir=indir,outdir=outdir,compress=compress,'
  print,'                  outfiles=outfiles,clobber=clobber,stp=stp'
  error = 'Not enough inputs'
  return
endif

npix = 2048L

; Get APOGEE directories if the environment variables exists
;data_dir = APGETDIR('APQLDATA_DIR',/exists,error=direrr)
;spectro_dir = APGETDIR('APQLSPECTRO_DIR',/exists,error=direrr)
;archive_dir = APGETDIR('APQLARCHIVE_DIR',/exists,error=direrr)

; ID8
frameid = string(long(frameid0),format='(I08)')

if n_elements(indir) eq 0 then begin
  daynum = strmid(frameid,0,4)
  mjd5 = string(long(daynum)+55562,format='(I5)')
  indir = data_dir+mjd5+'/'
endif
if n_elements(outdir) eq 0 then begin
  daynum = strmid(frameid,0,4)
  mjd5 = string(long(daynum)+55562,format='(I5)')
  outdir = archive_dir+mjd5+'/'
  ;outdir=indir
endif


; Find the files
files = file_search(indir+'/apRaw-'+frameid+'-???.fits',count=nfiles)
if nfiles eq 0 then begin
  print,'No files for ',frameid
stop
  return
endif
base = file_basename(files,'.fits')
num = strmid(base,15,3)
maxind = where(long(num) eq max(num),nmaxind)
print,'Nfiles = ',strtrim(nfiles,2),' for ',frameid

; Read in the FIRST frame and initialize the output files
;  this way the timestamp will be from the BEGINNING of
;  the exposure
minind = where(long(num) eq min(long(num)),nminind)
firstfile = files[minind[0]]
FITS_READ,firstfile,im,head,exten=0,/no_abort,message=message

chiptag = ['a','b','c']
outfiles = addslash(outdir)+'apR-'+chiptag+'-'+frameid+'.fits'
sz = size(im)

; Final Xsize
case sz[1] of
  ; Regular array output 3x2048
  6144: xsize = 2048L
  ; Extra reference 5th output, 4x2048
  8192: xsize = 2048L+512L
  else: begin
    error = 'Non-standard Xsize = '+strtrim(sz[1],2)+' not supported'
    print,error
    return
  end
endcase


; Check if the output files exist
test = file_test(outfiles)
if max(test) eq 1 and not keyword_set(clobber) then begin
  error = 'Output files exist already.  Set /clobber to overwrite'
  print,error
  return
endif

print,'Writing bundled file to ',addslash(outdir)+'apR-[abc]-'+frameid+'.fits'


; Initialize the three bundled files
ahead = head
sxaddpar,ahead,'SIMPLE','T',''
sxaddpar,ahead,'BITPIX',16,''
sxaddpar,ahead,'NAXIS',0
sxdelpar,ahead,'NAXIS1'
sxdelpar,ahead,'NAXIS2'
sxdelpar,ahead,'PCOUNT'
sxdelpar,ahead,'GCOUNT'
sxaddpar,ahead,'CHIP','a'
sxaddpar,ahead,'NREAD',nfiles
sxdelpar,ahead,'CHECKSUM'
sxdelpar,ahead,'DATASUM'
sxdelpar,ahead,'BZERO'
sxdelpar,ahead,'BSCALE'
fits_add_checksum, ahead, /no_timestamp
MWRFITS,0,outfiles[0],ahead,/silent,/create,/no_comment

bhead = ahead
sxaddpar,bhead,'CHIP','b'
sxdelpar,bhead,'CHECKSUM'
sxdelpar,bhead,'DATASUM'
fits_add_checksum, bhead, /no_timestamp
MWRFITS,0,outfiles[1],bhead,/silent,/create,/no_comment

chead = ahead
sxaddpar,chead,'CHIP','c'
sxdelpar,chead,'CHECKSUM'
sxdelpar,chead,'DATASUM'
fits_add_checksum, chead, /no_timestamp
MWRFITS,0,outfiles[2],chead,/silent,/create,/no_comment

; Read in the frames
For fpos=0,nfiles-1 do begin

  apgundef,message,im,head
  ; to use the fits_test_checksum, the data has to be read as a signed int first
  ; this is because the pyfits checksum is done on signed values
  FITS_READ,files[fpos],im,head,exten=0,/no_abort,/no_unsigned,/noscale,/no_pdu,message=message
  res=FITS_TEST_CHECKSUM(head, im, errmsg=errmsg)
  if res eq -1 then begin
      ; the checksum doesn't match -> send a warning and keep going
      message,/con,'BAD checksum for file '+files[fpos]
      return
  endif
  ; convert back to uint
  bzero=sxpar(head,'bzero',count=bzcount)
  if bzcount gt 0 and bzero ne 0 then im = uint(temporary(im)) - uint(bzero[0])
  bscale=sxpar(head,'bscale',count=bscount)
  if bscount gt 0 and bscale ne 1.0 then im = temporary(im) / bscale[0]

  ; The image is: red, green, blue
  ; if the 5th reference output is also there then
  ; this is: green-window (not used), blue-ref, green-ref, red-ref

  ; The RAW images:
  ; Basic orientation is long wavelengths on the left, to shorter to the
  ; right (higher column number).  readout direction is:
  ;         Red                       Green                 Blue           WG RefB RefG RefR
  ;|---->|<----|--->|<----||---->|<----|--->|<----||---->|<----|--->|<----||---->|<----|--->|<----|


  ; Chip a, RED
  ;--------------
  im1 = im[0:npix-1,*]
  ; Add the 5th output, 4th column
  ;  need to flip in X-direction so all three chips
  ;  have the reference reading out left-to-right (--->)
  if xsize eq 2560 then begin
    refim1 = im[6144+3*512:6144+4*512-1,*]
    refim1 = reverse(refim1,1)  ; flip X-direction
    im1 = [im1,refim1]
  endif
  ahead=head
  sxaddpar,ahead,'XTENSION','IMAGE','',before='SIMPLE'
  bitpix = sxpar(ahead,'BITPIX',count=bxcount)
  if bxcount gt 0 then sxaddpar,ahead,'BITPIX',bitpix,''
  sxaddpar,ahead,'NAXIS',2
  sxaddpar,ahead,'NAXIS1',xsize,'', after='NAXIS'
  sxaddpar,ahead,'NAXIS2',n_elements(im1[0,*]),'', after='NAXIS1'
  sxaddpar,ahead,'PCOUNT',0,'', after='NAXIS2'
  sxaddpar,ahead,'GCOUNT',1,'', after='PCOUNT'
  sxdelpar,ahead,'SIMPLE'      ; delete SIMPLE if present, only allowed in PDU
  sxaddpar,ahead,'CHIP','a'
  sxaddpar,ahead,'NREAD',nfiles
  sxdelpar,ahead,'CHECKSUM'
  sxdelpar,ahead,'DATASUM'
  ; re-write the bzero and bscale to add the '/' that mwrfits automatically adds
  ; so as to not mess up the checksum
  if bzcount gt 0 then sxaddpar,ahead,'BZERO',bzero
  if bscount gt 0 then sxaddpar,ahead,'BSCALE', bscale

  FITS_ADD_CHECKSUM, ahead, im1, /no_timestamp
  MWRFITS,im1,outfiles[0],ahead,/silent, /no_comment

  ; Chip b, GREEN
  ;---------------
  im2 = im[npix:2*npix-1,*]
  ; Add the 5th output, 3rd column
  ;  the reference readout direction is already left-to-right
  if xsize eq 2560 then begin
    refim2 = im[6144+2*512:6144+3*512-1,*]
    im2 = [im2,refim2]
  endif
  bhead=ahead
  sxaddpar,bhead,'NAXIS',2
  sxaddpar,bhead,'NAXIS1',xsize,'',after='NAXIS'
  sxaddpar,bhead,'NAXIS2',n_elements(im2[0,*]),'',after='NAXIS1'
  sxaddpar,bhead,'PCOUNT',0,'',after='NAXIS2'
  sxaddpar,bhead,'GCOUNT',1,'',after='PCOUNT'
  FITS_ADD_CHECKSUM, bhead, im2, /no_timestamp
  MWRFITS,im2,outfiles[1],bhead,/silent, /no_comment

  ; Chip c, BLUE
  ;--------------
  im3 = im[2*npix:3*npix-1,*]
  ; Add the 5th output, 2nd column
  ;  need to flip in X-direction so all three chips
  ;  have the reference reading out left-to-right (--->)
  if xsize eq 2560 then begin
    refim3 = im[6144+512:6144+2*512-1,*]
    refim3 = reverse(refim3,1)  ; flip X-direction
    im3 = [im3,refim3]
  endif
  chead=ahead
  sxaddpar,chead,'NAXIS',2
  sxaddpar,chead,'NAXIS1',xsize,'', after='NAXIS'
  sxaddpar,chead,'NAXIS2',n_elements(im3[0,*]),'', after='NAXIS1'
  sxaddpar,chead,'PCOUNT',0,'', after='NAXIS2'
  sxaddpar,chead,'GCOUNT',1,'', after='PCOUNT'
  FITS_ADD_CHECKSUM, chead, im3, /no_timestamp
  MWRFITS,im3,outfiles[2],chead,/silent, /no_comment

End ; read loop

; Compress
;----------
if keyword_set(compress) then begin

  print,'Compressing ',addslash(outdir)+'apR-[abc]-'+frameid+'.fits'
  APZIP,outfiles

  ; Return compressed filenames
  bundle_files = outfiles
  outfiles = file_dirname(bundle_files)+'/'+file_basename(bundle_files,'.fits')+'.apz'

end

dt = systime(1)-t0
print,'dt = ',strtrim(dt,2),' sec'

; 17 sec for a 60 read file.
; 204 sec = 3.4 min for 60 read file WITH compression

if keyword_set(stp) then stop

end
