;+
;
; AP2DPROC
;
; This program extracts a 2D APOGEE image.
; This is called from AP2D
;
; INPUTS:
;  inpfile      The name of the 2D APOGEE file.  This
;                 should be the directory and the ID8
;                 concatenated.
;  psffile      The name of the calibration PSF file to
;                 use. This should also be the directory
;                 and "base" name or ID concatenated.
;  extract_type The extraction method to use:
;                 1-Boxcar extraction (the default)
;                 2-PSF image extraction
;                 3-Gaussian PSF fitting extraction
;                 4-Jon's Empirical PSF extraction
;                 5-Full Gaussian-Hermite PSF fitting extraction
;  =outdir      The output directory.  By default the 1D extracted
;                 files are written to the same directory that
;                 the input 2D files are in.
;  =fluxcalfile The name of the relative flux calibration file to use.
;                 This should also be the directory and "base" name
;                 or ID concatenated.
;  =wavefile    The name of the wavelength calibration to use to
;                 add wavelengths to the output file.  This should
;                 also be the directory and "base" name or ID concatenated.
;  /skywave     To enable a pixel-shift to wavelength solution based on sky lines
;  =plugmap     To specify a plugmap for the sky-line wavelength solution;
;                 if plugmap is given, only use SKY fibers for the correction
;  =highrej     High rejection threshold for Gaussian PSF fitting
;                 The default is 7.
;  =lowrej      Low rejection threshold for Gaussian PSF fitting
;                 The default is 10.
;  =npolyback   The number of polynomial coeffiecients to use for
;                 the background.  Only for extract_type=3 for now.
;                 The default is npolyback=0.
;  /recenterfit Recenter the traces/PSF with one constant offset for
;                 each chip.  The shift is found empirically from
;                 the image itself.
;  /recenterln2 Recenter the traces/PSF with one constant offset for
;                 all chip.  The shift is determined by the LN2LEVEL
;                 header value in the image compared to the PSF file
;                 and an empirically derived relation between the
;                 LN2LEVEL and fiber trace shifts.  The LN2 level
;                 (due to gravity) slightly warps the instrument and
;                 shifts the traces (by fractions of a pixel).
;  /refpixzero  Set the image zeropoint using the reference pixels.
;  =fibers      Array of fibers to extract (0 for first fiber).
;  =chips       Array of chips to use (0 for first chip).
;  /fitsigma    Allow the sigma to vary for extract_type=3
;  /fixbadpix   Fix bad pixels using 2D interpolation of neighboring
;                 pixels.  This is the default for extract_type 1 and 2
;  /outlong     The output files should use LONG type intead of FLOAT.
;                 This actually takes up the same amount of space, but
;                 this can be losslessly compressed with FPACK.
;  /nowrite     Don't write the output to files.
;  /clobber     Overwrite existing files.
;  /verbose     Print a lot of information to the screen
;  /silent      Don't print anything to the screen
;  /stp         Stop at the end of the prrogram
;
; OUTPUTS:
;  1D extracted spectra are output.  One file for each frame.
;  =output      Structure with the extracted spectra for all three chips.
;  =outmodel    Structure with the model of the 2D image for all three chips.
;
; USAGE:
;  IDL>ap2dproc,inpfile,tracefile,outdir,1
;
; Written by D.Nidever  July 2010
;-

pro ap2dproc,inpfile,psffile,extract_type,outdir=outdir,clobber=clobber,fixbadpix=fixbadpix,$
             fluxcalfile=fluxcalfile,responsefile=responsefile,wavefile=wavefile,skywave=skywave,$
             plugmap=plugmap,highrej=highrej,lowrej=lowrej,$
             verbose=verbose,silent=silent,recenterfit=recenterfit,recenterln2=recenterln2,$
             fitsigma=fitsigma,refpixzero=refpixzero,outlong=outlong,output=output,outmodel=outmodel,$
             nowrite=nowrite,npolyback=npolyback,chips=chips,fibers=fibers,stp=stp,compress=compress

common savedepsf, savedepsffiles, epsfchip
if n_elements(savedepsfiles) eq 0 then savedepsffiles=[' ',' ',' ']  ; initialize if needed
if n_elements(epsfchip) eq 0 then epsfchip=0

if n_elements(verbose) eq 0 then verbose=0  ; NOT verbose by default

ninpfile = n_elements(inpfile)
npsffile = n_elements(psffile)

; Not enough inputs
if ninpfile eq 0 or npsffile eq 0 then begin
  if not keyword_set(silent) then begin
    print,'Syntax - ap2dproc,inpfile,psffile,extract_type,outdir=outdir,clobber=clobber,fixbadpix=fixbadpix,'
    print,'                  recenterfit=recenterfit,recenterln2=recenterln2,fluxcalfile=fluxcalfile,responsefile=responsefile,'
    print,'                  refpixzero=refpixzero,wavefile=wavefile,outlong=outlong,outmodel=outmodel,'
    print,'                  verbose=verbose,fitsigma=fitsigma,nowrite=nowrite,npolyback=npolyback,chips=chips,'
    print,'                  fibers=fibers,silent=silent,stp=stp'
  endif
  error = 'Not enough inputs'
  return
endif

; More than one file input
if ninpfile gt 1 then begin
  error = 'Only ONE FILE can be input at a time'
  if not keyword_set(silent) then print,error
  return
endif

; Default parameters
if n_elements(extract_type) eq 0 then extract_type=1   ; boxcar by default
;if n_elements(fixbadpix) eq 0 and (extract_type lt 3) then fixbadpix=1         ; fix bad pixels by default
if n_elements(fixbadpix) eq 0 then fixbadpix=0       ; don't fix bad pixels by default
if n_elements(outlong) eq 0 then outlong=0           ; use FLOAT by default
if n_elements(recenterfit) eq 0 then recenterfit=0   ; do NOT recenter by default
if n_elements(recenterln2) eq 0 then recenterln2=0   ; do NOT recenter by default
if n_elements(fitsigma) eq 0 then fitsigma=0         ; dot not fit sigma by default
if n_elements(refpixzero) eq 0 then refpixzero=0     ; no resetting of zeropoint by default
if n_elements(nowrite) eq 0 then nowrite=0           ; write to file by default
if n_elements(npolyback) eq 0 then npolyback=0       ; no background by default
if n_elements(chips) eq 0 then chips=indgen(3)       ; extract all chips by default
if not keyword_set(plugmap) then plugmap=0           ; default no plugmap file for ap1dwavecal

; Output directory
if n_elements(outdir) eq 0 then outdir=file_dirname(inpfile)+'/'
dirs=getdir()

chiptag = ['a','b','c']

;; OUTDIR must be a string
;if size(outdir,/type) ne 7 then begin
;  print,'OUTDIR must be a string'
;  return
;endif
; Does the output directory exist?
if file_test(outdir,/directory) eq 0 then begin
  if not keyword_set(silent) then print,'' & print,'Creating ',outdir
  FILE_MKDIR,outdir
endif

; No output requested
if keyword_set(nowrite) and arg_present(output) eq 0 and arg_present(outmodel) eq 0 then begin
  error = 'No output requested.'
  if not keyword_set(silent) then print,error
  return
endif

; Chips to extract
if n_elements(chips) gt 3 or min(chips) lt 0 or max(chips) gt 2 then begin
  error = 'CHIPS must have <=3 elements with values [0-2].'
  if not keyword_set(silent) then print,error
  return
endif

; Fibers to extract
if n_elements(fibers) gt 0 then $
  if n_elements(fibers) gt 300 or min(fibers) lt 0 or max(fibers) gt 299 then begin
  error = 'FIBERS must have <=300 elements with values [0-299].'
  if not keyword_set(silent) then print,error
  return
endif

; Make the filenames and check the files
dir = file_dirname(inpfile)
base = file_basename(inpfile)
if file_test(dir,/directory) eq 0 then begin
  error = 'DIRECTORY '+dir+' NOT FOUND'
  if not keyword_set(silent) then print,error
  return
endif
baseframeid = string(long(base),format='(I08)')
files = dir+'/'+dirs.prefix+'2D-'+chiptag+'-'+baseframeid+'.fits'
info = APFILEINFO(files,/silent)
framenum = info[0].fid8   ; the frame number
okay = (info.exists AND info.sp2dfmt AND info.allchips AND ((info.naxis eq 3) OR (info.exten eq 1)))
if min(okay) lt 1 then begin
  bd = where(okay eq 0,nbd)
  error = 'There is a problem with files: '+strjoin((files)(bd),' ')
  if not keyword_set(silent) then print,'halt: '+error
  stop
  return
endif

; Get PSF info
psf_dir = file_dirname(psffile)
psf_base = file_basename(psffile)
if file_test(psf_dir) eq 0 then begin
  error = 'PSF DIRECTORY '+psf_dir+' NOT FOUND'
  if not keyword_set(silent) then print,'halt: '+error
  stop
  return
endif
psfframeid = string(long(psf_base),format='(I08)')
psffiles = apogee_filename('PSF',num=psfframeid,chip=chiptag)
epsffiles = apogee_filename('EPSF',num=psfframeid,chip=chiptag)
pinfo = APFILEINFO(psffiles,/silent)
pokay = (pinfo.exists AND pinfo.allchips)
if min(pokay) lt 1 then begin
  pbd = where(pokay eq 0,nbd)
  error = 'There is a problem with PSF files: '+strjoin((psffiles)(pbd),' ')
  if not keyword_set(silent) then print,'halt: '+error
  stop
  return
endif

if not keyword_set(silent) then begin
  print,''
  print,'Extracting file ',inpfile
  print,'--------------------------------------------------'
endif

; Parameters
mjd5 = info[0].mjd5
nreads = info[0].nreads
if not keyword_set(silent) then print,'MJD5 = ',strtrim(mjd5,2)

; Check header
head = headfits(files[0],errmsg=errmsg)
if errmsg ne '' then begin
  error = 'There was an error loading the HEADER for '+file
  if not keyword_set(silent) then print,'halt: '+error
  stop
  return
endif

; Determine file TYPE
;----------------------
; dark - should be processed with 
; flat
; lamps
; object frame
;obstype = SXPAR(head,'OBSTYPE',count=nobs)
imagetyp = SXPAR(head,'IMAGETYP',count=nimagetyp)
if nimagetyp eq 0 then begin
  error = 'NO IMAGETYP keyword found for '+baseframeid
  if not keyword_set(silent) then print,error
  ;return
endif
;obstype = strlowcase(strtrim(obstype,2))
imagetyp = strlowcase(strtrim(imagetyp,2)) 

; Load the frame
APLOADFRAME,files,frame,error=loaderror
if n_elements(loaderror) ne 0 then return

; Double-check the flux calibration file
if n_elements(fluxcalfile) gt 0 then begin
  fluxcalfiles = file_dirname(fluxcalfile)+'/'+dirs.prefix+'Flux-'+chiptag+'-'+file_basename(fluxcalfile)+'.fits'
  ftest = file_test(fluxcalfiles)
  if total(ftest) lt 3 then begin
    error = 'Problems with FLUX CALIBRATION file '+fluxcalfile
    if not keyword_set(silent) then print,'halt: '+error
    stop
    return
  endif
endif

; Double-check the response calibration file
if n_elements(responsefile) gt 0 then begin
  responsefiles = file_dirname(responsefile)+'/'+dirs.prefix+'Response-'+chiptag+'-'+file_basename(responsefile)+'.fits'
  ftest = file_test(responsefiles)
  if total(ftest) lt 3 then begin
    error = 'Problems with RESPONSE CALIBRATION file '+responsefile
    if not keyword_set(silent) then print,'halt: '+error
    stop
    return
  endif
endif

; Double-check the wave calibration file
if keyword_set(wavefile) then begin
  wavefiles = file_dirname(wavefile)+'/'+dirs.prefix+'Wave-'+chiptag+'-'+file_basename(wavefile)+'.fits'
  wtest = file_test(wavefiles)
  if total(wtest) lt 3 then begin
    error = 'Problems with WAVELENGTH file '+wavefile
    if not keyword_set(silent) then print,'halt: '+error
    stop
    return
  endif

endif

; Wait if another process is working on this
lockfile = outdir+dirs.prefix+'1D-'+framenum ; lock file
if getlocaldir() then lockfile=getlocaldir()+'/'+dirs.prefix+'1D-'+framenum+'.lock' $
else lockfile=outdir+dirs.prefix+'1D-'+framenum+'.lock'

while file_test(lockfile) do apwait,lockfile,10

err=1
while err ne 0 do begin
  file_mkdir,file_dirname(lockfile)
  openw,lock,/get_lun,lockfile, error=err
  IF (err NE 0) then begin
     PRINTF, -2, !ERROR_STATE.MSG
    print,file_search('/scratch/local','*')
  endif
endwhile
free_lun,lock

; Since final ap1dwavecal requires simultaneous fit of all three chips, and
;  this required final output to be put off until after all chips are done,
;  all 3 need to be done here if any at all, so that data from all chips is loaded
; Output files
outfiles = outdir+dirs.prefix+'1D-'+chiptag+'-'+framenum+'.fits'  ; output file
outtest = FILE_TEST(outfiles)
if min(outtest) eq 0 then clobber=1 
if not keyword_set(clobber) then goto,apend

;----------------------------------
; Looping through the three chips
;----------------------------------
head_chip=strarr(n_elements(chips),5000) ;initialise an array to store the headers
apgundef,output,outmodel,outstr
ifirst=0
For i=0,n_elements(chips)-1 do begin

  t1=systime(/seconds)
  ichip = chips[i]   ; chip index, 0-first chip

  file = files[ichip]

  ; The chip structure
  chstr = frame.(ichip)

  ; Chip Trace filename
  ipsffile = psffiles[ichip]
  iepsffile = epsffiles[ichip]

  ; Output file
  outfile = outdir+dirs.prefix+'1D-'+chiptag[ichip]+'-'+framenum+'.fits'  ; output file

  if not keyword_set(silent) then begin
    if i gt 0 then print,''
    print,' Processing chip '+chiptag[ichip]+' - '+file_basename(file)
    print,'  PSF file = ',ipsffile
  endif

  ; Fix the Bad pixels and "Unfixable" pixels
  ;---------------------------------------------
  if keyword_set(fixbadpix) then $
    AP2DPROC_FIXPIX,chstr


  ;###############################################################
  ; NEED TO REMOVE THE LITTROW GHOST AND SECONDARY GHOST HERE!!!!!!!!
  ;###############################################################


  ; Restore the TRACE structure
  tracestr = MRDFITS(ipsffile,1,/silent)

  ; Fibers to extract
  ;if n_elements(fibers) eq 0 then fibers=indgen(n_elements(tracestr))
  if n_elements(fibers) gt 0 then $
    if max(fibers) gt n_elements(tracestr)-1 then begin
    error = 'MAX(Fibers) is larger than the number of fibers in PSF file.'
    if not keyword_set(silent) then print,'halt: '+error
    stop
    return
  endif


  ; Measuring the trace shift
  if keyword_set(recenterfit) then begin
    im = frame.(ichip).flux
    sz = size(im)
    npix = sz[1]
    nfibers = n_elements(tracestr)
    ; the red chip has problems on the left side,
    ;  so use columns farther to the right
    if ichip eq 0 then xmid=npix*0.75 else xmid=npix*0.5

    medspec = median(im[xmid-50:xmid+50,*],dim=1)
    gdpix = where(medspec gt 0.5*max(medspec),ngdpix)
    if ngdpix le 20 then begin
        ; we're probably trying to process a dark as object or flat
        ; i'm not sure if 20 is the right number but seems to work with darks
        if not keyword_set(silent) then print,'No signal was seen on any fibers for chip ',ichip
        xshift = 0.0d
    endif else begin
        medht = median(medspec[gdpix]) > 0.5*max(medspec)

        tpar = fltarr(nfibers*3)
        yfib = fltarr(nfibers)
        for l=0,nfibers-1 do yfib[l]=poly(xmid,tracestr[l].coef)
        tpar[0:3*nfibers-3:3] = medht
        ;tpar[1:3*nfibers-2:3] = xsol[xmid,*]
        tpar[1:3*nfibers-2:3] = yfib
        ;tpar[2:3*nfibers-1:3] = median(sigma2[xmid,*])
        tpar[2:3*nfibers-1:3] = 1.0 ;1.5
        x = findgen(npix)
        temp = gfunc(x,tpar)
        mask1d = long(medspec gt 0.5*max(medspec))
        ;XCORLB,temp,medspec,20,xsh,mask=mask1d

        lag = findgen(9)-4
        xc = c_correlate(temp,medspec*mask1d,lag)
        bestind = first_el(maxloc(xc))
        fitlo = (bestind-2) > 0
        fithi = (bestind+2) < 20
        estimates = [xc[bestind],lag[bestind],1,median(xc)]
        yfit = MPFITPEAK(lag[fitlo:fithi],xc[fitlo:fithi],par,nterms=4,/gaussian,/positive,estimates=estimates)
        xshift = par[1]
    endelse
    if not keyword_set(silent) then $
      print,'Recentering shift = ',stringize(xshift,ndec=3)

    ; this is an ADDITIVE OFFSET!

  endif  ; recenterfit

  ; Calculate the trace shift LN2LEVEL header values
  if keyword_set(recenterln2) then begin
    head_psf = headfits(ipsffile,exten=0)
    ln2level_psf = sxpar(head_psf,'LN2LEVEL',count=nln2level_psf)
    ln2level_im = sxpar(chstr.header,'LN2LEVEL',count=nln2level_im)

    if nln2level_psf gt 0 and nln2level_im gt 0 then begin

      ; The slope of trace shift vs. LN2LEVEL is (from green chip):  0.0117597
      ; Fits from check_traceshift.pro
      ; linear: coef=[ -1.02611, 0.0117597]
      ; quadratic:  coef=[-3.33460, 0.0613117, -0.000265449]
      ; A higher LN2LEVEL shifts the fiber DOWNWARDS
      xshift = (ln2level_im - ln2level_psf) * (-0.0117597)      
      if not keyword_set(silent) then $
        print,'Recentering shift = ',stringize(xshift,ndec=3)

      ; this is an ADDITIVE OFFSET!

    ; Don't have LN2LEVELs
    endif else begin
      if nln2level_psf eq 0 and not keyword_set(silent) then print,'Do NOT have header LN2LEVEL for PSF exposure'
      if nln2level_im eq 0 and not keyword_set(silent) then print,'Do NOT have header LN2LEVEL for this exposure'
      if not keyword_set(silent) then print,'CANNOT calculate fiber shift from LN2LEVEL in headers'
    endelse

  endif ; recenterln2

  ; Reset the zeropoint threshold using the reference pixels
  if keyword_set(refpixzero) then begin
    medref = median( [ chstr.flux[*,0:3], transpose(chstr.flux[0:3,*]), chstr.flux[*,2044:2047], transpose(chstr.flux[2044:2047,*]) ])
    if not keyword_set(silent) then print,'Setting image zeropoint using reference pixels.  Subtracting ',strtrim(medref,2)
    chstr.flux -= medref
  endif


  ; Initialize the Output header
  ;------------------------------
  head = chstr.header
  sxaddpar,head,'PSFFILE',ipsffile,' PSF file used'
  leadstr = 'AP2D: '
  sxaddpar,head,'V_APRED',getvers()
  sxaddhist,leadstr+systime(0),head
  info = GET_LOGIN_INFO()
  sxaddhist,leadstr+info.user_name+' on '+info.machine_name,head
  sxaddhist,leadstr+'IDL '+!version.release+' '+!version.os+' '+!version.arch,head
  ; Add Reduction pipeline version to the header
  sxaddhist,leadstr+' APOGEE Reduction Pipeline Version: '+getvers(),head
  sxaddhist,leadstr+'Output File:',head
  sxaddhist,leadstr+' HDU1 - image (ADU)',head
  sxaddhist,leadstr+' HDU2 - error (ADU)',head
  if (extract_type eq 1) then begin
    sxaddhist,leadstr+' HDU3 - flag mask (bitwise OR combined)',head
    sxaddhist,leadstr+'        1 - bad pixels',head
    sxaddhist,leadstr+'        2 - cosmic ray',head
    sxaddhist,leadstr+'        4 - saturated',head
    sxaddhist,leadstr+'        8 - unfixable',head
  endif else begin
    sxaddhist,leadstr+' HDU3 - flag mask',head
    sxaddhist,leadstr+'        0 - good pixels',head
    sxaddhist,leadstr+'        1 - bad pixels',head
  endelse
  if n_elements(wavefile) gt 0 then begin
    sxaddhist,leadstr+' HDU4 - wavelengths (Ang)',head
    sxaddhist,leadstr+' HDU5 - wavelength coefficients',head
  endif

  apgundef,outstr,ymodel

  ; Extraction type
  ;------------------
  CASE extract_type of


    ; BOXCAR Extraction
    ;-------------------
    1: begin
      if not keyword_set(silent) then print,'Using Boxcar Extraction'

      ; Update header
      sxaddhist,leadstr+'Extract_type=1 - Using Boxcar Extraction',head
      sxaddpar,head,'EXTRTYPE',1,'Extraction type'
      
      ; Recenter, shift the traces
      if keyword_set(recenterfit) or keyword_set(recenterln2) then begin
        tracestr.coef[0] += xshift
        tracestr.gaussy += xshift
        if keyword_set(recenterfit) and not keyword_set(recenterln2) then $
          sxaddhist,leadstr+' /RECENTERFIT set, shifting traces by '+stringize(xshift,ndec=3),head
        if keyword_set(recenterln2) then sxaddhist,leadstr+' /RECENTERLN2 set, shifting traces by '+stringize(xshift,ndec=3),head
      endif

      ; Extract the fibers
      APEXTRACT,chstr,tracestr,outstr,fibers=fibers

    end ; boxcar


    ; PSF Image Extraction
    ;----------------------
    2: begin
      if not keyword_set(silent) then print,'Using PSF Image Extraction'

      ; Load the PSF image
      FITS_READ,ipsffile,psfim,head_psfim,exten=2,/no_abort,message=message
      if message ne '' then begin
        if not keyword_set(silent) then print,'PSF file ',ipsffile,' does NOT contain a PSF image'
        goto,BOMB
      endif

      ; Update header
      sxaddhist,leadstr+'Extract_type=2 - Using PSF Image Extraction',head
      sxaddpar,head,'EXTRTYPE',2,'Extraction type'

      ; Recenter, shift the traces and the PSF image
      if keyword_set(recenterfit) or keyword_set(recenterln2) then begin
        tracestr.coef[0] += xshift
        tracestr.gaussy += xshift
        psfim0 = psfim
        psfim = IMDRIZZLE(psfim0,0.0,xshift)  ; shift the image with imdrizzle
        if keyword_set(recenterfit) and not keyword_set(recenterln2) then $
          sxaddhist,leadstr+' /RECENTERFIT set, shifting traces by '+stringize(xshift,ndec=3),head
        if keyword_set(recenterln2) then sxaddhist,leadstr+' /RECENTERLN2 set, shifting traces by '+stringize(xshift,ndec=3),head
      endif

      ; Extract the fibers
      APEXTRACTPSF,chstr,tracestr,psfim,outstr,model=ymodel,fibers=fibers

    end ; PSF image


    ; Gaussian PSF fitting
    ;----------------------
    ;   maybe use the idlspec2d extraction code for this
    3: begin
      if not keyword_set(silent) then print,'Using Gaussian PSF fitting Extraction'

      ; Update header
      sxaddhist,leadstr+'Extract_type=3 - Using Gaussian PSF fitting Extraction',head
      sxaddpar,head,'EXTRTYPE',3,'Extraction type'

      ; THE IDLSPEC2D PROGRAMS EXPECT THE FIBERS TO RUN ALONG THE Y
      ; TRANSPOSING THE ARRAYS FOR NOW

      ; Get the idlspec2d-style trace and widthset information
      FITS_READ,ipsffile,tset_coeff,tset_head,exten=3,message=message1,/no_abort
      tset = {func:strtrim(sxpar(tset_head,'FUNC'),2),xmin:sxpar(tset_head,'XMIN'),$
              xmax:sxpar(tset_head,'XMAX'),coeff:tset_coeff}
      FITS_READ,ipsffile,wset_coeff,wset_head,exten=4,message=message2,/no_abort
      widthset = {func:strtrim(sxpar(wset_head,'FUNC'),2),xmin:sxpar(wset_head,'XMIN'),$
                  xmax:sxpar(wset_head,'XMAX'),coeff:wset_coeff}
      proftype = sxpar(wset_head,'PROFTYPE')

      ; Get the trace and sigma arrays
      apgundef,ycen,xsol,xx,sigma2
      traceset2xy, tset, ycen, xsol
      traceset2xy, widthset, xx, sigma2

      ; Get the images ready
      img = transpose( float(frame.(ichip).flux) )
      ivar = transpose( 1/float(frame.(ichip).err)^2 )
      mask = transpose( frame.(ichip).mask )
      ;mask = 1-(mask eq 1 or mask eq 8)  ; 1-good, 0-bad
      ;mask = 1-(mask eq 1 or mask eq 4 or mask eq 8)  ; 1-good, 0-bad
      ;mask = 1-( ((mask and 1) eq 1) or ((mask and 8) eq 8) )
      mask = 1-( ((mask and 1) eq 1) or ((mask and 4) eq 4) or ((mask and 8) eq 8) )


      ; Recenter the traces
      if keyword_set(recenterfit) or keyword_set(recenterln2) then begin

        ; Need to ADD this to the traces
        xsol += xshift
        if keyword_set(recenterfit) and not keyword_set(recenterln2) then $
          sxaddhist,leadstr+' /RECENTERFIT set, shifting traces by '+stringize(xshift,ndec=3),head
        if keyword_set(recenterln2) then sxaddhist,leadstr+' /RECENTERLN2 set, shifting traces by '+stringize(xshift,ndec=3),head

      endif  ; recentering


      ;---------------------------------------------------------------------
      ; Extract the spectra
      ;---------------------------------------------------------------------

      if n_elements(highrej) eq 0 then highrej=7
      if n_elements(lowrej) eq 0 then lowrej=10
      ;highrej = 7 ; 15
      ;lowrej = 15
      ;npoly = 10
      ;npoly = 15  ; this seems to work better
      ; since the Gaussian is NOT a good fit use a lower
      ;  order background
      ;npoly = 1  ;3
      npoly = npolyback
      ;wfixed = [1,1] ; Just fit the first gaussian term
      wfixed = [1]   ; keep the sigmas fixed
      if keyword_set(fitsigma) then wfixed=[1,1]  ; fit sigma

      ; Only extract FIBERS
      if n_elements(fibers) gt 0 then begin
        xsol = xsol[*,fibers]
        sigma2 = sigma2[*,fibers]
      endif

      ;splog, 'Extracting arc'
      AP_EXTRACT_IMAGE, img, ivar, xsol, sigma2, $
          flux, fluxivar, proftype=proftype, wfixed=wfixed, $
          highrej=highrej, lowrej=lowrej, npoly=npoly, relative=1, $
          reject=[0.1, 0.6, 0.6],ymodel=ymodel,mask=mask,chisq=chisq

      ; Transpose the model
      ymodel = transpose(ymodel)

      ; Create outstr
      ;  bad pixels have fluxivar=0, they are given high ERR
      ;  MASK make it: 0-good, 1-bad
      outstr = {flux:flux, err:1/(sqrt(fluxivar>1d-12)), mask:fix(fluxivar*0)}
      outstr.mask = (fluxivar eq 0)  ; pixels with fluxivar=0 are bad
      ; negative pixels
      ;bd = where(outstr.flux lt 0,nbd)
      ;if nbd gt 0 then begin
      ;  outstr.flux[bd] = 0
      ;  outstr.err[bd] = 1e6
      ;  outstr.mask[bd] = 1  ; maybe give this a different value
      ;endif
      ; fix reference pixels
      outstr.flux[0:3,*] = 0
      outstr.flux[2040:2047,*] = 0
      outstr.err[0:3,*] = baderr()
      outstr.err[2040:2047,*] = baderr()
      outstr.mask[0:3,*] = 1
      outstr.mask[2040:2047,*] = 1

      ;stop

    end ; Gaussian fitting


    ; Empirical PSF Extraction
    ;--------------------------
    4: begin
      if not keyword_set(silent) then print,'Using Empirical PSF Extraction'

      ; Copied from holtz/approcess.pro

      if epsffiles[ichip] ne savedepsffiles[ichip] then begin

       ; Load Empirical PSF data
        tmp=mrdfits(iepsffile,0,phead,status=status_epsf,/silent)
        ntrace=sxpar(phead,'NTRACE')
        img=ptrarr(ntrace,/allocate_heap)
        for itrace=0,ntrace-1 do begin
          ptmp=mrdfits(iepsffile,itrace+1,/silent)
          *img[itrace] = ptmp.img
          p ={fiber: ptmp.fiber, lo: ptmp.lo, hi: ptmp.hi, img: img[itrace]}
          if itrace eq 0 then epsf=replicate(p,ntrace)

          epsf[itrace] = p
        endfor
        if ichip eq 0 then begin
          epsfchip=epsf 
          sz0=size(epsf,/dim)
          sz=sz0
        endif else begin
          sz=size(epsf,/dim)
          if sz eq sz0 then epsfchip=[[epsfchip],[epsf]]
        endelse
        if status_epsf ne 0 then begin
          if not keyword_set(silent) then print,'PSF file ',iepsffile,' does NOT contain an empirical PSF image'
          goto,BOMB
        endif
        if sz eq sz0 then savedepsffiles[ichip] = epsffiles[ichip]
      endif else epsf = reform(epsfchip[*,ichip])

      ; Update header
      sxaddhist,leadstr+'Extract_type=4 - Using Empirical PSF Extraction',head
      sxaddpar,head,'EXTRTYPE',4,'Extraction type'

      ; Recenter, shift the traces and the PSF image
      if keyword_set(recenterfit) or keyword_set(recenterln2) then begin
        ; shift the image with imdrizzle
        for l=0,n_elements(epsf)-1 do epsf[l].img = IMDRIZZLE(epsf[l].img,0.0,xshift)
        if keyword_set(recenterfit) and not keyword_set(recenterln2) then $
          sxaddhist,leadstr+' /RECENTERFIT set, shifting traces by '+stringize(xshift,ndec=3),head
        if keyword_set(recenterln2) then sxaddhist,leadstr+' /RECENTERLN2 set, shifting traces by '+stringize(xshift,ndec=3),head
      endif

      APEXTRACT_EPSF,frame.(ichip),epsf,outstr,model=ymodel,/scat ;,subonly=50*indgen(6)
      ;stop

    end ; EPSF


    ; Full Gaussian-Hermite PSF fitting
    ;-----------------------------------
    5: begin
      ;if not keyword_set(silent) then print,'Using Full Gaussian-Hermite PSF fitting Extraction'
      error = 'Full Gaussian-Hermite PSF fitting is not supported yet'
      if not keyword_set(silent) then print,error
      return
      ;; Update header
      ;sxaddhist,leadstr+'Extract_type=4 - Using Gaussian-Hermite PSF fitting Extraction',head
      ;sxaddpar,head,'EXTRTYPE',4,'Extraction type'
      return
    end


    ; Non-supported options
    else: begin
      error = 'Extraction Type not supported'
      if not keyword_set(silent) then print,error
      return
    end

  ENDCASE

  t2=systime(/seconds)
  ;stop


  ; Do the fiber-to-fiber throughput corrections and relative
  ; Flux calibration
  ;------------------------------------------------------------
  if n_elements(fluxcalfile) gt 0 then begin

    ; Restore the relative flux calibration correction file
    if not keyword_set(silent) then $
      print,'Flux calibrating with ',file_dirname(fluxcalfile)+'/'+dirs.prefix+'Flux-'+file_basename(fluxcalfile)
    fluxcalfiles = file_dirname(fluxcalfile)+'/'+dirs.prefix+'Flux-'+chiptag+'-'+file_basename(fluxcalfile)+'.fits'
    FITS_READ,fluxcalfiles[ichip],fluxcal,fluxcal_head,message=message,/no_abort
    outstr.flux /= fluxcal          ; correct flux
    bderr=where(outstr.err eq baderr(),nbd)
    outstr.err /= fluxcal           ; correct error
    if nbd gt 0 then outstr.err[bderr] = baderr()
    bd=where(finite(outstr.flux) eq 0,nbd)
    if nbd gt 0 then begin
       outstr.flux[bd]=0.
       outstr.err[bd]=baderr()
       outstr.mask[bd]=1
    endif

    ; Update header
    sxaddhist,leadstr+'Flux Calibrating the spectra with:',head
    sxaddhist,leadstr+fluxcalfiles[ichip],head
    sxaddpar,head,'FLUXFILE',fluxcalfile,' Flux Calibration file used'
  endif

  ; Response curve calibration
  ;--------------------
  if n_elements(responsefile) gt 0 then begin

    ; Restore the relative flux calibration correction file
    if not keyword_set(silent) then $
      print,'Response calibrating with ',file_dirname(responsefile)+'/'+dirs.prefix+'Flux-'+file_basename(responsefile)
    responsefiles = file_dirname(responsefile)+'/'+dirs.prefix+'Response-'+chiptag+'-'+file_basename(responsefile)+'.fits'
    FITS_READ,responsefiles[ichip],response,response_head,message=message,/no_abort

    sz=size(outstr.flux,/dim)
    outstr.flux *= response#replicate(1.,sz[1])         ; correct flux
    bderr=where(outstr.err eq baderr(),nbd)
    outstr.err *= response#replicate(1.,sz[1])           ; correct error
    if nbd gt 0 then outstr.err[bderr] = baderr()

    ; Update header
    sxaddhist,leadstr+'Applying response function:',head
    sxaddhist,leadstr+responsefiles[ichip],head
    sxaddpar,head,'RESPFILE',responsefile,' Response file used'
  endif

  ; Adding wavelengths
  ;--------------------
  if keyword_set(wavefile) then begin

    wavefiles = file_dirname(wavefile)+'/'+dirs.prefix+'Wave-'+chiptag+'-'+file_basename(wavefile)+'.fits'
    if not keyword_set(silent) then $
      print,'Adding wavelengths from ',file_dirname(wavefile)+'/'+dirs.prefix+'Wave-'+file_basename(wavefile)

    ; Get the Wavelength calibration data
    FITS_READ,wavefiles[ichip],wcoef,whead,exten=1
    FITS_READ,wavefiles[ichip],wim,whead2,exten=2
    ; This is now fixed in the apWave files
    ;wim = transpose(wim)  ; want it [Npix, Nfibers]

    sxaddhist,leadstr+'Adding wavelengths from',head
    sxaddhist,leadstr+wavefiles[ichip],head
    sxaddpar,head,'WAVEFILE',wavefile,' Wavelength Calibration file'
    sxaddpar,head,'WAVEHDU',5,' Wavelength coef HDU'

  endif

  ; Add header to structure
  head0 = head
  head = strarr(5000)
  nhead = n_elements(head0)
  head[0:nhead-1] = head0
  outstr=CREATE_STRUCT(outstr,'HEADER',head)

  head_chip[ichip,*]=head


  ; Add FIBERS to structure
  if n_elements(fibers) gt 0 then outstr=CREATE_STRUCT(outstr,'FIBERS',fibers)

  ; Output the 2D model spectrum
  if n_elements(ymodel) gt 0 then begin
    modelfile = outdir+dirs.prefix+'2Dmodel-'+chiptag[ichip]+'-'+framenum+'.fits'  ; model output file
    if not keyword_set(silent) then print,'Writing 2D model to: ',modelfile
    MWRFITS,ymodel,modelfile,/create
;    ; compress model and 2D image DONE IN AP2D
;    if keyword_set(compress) then begin
;      file_delete,modelfile+'.fz',/allow_nonexistent
;      SPAWN,'fpack -D -Y '+modelfile
;      origfile = outdir+dirs.prefix+'2D-'+chiptag[ichip]+'-'+framenum+'.fits'
;      if file_test(origfile) then begin
;        file_delete,origfile+'.fz',/allow_nonexistent
;        SPAWN,'fpack -D -Y '+origfile
;      endif
;    endif
  endif

  ; Add to output structure
  if ifirst eq 0 then begin
    output = CREATE_STRUCT('chip'+chiptag[ichip],outstr) 
    if n_elements(ymodel) gt 0 then outmodel=CREATE_STRUCT('chip'+chiptag[ichip],{model:ymodel}) 
    ifirst=1
  endif else begin
    output = CREATE_STRUCT(output,'chip'+chiptag[ichip],outstr) 
    if n_elements(ymodel) gt 0 then outmodel=CREATE_STRUCT(outmodel,'chip'+chiptag[ichip],{model:ymodel}) 
  endelse

  if keyword_set(logfile) then writeline,logfile,file_basename(outfile),string(format='(f8.3)',systime(/seconds)-t1)
  ;stop

  BOMB:

Endfor ; chip loop

; now we have output structure with three chips, each with tags header, flux, err, mask

; Add Wavelength information to the frame structure
;---------------------------------------------------------
; Loop through the chips
if n_elements(wavefiles) gt 0 then begin
 for k=0,2 do begin
  ; Get the Wavelength calibration data
  FITS_READ,wavefiles[k],wcoef,whead,exten=1
  ; Add to the chip structure
  chstr = output.(k)
  chstr = CREATE_STRUCT(temporary(chstr),'FILENAME',files[k],'WAVE_DIR',outdir,'WAVEFILE',wavefiles[k],'WCOEF',wcoef)
  ; Now add this to the final FRAME structure
  if k eq 0 then begin
    frame = CREATE_STRUCT('chip'+chiptag[k],chstr)
  endif else begin
    frame = CREATE_STRUCT(frame,'chip'+chiptag[k],chstr)
  endelse
 endfor
 apgundef,output   ; free up memory
 if file_test(outdir,/directory) eq 0 then file_mkdir,outdir
 plotfile = outdir+'/plots/pixshift-'+framenum
 if keyword_set(skywave) then $
   AP1DWAVECAL,frame,frame_wave,plugmap=plugmap,/verbose,/plot,pfile=plotfile else $
   AP1DWAVECAL,frame,frame_wave,/verbose,/noshift 
 apgundef,frame  ; free up memory
endif else frame_wave = output


; Write output file
;------------------
if not keyword_set(nowrite) then begin

  for i=0,n_elements(chips)-1 do begin
    ichip = chips[i]   ; chip index, 0-first chip
    ; Output file
    outfile = outdir+dirs.prefix+'1D-'+chiptag[ichip]+'-'+framenum+'.fits'  ; output file
    if not keyword_set(silent) then print,'Writing output to: ',outfile
    if keyword_set(outlong) and not keyword_set(silent) then print,'Saving FLUX/ERR as LONG instead of FLOAT'
    ; HDU0 - header only
    FITS_WRITE,outfile,0,reform(head_chip[ichip,*]),/no_abort,message=write_error    
    
    ; HDU1 - flux
    flux = frame_wave.(i).flux
    if keyword_set(outlong) then flux=round(flux)
    MKHDR,head1,flux,/image
    sxaddpar,head1,'CTYPE1','Pixel'
    sxaddpar,head1,'CTYPE2','Fiber'
    sxaddpar,head1,'BUNIT','Flux (ADU)'
    MWRFITS,flux,outfile,head1,/silent

    ; HDU2 - error
    err = errout(frame_wave.(i).err)
    if keyword_set(outlong) then err=round(err)
    MKHDR,head2,err,/image
    sxaddpar,head2,'CTYPE1','Pixel'
    sxaddpar,head2,'CTYPE2','Fiber'
    sxaddpar,head2,'BUNIT','Error (ADU)'
    MWRFITS,err,outfile,head2,/silent

    ; HDU3 - mask
    mask = frame_wave.(i).mask
    MKHDR,head3,mask,/image
    sxaddpar,head3,'CTYPE1','Pixel'
    sxaddpar,head3,'CTYPE2','Fiber'
    if (extract_type eq 1) then begin
      sxaddpar,head3,'BUNIT','Flag Mask (bitwise)'
      sxaddhist,'Explanation of BITWISE flag mask (OR combined)',head3
      sxaddhist,' 1 - bad pixels',head3
      sxaddhist,' 2 - cosmic ray',head3
      sxaddhist,' 4 - saturated',head3
      sxaddhist,' 8 - unfixable',head3
    endif else begin
      sxaddpar,head3,'BUNIT','Flag Mask'
      sxaddhist,'Explanation of flag mask',head3
      sxaddhist,' 0 - good pixels',head3
      sxaddhist,' 1 - bad pixels',head3
    endelse
    MWRFITS,mask,outfile,head3,/silent

    if n_elements(wavefiles) gt 0 then begin
      ; HDU4 - wavelengths
      MKHDR,head4,frame_wave.(i).wavelength,/image
      sxaddpar,head4,'CTYPE1','Pixel'
      sxaddpar,head4,'CTYPE2','Fiber'
      sxaddpar,head4,'BUNIT','Wavelength (Angstroms)'
      MWRFITS,frame_wave.(i).wavelength,outfile,head4,/silent
  
      ; HDU5 = Wavelength solution coefficients [DOUBLE]
      ;-----------------------------------------------------
      wcoef = double(frame_wave.(i).wcoef)
      MKHDR,head5,wcoef,/image
      sxaddpar,head5,'CTYPE1','Fiber'
      sxaddpar,head5,'CTYPE2','Parameters'
      sxaddpar,head5,'BUNIT','Wavelength Coefficients'
      sxaddhist,'Wavelength Coefficients to be used with PIX2WAVE.PRO:',head5
      sxaddhist,' 1 Global additive pixel offset',head5
      sxaddhist,' 4 Sine Parameters',head5
      sxaddhist,' 7 Polynomial parameters (first is a zero-point offset',head5
      sxaddhist,'                     in addition to the pixel offset)',head5
      MWRFITS,wcoef,outfile,head5,/silent
    endif

  endfor
endif

APEND:
file_delete,lockfile

if not keyword_set(silent) then print,'AP2PROC finished'

if keyword_set(stp) then stop

end
