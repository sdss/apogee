pro ap3dquick,filename,outfile,output=output,clobber=clobber,nfowler=nfowler,uptheramp=uptheramp,$
                 saturation=saturation0,detfile=detfile,bpmfile=bpmfile,satfix=satfix,error=error,silent=silent,$
                 outlong=outlong,stp=stp

;+
;
; AP3DQUICK  (note version in apgquicklook vs version in apogeereduce, which is now obsolete)
;
; APOGEE 3D->2D Quick Reduction
;
; INPUTS:
;  filename     Filename for the bundled APOGEE datacube to process
;  outfile      Output filename
;  =nfowler     The number of Fowler samples to use at beginning
;                 and end.
;  /uptheramp   Use Up-the-ramp method to collapse the datacube.
;                 This is the default.
;  /satfix      Fix saturated pixels.  This only works with /uptheramp.
;                 This will increase the runtime by ~30%.
;                 The default is satfix=0.
;  =bpmfile     Bad pixel mask file.
;  =detfile     Detector file with gain and readout noise
;  =saturation  The saturation level.  The default is 65000.
;  /outlong     The output files should use LONG type intead of FLOAT.
;                 This actually takes up the same amount of space, but
;                 this can be losslessly compressed with FPACK.
;  /clobber     Overwrite existing output file.
;  /silent      Don't print anything to the screen
;  /stp         Stop at the end of the program.
;
; OUTPUTS:
;  The collapsed 2D images are saved.
;  =output      The final output data [Nx,Ny,3].
;  =error       The error message if one occurred.
;
; USAGE:
;  IDL>apquickred,'apR-a-00000111.fits',outdir='/net/stream/apogee/data/20110117/'
;
; By D.Nidever  Dec 2010
;-

t0 = systime(1)


; Not enough inputs
nfilename = n_elements(filename)
if nfilename eq 0 then begin
  print,'Syntax - ap3dquick,filename,outfile,output=output,clobber=clobber,nfowler=nfowler,'
  print,'                      uptheramp=uptheramp,saturation=saturation,bpmfile=bpmfile,detfile=detfile,'
  print,'                      outlong=outlong,satfix=satfix,error=error,silent=silent,stp=stp'
  error = 'Not enough inputs'
  return
endif

; No output requested
if n_elements(outfile) eq 0 and not arg_present(output) then begin
  error = 'No output requested'
  if not keyword_set(silent) then print,error
  return
endif

; Test that the file exists
if file_test(filename) eq 0 then begin
  error = filename+' NOT FOUND'
  if not keyword_set(silent) then print,error
  return
endif

; Check that it's FITS
len = strlen(filename)
if strmid(filename,len-5,5) ne '.fits' then begin
  error = filename+' is NOT a FITS file'
  if not keyword_set(silent) then print,error
  return
endif

; Check that the header can be read
head = headfits(filename,exten=0,errmsg=errmsg)
if errmsg ne '' then begin
  error = 'Error reading header for '+filename
  if not keyword_set(silent) then print,error
  return
endif

; Test if the output file already exists
if n_elements(outfile) gt 0 then begin
  test = file_test(outfile)
  if test eq 1 and not keyword_set(clobber) then begin
    error = 'OUTFILE = '+outfile+' ALREADY EXISTS.  Set /clobber to overwrite.'
    if not keyword_set(silent) then print,error
    return
  endif
endif

; Getting Nreads
nreads = sxpar(head,'NREAD')

; Need at least two reads
if nreads lt 2 then begin
  error = 'Nreads = '+strtrim(nreads,2)+'.  NEED AT LEAST 2 READS'
  if not keyword_set(silent) then print,error
  return
endif

if not keyword_set(silent) then $
  print,'Processing ',filename,'  Nreads=',strtrim(nreads,2)

; Default parameters
if n_elements(nfowler) eq 0 and n_elements(uptheramp) eq 0 then uptheramp=1
if n_elements(nfowler) gt 0 then use_nfowler=nfowler < (nreads/2)
if keyword_set(uptheramp) then ap3d_method='uptheramp' else ap3d_method='fowler'
if n_elements(saturation0) eq 0 then saturation = 65000L else saturation=saturation>1  ; must be positive
if n_elements(satfix) eq 0 then satfix=0     ; don't fix saturated pixels by default

; Get dimensions
FITS_READ,filename,im1,head1,message=message,exten_no=1,/no_abort
if message ne '' then begin
  error = 'Error reading extension 1 of '+filename
  if not keyword_set(silent) then print,error
  return
endif
sz = size(im1)
nx = sz[1]
ny = sz[2]
npix = 2048L


; Load the bad pixel mask file
if n_elements(bpmfile) gt 0 then begin
  FITS_READ,bpmfile[0],bpmim,bpmhead,message=bpm_message,/no_abort
  if bpm_message ne '' then begin
    error = 'Error reading BPM file '+bpmfile
    if not keyword_set(silent) then print,error
    return
  endif
  szbpm = size(bpmim)
  if not keyword_set(silent) then print,'Using BPM file = ',bpmfile
endif

; Load the detector file
if n_elements(detfile) gt 0 then begin
  dethead = headfits(detfile)
  FITS_READ,detcorr,rdnoiseim,noisehead,message=message2,/no_abort,exten=1
  FITS_READ,detcorr,gainim,gainhead,message=message1,/no_abort,exten=2
  ;  This should be 2048x2048x3 (each pixel) or 4x3 (each output),
  ;  where the 2 is for a quadratic polynomial
  FITS_READ,detcorr,lindata,linhead,message=message3,/no_abort,exten=3
  if not keyword_set(silent) then print,'Using Detector file = ',detfile
endif else begin
  print,'No detector file found: ',detfile
  dirs=getdir()
  if dirs.instrument eq 'apogee-n' then begin
    gain=1.9
    rdnoise=8
  endif
  if dirs.instrument eq 'apogee-s' then begin
    gain=3.0
    rdnoise=12
  endif
  print,'Using default gain: ',gain
endelse

; Calculate the background/read noise using the reference output
;  this is an underestimate
if nx gt 2048 then noise = MAD(im1[2048:*,*]) > rdnoise else noise=15


; Reference Pixel Correction
;---------------------------
; do something simple and fast
; just subtract median of reference pixels (perimeter)


; Compress the datacube
;------------------------
apgundef,im1,output
CASE ap3d_method of

  ;------------
  ; UP-THE-RAMP
  ;------------
  'uptheramp':  begin

    if not keyword_set(silent) then print,'Using UP-THE-RAMP sampling'
    if keyword_set(satfix) and not keyword_set(silent) then print,'Fixing saturated pixels'

    ; Initializing the cube
    cube = lonarr(nx,ny,nreads)    ; long is big enough and takes up less memory than float
  
    ; Read in the extensions
    for k=1,nreads do begin
      FITS_READ,filename,extim,exthead,exten_no=k,message=message,/no_abort   ; UINT
      cube[*,*,k-1] = extim
    end
    head=headfits(filename)

    ;reference pixel correction
    tmp=APREFCORR(cube,head,mask,readmask=readmask)
    cube=tmp

    bdreads2 = where(readmask eq 1,nbdreads2)
    if nbdreads2 gt 0 then PUSH,bdreads,bdreads2
    nbdreads = n_elements(uniq(bdreads,sort(bdreads)))
  
    if nbdreads gt (nreads-2) then begin
      error = 'Not enough good reads'
      if not keyword_set(silent) then print,'halt: '+error
      stop
      return
    endif

    gdreads = lindgen(nreads)
    REMOVE,bdreads,gdreads
    ngdreads = n_elements(gdreads)
    sz = size(cube)

    sumts = fltarr(sz[1],sz[2])   ; SUM t*s
    sums = fltarr(sz[1],sz[2])    ; SUM s
    sum = intarr(sz[1],sz[2])    ; SUM 
    sumt = fltarr(sz[1],sz[2])   ; SUM t
    sumt2 = fltarr(sz[1],sz[2])   ; SUM t^2
    for k=0L,ngdreads-1 do begin
      slice = cube[*,*,gdreads[k]]
      good = where(finite(slice))   ;might need to include saturation mask here?  See ap3dproc.
      sumts[good] += gdreads[k]*reform(slice[good])
      sums[good] += reform(slice[good])
      sum[good] += 1
      sumt[good] += gdreads[k]
      sumt2[good] += gdreads[k]^2
    end

    ; The slope in Counts per read, similar to med_dCounts_im
    slope_im = (sum*sumts - sumt*sums)/(sum*sumt2 - sumt^2)
    ; To get the total counts just multiply by nread
    im = slope_im * (ngdreads-1L)

    satmask = (sum lt ngdreads)

    ; No saturation fixing, set saturated pixels to the saturation level
    if not keyword_set(satfix) then begin
      im = im*(1-satmask) + satmask*saturation

    ; Saturation fixing, set unfixable pixels to saturation level
    ;  calculate extrapolation error
    endif else begin
      ; set unfixable pixels to saturation level
      unfmask = (sum lt 2)
      im = im*(1-unfmask) + unfmask*saturation

      ; Get approximate error for extrapolation of saturated pixels
      ;  multiply sigma_slope by number of extrapolated reads
      ;  SUM is the number of "good" reads
      sat_extrap_error = ( sigma_slope * (nreads-sum) )*satmask
      ; set to zero for unfixable pixels
      sat_extrap_error *= (1-unfmask)

    endelse

    ; the first read doesn't really add any signal, just a zero-point
  
    ; See Equation 1 in Rauscher et al.(2007), SPIE
    ;  with m=1
    ;  noise and image/flux should be in electrons, sample_noise is in electrons
    sample_noise = sqrt( 12*(ngdreads-1.)/(nreads*(ngdreads+1.))*noise^2 + 6.*(ngdreads^2+1)/(5.*ngdreads*(ngdreads+1))*im*gain )
    sample_noise /= gain      ; convert to ADU
  end ; up-the-ramp


  ;---------
  ; FOWLER
  ;---------
  'fowler': begin


    ; Make sure that Nfowler isn't too large
    Nfowler_used = Nfowler
    if Nfowler gt Nreads/2 then Nfowler_used=Nreads/2

    if not keyword_set(silent) then $
      print,'Using FOWLER sampling. NFOWLER=',strtrim(Nfowler_used,2)

    ; Use the mean of Nfowler reads

    ; Go through the caes
    CASE 1 of

      ; Enough reads to average beg/end
      (nreads gt 3): begin
          ; Beginning sample
          lo1 = 1
          hi1 = nfowler_used
          ; End sample
          lo2 = nreads-nfowler_used+1
          hi2 = nreads
        end

      ; Use the middle read twice
      (nreads eq 3): begin
          Nfowler_used = 2
          ; Beginning sample
          lo1 = 1
          hi1 = 2
          ; End sample
          lo2 = 2
          hi2 = 3
        end

      ; No averaging, simple correlated double sampling
      (nreads eq 2): begin
          ; Beginning sample
          lo1 = 1
          hi1 = 1
          ; End sample
          lo2 = 2
          hi2 = 2
        end

      else: stop

    ENDCASE
 
    ; Get first reads
    totim1 = fltarr(nx,ny)
    numgood1 = lonarr(nx,ny)
    For i=lo1,hi1 do begin
      FITS_READ,filename,im1,head1,exten_no=i,message=message,/no_abort
      if message ne '' then begin
        if not keyword_set(silent) then error = 'Error reading extension '+strtrim(i,2)+' of '+filename
        print,error
        return
      endif
      ;im1 = long(im1[0:npix-1,*])
      im1 = long(im1)
      numgood1 += (im1 lt saturation)
      totim1 += im1
    end
    im_beg = totim1/(hi1-lo1+1)

    ; Get last reads
    totim2 = fltarr(nx,ny)
    numgood2 = lonarr(nx,ny)
    For i=lo2,hi2 do begin
      FITS_READ,filename,im1,head1,exten_no=i,message=message,/no_abort
      if message ne '' then begin
        if not keyword_set(silent) then error = 'Error reading extension '+strtrim(i,2)+' of '+filename
        print,error
        return
      endif
      ;im1 = long(im1[0:npix-1,*])
      im1 = long(im1)
      numgood2 += (im1 lt saturation)
      totim2 += im1
    end
    im_end = totim2/(hi2-lo2+1)


    ; Subtract beginning from end
    im = im_end - im_beg

    ; Make slope image
    mnread1 = mean( lindgen(hi1-lo1+1)+lo1 )  ; mean X of the beg group
    mnread2 = mean( lindgen(hi2-lo2+1)+lo2 )  ; mean X of the end group
    delta_reads = mnread2-mnread1  ; read difference between Fowler averaged images
    slope_im = im/float(delta_reads)

    ; Set saturated pixels to the saturation value
    satmask = (numgood1 lt (hi1-lo1+1)) or (numgood2 lt (hi2-lo2+1))
    im = im*(1-satmask) + satmask*saturation

    ; Correct for the "missing" reads
    ;im *= float(nreads-1)/float(delta_reads)

    ; Noise contribution to the variance
    ;   NEED TO ACCOUNT FOR BAD READS!!!
    sample_noise = noise * sqrt(2.0/Nfowler_used)
    ;sample_noise *= float(nreads-1)/float(delta_reads)  ; Correct for the "missing" reads

  end

  else:  begin
    error = 'This sampling method not supported'
    if not keyword_set(silent) then print,error
    return
  end

ENDCASE


; Trim off the reference output
im = im[0:npix-1,*]
slope_im = slope_im[0:npix-1,*]
satmask = satmask[0:npix-1,*]
if n_elements(unfmask) gt 0 then unfmask=unfmask[0:npix-1,*]

; Mask bad pixels
if n_elements(bpmfile) gt 0 then begin
  ; 0-good pixels, 1-bad pixels
  im = im*(1-bpmim)  ; set bad pixels to zero
endif


;------------------------
; Calculate the Variance
;------------------------
; the total variance is the sum of the rdnoise and poisson noise
; poisson noise = sqrt(Counts) so variance is just Counts
;   need to divide by gain^2 ??
; gain is normally in e/count
; is rdnoise normally in counts or e-??
; do "help apvariance" in IRAF and look for the noise model
;
; variance(ADU) = N(ADU)/Gain + rdnoise(ADU)^2

; Initialize varim
varim = fltarr(npix,npix)         ; variance in ADU

; 1. Poisson Noise
if not keyword_set(uptheramp) then varim = im/gain > 0

; 2. Readnoise
varim += sample_noise^2

; 3. Saturation error
;      We used median(dCounts) to extrapolate the saturated pixels
;      Use the variability in dCounts to estimate the error of doing this
if keyword_set(satfix) then begin
  ; add saturation extrapolation error
  varim += sat_extrap_error     ; zero for non-saturated and unfixable pixels
endif else begin
  varim = varim*(1-satmask) + satmask*1e12   ; saturated pixels are bad!
endelse
; Unfixable pixels
if n_elements(unfmask) gt 0 then $
  varim = varim*(1-unfmask) + unfmask*1e12         ; unfixable pixels are bad!

; Bad pixels
if n_elements(bpmim) gt 0 then begin
  varim = varim*(1-bpmim) + bpmim*1e12               ; bad pixels are bad!
endif


; Construct final mask
;-----------------------
; 1 - bad pixels
; 2 - cosmic rays
; 4 - saturated
; 8 - unfixable
; Saturated pixels
mask = satmask*4L
; Unfixable pixels
if n_elements(unfmask) gt 0 then mask+=unfmask*8L
; Bad pixels
if n_elements(bpmim) gt 0 then begin
  mask += bpmim
endif

;----------------------------
; Construct output datacube
;  [image, variance, mask]
;----------------------------
output = fltarr(npix,npix,3)   ; should this be LONG?
output[*,*,0] = im
output[*,*,1] = sqrt(varim) > 1 ; must be greater than zero
output[*,*,2] = mask


;-----------------------------
; Update header
;-----------------------------
leadstr = 'AP3DQUICK: '
sxaddhist,leadstr+systime(0),head
info = GET_LOGIN_INFO()
sxaddhist,leadstr+info.user_name+' on '+info.machine_name,head
sxaddhist,leadstr+'IDL '+!version.release+' '+!version.os+' '+!version.arch,head
sxaddhist,leadstr+'Output File:',head
sxaddhist,leadstr+' HDU1 - image (ADU)',head
sxaddhist,leadstr+' HDU2 - error (ADU)',head
sxaddhist,leadstr+' HDU3 - flag mask',head
sxaddhist,leadstr+'        1 - bad pixels',head
sxaddhist,leadstr+'        2 - cosmic ray',head
sxaddhist,leadstr+'        4 - saturated',head
sxaddhist,leadstr+'        8 - unfixable',head
maxlen = 72-strlen(leadstr)
; Bad pixel mask file
if n_elements(bpmim) gt 0 then begin
  line = 'BAD PIXEL MASK file="'+bpmfile+'"'
  if strlen(line) gt maxlen then begin
    line1 = strmid(line,0,maxlen)
    line2 = strmid(line,maxlen,100)
    sxaddhist,leadstr+line1,head
    sxaddhist,leadstr+line2,head
  endif else sxaddhist,leadstr+line,head
end
; Bad pixels
if n_elements(bpmim) gt 0 then begin
  totbpm = total(bpmim)
  sxaddhist,leadstr+strtrim(long(totbpm),2)+' pixels are bad',head
endif
; Saturated pixels
;satmask = long((long(mask) AND 4) eq 4)  ; saturated
totsat = total(satmask)
sxaddhist,leadstr+strtrim(long(totsat),2)+' pixels are saturated',head
; Unfixable pixels
if n_elements(unfmask) gt 0 then begin
  totunf = total(unfmask)
  sxaddhist,leadstr+strtrim(long(totunf),2)+' pixels are unfixable',head
endif
; Sampling
if keyword_set(uptheramp) then sxaddhist,leadstr+'UP-THE-RAMP Sampling',head else $
  sxaddhist,leadstr+'Fowler Sampling, Nfowler='+strtrim(long(Nfowler_used),2),head 


;-----------------------
; Output the final image
;-----------------------

; Does the output directory exist?
if file_test(file_dirname(outfile),/directory) eq 0 then begin
  if not keyword_set(silent) then print,'Creating ',file_dirname(outfile)
  FILE_MKDIR,file_dirname(outfile)
endif

; Test if the output file already exists
test = file_test(outfile)

if test eq 1 and keyword_set(clobber) and not keyword_set(silent) then $
  print,'OUTFILE = ',outfile,' ALREADY EXISTS.  OVERWRITING'
if test eq 1 and not keyword_set(clobber) and not keyword_set(silent) then $
  print,'OUTFILE = ',outfile,' ALREADY EXISTS. '
  
  ; Writing file
if test eq 0 or keyword_set(clobber) then begin

  if not keyword_set(silent) then print,'Writing output to: ',outfile
  if keyword_set(outlong) then print,'Saving FLUX/ERR as LONG instead of FLOAT'
  ; HDU0 - header only
  FITS_WRITE,outfile,0,head,/no_abort,message=write_error    
  ; HDU1 - flux
  flux = reform(output[*,*,0])
  if keyword_set(outlong) then flux=round(flux)
  MKHDR,head1,flux,/image
  sxaddpar,head1,'CTYPE1','Pixel'
  sxaddpar,head1,'CTYPE2','Pixel'
  sxaddpar,head1,'BUNIT','Flux (ADU)'
  MWRFITS,flux,outfile,head1,/silent
  ; HDU2 - error
  err = reform(output[*,*,1]) > 1  ; must be greater than zero
  if keyword_set(outlong) then err=round(err)
  MKHDR,head2,err,/image
  sxaddpar,head2,'CTYPE1','Pixel'
  sxaddpar,head2,'CTYPE2','Pixel'
  sxaddpar,head2,'BUNIT','Error (ADU)'
  MWRFITS,err,outfile,head2,/silent
  ; HDU3 - mask
  flagmask = fix(reform(output[*,*,2]))
  MKHDR,head3,flagmask,/image
  sxaddpar,head3,'CTYPE1','Pixel'
  sxaddpar,head3,'CTYPE2','Pixel'
  sxaddpar,head3,'BUNIT','Flag Mask (bitwise)'
  sxaddhist,'Explanation of BITWISE flag mask',head3
  sxaddhist,' 1 - bad pixels',head3
  sxaddhist,' 2 - cosmic ray',head3
  sxaddhist,' 4 - saturated',head3
  sxaddhist,' 8 - unfixable',head3
  MWRFITS,flagmask,outfile,head3,/silent

endif

dt = systime(1)-t0
if not keyword_set(silent) then print,'dt = ',strtrim(dt,2),' sec'

if keyword_set(stp) then stop

end
