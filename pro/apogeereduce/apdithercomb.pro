pro apdithercomb,allframes,shiftstr,pairstr,plugmap,outframe,noscale=noscale,$
                 globalwt=globalwt,nodither=nodither,verbose=verbose,stp=stp,error=error,$
                 newerr=newerr,npad=npad,median=median

;+
;
; APDITHERCOMB
;
; This combines a number of dithered APOGEE frames.
;
; INPUTS:
;  allframes    An array of "frame" structures with the three chip headers and data for
;                 all of the frames.
;  shiftstr     A structure that gives the shifts for all of the frames
;                 relative to the first. A positive shift means that
;                 frame2 is to the RIGHT of frame1 (the reference frame).
;                 This is normally measured by APDITHERSHIFT.PRO
;  plugmap      The Plug Map structure for this plate
;  /noscale     Do NOT scale the two frames.  This would be used for
;                 ThAr or other frames where there are no throughput
;                 variations.
;  /globalwt    When combining the dither-combined spectra use a
;                 single weight for each spectrum (NOT pixel-by-pixel).
;  /nodither    No dithers were performed. Just combine the spectra.
;  /verbose     Print a lot to the screen.
;  /stp         Stop at the end of the program.
;
; OUTPUTS:
;  outframe      A structure that contains the combined images
;                 and headers of the three chips for the dither pair.
;
; USAGE:
;  IDL>apdithercomb,allframes,shiftstr,plugmap,outframe
;
; By D.Nidever  March 2010
; Significant revisions: Jon Holtzman Feb 2012
;-


apgundef,outframe

nallframes = n_elements(allframes)
nshiftstr = n_elements(shiftstr)
nplugmap = n_elements(plugmap)

; Not enough inputs
if nallframes eq 0 or nshiftstr eq 0 or nplugmap eq 0 then begin
  print,'Syntax - apdithercomb,allframes,shiftstr,plugmap,outframe,noscale=noscale,'
  print,'                      globalwt=globalwt,nodither=nodither,verbose=verbose,stp=stp'
  return
endif

chiptag = ['a','b','c']

; Checking the tags of the input structure
for f=0,nallframes-1 do begin
  tags = tag_names(allframes[f])
  needtags1 = ['CHIPA','CHIPB','CHIPC']
  for i=0,n_elements(needtags1)-1 do begin
    if (where(tags eq needtags1[i]))[0] eq -1 then begin
      print,'TAG ',needtags1[i],' NOT FOUND in input structure'
      return
    end
  end
  needtags2 = ['HEADER','FLUX','ERR','MASK']
  for i=0,2 do begin
    tags2 = tag_names(allframes[f].(i))
    for j=0,n_elements(needtags2)-1 do begin
      if (where(tags2 eq needtags2[j]))[0] eq -1 then begin
        print,'TAG ',needtags2[j],' NOT FOUND in input structure'
        return
      end
    end
  end
end ; allframes loop

; Only one exposure
if nallframes eq 1 then begin
  print,'Only one exposure.  Nothing to combine.'
  outframe = allframes
  return
endif

; No dither
if keyword_set(nodither) then begin
  print,'NO Dithers, just combining'
  allcombframes = allframes
  npix = n_elements(allframes[0].(0).flux[*,0])
  nfibers = n_elements(allframes[0].(0).flux[0,*])
  npairs = nallframes  ; Npairs=Nexposures
  y2 = findgen(npix)
  goto,combine
endif

;-----------------------------
; PART I - PAIR UP THE FRAMES
;-----------------------------
APDITHERPAIRS,shiftstr,pairstr,error=error,verbose=verbose,/snsort
npairs = n_elements(pairstr)

if npairs eq 0 and not keyword_set(nodither) then begin
  print,'Error: no dither pairs.'  
  return
endif

if npairs eq 0 then begin
  print,'No dither pairs.  Assuming /nodither and just combining'
  nodither = 1
  allcombframes = allframes
  npix = n_elements(allframes[0].(0).flux[*,0])
  nfibers = n_elements(allframes[0].(0).flux[0,*])
  npairs = nallframes  ; Npairs=Nexposures
  y2 = findgen(npix)
  goto,combine
  return
endif
if n_elements(error) gt 0 then begin
  error = 'There was a problem with the pairing'
  print,error
  return
endif

; Use the first frame of the first pair as the reference frame

refind = pairstr[0].index[0]
refframe = allframes[refind]
sz = size(refframe.chipa.flux)
npix = sz[1]
nfibers = sz[2]
ncol = sz[3]

;; Size not the same
;if total(abs(sz1-sz2)) ne 0 then begin
;  print,'Sizes are not the same'
;  return
;endif



; THINGS TO ADD:
; -NEED TO SINC INTERPOLATE THE DITHER PAIRS ONTO THE SAME FINAL PIXEL
;   ARRAY
;  --> I believe this is being done
; -NEED TO DEAL WITH "MISSING" PIXELS AND "EXTRA" PIXELS AT THE ENDS
; -NEED TO DEAL WITH NOT COUNTING FRAMES MULTIPLE TIMES IN THE
;   VARIANCE ARRAY
;  --> don't use frames multiple times!
; -WHEN COMBINING THE WELL-SAMPLED SPECTRA I NEED TO RESCALE THEM
;   AGAIN.
;  --> I don't think so, they should be weighted by errors
; -PUT ALL THE INFORMATION IN THE HEADER

; I think not counting frames multiple times needs to be done in
; part II just before doing sincinterlace for the variance.
; need to inflate the variance for the spectrum/frame that is being
; multiple times.

; The first spectrum of the first pair has the largest shift to the
; left.  Use this as the reference frame.


;---------------------------------------------
; PART II - SINC INTERLACE THE DITHER PAIRS
;---------------------------------------------
; Loop through the pairs
FOR p=0,npairs-1 do begin

  ipairstr = pairstr[p]
  frame1 = allframes[ipairstr.index[0]]
  frame2 = allframes[ipairstr.index[1]]
  shift = ipairstr.relshift
  ; reference frame is first frame of first pari
  i0=pairstr[0].index[0]
  ; indices of the frames for this pair
  i1=ipairstr.index[0]
  i2=ipairstr.index[1]
  print,'Combining Pair ',strtrim(p+1,2),' - ',ipairstr.framename[0],' + ',ipairstr.framename[1]

  ; Initialize combframe structure. Need 2x as many pixels
  For i=0,2 do begin
    chstr0 = frame1.(i)
    tags = tag_names(chstr0)
    apgundef,chstr

    ; Make the chip structure which will go into the combframe structure
    for j=0,n_elements(tags)-1 do begin
      arr = chstr0.(j)
      type = size(arr,/type)
      ; Data arrays. These are [NPix,Nfiber]
      dum = where(stregex(['FLUX','ERR','MASK','WAVELENGTH','SKY','SKYERR','TELLURIC','TELLURICERR'],tags[j],/boolean) eq 1,ndata)
      if ndata gt 0 then arr=make_array(2*npix,nfibers,type=type)

      ; Add normal tags/data
      if tags[j] ne 'FILENAME' then begin
        if n_elements(chstr) eq 0 then begin
          chstr = CREATE_STRUCT(tags[j],arr)
        endif else begin
          chstr = CREATE_STRUCT(chstr,tags[j],arr)
        endelse

      ; Add FILENAME1 and FILENAME2
      endif else begin
        if n_elements(chstr) eq 0 then begin
          chstr = CREATE_STRUCT('FILENAME1',frame1.(0).filename,'FILENAME2',frame2.(0).filename)
        endif else begin
          chstr = CREATE_STRUCT(chstr,'FILENAME1',frame1.(0).filename,'FILENAME2',frame2.(0).filename)
        endelse
      endelse
    endfor ; tag loop

    ; Add to the final COMBFRAME
    if i eq 0 then begin
      combframe = CREATE_STRUCT('chip'+chiptag[i],chstr)
    endif else begin
      combframe = CREATE_STRUCT(combframe,'chip'+chiptag[i],chstr)
    endelse

  Endfor ; chip loop
  combtags = tag_names(combframe.(0))

  ;----------------------
  ; COMBINE THE FRAMES
  ;----------------------

;  shift_rnd = round(shift*1000.0)/1000.0  ; round to nearest 1/1000th of a pixel
;  shift=shift_rnd

  ; Combine the data with SINC interlacing
  ;------------------------------------------
  y = findgen(npix)

  ;print,'Combining the dither images with SINC interlacing'

  ; Make dummy chip structure with quantities we need to combine
  apgundef,usetags,usetagsnum,f0ch
  tags = tag_names(frame1.(0))
  postags = ['FLUX','ERR','MASK','WAVELENGTH','SKY','SKYERR','TELLURIC','TELLURICERR']  ; use these if they exist
  for k=0,n_elements(tags)-1 do begin
    gd = where(postags eq tags[k],ngd)
    if ngd gt 0 then begin
      PUSH,usetags,tags[k]
      PUSH,usetagsnum,k
      if n_elements(f0ch) eq 0 then begin
        f0ch = CREATE_STRUCT(tags[k],(frame1.(0).(k))[*,0])
      endif else begin
        f0ch = CREATE_STRUCT(f0ch,tags[k],(frame1.(0).(k))[*,0])
      endelse
    end
  end
  PUSH,usetags,'SCALE'
  f0ch = CREATE_STRUCT(f0ch,'SCALE',fltarr(npix))
  nusetags = n_elements(usetags)

  ; we will "pad" the spectra at either end with masked pixels
  ;   to determine which of the edge pixels should be declared
  ;   bad
  if not keyword_set(npad) then npad=0
  if npad gt 0 then begin
    lftpad = fltarr(npad)
    rtpad = fltarr(npad)
    maskpad=fltarr(npad)
    maskpad[*]=1
  endif

  ; Loop through the fibers
  ;-------------------------
  for j=0,nfibers-1 do begin

    if (j+1) mod 50 eq 0 then print,strtrim(j+1,2),'/',strtrim(nfibers,2)

    ; The plugmap index for this fiber
    ; fiberid=1 is at the top of the detector or index=299
    ; index = 300-fiberid
    iplugind = where(plugmap.fiberdata.spectrographid eq 2 and $
                     plugmap.fiberdata.fiberid eq 300-j,niplugind)
    ; No information for this fiber
    if niplugind eq 0 then begin
       print,'No information for Fiber=',strtrim(300-j,2),' in the plugmap file'
       goto,BOMB
    endif

    ; Fiber type from the plugmap structure
    fiberobjtype = plugmap.fiberdata[iplugind].objtype
    fiberobjid = plugmap.fiberdata[iplugind].tmass_style

    ; Loop through the chips
    for ichip=0,2 do begin

      ; Insert the quantities that we want to combine
      data1 = f0ch
      for k=0,nusetags-2 do data1.(k) = (frame1.(ichip).(usetagsnum[k]))[*,j]
      data2 = f0ch
      for k=0,nusetags-2 do data2.(k) = (frame2.(ichip).(usetagsnum[k]))[*,j]

      ; replace bad pixels and errors with smoothed version to minimize their impact on adjacent pixels
      ; we will extend masks if these pixels are important
      if keyword_set(median) then begin
        bd=where(data1.mask AND badmask(),nbd)
        if nbd gt 0 then data1.flux[bd]=!values.f_nan
	scale1=smooth(medfilt1d(data1.flux,501,edge=2),100,/nan)
	tmperr=smooth(medfilt1d(data1.err,501,edge=2),100,/nan)
        if nbd gt 0 then begin
          data1.flux[bd]=scale1[bd]
          data1.err[bd]=tmperr[bd]
        endif
        bd=where(data2.mask and badmask(),nbd)
        if nbd gt 0 then data2.flux[bd]=!values.f_nan
	scale2=smooth(medfilt1d(data2.flux,501,edge=2),100,/nan)
	temperr=smooth(medfilt1d(data2.err,501,edge=2),100,/nan)
        if nbd gt 0 then begin
          data2.flux[bd]=scale2[bd]
          data2.err[bd]=tmperr[bd]
        endif
      endif else begin
        ; Fit a low-order polynomial to the data
        scalecoef1 = ROBUST_POLY_FIT(y,data1.flux,5)
        scale1 = POLY(y,scalecoef1)
        scalecoef2 = ROBUST_POLY_FIT(y,data2.flux,5)
        scale2 = POLY(y,scalecoef2)
      endelse

      ; Need to normalize them if they are object spectra
      ; Don't do this for sky spectra, since they may have real variations
      ;   of course, this means you can't really dither-combine the sky spectra, either...
      ;   but we will just to put something in the dither-combined output
      ; Use the plugmap information to see which fibers are sky
      ; Only normalize the flux and errors, not sky or telluric!
      if (fiberobjtype ne 'SKY') and not keyword_set(noscale) then begin

        ; "Normalize" the spectra in case there's been some variation in response
        data1.flux = data1.flux/scale1        ; normalize spectrum
        data1.err = data1.err/scale1          ; normalize error
        data1.scale = scale1

        data2.flux = data2.flux/scale2        ; normalize spectrum
        data2.err = data2.err/scale2          ; normalize error
        data2.scale = scale2
      endif

      ; We want to interpolate onto a 1/2 pixel scale.  So we only
      ; really need to interpolate half of the pixels, the other half
      ; just stay the way they are.
      ;-------------------
      ; FRAME TWO is FIRST
      ;-------------------
      ; Frame2 is a fraction of a pixel to the RIGHT of Frame1
      ;  and so for Frame2 the detector moved to the LEFT
      ;  and Frame2 should be interleaved FIRST, i.e.
      ;  [Frame2, Frame1]

      ;  Now the planes are: [spec, wave, error, flag, sky, errsky,
      ;    telluric, error_telluric]

      ; What is the relative shift between this pair and the ABSOLUTE
      ; frame?
      abs_shift = ipairstr.refshift
      ;print,ichip,j,shift,abs_shift

      ; Dec 2018: use chip and fiber dependent shifts
      ; shifts are all measured in the reverse direction from the maximum shift (which is put in first pair, see apditherpairs)
      ;new_abs_shift = allframes[i0].shift.shiftfit[0] - allframes[i1].shift.shiftfit[0]
      ;new_shift = allframes[i0].shift.shiftfit[0] - allframes[i2].shift.shiftfit[0] - new_abs_shift
      ;shift_i0 = allframes[i0].shift.shiftfit[0]
      ;shift_i1 = allframes[i1].shift.shiftfit[0]
      ;shift_i2 = allframes[i2].shift.shiftfit[0]
      shift_i0 = allframes[i0].shift.chipfit[0]*j+allframes[i0].shift.chipfit[ichip+1] 
      shift_i1 = allframes[i1].shift.chipfit[0]*j+allframes[i1].shift.chipfit[ichip+1] 
      shift_i2 = allframes[i2].shift.chipfit[0]*j+allframes[i2].shift.chipfit[ichip+1] 
      new_abs_shift = shift_i0 - shift_i1
      new_shift = shift_i1 - shift_i2
      ;print,new_shift,new_abs_shift
      shift = new_shift
      abs_shift = new_abs_shift

      ; if abs_shift is positive then we are to the RIGHT of the
      ; absolute frame and we want pixels to the LEFT

      ; If the shift is large, then SINCINTERLACED.PRO will return
      ; garbage for the "overhanging/missing" pixels.  We just need
      ; modify these to show they are bad.
      ; This is done at the end.

      ; Combine the SPECTRA
      ;---------------------
      ; pad the data and error arrays
      if npad eq 0 then begin
        spec1=data1.flux
        spec2=data2.flux
        err1=data1.err
        err2=data2.err
      endif else begin
        lftpad[*] = data1.flux[0]
        rtpad[*] = data1.flux[npix-1]
        spec1 = [lftpad,data1.flux,rtpad]   ; the leftmost (in wavelength)
        lftpad[*] = data1.err[0]
        rtpad[*] = data1.err[npix-1]
        err1 = [lftpad,data1.err,rtpad]
        lftpad[*] = data2.flux[0]
        rtpad[*] = data2.flux[npix-1]
        spec2 = [lftpad,data2.flux,rtpad]   ; starting to the right
        lftpad[*] = data2.err[0]
        rtpad[*] = data2.err[npix-1]
        err2 = [lftpad,data2.err,rtpad]
      endelse

      ; Do the interlaced sinc interpolation
      ;;  If shift = 0.5 then it will output the appropriate spectrum
      ;;  and NOT do the interpolation

      ; Need to get the interpolation on the right absolute scale
      ;  spec1 is the leftmost
      ;  spec2 starts at a higher value, but goes first
      spec_lfthalf = SINCINTERLACED(spec2,spec1,shift, 0.0-abs_shift,$
                        err1=err2,err2=err1,errout=err_lfthalf)
      spec_rthalf = SINCINTERLACED(spec2,spec1,shift, 0.5-abs_shift,$
                        err1=err2,err2=err1,errout=err_rthalf)

      ; Now combine the spectra, without the padded pixels
      combspec = dblarr(npix*2)
      combspec[0:npix*2-2:2] = spec_lfthalf[npad:npad+npix-1]
      combspec[1:npix*2-1:2] = spec_rthalf[npad:npad+npix-1]
      comberr = dblarr(npix*2)
      comberr[0:npix*2-2:2] = sqrt(err_lfthalf[npad:npad+npix-1])
      comberr[1:npix*2-1:2] = sqrt(err_rthalf[npad:npad+npix-1])
      ; RENORMALIZE THE DATA back up
      ;  use the average of scale1 + scale2
      if fiberobjtype ne 'SKY' and not keyword_set(noscale) then begin
        ;y2 = findgen(npix*2)/2
        ;; Using continuum poly fit for Frame1
        ;;if data1[npix/2,3] gt data2[npix/2,3] then begin
        ;if scale1[npix/2] gt scale2[npix/2] then begin
        ;  rescale = POLY(y2,scalecoef1) / scale1[npix/2] * MEAN([ scale1[npix/2], scale2[npix/2] ])
        ;; Using continuum poly fit for Frame2
        ;endif else begin
        ;  rescale = POLY(y2,scalecoef2) / scale2[npix/2] * MEAN([ scale1[npix/2], scale2[npix/2] ])
        ;endelse
        samp=indgen(npix*2)/2.
        rescale=interpolate((scale1+scale2)/2.,samp)
        combspec = combspec * rescale                       ; rescale spectrum
        comberr = comberr * rescale                         ; rescale error
      endif

      ; old way of doing the error is to interlace the scaled error arrays,
      ;  but since there may be a mismatch in S/N, this doesn't seem
      ;  like a good idea
      if not keyword_set(newerr) then begin
        ; Combine the ERROR
        ;----------------------
        err1 = data1.err
        err2 = data2.err
        ; Do the interlaced sinc interpolation
        ;;  If shift = 0.5 then it will output the appropriate fiber
        ;;  and NOT do the interpolation
        ; Need to get the interpolation on the right absolute scale
        ;  err1 is the leftmost
        ;  err2 starts at a higher value, but goes first
        err_lfthalf = SINCINTERLACED(err2,err1,shift, 0.0-abs_shift)
        err_rthalf = SINCINTERLACED(err2,err1,shift, 0.5-abs_shift)

        ; Now combine the error
        comberr = dblarr(npix*2)
        comberr[0:npix*2-2:2] = err_lfthalf
        comberr[1:npix*2-1:2] = err_rthalf
        if fiberobjtype ne 'SKY' and not keyword_set(noscale) then $
          comberr = comberr * rescale                         ; rescale error
        ; Make sure we don't get negative or zero values
        bderr = where(comberr le 0.0,nbderr)
        if nbderr gt 0 then begin
          gderr = where(comberr gt 0.0,ngderr)
          comberr[bderr] = min(comberr[gderr])  ; use the minimum "good" error
          ; should these pixels be considered "bad" and put
          ; to a high value???
        endif
      endif

      ; put results in output structure
      combframe.(ichip).flux[*,j] = combspec
      combframe.(ichip).err[*,j] = comberr

      ; get the contribution of masked pixels to the output. Do this separately
      ;   for each bit to keep track of the flags
      combframe.(ichip).mask[*,j]=0
      getmaskvals,flag,badflag,maskcontrib
      for ibit=0,n_elements(flag)-1 do begin
        vmask=2^ibit
        mask1=(data1.mask and vmask)<1
        mask2=(data2.mask and vmask)<1
        junk=where(mask1 gt 0,nmask1)
        junk=where(mask2 gt 0,nmask2)
        ; only continue if we have any of these bits set!
        if nmask1 gt 0 or nmask2 gt 0 then begin
          if npad eq 0 then begin
            mask1=float(mask1)
            mask2=float(mask2)
          endif else begin
            mask1 = [maskpad,float(mask1),maskpad]   ; the leftmost (in wavelength)
            mask2 = [maskpad,float(mask2),maskpad]   ; starting to the right
          endelse
          spec_lfthalf = SINCINTERLACED(mask2,mask1,shift, 0.0-abs_shift)
          spec_rthalf = SINCINTERLACED(mask2,mask1,shift, 0.5-abs_shift)
          combmask = dblarr(npix*2)
          combmask[0:npix*2-2:2] = spec_lfthalf[npad:npad+npix-1]
          combmask[1:npix*2-1:2] = spec_rthalf[npad:npad+npix-1]
          ; any pixel that has more than allowed contribution from a bad pixel gets this mask set
          bd=where(abs(combmask) gt maskcontrib[ibit],nbd)
          if nbd gt 0 then combframe.(ichip).mask[bd,j]=combframe.(ichip).mask[bd,j] OR vmask
        endif
      endfor
      ; if this maskval corresponds to a bad pixel, inflate the error
      bd = where(combframe.(ichip).mask[*,j] and badmask(),nbd)
      if nbd gt 0 then begin
        combframe.(ichip).err[bd,j] *= 10.
      ;  combframe.(ichip).flux[bd,j]=!values.f_nan
      endif

      ; Flag "bad/missing" pixels at the ends
      ;--------------------------------------
      if npad eq 0 and abs(abs_shift) gt 0.6 then begin
        ; how many pixels are bad,  we can interpolate on pixel at ends
        ;  the shift is in original pixels
        ;  we are flagging dither combined pixels
        nbadpix = floor(abs(abs_shift*2)) - 1

        ; pixels at beginning are bad
        if abs_shift gt 0.0 then begin
          ; just set the flag value to NAN
          ;  this might already happen above when the non-interpolated
          ;  values are shifted
          ;combframe.(ichip).data[0:nbadpix-1,j,*] = !values.f_nan
          ;combframe.(ichip).mask[0:nbadpix-1,j] = !values.f_nan
          combframe.(ichip).mask[0:nbadpix-1,j] = maskval('BADPIX')

        ; pixels at the end are bad
        endif else begin
          ; just set the flag value to NAN
          ;combframe.(ichip).data[npix*2-nbadpix:npix*2-1,j,*] = !values.f_nan
          ;combframe.(ichip).mask[npix*2-nbadpix:npix*2-1,j] = !values.f_nan
          combframe.(ichip).mask[npix*2-nbadpix:npix*2-1,j] = maskval('BADPIX')
        endelse

      end

      ; Make the combined wavelength array
      ;-------------------------------------
      ;  Just use the wavelength coefficients of Frame 1
      ;    since it is on the left (lowest wavelength)
      dum = where(usetags eq 'WAVELENGTH',Nwavelength)
      if Nwavelength gt 0 then begin
        w1 = data1.wavelength
        w2 = data2.wavelength
        wcoef1 = reform(frame1.(ichip).wcoef[j,*])

        ; We need wavelengths on the absolute scale
        ;  same as for the sinc interpolation
        wave_lfthalf = PIX2WAVE(y+0.0-abs_shift,wcoef1)
        wave_rthalf = PIX2WAVE(y+0.5-abs_shift,wcoef1)

        combwave = dblarr(npix*2)
        combwave[0:npix*2-2:2] = wave_lfthalf
        combwave[1:npix*2-1:2] = wave_rthalf
        ; Copy the Wavelength coefficients to the output frame
        ;  modify wcoef2 for the absolute shift
        wcoef1[0] = wcoef1[0]-abs_shift
        combframe.(ichip).wcoef[j,*] = wcoef1

        ; Stuff it in the output structure
        combframe.(ichip).wavelength[*,j] = combwave
      endif

      ; Combine flags, sky, errsky, telluric, error_telluric
      ;-------------------------------------------------------
      ;   the sky and telluric come from different dithered
      ;   exposures and will be on different levels, so we can't sinc interlace them
      ;   just do a simple interlace (really, we should shift these!)
      extratags = ['SKY','SKYERR','TELLURIC','TELLURICERR']  ; all possible extras
      nextra = n_elements(extratags)

      ; Loop through the extras
      for k=0,nextra-1 do begin

        xtraind = where(usetags eq extratags[k],nxtraind)
        combind = where(combtags eq extratags[k],ncombind)

        ; We have this extra
        if nxtraind gt 0 then begin

          ; Combine the data
          combextra = dblarr(npix*2)
          combextra[0:npix*2-2:2] = data2.(xtraind)
          combextra[1:npix*2-1:2] = data1.(xtraind)

          ; SHIFT NON-INTERPOLATED VALUES
          ;  if abs_shift is greater than 1/2 dither combined pixel
          ;  then we need to shift these values
          if abs(abs_shift*2) gt 0.5 then begin
            combextra0 = combextra
            xtra_shift = ceil(abs(abs_shift*2))

            ; if abs_shift is positive then extra pixels will be "added"
            ; on the left side and our arrays need to shifted to the right

            ; shift to the right
            if abs_shift gt 0 then begin
              combextra = shift(combextra,xtra_shift)
              combextra[0:xtra_shift-1] = !values.f_nan  ; these values are garbage

            ; shift to the left
            endif else begin
               combextra = shift(combextra,-xtra_shift)
               combextra[npix*2-xtra_shift:npix*2-1] = !values.f_nan   ; these values are garbage
            endelse

          endif ; shifting

          ; Stuff in the output structure
          combframe.(ichip).(combind)[*,j] = combextra

        endif ; we have this extra


      endfor    ; extras loop

      ;stop


    End ; chip loop

    ; Plotting
    ;----------
    pl = 0 ;1
    if keyword_set(pl) then begin

      xx = [findgen(npix*2), findgen(npix*2)+npix*2+300, findgen(npix*2)+npix*4+2*300]
      yy = [combframe.(0).flux[*,j], combframe.(1).flux[*,j], combframe.(2).flux[*,j] ]

      ;xr = [480,490]
      ;yr = minmax(combspec[xr[0]*2:xr[1]*2])
      ;xr = [0,npix*2]
      ;yr = [min(combspec),max(combspec)]
      xr = minmax(xx)
      yr = minmax(yy)
 
      plot,[0],[0],/nodata,xr=xr,yr=yr,xs=1,ys=1,xtit='Pixel',ytit='Counts',$
           tit=strtrim(fiberobjid,2)+' (OBJTYPE='+fiberobjtype+')'
      for k=0,2 do oplot,findgen(npix*2)+npix*2*k+k*300,combframe.(k).flux[*,j]
      ;oplot,combspec
      ;oplot,sqrt(combvar),co=250
      ;oplot,xx,yy

      wait,1
      ;stop

    end

    ;stop

    BOMB:

  End ; fiber loop

  ;-------------------------------------------------------------------------
  ; Convert the WAVELENGTH solution coefficients to dither combined pixels
  ;-------------------------------------------------------------------------
  ; The equation is:
  ; wave = P[1]*( SIN( (Y+P[0]+P[2])/P[3]/radeg ) + P[4]) + POLY(Y+P[0]+P[5],P[6:*])

  ; Since we changing Y->2*Y we need to do:
  ; P[0] -> P[0]*2
  ; P[2] -> P[1]*2
  ; P[3] -> P[3]*2
  ; P[5] -> P[5]*2
  ; P[6] -> P[6]  constant term
  ; P[7] -> P[7]/2
  ; P[8] -> P[8]/2^2
  ; P[9] -> P[9]/2^3  and so on
  dum = where(combtags eq 'WCOEF',nwcoef)
  if nwcoef gt 0 then begin

    ; Loop through the chips
    for i=0,2 do begin
      wcoef_old = combframe.(i).wcoef
      wcoef_new = wcoef_old
      nwcoef = n_elements(wcoef_new[0,*])
      npolycoef = nwcoef-6
      wcoef_new[*,0] *= 2
      wcoef_new[*,2] *= 2
      wcoef_new[*,3] *= 2
      wcoef_new[*,5] *= 2
      for j=0,npolycoef-1 do wcoef_new[*,j+6] /= 2^j

      combframe.(i).wcoef = wcoef_new
    endfor

  endif ; WCOEF is there

  ;-------------------------------------------------------
  ; Convert the LSF parameters to dither combined pixels
  ;-------------------------------------------------------

  ; Since we're changing Y->2*Y we need to do:
  ;  polynomial coefficients change depending on the power
  ; P[0] -> P[0]      constant term
  ; P[1] -> P[1]/2    linear term
  ; P[2] -> P[2]/2^2  quadratic term
  ; P[3] -> P[3]/2^3  and so on

  ; Only convert if the LSF parameters were created from non-dithered
  ;  exposures and have binsize=1
  dum = where(combtags eq 'LSFCOEF',nlsfcoef)
  if nlsfcoef gt 0 then if combframe.(0).lsfcoef[0,0] eq 1 then begin

    ; Loop through the chips
    for i=0,2 do begin
      ; Loop through the fibers
      for j=0,nfibers-1 do begin

        lsfpar = reform(combframe.(i).lsfcoef[j,*])
        npar = n_elements(lsfpar)

        ; Breaking up the parameters
        binsize = lsfpar[0]
        Xoffset = lsfpar[1]   ; Additive Xoffset
        Horder = lsfpar[2]
        Porder = lsfpar[3:Horder+3]   ; Horder+1 array
        nGHcoefs = total(Porder+1)

        ; Getting the GH parameters that vary globally
        cpar = lsfpar[Horder+4:Horder+4+nGHcoefs-1]
        coefarr = dblarr(Horder+1,max(Porder)+1)
        cstart = [0,TOTAL(Porder+1,/cum)]  ; extra one at the end
        ; Coefarr might have extra zeros at the end, but it shouldn't
        ;  make a difference.
        for k=0,Horder do $
          coefarr[k,0:Porder[k]] = cpar[cstart[k]:cstart[k]+Porder[k]]

        lsfpar_new = lsfpar ; initialize the new lsf parameter array
        lsfpar_new[0] *= 2  ; update binsize
        lsfpar_new[1] *= 2  ; update Xoffset

        ; Update polynomial coefficients for factor of 2x change
        cpar_new = cpar
        coefarr_new = coefarr
        for k=1,max(Porder) do coefarr_new[*,k] /= 2^k  ; correct for 2 starting w linear term

        ; Now we need to update the GH parameters themselves to
        ;  account for the change in X.
        ;  Need to multiply all polynomial coefficients by the
        ;   appropriate factor
        ; Just need to scale SIGMA
        coefarr_new[0,*] *= 2  ; sigma -> sigma*2

        ; Stuff back in
        for k=0,Horder do $
          cpar_new[cstart[k]:cstart[k]+Porder[k]] = coefarr_new[k,0:Porder[k]]
        lsfpar_new[Horder+4:Horder+4+nGHcoefs-1] = cpar_new  ; stuff it back in


        ; Wing parameters
        if npar gt (3+Horder+1+nGHcoefs) then begin
          wpar = lsfpar[3+Horder+1+nGHcoefs:*]

          ; Nwpar     number of W parameters
          ; WPorder   the polynomial order for each
          ; Wing coefficients
          wproftype = wpar[0]
          nWpar = wpar[1]
          wPorder = wpar[2:2+nWpar-1]
          nWcoefs = total(wPorder+1)

          ; Getting the Wing parameters that vary globally
          wcoef = wpar[nWpar+2:*]
          wcoefarr = dblarr(nWpar,max(wPorder)+1)
          wcstart = [0,TOTAL(wPorder+1,/cum)]  ; extra one at the end
          ; wcoefarr might have extra zeros at the end, but it shouldn't
          ;  make a difference.
          for k=0,nWpar-1 do $
            wcoefarr[k,0:wPorder[k]] = wcoef[wcstart[k]:wcstart[k]+wPorder[k]]

          ; Wing input parameters
          ; 1st par - area under curve
          ; 2nd par - center

          ; Update polynomial coefficients
          wcoef_new = wcoef
          wcoefarr_new = wcoefarr
          for k=1,max(wPorder) do wcoefarr_new[*,k] /= 2^k  ; correct for 2 starting w linear term

          ; Now we need to update the GH parameters themselves to
          ;  account for the change in X.
          ;  Need to multiply all polynomial coefficients by the
          ;   appropriate factor
          ; Just need to scale SIGMA
          wcoefarr_new[1,*] *= 2  ; sigma -> sigma*2

          ; Stuff it back in
          for k=0,nWpar-1 do $
            wcoef_new[wcstart[k]:wcstart[k]+wPorder[k]] = wcoefarr_new[k,0:wPorder[k]]
          wpar_new = wpar
          wpar_new[nWpar+2:*] = wcoef_new
          lsfpar_new[3+Horder+1+nGHcoefs:*] = wpar_new  ; stuff it back in

        endif ; wing parameters

        ; Stuff the new LSF parameters into the combframe structure
        combframe.(i).lsfcoef[j,*] = lsfpar_new

      endfor ; fiber loop
    endfor  ; chip loop

  endif ; LSFCOEF is there

  ; Update the headers
  ;----------------------
  leadstr = 'APDITHERCOMB: '
  maxlen = 72-strlen(leadstr)
  line = 'Dither combined frames '+combframe.chipa.filename1+' and '+combframe.chipa.filename2
  ncuts = ceil(strlen(line)/float(maxlen))
  for l=0,ncuts-1 do apaddpar,combframe,leadstr+strmid(line,l*maxlen,maxlen),/history


  ; Add to structure of all dither combined frames
  ;------------------------------------------------
  if p eq 0 then allcombframes=combframe else allcombframes=[allcombframes,combframe]

ENDFOR  ; pair loop

;stop

;---------------------------------------------------------------
; PART III - COMBINE all fully-sampled (dither combined) frames
;---------------------------------------------------------------
COMBINE:

npix2 = n_elements(allcombframes[0].(0).flux[*,0])

; Initialize OUTFRAME structure. Need 2x as many pixels
For i=0,2 do begin
  ;chstr0 = allframes[0].(i)
  chstr0 = allcombframes[0].(i)
  tags = tag_names(chstr0)
  apgundef,chstr

  ; Make the new chip structure
  for j=0,n_elements(tags)-1 do begin
    arr = chstr0.(j)
    sz = size(arr)
    type = size(arr,/type)
    ; Data arrays. These are [NPix,Nfibers]
    dum = where(stregex(['FLUX','ERR','MASK','WAVELENGTH','SKY','SKYERR','TELLURIC','TELLURICERR'],tags[j],/boolean) eq 1,ndata)
    if ndata gt 0 then arr=make_array(npix2,nfibers,type=type)
    ; I think this is redundant now that we are starting with ALLCOMBFRAME
    ; instead of ALLFRAMES

    ; Skip FILENAMEs
    if strmid(tags[j],0,8) ne 'FILENAME' then begin
      if n_elements(chstr) eq 0 then begin
        chstr = CREATE_STRUCT(tags[j],arr)
      endif else begin
        chstr = CREATE_STRUCT(chstr,tags[j],arr)
      endelse
    endif
  endfor ; tag loop

  ; Reset to start with blank header
  mkhdr,header,0
  chstr.header=header

  ; Add to the final OUTFRAME
  if i eq 0 then begin
    outframe = CREATE_STRUCT('chip'+chiptag[i],chstr)
  endif else begin
    outframe = CREATE_STRUCT(outframe,'chip'+chiptag[i],chstr)
  endelse

Endfor ; chip loop

; Put information in the header
leadstr = 'APDITHERCOMB: '
maxlen = 72-strlen(leadstr)
apaddpar,outframe,leadstr+'Combining '+strtrim(npairs,2)+' dither pairs',/history
if keyword_set(globalwt) then wtmethod = 'Using global spectrum weighting' else $ ; weight method
  wtmethod = 'Using pixel-by-pixel weighting'
apaddpar,outframe,leadstr+wtmethod,/history
if not keyword_set(nodither) then begin

 ; Add dither information to headers
  for i=0,npairs-1 do begin
    apaddpar,outframe,leadstr+'Pair '+strtrim(i+1,2),/history
    apaddpar,outframe,leadstr+file_basename(allcombframes[i].(0).filename1),/history
    apaddpar,outframe,leadstr+'Shift='+strtrim(pairstr[i].shift[0],2)+' at row 0',/history
    apaddpar,outframe,leadstr+file_basename(allcombframes[i].(0).filename2),/history
    apaddpar,outframe,leadstr+'Shift='+strtrim(pairstr[i].shift[1],2)+' at row 0',/history
  end
  ; Add frame names to header
  apgundef,framenames
  apaddpar,outframe,'NPAIRS',strtrim(npairs,2)
  apaddpar,outframe,'NCOMBINE',strtrim(npairs,2)*2
  for i=0,npairs-1 do PUSH,framenames,pairstr[i].framename
  nframes = n_elements(framenames)
  ;for i=0,nframes-1 do apaddpar,outframe,'FRAME'+strtrim(i+1,2),framenames[i]
  ii=0
  for i=0,npairs-1 do begin
    apaddpar,outframe,'SHIFT'+strtrim(ii+1,2),pairstr[i].shift[0]
    apaddpar,outframe,'FRAME'+strtrim(ii+1,2),pairstr[i].framenum[0]
    ii+=1
    apaddpar,outframe,'SHIFT'+strtrim(ii+1,2),pairstr[i].shift[1]
    apaddpar,outframe,'FRAME'+strtrim(ii+1,2),pairstr[i].framenum[1]
    ii+=1
 endfor
endif else begin
  ; No dither combining
  ; Add frame names to header
  apgundef,framenames
  for i=0,npairs-1 do PUSH,framenames,shiftstr[i].framenum
  nframes = n_elements(framenames)
  apaddpar,outframe,'NCOMBINE',nframes
  for i=0,nframes-1 do apaddpar,outframe,'FRAME'+strtrim(i+1,2),framenames[i]

endelse

; Update the exposure time
;  We sum the exposures so we should sum the exposure times
totexptime = 0.0
for i=0,npairs-1 do totexptime+=sxpar(allcombframes[i].(0).header,'EXPTIME')
apaddpar,outframe,'EXPTIME',totexptime,' Total visit exposure time per dither pos'

; Only one pair
combtags = tag_names(allcombframes[0].(0))
outtags = tag_names(outframe.(0))
if npairs eq 1 then begin
  ; Loop through the tags
  for i=0,n_elements(outtags)-1 do begin
    ind = where(combtags eq outtags[i],nind)
    for k=0,2 do outframe.(k).(i) = combframe.(k).(ind)
  endfor
  return
endif

; Add UT-MID, JD-MID and EXPTIME to the header
;---------------------------------------------
; get JD-MID and EXPTIME for all exposures
jdmid = dblarr(nallframes)
exptime = fltarr(nallframes)
for i=0,nallframes-1 do begin
  jdmid1 = sxpar(allframes[i].(0).header,'JD-MID',count=njdmid)
  exptime1 = sxpar(allframes[i].(0).header,'EXPTIME')
  if njdmid eq 0 then begin
    dateobs = sxpar(allframes[i].(0).header,'DATE-OBS')
    jd = date2jd(dateobs)
    jdmid1 = jd + (0.5*exptime1)/24./3600.d0
  endif
  jdmid[i] = jdmid1-2400000.
  exptime[i] = exptime1
endfor
; calculate exptime-weighted mean JD-MID
comb_jdmid = total(exptime*jdmid)/total(double(exptime))+2400000.
comb_utmid = jd2date(comb_jdmid)
apaddpar,outframe,'UT-MID',comb_utmid,' Date at midpoint of visit'
apaddpar,outframe,'JD-MID',comb_jdmid,' JD at midpoint of visit'

; Use DATE-OBS of the first exposure for the combined frame
minind = first_el(minloc(jdmid))
apaddpar,outframe,'DATE-OBS',sxpar(allframes[minind].(0).header,'DATE-OBS')


print,'Combining ',strtrim(npairs,2),' dither pairs'

; Loop through the fibers
;------------------------
For j=0,nfibers-1 do begin

  if (j+1) mod 50 eq 0 then print,strtrim(j+1,2),'/',strtrim(nfibers,2)

  ; Get object type
  ; The plugmap index for this fiber
  iplugind = where(plugmap.fiberdata.spectrographid eq 2 and $
                   plugmap.fiberdata.holetype eq 'OBJECT' and $
                   plugmap.fiberdata.fiberid eq 300-j,niplugind)
  ; Getting information on the object
  if niplugind gt 0 then begin
    ; Fiber type from the plugmap structure
    fiberobjtype = plugmap.fiberdata[iplugind].objtype
    fiberobjid = plugmap.fiberdata[iplugind].tmass_style
  endif else begin
    ; No information for this fiber
     print,'No information for FiberID=',strtrim(300-j,2),' in the plugmap file'
     fiberobjtype = ''
     fiberobjid = -1
  endelse

  ; Measure the average wavelength zeropoint for all pairs
  ;--------------------------------------------------------
  ;  The accuracy of each wavelength solution is about the same,
  ;  and so it's probably best to give them all the same weight
  wpix0 = allcombframes.(0).wcoef[j,0]
  wpix0 = wpix0-wpix0[0]   ; relative to the first one
  wpix0_offset = MEAN(wpix0)

  ; Loop through the chips
  ;-----------------------
  For ichip=0,2 do begin

    ; Initialize the "data" structure for all spectra
    dumstr = {flux:fltarr(npix2),err:fltarr(npix2),mask:intarr(npix2),$
             wavelength:dblarr(npix2),sky:fltarr(npix2),skyerr:fltarr(npix2),$
             telluric:fltarr(npix2),telluricerr:fltarr(npix2),scale:fltarr(npix2)}
    data = REPLICATE(dumstr,npairs)

    ; Loop through the spectra/pairs
    ;-------------------------------
    For k=0,npairs-1 do begin

      ; Load the data
      data[k].flux = allcombframes[k].(ichip).flux[*,j]
      data[k].err = allcombframes[k].(ichip).err[*,j]
      data[k].mask = allcombframes[k].(ichip).mask[*,j]
      data[k].wavelength = allcombframes[k].(ichip).wavelength[*,j]
      data[k].sky = allcombframes[k].(ichip).sky[*,j]
      data[k].skyerr = allcombframes[k].(ichip).skyerr[*,j]
      data[k].telluric = allcombframes[k].(ichip).telluric[*,j]
      data[k].telluricerr = allcombframes[k].(ichip).telluricerr[*,j]

      ; Need to normalize them if they are not sky spectra
      if (fiberobjtype ne 'SKY') and not keyword_set(noscale) then begin

        ; Fit a low-order polynomial to the data
        ;  bad/missing values at the ends have flag=NAN
        if keyword_set(median) then begin
          bd=where(data[k].mask and badmask(),nbd)
          if nbd gt 0 then data[k].flux[bd]=!values.f_nan
	  scale=smooth(medfilt1d(data[k].flux,501,edge=2),100,/nan)
          if nbd gt 0 then data[k].flux[bd]=scale[bd]
        endif else begin
          gd = where(finite(data[k].flux) eq 1 and $
                   data[k].sky lt (2*median(data[k].flux)>2*median(data[k].sky)) and $
                   (data[k].mask and badmask()) eq 0,ngd)    ; only want good values
          if ngd gt 100 then  begin
            if keyword_set(nointerp) then y2=findgen(npix) else y2 = findgen(npix2)/2
            scalecoef = AP_ROBUST_POLY_FIT(y2[gd],data[k].flux[gd],5,status=status)
            scale = POLY(y2,scalecoef)
          endif else scale = dblarr(npix2)+1.0
        endelse
        ;if status eq 0 then stop
      endif else begin
        scale = dblarr(npix2)+1.0
      endelse

      data[k].flux = data[k].flux/scale      ; normalize spectrum
      data[k].err = data[k].err/scale        ; normalize error
      data[k].scale = scale

      ;stop

    Endfor ; pairs/spectra loop
    scales = data.scale
    ;scales = reform(data[*,*,ncol])
    sumscales = TOTAL(scales,2)

    ; Now combine the spectra and errors
    ;-------------------------------------
    dataspec = reform(data.flux)
    dataerr = reform(data.err)
    datamask = reform(data.mask)

    ; Combine the errors for the mean
    ;  this will be redone below
    comberr = sqrt( TOTAL(dataerr^2,2)/npairs^2 )    ; prop. of errors for mean

    nbdpix=0
    bdpix=where(datamask and badmask(),nbdpix)

    ; Do outlier rejection: skip, can fail in case of inhomogeneous S/N, etc.
;    medspec = MEDIAN(dataspec,dim=2)                        ; median spectrum
;    diffspec = dataspec - medspec#replicate(1.0,npairs)     ; difference spectrum
;    bdpix = where(abs(diffspec/dataerr) gt 5 OR $
;                  data.mask and badmask(),nbdpix)  ; 5 sigma outliers
    maskspec = dataspec
    maskerr = dataerr
    ;if nbdpix gt 0 then maskspec[bdpix] = !values.f_nan     ; set bad pix to NAN
    ;if nbdpix gt 0 then maskerr[bdpix] = !values.f_nan      ; set bad pix to NAN
    ;if nbdpix gt 0 then maskerr[bdpix] *= 10. 
    ; Number of "good" points for each pixel
    ;masknpairs = replicate(1L,npix2,npairs)
    ;if nbdpix gt 0 then masknpairs[bdpix] = 0
    ;npairs2d = TOTAL(masknpairs,2)               ; number of spectra/pairs for each pixel

    ; Take a weighted mean of the spectrum
    ;  wt = 1/error^2
    ;  weighted mean = Sum(wt*X)/Sum(wt)

    ; Global weight per spectrum
    ;---------------------------
    if keyword_set(globalwt) then begin

      ; Create the error array using same values for all
      ;   pixels of a spectrum
      mederr = MEDIAN(dataerr,dim=1)              ; median err per spectrum
      meddataerr = replicate(1.0,npix2)#mederr   ; make 2d array
      ; mask bad pixels
      if nbdpix gt 0 then meddataerr[bdpix] = !values.f_nan

      ; Now do a weighted mean (global) while ignoring the outliers (NANs)
      combspec = TOTAL(maskspec/meddataerr^2,2,/NAN) / TOTAL(1.0/meddataerr^2,2,/NAN)
      ;comberr = sqrt( TOTAL(maskerr^2,2,/NAN)/npairs2d^2 )  ; ignore NANs
      ;comberr = sqrt( 1./TOTAL(1./maskerr^2,1,/NAN) )   ; ignore NANs
      comberr = sqrt( 1./TOTAL(1./meddataerr^2,2,/NAN) )   ; ignore NANs

    ; Pixel-by-pixel weighting
    ;-------------------------
    endif else begin

      ; Now do a weighted mean (pixel-by-pixel) while ignoring the outliers (NANs)
      combspec = TOTAL(maskspec/maskerr^2,2,/NAN) / TOTAL(1.0/maskerr^2,2,/NAN)
      ; Combine the errors for mean and ignore outlier points
      ;comberr = sqrt( TOTAL(maskerr^2,2,/NAN)/npairs2d^2 )  ; ignore NANs
      comberr = sqrt( 1./TOTAL(1./maskerr^2,2,/NAN) )   ; ignore NANs

    endelse

    ; Initialize combined mask
    combmask = fix(combspec*0)

    ; Set bad pixels to zero
    bdpix = where(finite(combspec) eq 0,nbdpix)
    if nbdpix gt 0 then begin
      combspec[bdpix] = 0.0
      comberr[bdpix] = baderr()
      combmask[bdpix] = combmask[bdpix] or maskval('BADPIX')
    endif

    ; Rescale spec/error using the SUM of the continua
    combspec = combspec * sumscales
    comberr = comberr * sumscales

    ; Wavelength array
    ;  They are going to be slightly different because the measured dither
    ;  shift and wavelength zeropoint offset of done independently.
    ;  Use a mean of the wavelength zeropoint offset all other wavelength
    ;  coefficients should be the same (since they used the same apWave file).
    ; The offset is relative to the wavelength zeropoint of the first pair
    wcoef = reform(allcombframes[0].(ichip).wcoef[j,*])
    wcoef[0] = wcoef[0] + wpix0_offset   ; relative to the first
    combwave = PIX2WAVE(dindgen(npix2),wcoef)

    ; flags, sky, skyerr, telluric, telluric_err
    ;combflags = TOTAL(data.mask,2)            ; flags????
    for k=0,npairs-1 do combmask = combmask OR data[k].mask
    combsky = TOTAL(data.sky,2)               ; sum sky
    combskyerr = TOTAL(data.skyerr,2)         ; sum sky error
    combtel = TOTAL(data.telluric,2)/npairs   ; average telluric
    combtelerr = TOTAL(data.telluricerr,2)    ; sum telluric error

    ; Now stuff it in the output structure
    ;-------------------------------------
    outframe.(ichip).flux[*,j] = combspec
    outframe.(ichip).wavelength[*,j] = combwave
    outframe.(ichip).err[*,j] = comberr
    outframe.(ichip).mask[*,j] = combmask  ; combflags
    outframe.(ichip).sky[*,j] = combsky
    outframe.(ichip).skyerr[*,j] = combskyerr
    outframe.(ichip).telluric[*,j] = combtel
    outframe.(ichip).telluricerr[*,j] = combtelerr

    ; Update the wavelength coefficients
    outframe.(ichip).wcoef[j,*] = wcoef

    ; DO I NEED TO MAKE A NEW LSF ARRAY???????
    ; They are all using the SAME LSF calibration frame and it's
    ; only going to be slightly shifted.  So I think we can just use
    ; the LSF coefficients of the first frame

    ;if j eq 77 then stop

  Endfor ; chip loop

  ; Plotting
  ;----------
  ;pl = 1 ; 0 ;1
  if keyword_set(pl) then begin

    xx = [findgen(npix2), findgen(npix2)+npix2+300, findgen(npix2)+2*npix2+2*300]
    yy = [outframe.(0).flux[*,j], outframe.(1).flux[*,j], outframe.(2).flux[*,j] ]

    ;xr = [480,490]
    ;yr = minmax(combspec[xr[0]*2:xr[1]*2])
    ;xr = [0,npix*2]
    ;yr = [min(combspec),max(combspec)]
    xr = minmax(xx)
    ;yr = minmax(yy)
    yr = [0.0, median(yy)*2] 

    plot,[0],[0],/nodata,xr=xr,yr=yr,xs=1,ys=1,xtit='Pixel',ytit='Counts',$
         tit=strtrim(fiberobjid,2)+' (OBJTYPE='+fiberobjtype+')'
    for k=0,2 do oplot,findgen(npix2)+npix2*k+k*300,outframe.(k).flux[*,j]
    ;oplot,combspec
    ;oplot,sqrt(combvar),co=250
    ;oplot,xx,yy

    wait,1
    ;stop

  endif

  ;stop

  BOMB2:

End ; fiber loop

;stop

if keyword_set(stp) then stop

end
