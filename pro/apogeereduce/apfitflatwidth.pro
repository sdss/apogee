;+
; NAME:
;   apfitflatwidth
;
; PURPOSE:
;   Fit a traceset to the first-order corrected width of the flat field
;
; CALLING SEQUENCE:
;   widthset = fitflatwidth(flux, fluxivar, ansimage, [ fibermask, $
;    ncoeff=, sigma=, medwidth= ])
;
; INPUTS:
;   flux       - flat-field extracted flux
;   fluxivar   - corresponding inverse variance
;   ansimage   - output from extract image which contains parameter values
;
; OPTIONAL INPUTS:
;   fibermask  - nTrace bit mask, which marks bad fibers
;   ncoeff     - Order of legendre polynomial to apply to width vs. row;
;                default to 5.
;   sigma      - The SIGMA input to EXTRACT_IMAGE when determining ANSIMAGE;
;                default to 1.0 pix.  This can be a scalar, an [NFIBER] vector,
;                or an [NROW,NFIBER] array.
;
; OUTPUTS:
;   widthset   - Traceset structure containing fitted coefficients
;
; OPTIONAL OUTPUTS:
;   medwidth  - Median dispersion widths in each of the 4 quadrants
;               of the CCD, ordered LL,LR,UL,UR.
;
; COMMENTS:
;   The widths are forced to be the same as a function of row number
;   for all 16 fibers in each fiber bundle.
;
;   Used to fill flatstruct.widthset, which can then be applied
;   to object extraction (known profile widths).
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   xy2traceset
;
; REVISION HISTORY:
;   01-Mar-2000  Written by S. Burles, FNAL
;-
;------------------------------------------------------------------------------
function apfitflatwidth, flux, fluxivar, ansimage, fibermask, $
 ncoeff=ncoeff, sigma=sigma, medwidth=medwidth

   if (NOT keyword_set(ncoeff)) then ncoeff = 5
   if (NOT keyword_set(sigma)) then sigma = 1.0

   ntrace = (size(flux,/dimen))[1]
   nrow = (size(flux,/dimen))[0]
   ;if (ntrace NE 300) then $
   ; message, 'Must have 300 traces!'

   ;----------
   ; Generate a mask of good measurements based only upon the fibermask.

   mask = (flux GT 0) * (fluxivar GT 0) 

   if (keyword_set(fibermask)) then begin
      badflats = where(fibermask NE 0)
      if (badflats[0] NE -1) then mask[*,badflats] = 0
   endif

   ;----------
   ; Determine the widths from the output array from EXTRACT_IMAGE.

   igood = where(mask)
   widthterm = transpose(ansimage[lindgen(ntrace)*2+1,*])
   width = make_array(size=size(flux), /float)
   if (igood[0] NE -1) then $
    width[igood] = (1 + widthterm[igood] / flux[igood])

   ndim = size(sigma, /n_dimen)
   if (n_elements(sigma) EQ 1) then begin
      width = width * sigma[0]
   endif else if (ndim EQ 1) then begin
      for itrace=0, ntrace-1 do $
       width[*,itrace] = width[*,itrace] * sigma[itrace]
   endif else if (ndim EQ 2) then begin
      if (n_elements(sigma) NE n_elements(width)) then $
       message, 'Dimensions of SIGMA and WIDTH do not agree'
      width = width * sigma
   endif else begin
      message, 'Unsupported number of elements for SIGMA'
   endelse

   ;----------
   ; Compute the widths in each of 4 quandrants on the CCD

;   medwidth = [ median(width[0:nrow/2-1,0:ntrace/2-1]), $
;                median(width[0:nrow/2-1,ntrace/2:ntrace-1]), $
;                median(width[nrow/2:nrow-1,0:ntrace/2-1]), $
;                median(width[nrow/2:nrow-1,ntrace/2:ntrace-1]) ]
;
;   splog, 'Median spatial widths = ' $
;    + string(medwidth,format='(4f5.2)') + ' pix (LL LR UL UR)'


; APOGEE 
;
;   ;----------
;   ; Perform median across bundles on good arclines only
;   ; somewhat tedious, but it works
;
;   width = reform(width,nrow,20,16)
;   mask = reform(mask,nrow,20,16)
;   width_bundle = fltarr(nrow,16)
;
;   for irow=0, nrow-1 do begin
;      for j=0, 15 do begin
;         ss = where(mask[irow,*,j])
;         if (ss[0] NE -1) then $
;          width_bundle[irow,j] = djs_median(width[irow,ss,j])
;      endfor
;   endfor
;
;   width_final = rebin(width_bundle, nrow, ntrace, /sample)
width_final = width

   ; Make sure that all of the widths are positive
   ;------------------------------------------------
   bad = where(width_final le 0.0,nbad)
   if nbad gt 0 then begin

     sm_width_final0 = MEDFILT2D(width_final,21,dim=1,/edge)
     minfiber = min(sm_width_final0,dim=1)
     bdfiber = where(minfiber le 0.0,nbdfiber)

     ; Replace bad points with polynomial interpolation
     x = findgen(nrow)
     for i=0,nbdfiber-1 do begin
       bdpix = where(width_final[*,bdfiber[i]] le 0.0,nbdpix,complement=gdpix)
       ;coef = poly_fit(x[gdpix],width_final[gdpix,bdfiber[i]],8) ;ncoeff)
       if n_elements(gdpix) gt nrow/2 or i eq 0 then begin
         coef = poly_fit(x[gdpix],width_final[gdpix,bdfiber[i]],2)
         width_final[bdpix,bdfiber[i]] = poly(x[bdpix],coef)
       endif else width_final[bdpix,bdfiber[i]] = width_final[bdpix,bdfiber[i-1]]
     end

     sm_width_final = MEDFILT2D(width_final,51,dim=1,/edge)

     ;med_width = median(width_final,dim=1)
     ;stop

   end


   ;----------
   ; Turn the widths back into a traceset.

   apxy2traceset, findgen(nrow) # replicate(1,ntrace), $
    width_final, widthset, ncoeff=ncoeff, xmin=xmin, xmax=xmax


   ; Identify fibers with negative widths and fix them
   ;--------------------------------------------------
   traceset2xy, widthset, xx, width_fit0
   minfiber = min(width_fit0,dim=1)
   bdfiber = where(minfiber le 0.0,nbdfiber)

   ; use median smoothed coeffs for the bad fibers
   sm_coeff = MEDFILT2D(widthset.coeff,11,dim=2,/edge)  ; smooth the coeffs
   for i=0,nbdfiber-1 do begin
     widthset.coeff[*,bdfiber[i]] = sm_coeff[*,bdfiber[i]]
   end

   ;----------
   ; Compute the widths in each of 4 quandrants on the CCD.

   traceset2xy, widthset, xx, width_fit
   medwidth = [ median(width_fit[0:nrow/2-1,0:ntrace/2-1]), $
                median(width_fit[0:nrow/2-1,ntrace/2:ntrace-1]), $
                median(width_fit[nrow/2:nrow-1,0:ntrace/2-1]), $
                median(width_fit[nrow/2:nrow-1,ntrace/2:ntrace-1]) ]
   splog, 'Median spatial widths = ' $
    + string(medwidth,format='(4f5.2)') + ' pix (LL LR UL UR)'

   rms = mad(width-width_fit,/zero)
   print,'  RMS of width fit = ',string(rms,format='(F5.3)'),' pixels'

   ;stop

   return, widthset
end
;------------------------------------------------------------------------------
