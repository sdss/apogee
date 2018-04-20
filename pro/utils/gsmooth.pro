function gsmooth,arr,fwhm,widfwhm=widfwhm0,stp=stp,squared=squared,error=error

;+
;
; GSMOOTH
;
; This does Gaussian smoothing on a 1D or 2D array.
;
; INPUTS:
;  arr       The array to smooth
;  fwhm      The Gaussian FWHM to use.  Can be a scalar
;              or two elements if an elliptical Gaussian
;              is desired for a 2D image.  To smooth a 2D
;              image in only one dimension set the FWHM
;              of the other dimension to 0 (i.e. FWHM=[0,10]).
;  =widfwhm  Width in FWHM.  Default is 6.0
;  /squared  Convolve with the SQUARE of the Gaussian.  This is
;              useful for propagating errors properly.  I.e. for
;              Gaussian smoothing you should do,
;                smerr = sqrt(gsmooth(err^2,fwhm,/squared))
;  =error    The error message if one occurred.
;  /stp      Stop at the end of the program
;
; OUTPUTS:
;  sm      The smoothed array
;
; USAGE:
;  IDL>sm = gsmooth(arr,5.0)
;
; By D.Nidever   Oct 20008
;-

sz = size(arr)
narr = n_elements(arr)
nfwhm = n_elements(fwhm)

; Not enough inputs
if narr eq 0 or nfwhm eq 0 then begin
  print,'Syntax -  sm = gsmooth(arr,fwhm,widfwhm=widfwhm,error=error,/stp)'
  error = 'Not enough inputs'
  return,-1
endif

; Dimensionality too large
if (sz[0] gt 2) then begin
  error = 'Only 1D or 2D arrays are allowed'
  print,error
  return,-1
endif

;; Zero or negative FWHM
;if min(fwhm) le 0.0 then begin
;  error = 'FWHM must be positive'
;  print,error
;  return,arr
;endif

; ONE-DIMENSIONAL
;----------------
if (sz[0] eq 1) then begin

  ; Zero or negative FWHM
  if min(fwhm) le 0.0 then begin
    error = 'FWHM must be positive'
    print,error
    return,arr
  endif

  if n_elements(widfwhm0) gt 0 then widfwhm=widfwhm0[0] else widfwhm=6.0

  npix = round(widfwhm*fwhm[0]) < narr
  psf = psf_gaussian(npixel=npix,fwhm=fwhm[0],ndim=1,/norm)
  ;psf = psf_gaussian(npixel=sz[1],fwhm=fwhm[0],ndim=1,/norm)
  if not keyword_set(squared) then begin
    sm = CONVOL(arr,psf,/center,/edge_truncate,/nan,/normalize)
    ;sm = CONVOL(arr,psf,/center,/edge_wrap,/nan,/normalize)
  endif else begin
    sm = CONVOL(arr,psf^2,/center,/edge_truncate,/nan)
  endelse

endif

; TWO-DIMENSIONAL
;----------------
if (sz[0] eq 2) then begin

  if n_elements(widfwhm0) gt 0 then widfwhm=widfwhm0[0] else widfwhm=6.0

  if nfwhm eq 2 then fw = fwhm else fw=[1,1]*fwhm[0]

  ; Only smooth in one dimension - X
  if fw[0] gt 0 and fw[1] eq 0 then begin

    npix = round(widfwhm*fw[0]) < sz[1]
    psf = psf_gaussian(npixel=[npix,1],fwhm=[fw[0],1d10],ndim=2,/norm)
    ;psf = psf_gaussian(npixel=sz[1],fwhm=fwhm[0],ndim=1,/norm)
    if not keyword_set(squared) then begin
      sm = CONVOL(arr,psf,/center,/edge_truncate,/nan,/normalize)
    endif else begin
      sm = CONVOL(arr,psf^2,/center,/edge_truncate,/nan)
    endelse

    return,sm
  end
  ; Only smooth in one dimension - Y
  if fw[0] eq 0 and fw[1] gt 0 then begin

    npix = round(widfwhm*fw[1]) < sz[2]
    psf = psf_gaussian(npixel=[1,npix],fwhm=[1d10,fw[1]],ndim=2,/norm)
    ;psf = psf_gaussian(npixel=sz[1],fwhm=fwhm[0],ndim=1,/norm)
    if not keyword_set(squared) then begin
      sm = CONVOL(arr,psf,/center,/edge_truncate,/nan,/normalize)
    endif else begin
      sm = CONVOL(arr,psf^2,/center,/edge_truncate,/nan)
    endelse

    return,sm
  end


  ; Zero or negative FWHM
  if min(fwhm) le 0.0 then begin
    error = 'FWHM must be positive'
    print,error
    return,arr
  endif

  sm = filter_image(arr,fwhm=fw,/all_pixels,psf=psf)

  if keyword_set(squared) then $
    sm = filter_image(arr,fwhm=fw,/all_pixels,psf=psf^2)

  ;npix1 = round(6.0*fw[0])
  ;npix2 = round(6.0*fw[1])
  ;psf = psf_gaussian(npixel=[npix1,npix2],fwhm=fw,ndim=2,/norm)
  ;;psf = psf_gaussian(npixel=[sz[1],sz[2]],fwhm=fw,ndim=2,/norm)
  ;sm = CONVOL(arr,psf,/center,/edge_truncate)

endif

if keyword_set(stp) then stop

return,sm

end
