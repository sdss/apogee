;+
; NAME:
;   sn
;
; PURPOSE:
;   Compute spectroscopic S/N (median and at 1.6 um)
;
; CALLING SEQUENCE:
;   snmed = sn(loglam, objflux, objivar)
;
; INPUTS:
;   loglam     - Log10 wavelength [NPIX]
;   objflux    - Flux [NPIX,NOBJ]
;   objivar    - Inverse variance in same units as OBJFLUX [NPIX,NOBJ]
;
; OPTIONAL KEYWORDS:
;
; OUTPUTS:
;   sn       - Median S/N in each of the SDSS filters per pix [5,NOBJ]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The median S/N is computed for only good data points that fall
;   within some hard-wired wavelength ranges defined by those
;   The routine is a copy of sn_median by Schlegel
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   16-Feb-2011  Written by D. Schlegel, LBL
;-  18-May-2011  Modified by A.E. Garcia Perez, UVA to calculate the
;  SN_median for APOGEE spectra
;------------------------------------------------------------------------------
function aspcap_snmedian, loglam, objflux, objivar

   w1 = [1514.0,1585.5,1647.4,1599.9]*1d1
   w2 = [1580.7,1643.3,1695.6,1600.1]*1d1

   ndim = size(objflux, /n_dimen)
   dims = size(objflux, /dimens)

   if (ndim EQ 1) then nobj = 1 $
    else nobj = dims[1]
   snmed = fltarr(4,nobj)
   sn=fltarr(nobj,2)
   for iobj=0L, nobj-1L do begin
      for j=0, 3 do begin
         indx = where(loglam GE alog10(w1[j]) AND loglam LT alog10(w2[j]) and objivar[*,iobj] gt 1e-35 $
          , ncts)
         if (ncts GT 1) then $
          snmed[j,iobj] = median( sqrt(objivar[indx,iobj])*objflux[indx,iobj] )
     endfor

      ;    sn[iobj,*]=[snmed[0,iobj],median(snmed[1:3,iobj])]
     sn[iobj,*]=[snmed[3,iobj],median(snmed[0:2,iobj])] ; to check
 
  endfor
   return, sn
end
;------------------------------------------------------------------------------
