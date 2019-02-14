pro ap1dfluxing,frame,plugmap,outframe,silent=silent,verbose=verbose,pl=pl,stp=stp
;+
;
; AP1DFLUXING
;
; This does a rough flux calibration
;
; INPUTS:
;  frame     A structure with the header/data information for a
;                dither combined frame that has also been wavelength
;                calibrated and airglow line subtracted
;  plugmap   The Plug Map structure for this plate
;  /silent   Don't print anything to the screen.
;  /verbose  Print lots of information to the screen
;  /pl       Make some plots
;  /stp      Stop at the end of the program.
;
; OUTPUTS:
;  outframe  The same frame but with the spectra/errors flux calibrated.
;
; USAGE:
;  IDL>ap1dfluxing,frame,plugmap,outframe
;
; By D. Nidever  May 2010
; mods, J. Holtzman 2013 and 2019
;-

apgundef,outframe

nframe = n_elements(frame)
nplugmap = n_elements(plugmap)
sz = size(frame.chipa.flux)
npix = sz[1]
nfibers = sz[2]

; Not enough inputs
if nframe eq 0 or nplugmap eq 0 then begin
  print,'Syntax - ap1dfluxing,frame,plugmap,outframe,silent=silent,verbose=verbose,pl=pl,stp=stp'
  stop
  return
endif

outframe = frame   ; the format is the same

; first get relative flux curve using tellurics
j = where(plugmap.fiberdata.spectrographid eq 2 and $
             plugmap.fiberdata.holetype eq 'OBJECT' and $
             plugmap.fiberdata.objtype eq 'HOT_STD',nstars)
tell=plugmap.fiberdata[j].fiberid

; do polynomial fit to log(flux), with 4th order plus offset fo each star,
; using every 10th pixel in each chip, so we have 190 pixels * 3 chips * ntelluric data points
; and 4 + ntellurics parameters
npix=190
design=fltarr(3*npix*nstars,4+nstars)
y=fltarr(3*npix*nstars)
for ichip=0,2 do begin
  x=outframe.(ichip).wavelength - 16000.
  for irow=0,nstars-1 do begin
    row=300-tell[irow]
    design[irow*3*npix+ichip*npix:irow*3*npix+ichip*npix+npix-1,0] = x[100:1990:10,row]^4
    design[irow*3*npix+ichip*npix:irow*3*npix+ichip*npix+npix-1,1] = x[100:1990:10,row]^3
    design[irow*3*npix+ichip*npix:irow*3*npix+ichip*npix+npix-1,2] = x[100:1990:10,row]^2
    design[irow*3*npix+ichip*npix:irow*3*npix+ichip*npix+npix-1,3] = x[100:1990:10,row]
    design[irow*3*npix+ichip*npix:irow*3*npix+ichip*npix+npix-1,4+irow] = 1.
    y[irow*3*npix+ichip*npix:irow*3*npix+ichip*npix+npix-1] = alog10(outframe.(ichip).flux[100:1990:10,row])
  endfor
endfor
gd=where(finite(y) eq 1)
design=design[gd,*]
y=y[gd]
; do the fit
a=matrix_multiply(design,design,/atranspose)
b=matrix_multiply(design,y,/atranspose)
pars=invert(a)#b

; with plot option, show results for tellurics
if keyword_set(pl) then begin
 for irow=0,nstars-1 do begin
  row=300-tell[irow]
  for ichip=0,2 do begin
    w=outframe.(ichip).wavelength[*,row]
    spec=outframe.(ichip).flux[*,row]
    x = w-16000.
    logflux = pars[0]*x^4 + pars[1]*x^3 + pars[2]*x^2 + pars[3]*x
    logflux += 2*alog10(w/16000.)
    resp= 10.^logflux
    if ichip eq 0 then plot,w,spec,xr=[15100,17000] else oplot,w,spec
    oplot,w,spec/resp,color=255

  endfor
  stop
 endfor
endif

; apply the fit. Note that a term is added so that response gives 1/lambda**-4 shape
for ichip=0,2 do begin
  for irow=0,299 do begin 
    w=outframe.(ichip).wavelength[*,row]
    spec=outframe.(ichip).flux[*,row]
    x = w-16000.
    logflux = pars[0]*x^4 + pars[1]*x^3 + pars[2]*x^2 + pars[3]*x
    logflux += 4*alog10(w/16000.)
    resp= 10.^logflux
    outframe.(ichip).flux[*,irow] /= resp
    bderr=where(outframe.(ichip).err[*,irow] eq baderr(),nbd)
    outframe.(ichip).err[*,irow] /= resp
    if nbd gt 0 then outframe.(ichip).err[bderr,irow] = baderr()
    outframe.(ichip).sky[*,irow] /= resp
    outframe.(ichip).skyerr[*,irow] /= resp
  endfor
endfor

; simple absolute normalization based on H magnitude, since conversion to F_lambda has already been done with response curve
skyind = where(plugmap.fiberdata.spectrographid eq 2 and $
               plugmap.fiberdata.holetype eq 'OBJECT' and $
               plugmap.fiberdata.objtype eq 'SKY',nskyind)
if nskyind gt 0 then begin
  fiber=plugmap.fiberdata[skyind].fiberid
  medsky=median(frame.chipb.flux[*,300-fiber])
endif else medsky=0
objind = where(plugmap.fiberdata.spectrographid eq 2 and $
               plugmap.fiberdata.holetype eq 'OBJECT' and $
               plugmap.fiberdata.objtype ne 'SKY',nobjind)
zero=fltarr(n_elements(objind))
norm=fltarr(n_elements(objind))
print,'Fluxing: '
print,'FIBER   H   Zero   Norm'
fluxcorr = fltarr(nfibers)
for istar=0,n_elements(objind)-1 do begin
 fiber=plugmap.fiberdata[objind[istar]].fiberid
 if fiber gt 0 and fiber le 300 then begin  
  fiber_mag = plugmap.fiberdata[objind[istar]].mag
  hmag=fiber_mag[1]
  fiber_ra = plugmap.fiberdata[objind[istar]].ra
  tmass_name = plugmap.fiberdata[objind[istar]].tmass_style
  ; try to look up H mag from catalog, since plugmaps can sometimes have zero!
  ; now that we get "plugmap" info from plateHoles, should be OK?
  ;cat=getcat(strtrim(tmass_name,2),ra=fiber_ra)
  ;if size(cat,/type) eq 8 then hmag=cat.h
  medflux=median(frame.chipb.flux[*,300-fiber])
  if hmag ne 0 and hmag lt 30 and medflux-medsky gt 100 then begin
    zero[istar]=hmag+2.5*alog10((medflux-medsky)>1.)

  ;  http://www.ipac.caltech.edu/2mass/releases/allsky/doc/sec6_4a.html
  ; Table 1 - 2MASS Isophotal Bandpasses and Fluxes-for-0-magnitude from Cohen et al. (2003)
  ; Band   Lambda (um) 	   Bandwidth (um)   Fnu - 0 mag (Jy) Flambda - 0 mag (W cm-2 µm-1)
  ; H 	1.662 +/- 0.009   0.251 +/- 0.002   1024 +/- 20.0   1.133E-13 +/- 2.212E-15
    norm[istar]=10^(-0.4*hmag)*1.133e-10/((medflux-medsky)>1.)   ; ers/cm^2/s/A
  endif else begin
    zero[istar]=!values.f_nan
    norm[istar]=!values.f_nan
  endelse
 endif
endfor
mednorm=median(norm)
if not finite(mednorm) then mednorm=1.

for istar=0,n_elements(objind)-1 do begin
 fiber=plugmap.fiberdata[objind[istar]].fiberid
 if fiber gt 0 and fiber le 300 then begin  
  ; deal with objects that have no hmag
  if finite(norm[istar]) eq 0 then norm[istar]=mednorm
  fluxcorr[300-fiber]=norm[istar]
  print,fiber,zero[istar],norm[istar]
  for ichip=0,2 do begin
    ; Calibrate the flux,error,sky,and skyerror spectra
    outframe.(ichip).flux[*,300-fiber] *= norm[istar]
    bderr=where(outframe.(ichip).err[*,300-fiber] eq baderr(),nbd)
    outframe.(ichip).err[*,300-fiber] *= norm[istar]
    if nbd gt 0 then outframe.(ichip).err[bderr,300-fiber] = baderr()
    outframe.(ichip).sky[*,300-fiber] *= norm[istar]
    outframe.(ichip).skyerr[*,300-fiber] *= norm[istar]
  endfor
 endif
endfor
; correct the sky fibers using the median normalization
for istar=0,nskyind-1 do begin
 fiber=plugmap.fiberdata[skyind[istar]].fiberid
 if fiber gt 0 and fiber le 300 then begin  
  fluxcorr[300-fiber]=mednorm
  for ichip=0,2 do begin
    outframe.(ichip).flux[*,300-fiber] *= mednorm
    bderr=where(outframe.(ichip).err[*,300-fiber] eq baderr(),nbd)
    outframe.(ichip).err[*,300-fiber] *= mednorm
    if nbd gt 0 then outframe.(ichip).err[bderr,300-fiber] = baderr()
    outframe.(ichip).sky[*,300-fiber] *= mednorm
    outframe.(ichip).skyerr[*,300-fiber] *= mednorm
  endfor
 endif
endfor
; Add flux correction factors to structure
outframe = create_struct(outframe,'FLUXCORR',fluxcorr)
return

;==========================
; remainder is unused
;=========================
; Constants
;h = 6.626076d-34  ; planck's constant (J s)
hconst_cgs = 6.626076d-27          ; planck's constant (erg s)
cspeed_cgs = 2.99792458d10         ; speed of light in cm/s
hc_cgs = hconst_cgs * cspeed_cgs   ; erg cm 

; Get APOGEE directories
dirs=getdir(apogee_dir,cal_dir,spectro_dir,apogee_vers,lib_dir)
fluxlib_dir = lib_dir+'flux/'
if FILE_TEST(fluxlib_dir,/directory) eq 0 then begin
  print,'FLUX LIBRARY Directory ',flux_dir,' NOT FOUND'
  return
endif

; Checking the tags of the input structure
tags = tag_names(frame)
needtags1 = ['CHIPA','CHIPB','CHIPC']
for i=0,n_elements(needtags1)-1 do begin
  if (where(tags eq needtags1[i]))[0] eq -1 then begin
    print,'TAG ',needtags1[i],' NOT FOUND in input structure'
    return
  end
end
needtags2 = ['HEADER','FLUX','ERR','MASK','WAVELENGTH','SKY','SKYERR',$
             'TELLURIC','TELLURICERR','LSFCOEF','WCOEF']
for i=0,2 do begin
  tags2 = tag_names(frame.(i))
  for j=0,n_elements(needtags2)-1 do begin
    if (where(tags2 eq needtags2[j]))[0] eq -1 then begin
      print,'TAG ',needtags2[j],' NOT FOUND in input structure'
      return
    end
  end
end

sz = size(frame.chipa.flux)
npix = sz[1]
pix = findgen(npix)
nfibers = sz[2]
; Is this a dither-combined spectrum?
xscale = 1    ; assume original non-dither combined spectrum
if npix eq 4096 then xscale = 2

; Initialize outframe, add telluric and error_telluric
outframe = frame   ; the format is the same


; Load the H-band passband curve
;  http://www.ipac.caltech.edu/2mass/releases/allsky/doc/sec6_4a.html
Hpassband0 = IMPORTASCII(fluxlib_dir+'hband_passband.txt',fieldnames=['WAVE','FLUX'],$
                        fieldtypes=[5,5],/silent)
; wavelengths in microns, convert to Angstroms
Hpassband0.wave *= 1e4

; Spline onto a new grid with a similar disperion to the observed spectrum
nHflux = range(Hpassband0.wave)/0.14 + 1
Hpassband = REPLICATE({wave:0.0d0,flux:0.0d0},nHflux)
Hpassband.wave = dindgen(nHflux)*0.14d0 + min(Hpassband0.wave)
Hpassband.flux = SPLINE(Hpassband0.wave,Hpassband0.flux,Hpassband.wave,/double)
; Compute TOTALS to be use to correct for the fact that we don't cover
; the entire H-band.
dHwave = slope(Hpassband0.wave)
dHwave = [dHwave,first_el(dHwave,/last)]
;Htot0 = TOTAL(Hpassband0.wave*Hpassband0.flux*dHwave)
;Htot = TOTAL(Hpassband.wave*Hpassband.flux*0.14)
Htot0 = TOTAL(Hpassband0.flux*dHwave)
Htot = TOTAL(Hpassband.flux*0.14)


; Get global observation parameters
gain = sxpar(frame.chipa.header,'GAIN',count=ngain)
if ngain eq 0 then gain=2.0  ; 1.0
exptime = sxpar(frame.chipa.header,'EXPTIME')
;airmass = sxpar(frame.chipa.header,'AIRMASS')
altitude = sxpar(frame.chipa.header,'ALT',count=nalt)
if nalt eq 0 or strtrim(altitude,2) eq 'NAN' then airmass=1.0
if nalt gt 0 then airmass=1.0/cos((90.0-altitude)/!radeg)


; SDSS telescope collecting area
;  2.5m diameter primary mirror
;  1.08m diamter secondary mirror
telarea = !dpi*(250.0/2.)^2 - !dpi*(108.0/2.)^2   ; in cm^2


;--------------------------
; Flux Correct all fibers
;--------------------------

; Loop through the fibers
fluxcorr_array = fltarr(npix,3,nfibers)
For i=0,nfibers-1 do begin

  ifiber = i

  if keyword_set(verbose) then $
    print,'Fiber ',strtrim(ifiber+1,2)


  ;-------------------------------------------------------
  ; Set on the absolute flux scale using 2MASS photometry
  ;-------------------------------------------------------

  ; Query the database to get the 2MASS photometry for this star
  phot = {j:99.0,j_err:0.02,h:99.0,h_err:0.02,K:99.0,K_err:0.02}
  plugind = where(plugmap.fiberdata.spectrographid eq 2 and $
                  plugmap.fiberdata.fiberid eq 300-i,nplugind)
  if nplugind gt 0 then begin
    objtype = plugmap.fiberdata[plugind].objtype
    ; Real object, use Hmag
    if objtype ne 'SKY' then begin
      phot.j = plugmap.fiberdata[plugind].mag[0]
      phot.h = plugmap.fiberdata[plugind].mag[1]
      phot.k = plugmap.fiberdata[plugind].mag[2]
    ; Sky fiber, use average flux calibration factor
    endif else begin
      ;print,'Fiberid=',strtrim(300-i,2),' is a sky fiber. Using average flux calibration factor'
    endelse
  endif else begin
    print,'No information for FiberID=',strtrim(300-i,2),' in the plugmap file. Using average flux calibration factor'
  endelse


  ; We need to see what the flux is in the H-band for this spectrum.
  
  ;  http://www.ipac.caltech.edu/2mass/releases/allsky/doc/sec6_4a.html
  ; Table 1 - 2MASS Isophotal Bandpasses and Fluxes-for-0-magnitude from Cohen et al. (2003)
  ; Band   Lambda (um) 	   Bandwidth (um)   Fnu - 0 mag (Jy) Flambda - 0 mag (W cm-2 µm-1)
  ; J 	1.235 +/- 0.006   0.162 +/- 0.001   1594 +/- 27.8   3.129E-13 +/- 5.464E-15
  ; H 	1.662 +/- 0.009   0.251 +/- 0.002   1024 +/- 20.0   1.133E-13 +/- 2.212E-15
  ; Ks 	2.159 +/- 0.011   0.262 +/- 0.002   666.7 +/- 12.6  4.283E-14 +/- 8.053E-16
  ;
  ; Figures 1, 2, and 3 present the 2MASS J, H and Ks relative spectral response curves (RSRs),
  ; peak-normalized to unity, derived by Cohen et al. (2003). As stated by these authors, these
  ; curves "are designed to be integrated directly over stellar spectra in Flambda form, in order 
  ; to calculate synthetic photometric magnitudes. The QE-based component was converted to yield 
  ; photon-counting RSRs by multiplying by wavelength and renormalized, as described by Bessel
  ; (2000)." These RSRs are consistent with the absolute calibration of 2MASS given in Table 1 above. 
  ; Bessel (2000), PASP, 112, 961
  ; Cohen et al. (2003), AJ, 126, 1090

  ; The equation from Bessel (2000) is:
  ;  Integral( [f(lambda)/(h*nu)] * Rx(lambda) * dlambda  )
  ;  = (1/hc) * Integral( f(lambda) * [lambda * Rx(lambda)] * dlambda  )
  ;  where Rx(lambda) is the response function (our Hpassband file).



  ; ABSOLUTE CALIBRATION STEPS:
  ; Need to correct for exptime, collecting area, airmass
  ; gain, QA (?) in order to get the flux in physical units.
  ; Also need to convert Counts into energy units by using h*nu

  ; I think it's best to convert the spectrum into flux units
  ; using the steps outline above, then filter with the passband
  ; and integrate over all pixels.
  ; Then compare it to the expected total integrated flux based
  ; on the 2MASS zero-points rescaled to the observed H-band
  ; magnitude of this star.

  ; I probably should keep a running array of the conversion
  ; factor so I can then apply it to the other arrays (variance,
  ; sky, etc.).

  ; Computing "global" conversion factor, multiplicative
  ;  to convert from ADU to flux in ergs/s/cm^2/Ang
  convfactor = 1.0d0      ; initialize
  convfactor *= gain      ; multiply by gain to convert ADU->electrons
 ; DON'T NEED THE GAIN CORRECTION, since we already do this in ap3dproc.pro!!!!
  convfactor /= exptime   ; divide by exptime in seconds
  convfactor /= telarea   ; divide by telescope area in cm^2
  ; Airmass/extinction correction
  ;   http://www.ipac.caltech.edu/2mass/releases/allsky/doc/sec4_8.html#zeropoint
  ;   http://www.ipac.caltech.edu/2mass/releases/allsky/doc/sec4_8.tbl2.html
  ;  but 2MASS used c2*(X-1.0), not sure why they subtracted 1.
  Hext_term = 0.03    ; mag/airmass
  convfactor /= 10.^(-Hext_term * airmass/2.5)   ; correct for airmass
  ; Convert to energy units, h*nu = h*c/lambda, units of 
  convfactor *= hc_cgs
  ; This multiplicative factor converts from ADU -> ergs cm/s/cm^2
  ; Still need to divide by lambda in cm, and dlambda in Ang to
  ; get the final correct flux units of ergs/s/cm^2/Ang



  ; Convert Spectrum to flux units, filter, and integrate.
  spec1 = outframe.(0).flux[*,ifiber]
  spec2 = outframe.(1).flux[*,ifiber]
  spec3 = outframe.(2).flux[*,ifiber]
  spec = [spec1,spec2,spec3]
  wave1 = outframe.(0).wavelength[*,ifiber]  ; Wavelength in Angstroms
  wave2 = outframe.(1).wavelength[*,ifiber]
  wave3 = outframe.(2).wavelength[*,ifiber]
  wave = [wave1,wave2,wave3]
  sky = [outframe.(0).sky[*,ifiber], outframe.(1).sky[*,ifiber], outframe.(2).sky[*,ifiber]]

  ; Our pixel width is twice the dispersion since we have dither combined.
  dwave1 = abs(slope(wave1))                 ; dWavelength in Angstroms
  dwave1 = [dwave1,first_el(dwave1,/last)]
  dwave2 = abs(slope(wave2))
  dwave2 = [dwave2,first_el(dwave2,/last)]
  dwave3 = abs(slope(wave3))
  dwave3 = [dwave3,first_el(dwave3,/last)]
  dwave = [dwave1,dwave2,dwave3] * xscale    ; dither-combined pixels?

  ; Convert the spectrum from ADU to flux units
  ;---------------------------------------------
  spec_flux = spec * convfactor
  spec_flux /= wave*1d-8  ; divide by lambda in cm
  spec_flux /= dwave      ; divide by dlambda in Ang
  ; the units are now ergs/s/cm^2/Ang


  ; Filter the fluxed spectrum and integrate
  ;--------------------------------------------

  ; Spline the passband to our wavelengths
  Hpassobs = wave*0.0d0
  wsi = sort(wave)
  Hpassobs[wsi] = spline(Hpassband0.wave,Hpassband0.flux,wave[wsi],/double)

  ; Do the filtering and integrate the spectrum
  ;  need to multiply by dWave to get actual flux
  gdpix = where(spec gt 1 and finite(spec) eq 1,ngdpix)
  if ngdpix gt 0 then begin
    ;spectot_cgs = TOTAL(spec_flux[gdpix] * Hpassobs[gdpix] * dwave[gdpix])

    ; Use a polynomial fit to the continuum, robust to sky lines
    pcoef = ROBUST_POLY_FITQ( (wave-16000.)/900., spec_flux, 4)
    spec_flux_cont = poly( (wave-16000.)/900., pcoef)
    spectot_cgs = TOTAL(spec_flux_cont[gdpix] * Hpassobs[gdpix] * dwave[gdpix])
  endif else begin
    ; This is at least ballpark
    spectot_cgs = TOTAL( 1 * Hpassobs * dwave)
  endelse

  ; We need to correct for the fact that we don't observe the entire H-band
  ;--------------------------------------------------------------------------
  ; Figure out how much of the H-band "flux" we cover and correct for that.
  wr = [minmax(outframe.(0).wavelength[*,ifiber]),$  ; wavelength coverage
        minmax(outframe.(1).wavelength[*,ifiber]),$
        minmax(outframe.(2).wavelength[*,ifiber])]
  ; These are average values
  if max(wr) lt 1e4 then $
    wr = [16472.948d, 16954.926d, 15856.186d, 16434.144d, 15142.627d, 15809.923d]

  ; Get interpolated Hpassband pixels that are within the windows we have
  ;  observed
  gd = where( ( Hpassband.wave ge min(wr[0:1]) and Hpassband.wave le max(wr[0:1]) ) OR $
              ( Hpassband.wave ge min(wr[2:3]) and Hpassband.wave le max(wr[2:3]) ) OR $
              ( Hpassband.wave ge min(wr[4:5]) and Hpassband.wave le max(wr[4:5]) ), ngd)
  ; Compute fraction of H-band we have observed
  ;Htotobs = TOTAL(Hpassband[gd].wave*Hpassband[gd].flux*0.14)
  ;Hbandfrac = Htotobs/Htot0
  ;Htotobs = TOTAL(Hpassband[gd].flux*0.14)
  Hbandfrac = TOTAL(Hpassband[gd].flux*0.14)/TOTAL(Hpassband.flux*0.14)
  ; this is around 0.61

  ; Correct the integrated spectrum
  spectot_cgs_corr = spectot_cgs / Hbandfrac


  ; Figure out what total integrated flux we expect from the 2MASS H-band
  ; zeropoints and the magnitude of the star
  ;-----------------------------------------------------------------------
  ;  Flambda = 1.33E-13 for a 0th magnitude star in the H-band
  ;   in W/cm^2/micron
  spectot_expected = 1.33d-13
  spectot_expected *= 10.^(-phot.h/2.5)  ; correct for the star's magnitude
  ; convert to ergs/s/cm^2
  ;  1 ergs/s = 1E-7 W
  ;  1 Ang = 1E-4 micron
  spectot_expected *= 1e7
  ; This is still a flux density in microns
  ; Multiply by filter width to get actual flux
  ; The H-band filter width is 0.251 microns (see table above)
  spectot_expected *= 0.251   ; ergs/s/cm^2

  ; Calculate the absolute flux correction
  ;----------------------------------------
  ; This should account for the efficiency of the telecope/fiber/instrument
  ; This is really 1/efficiency of the system which is ~10%
  ; So should be close to ~10
  ;  multiplicative,  spec_absflux = spec * absflux_correction
  absflux_correction = spectot_expected / spectot_cgs_corr
  
  ; Constrain this to a reasonable range
  absflux_correction = ( absflux_correction > 0.5 ) < 50.0
  if finite(absflux_correction) eq 0 then absflux_correction = 8.0

  ; Use average value for stars with no plugmap/magnitude information
  ;  and sky fibers.  They have Hmag=99.0
  if phot.h gt 50. then absflux_correction = 8.0


  ;---------------------------------------------------
  ; Absolute Flux correct all of the flux quantities
  ;---------------------------------------------------
  ; Three steps:
  ; 1.) multiply by "convfactor"
  ; 2.) divide by lambda in cm and dlambda in Ang (to get units of ergs/s/cm^2/Ang)
  ; 3.) multiply by "absflux_correction"

  ; Loop through the chips
  For j=0,2 do begin

    ; Get lambda and dlambda in the correct units
    lambda = reform( outframe.(j).wavelength[*,ifiber] )
    dlambda = abs(slope(lambda))
    dlambda = [dlambda, first_el(dlambda,/last)] * xscale  ;in Ang
    lambda *= 1d-8  ; convert from Ang to cm
    
    ; Make a correction array
    correction_array = convfactor
    correction_array /= lambda
    correction_array /= dlambda
    correction_array *= absflux_correction

    ; Add to fluxcorr_array
    fluxcorr_array[*,j,i] = correction_array

    ;  The planes are: [spec, wave, error, flag, sky, errsky,
    ;       telluric, error_telluric]

    ; Calibrate the spectrum
    outframe.(j).flux[*,ifiber] *= correction_array
    ; Calibrate the error spectrum
    bderr=where(outframe.(j).err[*,ifiber] eq baderr(),nbd)
    outframe.(j).err[*,ifiber] *= correction_array
    if nbd gt 0 then outframe.(j).err[bderr,ifiber] = baderr()
    ; Calibrate the sky spectrum
    outframe.(j).sky[*,ifiber] *= correction_array
    ; Calibrate the sky error spectrum
    outframe.(j).skyerr[*,ifiber] *= correction_array

    ;stop

  End

  BOMB:

  ;stop

End

; Add flux correction factors to structure
outframe = create_struct(outframe,'FLUXCORR',fluxcorr_array)

;stop

if keyword_set(stp) then stop

end
