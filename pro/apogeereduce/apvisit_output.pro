pro apvisit_output,frame,plugmap,shiftstr,pairstr,silent=silent,stp=stp,single=single,iter=iter,mjdfrac=mjdfrac,survey=survey

;+
;
; APVISIT_OUTPUT
;
; This outputs the final dither-combined and sky-corrected
; plate frame and individual visit spectra.
; This is called from ap1dvisit.pro
;
; INPUTS:
;  frame         The structure that contains the dither-combined
;                 and sky-corrected frame with all three chips.
;  plugmap       The Plug Map structure for this plate
;  /silent       Don't print anything to the screen.
;  /stp          Stop at the end of the program.
;
; OUTPUTS:
;  The frame is written to the "plate" directory, and the
;  individual spectra are written to the "visit" directory.
;
; USAGE:
;  IDL>apvisit_output,frame
;
; By D.Nidever  May 2010
; Modifications J. Holtzman 2011+
;-

dirs=getdir(apred_vers=apred_vers)
nframe = n_elements(frame)
nplugmap = n_elements(plugmap)

; Not enough inputs
if nframe eq 0 or nplugmap eq 0 then begin
  print,'Syntax - apvisit_output,frame,plugmap,silent=silent,stp=stp'
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


; Get the frame information
plate = plugmap.plateid
mjd = plugmap.mjd
platemjd5 = strtrim(plate,2)+'-'+strtrim(mjd,2)
locid = plugmap.locationid

sz = size(frame.chipa.flux)
npix = sz[1]
nfibers = sz[2]

svn_vers = getvers()


;##############################################################
;  OUTPUT THE PLATE FILE
;##############################################################

; apPlate files contain the following:
;
;    * HDU #0 = Header only
;    * HDU #1 = Flux in units of 10^(-17) erg/s/cm^2/Ang [FLOAT]
;    * HDU #2 = Error for the above [FLOAT]
;    * HDU #3 = Mask [16-bit INT]
;    * HDU #4 = Wavelength in units of Ang [DOUBLE]
;    * HDU #5 = Average sky flux in units of 10^(-17) erg/s/cm^2/Ang [FLOAT]
;    * HDU #6 = Sky error in units of 10^(-17) erg/s/cm^2/Ang [FLOAT]
;    * HDU #7 = Telluric absorption flux in units of [FLOAT]
;    * HDU #8 = Telluric error in units of [FLOAT]
;    * HDU #9 = Wavelength solution coefficients [BINARY FITS TABLE or DOUBLE]
;    * HDU #10 = LSF coefficients [BINARY FITS TABLE or DOUBLE]
;    * HDU #11 = Plug-map structure from plPlugMapM file [BINARY FITS TABLE]
;    * HDU #12 = Plug-map header from plPlugMapM file [BINARY FITS TABLE]
;    * HDU #13 = Shift structure from ap1Dvisit
;    * HDU #14 = Pair structure from ap1dvisit
;    * HDU #15 = Counts to flux units conversion factors

; There is a separate file for each chip - [abc]


chiptag = ['a','b','c']

if ~keyword_set(single) then begin

; Loop through the three chips
For i=0,2 do begin

  ; Update the header:
  ;-------------------
  header = frame.(i).header

  ; Remove HISTORY tags from AP3D and AP2D
  remhead = where(stregex(header,'HISTORY AP3D: ',/boolean) eq 1 or $
                  stregex(header,'HISTORY AP2D: ',/boolean) eq 1,nremhead)
  if nremhead gt 0 then REMOVE,remhead,header

  ; Remove the trailing blank lines
  indend = where(stregex(header,'^END',/boolean) eq 1,nindend)
  if indend[0] eq -1 then indend=n_elements(header)-1
  header = header[0:indend[0]]

  ; Add extension explanations
  ;----------------------------
  leadstr = 'AP1DVISIT: '
  sxaddpar,header,'V_APRED',svn_vers,'apogeereduce software version'
  sxaddhist,leadstr+systime(0),header
  info = GET_LOGIN_INFO()
  sxaddhist,leadstr+info.user_name+' on '+info.machine_name,header
  ; Add Reduction pipeline version to the header
  sxaddhist,leadstr+' APOGEE Reduction Pipeline Version: '+svn_vers,header
  sxaddhist,leadstr+'IDL '+!version.release+' '+!version.os+' '+!version.arch,header
  sxaddhist,leadstr+'Output File:',header
  sxaddhist,leadstr+' HDU0 - Header only',header
  sxaddhist,leadstr+' All image extensions HDU1-HDU8 have 3 chips in 3 rows',header
  sxaddhist,leadstr+' HDU1 - Flux (10^-17 ergs/s/cm^2/Ang)',header
  sxaddhist,leadstr+' HDU2 - Error (10^-17 ergs/s/cm^2/Ang)',header
  sxaddhist,leadstr+' HDU3 - flag mask (bitwise OR combined)',header
  sxaddhist,leadstr+' HDU4 - Wavelength (Ang)',header
  sxaddhist,leadstr+' HDU5 - Sky (10^-17 ergs/s/cm^2/Ang)',header
  sxaddhist,leadstr+' HDU6 - Sky Error (10^-17 ergs/s/cm^2/Ang)',header
  sxaddhist,leadstr+' HDU7 - Telluric',header
  sxaddhist,leadstr+' HDU8 - Telluric Error',header
  sxaddhist,leadstr+' HDU9 - Wavelength coefficients',header
  sxaddhist,leadstr+' HDU10 - LSF coefficients',header
  sxaddhist,leadstr+' HDU11 - Plugmap structure',header
  sxaddhist,leadstr+' HDU12 - Plugmap header',header
  sxaddhist,leadstr+' HDU13 - Shift structure',header
  sxaddhist,leadstr+' HDU14 - Pair structure',header
  if tag_exist(frame,'FLUXCORR') then $
    sxaddhist,leadstr+' HDU15 - Counts to flux units conversion factors',header

  ; Add plate information
  sxaddpar,header,'PLATE',plate,'Plate iD'
  sxaddpar,header,'MJD5',mjd,'MJD of observation'
  sxaddpar,header,'LOCID',locid,'Location ID of field'
  if keyword_set(single) then $
    sxaddpar,header,'TELESCOP','apo1m','Telescope' else $
    sxaddpar,header,'TELESCOP','apo25m','Telescope'

  ; Create filename
  ;   apPlate-[abc]-PLATE4-MJD5.fits 
  outfile = apogee_filename('Plate',chip=chiptag[i],plate=plate,mjd=mjd)
  if not keyword_set(silent) then $
    print,'Writing Plate frame to ',outfile

  ; HDU #0 = Header only
  ;----------------------
  FITS_WRITE,outfile,0,header
  
  ; HDU #1 = Flux in units of 10^(-17) erg/s/cm^2/Ang [FLOAT]
  ;-----------------------------------------------------------
  flux = float(frame.(i).flux) / 1e-17
  bad = where(finite(flux) eq 0,nbad)
  if nbad gt 0 then flux[bad] = 0.
  MKHDR,header1,flux,/image
  sxaddpar,header1,'CTYPE1','Pixel'
  sxaddpar,header1,'CTYPE2','Fiber'
  sxaddpar,header1,'BUNIT','Flux (10^-17 ergs/s/cm^2/Ang)'
  MWRFITS,flux,outfile,header1,/silent

  ; HDU #2 = Flux Error 10^(-17) erg/s/cm^2/Ang [FLOAT]
  ;------------------------------------------------------
  err = float(frame.(i).err) 
  bderr=where(err eq baderr() or finite(err) eq 0,nbderr)
  err /= 1e-17
  if nbderr gt 0 then err[bderr] = baderr()
  MKHDR,header2,err,/image
  sxaddpar,header2,'CTYPE1','Pixel'
  sxaddpar,header2,'CTYPE2','Fiber'
  sxaddpar,header2,'BUNIT','Flux Error (10^-17 ergs/s/cm^2/Ang)'
  MWRFITS,errout(err),outfile,header2,/silent

  ; HDU #3 = Pixel mask [16-bit INT]
  ;---------------------------------
  mask = fix(frame.(i).mask)
  if nbad gt 0 then mask[bad] = mask[bad] or maskval('BADPIX')
  MKHDR,header3,mask,/image
  sxaddpar,header3,'CTYPE1','Pixel'
  sxaddpar,header3,'CTYPE2','Fiber'
  sxaddpar,header3,'BUNIT','Flag Mask (bitwise)'
  ;sxaddhist,'Explanation of BITWISE flag mask (OR combined)',header3
  ;sxaddhist,' 1 - bad pixels',header3
  ;sxaddhist,' 2 - cosmic ray',header3
  ;sxaddhist,' 4 - saturated',header3
  ;sxaddhist,' 8 - unfixable',header3
  MWRFITS,mask,outfile,header3,/silent

  ; HDU #4 = Wavelength in units of A [DOUBLE]
  ;---------------------------------------------
  wave = double(frame.(i).wavelength)
  MKHDR,header4,wave,/image
  sxaddpar,header4,'CTYPE1','Pixel'
  sxaddpar,header4,'CTYPE2','Fiber'
  sxaddpar,header4,'BUNIT','Wavelength (Ang)'
  MWRFITS,wave,outfile,header4,/silent

  ; HDU #5 = Sky flux in units of 10^(-17) erg/s/cm^2/Ang [FLOAT]
  ;---------------------------------------------------------------
  sky = float(frame.(i).sky) / 1e-17
  MKHDR,header5,sky,/image
  sxaddpar,header5,'CTYPE1','Pixel'
  sxaddpar,header5,'CTYPE2','Fiber'
  sxaddpar,header5,'BUNIT','Sky (10^-17 ergs/s/cm^2/Ang)'
  MWRFITS,sky,outfile,header5,/silent

  ; HDU #6 = Sky error in units of 10^(-17) erg/s/cm^2/Ang [FLOAT]
  ;----------------------------------------------------------------
  bderr=where(frame.(i).skyerr eq baderr(),nbd)
  skyerr = float(frame.(i).skyerr) / 1e-17
  if nbd gt 0 then skyerr[bderr] = baderr()
  MKHDR,header6,skyerr,/image
  sxaddpar,header6,'CTYPE1','Pixel'
  sxaddpar,header6,'CTYPE2','Fiber'
  sxaddpar,header6,'BUNIT','Sky Error (10^-17 ergs/s/cm^2/Ang)'
  MWRFITS,errout(skyerr),outfile,header6,/silent

  ; HDU #7 = Telluric absorption flux in units of [FLOAT]
  ;--------------------------------------------------------
  telluric = float(frame.(i).telluric)
  MKHDR,header7,telluric,/image
  sxaddpar,header7,'CTYPE1','Pixel'
  sxaddpar,header7,'CTYPE2','Fiber'
  sxaddpar,header7,'BUNIT','Telluric'
  MWRFITS,telluric,outfile,header7,/silent

  ; HDU #8 = Telluric error [FLOAT]
  ;-------------------------------------
  telerr = float(frame.(i).telluricerr)
  MKHDR,header8,telerr,/image
  sxaddpar,header8,'CTYPE1','Pixel'
  sxaddpar,header8,'CTYPE2','Fiber'
  sxaddpar,header8,'BUNIT','Telluric Error'
  MWRFITS,telerr,outfile,header8,/silent

  ; HDU #9 = Wavelength solution coefficients [DOUBLE]
  ;-----------------------------------------------------
  wcoef = double(frame.(i).wcoef)
  MKHDR,header9,wcoef,/image
  sxaddpar,header9,'CTYPE1','Fiber'
  sxaddpar,header9,'CTYPE2','Parameters'
  sxaddpar,header9,'BUNIT','Wavelength Coefficients'
  sxaddhist,'Wavelength Coefficients to be used with PIX2WAVE.PRO:',header9
  sxaddhist,' 1 Global additive pixel offset',header9
  sxaddhist,' 4 Sine Parameters',header9
  sxaddhist,' 7 Polynomial parameters (first is a zero-point offset',header9
  sxaddhist,'                     in addition to the pixel offset)',header9
  MWRFITS,wcoef,outfile,header9,/silent

  ; HDU #10 = LSF coefficients [DOUBLE]
  ;-------------------------------------
  lsfcoef = double(frame.(i).lsfcoef)
  MKHDR,header10,lsfcoef,/image
  sxaddpar,header10,'CTYPE1','Fiber'
  sxaddpar,header10,'CTYPE2','Parameters'
  sxaddpar,header10,'BUNIT','LSF Coefficients'
  sxaddhist,'LSF Coefficients to be used with LSF_GH.PRO:',header10
  sxaddhist,'  binsize  The width of a pixel in X-units.  If this is non-zero',header10
  sxaddhist,'             then a "binned" Gauss-Hermite function is used.  If',header10
  sxaddhist,'             binsize=0 then a "normal, unbinned" Gauss-Hermite',header10
  sxaddhist,'             function is used.',header10
  sxaddhist,'  X0       An additive x-offset.  This is only used to',header10
  sxaddhist,'             evaluate the GH parameters that vary globally',header10
  sxaddhist,'             with X.',header10
  sxaddhist,'  Horder   The highest Hermite order, Horder=0 means',header10
  sxaddhist,'             only a constant term (i.e. only Gaussian).',header10
  sxaddhist,'             There are Horder Hermite coefficients (since we fix H0=1).',header10
  sxaddhist,'  Porder   This array gives the polynomial order for the',header10
  sxaddhist,'             global variation (in X) of each LSF parameter.',header10
  sxaddhist,'             That includes sigma and the Horder Hermite',header10
  sxaddhist,'             coefficients (starting with H1 because we fix H0=1)',header10
  sxaddhist,'             There will be Porder[i]+1 coefficients for',header10
  sxaddhist,'             parameter i.',header10
  sxaddhist,'  GHcoefs  The polynomial coefficients for sigma and the',header10
  sxaddhist,'             Horder Hermite parameters.  There are Porder[i]+1',header10
  sxaddhist,'             coefficients for parameter i.  The Hermite parameters',header10
  sxaddhist,'             start with H1 since we fix H0=1.',header10
  MWRFITS,lsfcoef,outfile,header10,/silent

  ; HDU #11 = Plug-map structure from plPlugMapM file [BINARY FITS TABLE]
  ;----------------------------------------------------------------------
  plugdata = plugmap.fiberdata
  MWRFITS,plugdata,outfile,/silent ; first write the data with no header
                                   ; MWRFITS will add the necessary info

  ; HDU # 12 = Plug-map header values
  ;-------------------------------------
  ; remove FIBERDATA
  pltags = tag_names(plugmap)
  plind = indgen(n_elements(pltags))
  bd = where(pltags eq 'FIBERDATA',nbd)
  tmp = where(pltags eq 'GUIDEDATA',nbd)
  if nbd gt 0 then bd=[bd,tmp]
  REMOVE,bd,plind
  newplug = CREATE_STRUCT(pltags[plind[0]],plugmap.(plind[0]))
  for k=1,n_elements(plind)-1 do newplug = CREATE_STRUCT(newplug,pltags[plind[k]],plugmap.(plind[k]))
  MWRFITS,newplug,outfile,/silent

  ; HDU # 13 = Shift structure
  ;-------------------------------------
  MWRFITS,shiftstr,outfile,/silent

  ; HDU # 14 = Pair structure
  ;-------------------------------------
  MWRFITS,pairstr,outfile,/silent

  ; HDU # 15 = Flux Correction Factors
  ;-------------------------------------
  if tag_exist(frame,'FLUXCORR') then $
    MWRFITS,frame.fluxcorr,outfile,/silent

  ;; Now modify the header
  ;header11 = HEADFITS(outfile,exten=11)
  ;pltags = tag_names(plugmap)
  ;; Add plate information to output header
  ;for j=0,n_elements(pltags)-1 do begin
  ;  sz = size(plugmap.(j))
  ;  type = size(plugmap.(j),/type)
  ;  if pltags[j] ne 'HDR' and pltags[j] ne 'FIBERDATA' and sz[0] lt 2 then begin
  ;    ;sxaddpar,header11,pltags[j],plugmap.(j),' PLUGMAPTAG'
  ;    line = 'PLMAPTG: '+strtrim(pltags[j],2)+'='
  ;    if type eq 7 then line=line+"'"+strtrim(plugmap.(j),2)+"'" else $
  ;      line=line+strtrim(plugmap.(j),2)
  ;    sxaddhist,line,header11
  ;  endif
  ;end
  ;; Add original plPlugMap header as history lines
  ;for j=0,n_elements(plugmap.hdr)-1 do begin
  ;  if strtrim(plugmap.hdr[j],2) ne '' then $
  ;    sxaddhist,'PLMAPHD: '+plugmap.hdr[j],header11
  ;  end
  ;; Now modify the header in the file
  ;MODFITS,outfile,0,header11,exten_no=11

  ;stop

endfor

endif


;##############################################################
;  OUTPUT THE INDIVIDUAL VISIT SPECTRA
;##############################################################

; Fiber Loop
For i=0,nfibers-1 do begin

  ifiber = i
  ifiberid = 300-ifiber

  ; The plugmap index for this fiber
  ; fiberid=1 is at the top of the detector or index=299
  ; index = 300-fiberid
  ;iplugind = where(plugmap.fiberdata.fiberid eq ifiber+1,niplugind)
  iplugind = where(plugmap.fiberdata.spectrographid eq 2 and $
                   plugmap.fiberdata.fiberid eq ifiberid,niplugind)

  ; No information for this fiber
  if niplugind eq 0 then begin
    print,'No information for Fiber=',strtrim(ifiberid,2),' in the plugmap file'
    goto,BOMB
  endif

  ; Only output object and telluric spectra
  ;  type: 0=target, 1=sky, 2=telluric
  ;fibertype = plugmap.fiberdata[i].type
  ;if fibertype eq 0 or fibertype eq 2 then begin
  fiberholetype = plugmap.fiberdata[iplugind].holetype
  fiberobjtype = plugmap.fiberdata[iplugind].objtype
  ;fiberobjid = plugmap.fiberdata[iplugind].objid
  fiber_ra = plugmap.fiberdata[iplugind].ra
  fiber_dec = plugmap.fiberdata[iplugind].dec
  fiber_mag = plugmap.fiberdata[iplugind].mag
  tmass_name = plugmap.fiberdata[iplugind].tmass_style
  targ1 = plugmap.fiberdata[iplugind].target1
  targ2 = plugmap.fiberdata[iplugind].target2
  targ3 = plugmap.fiberdata[iplugind].target3
  ak_targ = plugmap.fiberdata[iplugind].ak_targ
  ak_targ_method = plugmap.fiberdata[iplugind].ak_targ_method
  ak_wise = plugmap.fiberdata[iplugind].ak_wise
  sfd_ebv = plugmap.fiberdata[iplugind].sfd_ebv
  ; neighboring fibers
  iplugind1 = where(plugmap.fiberdata.spectrographid eq 2 and $
                   plugmap.fiberdata.fiberid eq ifiberid+1,niplugind)
  if niplugind gt 0 then begin
    if plugmap.fiberdata[iplugind1].mag[1] lt 0.01 then hplus=99.99 else $
       hplus = plugmap.fiberdata[iplugind1].mag[1]-fiber_mag[1] 
  endif else hplus=99.99

  iplugind1 = where(plugmap.fiberdata.spectrographid eq 2 and $
                   plugmap.fiberdata.fiberid eq ifiberid-1,niplugind)
  if niplugind gt 0 then begin
    if plugmap.fiberdata[iplugind1].mag[1] lt 0.01 then hminus=99.99 else $
       hminus = plugmap.fiberdata[iplugind1].mag[1]-fiber_mag[1] 
  endif else hminus=99.99

  if fiberholetype eq 'OBJECT' and (fiberobjtype ne 'SKY') and $
     tmass_name ne '-' then begin

    ; Each extension contains a different data type and is a 2D
    ; array of length [3,4096], although the wavelength coefficients
    ; and lsf coefficents will have different dimensions

    ;  Each file contains the following HDUs:
    ;
    ; * HDU 0 = Header only
    ; * HDU 1 = Flux: each chip in a separate row (i.e. 4096x3)
    ; * HDU 2 = Error
    ; * HDU 3 = Mask [16-bit INT]
    ; * HDU 4 = Wavelength in units of Ang [DOUBLE]
    ; * HDU 5 = Sky spectrum for each chip
    ; * HDU 6 = Sky error
    ; * HDU 7 = Telluric correction for each chip
    ; * HDU 8 = Telluric error
    ; * HDU 9 = Wavelength coefficients
    ; * HDU 10 = LSF coefficients
    ; * HDU 11 = ADU to flux conversion factor array

    ; Create filename
    ;   apVisit-VERS-PLATE4-MJD5-FIBER3.fits 
    outfile = apogee_filename('Visit',plate=plate,mjd=mjd,fiber=ifiberid,reduction=strtrim(tmass_name,2))
    if keyword_set(mjdfrac) then begin
      cmjd=strtrim(string(mjd),2)
      s=strsplit(outfile,cmjd,/extract,/regex)
      outfile=s[0]+cmjd+s[1]+string(format='(f8.2)',mjdfrac)+s[2]
    endif
    if not keyword_set(silent) then $
      print,'Writing Individual Visit Spectrum to ',outfile

    ; Create main header:
    ;---------------------
    MKHDR,header,0

    ; Add plate/observation information
    sxaddpar,header,'LOCID',locid,'Location ID of field'
    sxaddpar,header,'PLATE',plate,'Plate ID'
    if keyword_set(single) then $
      sxaddpar,header,'TELESCOP','apo1m','Telescope' else $
      sxaddpar,header,'TELESCOP','apo25m','Telescope'
    sxaddpar,header,'MJD5',mjd,'MJD of observation'
    sxaddpar,header,'FIBERID',ifiberid,' APOGEE Fiber ID 1-300'
    sxaddpar,header,'DATE-OBS',sxpar(frame.(0).header,'DATE-OBS')
    ; add JD, HJD ??
    sxaddpar,header,'EXPTIME',sxpar(frame.(0).header,'EXPTIME'),'Total visit exptime per dither pos'
    sxaddpar,header,'JD-MID',sxpar(frame.(0).header,'JD-MID'),' JD at midpoint of visit'
    sxaddpar,header,'UT-MID',sxpar(frame.(0).header,'UT-MID'),' Date at midpoint of visit'
    ncombine = sxpar(frame.(0).header,'NCOMBINE',count=num_ncombine)
    if num_ncombine eq 0 then ncombine=1
    sxaddpar,header,'NCOMBINE',ncombine
    for j=0,ncombine-1 do sxaddpar,header,'FRAME'+strtrim(j+1,2),sxpar(frame.(0).header,'FRAME'+strtrim(j+1,2)),'Constituent frame'
    sxaddpar,header,'NPAIRS',sxpar(frame.(0).header,'NPAIRS'),' Number of dither pairs combined'

    ; Add star information
    sxaddpar,header,'OBJID',strtrim(tmass_name,2),' Object ID'
    sxaddpar,header,'OBJTYPE',fiberobjtype,' Object type'
    ; get catalog information downloaded from database
    ; with plateHoles info, don't bother with this
    ;cat=getcat(strtrim(tmass_name,2),ra=fiber_ra)
    ;if size(cat,/type) eq 8 then begin
    ;  sxaddpar,header,'RA',cat.ra,' catalog right ascension, deg, J2000'
    ;  sxaddpar,header,'DEC',cat.dec,' catalog declination, deg, J2000'
    ;  sxaddpar,header,'RA_TARG',fiber_ra,' targeting right ascension, deg, J2000'
    ;  sxaddpar,header,'DEC_TARG',fiber_dec,' targeting declination, deg, J2000'
    ;  sxaddpar,header,'AKTARG',cat.ak,' Extinction used for targeting'
    ;  sxaddpar,header,'AKMETHOD',cat.method,' Extinction method used for targeting'
    ;  sxaddpar,header,'AKWISE',cat.ak_wise,' WISE all-sky extinction'
    ;  sxaddpar,header,'SFD_EBV',cat.sfd_ebv,' SFD E(B-V)'
    ;  sxaddpar,header,'J',cat.j>0.?cat.j:99.999,' 2MASS J magnitude'
    ;  sxaddpar,header,'J_ERR',cat.j_err,' 2MASS J magnitude uncertainty'
    ;  sxaddpar,header,'H',cat.h>0.?cat.h:99.999,' 2MASS H magnitude'
    ;  sxaddpar,header,'H_ERR',cat.h_err,' 2MASS H magnitude uncertainty'
    ;  sxaddpar,header,'K',cat.ks>0.?cat.ks:99.999,' 2MASS Ks magnitude'
    ;  sxaddpar,header,'K_ERR',cat.ks_err,' 2MASS Ks magnitude uncertainty'
    ;endif else begin
      print,'Object not found in catalog file. You can continue, but no extinction in header'
      sxaddpar,header,'AKMETHOD','NONE',' Extinction method'
      sxaddpar,header,'RA',fiber_ra,' targeting right ascension, deg, J2000'
      sxaddpar,header,'DEC',fiber_dec,' targeting declination, deg, J2000'
      sxaddpar,header,'J',fiber_mag[0]>0.?fiber_mag[0]:99.999,' 2MASS J magnitude'
      sxaddpar,header,'H',fiber_mag[1]>0.?fiber_mag[1]:99.999,' 2MASS H magnitude'
      sxaddpar,header,'K',fiber_mag[2]>0.?fiber_mag[2]:99.999,' 2MASS Ks magnitude'
    ;endelse
    snr=median(frame.(1).flux[*,ifiber]/frame.(1).err[*,ifiber])
    sxaddpar,header,'SNR',snr,' median S/N, middle chip'
    if keyword_set(survey) then sxaddpar,header,'SURVEY',survey,' Survey definition (for targeting flags)'
    sxaddpar,header,'TARG1',targ1,' First APOGEE targeting flag (bitwise, see docs)'
    sxaddpar,header,'TARG2',targ2,' Second APOGEE targeting flag (bitwise, see docs)'
    sxaddpar,header,'TARG3',targ3,' Third APOGEE targeting flag (bitwise, see docs)'
    sxaddpar,header,'AKTARG',ak_targ,'Extinction used in targeting'
    sxaddpar,header,'AKMETHOD',ak_targ_method,'Extinction method'
    sxaddpar,header,'AKWISE',ak_wise,'WISE extinction'
    sxaddpar,header,'SFD_EBV',sfd_ebv,'SFD E(B-V)'
    sxaddpar,header,'HPLUS',hplus,' Delta H (neigh-obj) mag of neighboring (+1) fiber'
    sxaddpar,header,'HMINUS',hminus,' Delta H (neigh-obj) mag of neighboring (-1) fiber'

    ; Start a FLAG for this object
    flag=0L
    ; is this low S/N?
    if snr lt 5 then flag=flag or starflagval('LOW_SNR')

    ; is star from comissioning data?
    if mjd le 55761L then flag = flag OR starflagval('COMMISSIONING')

    ; does star have a bright neighbor?
    if hplus lt -5 or hminus lt -5 then flag=flag or starflagval('VERY_BRIGHT_NEIGHBOR') else $
    if hplus lt -2.5 or hminus lt -2.5 then flag=flag or starflagval('BRIGHT_NEIGHBOR')

    ; get the pixel mask
    mask = intarr(npix,3)
    flux = fltarr(npix,3)
    for j=0,2 do mask[*,j]=fix(frame.(j).mask[*,ifiber])
    for j=0,2 do flux[*,j]=float(frame.(j).flux[*,ifiber])

    ; does star have a significant number of persistence pixels?
    junk=where(mask and maskval('PERSIST_HIGH'),nbad)
    if float(nbad)/n_elements(mask) gt 0.2 then flag = flag OR starflagval('PERSIST_HIGH') else begin
      junk=where((mask and maskval('PERSIST_HIGH')) or (mask and maskval('PERSIST_MED')),nbad)
      if float(nbad)/n_elements(mask) gt 0.2 then flag = flag OR starflagval('PERSIST_MED') else begin
        junk=where((mask and maskval('PERSIST_HIGH')) or (mask and maskval('PERSIST_MED')) or $
                   (mask and maskval('PERSIST_LOW')),nbad)
        if float(nbad)/n_elements(mask) gt 0.2 then flag = flag OR starflagval('PERSIST_LOW')
      endelse
    endelse

    ; does star show evidence of significant persistence?
    bmed=median(frame.(2).flux[*,ifiber])
    gmed=median(frame.(1).flux[*,ifiber])
    if bmed gt 1.5*gmed*(1.55/1.62)^(-4) then flag = flag OR starflagval('PERSIST_JUMP_POS')
    if bmed lt 0.667*gmed*(1.55/1.62)^(-4) then flag = flag OR starflagval('PERSIST_JUMP_NEG')
    ; if persistence in evident, flag persistence pixels as bad 
    if (flag AND  starflagval('PERSIST_JUMP_POS')) gt 0 then begin
      junk=where((mask and maskval('PERSIST_LOW') gt 0) or (mask and maskval('PERSIST_MED') gt 0) or (mask and maskval('PERSIST_HIGH') gt 0),nbad)
      if nbad gt 0 then mask[junk] = mask[junk] or maskval('BADPIX')
    endif

    ; any unaccounted for NaNs? make sure BADPIX is set!
    junk=where(finite(flux) eq 0,nbad)
    if nbad gt 0 then mask[junk] = mask[junk] or maskval('BADPIX')

    ; does star have a significant number of bad pixels?
    junk=where(mask and badmask(),nbad)
    if float(nbad)/n_elements(mask) gt 0.2 then flag = flag OR starflagval('BAD_PIXELS')

    ; add to header
    sxaddpar,header,'STARFLAG',flag,' Star data quality flag'

    ; Add flux correction factor
    if tag_exist(frame,'FLUXCORR') then begin
      szflux=size(frame.fluxcorr)
      if szflux[0] eq 1 then $
        sxaddpar,header,'FLUXFLAM',frame.fluxcorr[ifiber],' ADU to flux units conv factor (ergs/s/cm^2/A)' $
      else $
        sxaddpar,header,'FLUXFLAM',median(frame.fluxcorr[*,*,ifiber]),' ADU to flux units conv factor (ergs/s/cm^2/A)'
    endif

    ; Description of the extensions
    leadstr = 'AP1DVISIT: '
    sxaddpar,header,'V_APRED',svn_vers,'apogeereduce software version'
    sxaddhist,leadstr+systime(0),header
    info = GET_LOGIN_INFO()
    sxaddhist,leadstr+info.user_name+' on '+info.machine_name,header
    sxaddhist,leadstr+' APOGEE Reduction Pipeline Version: '+svn_vers,header
    sxaddhist,leadstr+'IDL '+!version.release+' '+!version.os+' '+!version.arch,header
    sxaddhist,leadstr+'The data for the individual visit spectrum are in',header
    sxaddhist,leadstr+'  separate extensions.  All HDUs have data for each chip in',header
    sxaddhist,leadstr+'  a separate row (i.e. 4096x3) except for Wavelength and LSF coefficients',header
    sxaddhist,leadstr+' HDU0 = Header only',header
    sxaddhist,leadstr+' HDU1 - Flux (10^-17 ergs/s/cm^2/Ang)',header
    sxaddhist,leadstr+' HDU2 - Error (10^-17 ergs/s/cm^2/Ang)',header
    sxaddhist,leadstr+' HDU3 - Flag mask (bitwise OR combined)',header
    sxaddhist,leadstr+' HDU4 - Wavelength (Ang)',header
    sxaddhist,leadstr+' HDU5 - Sky (10^-17 ergs/s/cm^2/Ang)',header
    sxaddhist,leadstr+' HDU6 - Sky Error (10^-17 ergs/s/cm^2/Ang)',header
    sxaddhist,leadstr+' HDU7 - Telluric',header
    sxaddhist,leadstr+' HDU8 - Telluric Error',header
    sxaddhist,leadstr+' HDU9 - Wavelength coefficients',header
    sxaddhist,leadstr+' HDU10 - LSF coefficients',header
    if tag_exist(frame,'FLUXCORR') and szflux[0] gt 1 then $
      sxaddhist,leadstr+' HDU11 - ADU to flux conversion factor array',header

    ; Do we need to add information from the plugmap structure??

    ; HDU #0 = Header only
    ;----------------------
    FITS_WRITE,outfile,0,header

    ; HDU #1 = Flux in units of 10^(-17) erg/s/cm^2/Ang [FLOAT]
    ;------------------------------------------------------------
    flux = fltarr(npix,3)
    for j=0,2 do flux[*,j]=float(frame.(j).flux[*,ifiber])
    bad = where(finite(flux) eq 0,nbad)
    if nbad gt 0 then flux[bad] = 0.
    flux /= 1e-17
    MKHDR,header1,flux,/image
    sxaddpar,header1,'CTYPE1','Pixel'
    sxaddpar,header1,'CTYPE2','Chip'
    sxaddpar,header1,'BUNIT','Flux (10^-17 erg/s/cm^2/Ang)'
    MWRFITS,flux,outfile,header1,/silent

    ; HDU #2 = Error in units of 10^(-17) erg/s/cm^2/Ang [FLOAT]
    ;------------------------------------------------------------
    err = fltarr(npix,3)
    for j=0,2 do err[*,j]=float(frame.(j).err[*,ifiber])
    bderr=where(err eq baderr() or finite(err) eq 0,nbd)
    err /= 1e-17
    if nbd gt 0 then err[bderr] = baderr()
    MKHDR,header2,err,/image
    sxaddpar,header2,'CTYPE1','Pixel'
    sxaddpar,header2,'CTYPE2','Chip'
    sxaddpar,header2,'BUNIT','Flux Error (10^-17 erg/s/cm^2/Ang)'
    MWRFITS,errout(err),outfile,header2,/silent

    ; HDU #3 = Pixel mask [16-bit INT]
    ;---------------------------------
    MKHDR,header3,mask,/image
    sxaddpar,header3,'CTYPE1','Pixel'
    sxaddpar,header3,'CTYPE2','Chip'
    sxaddpar,header3,'BUNIT','Flag Mask (bitwise)'
    ;sxaddhist,'Explanation of BITWISE flag mask (OR combined)',header3
    ;sxaddhist,' 1 - bad pixels',header3
    ;sxaddhist,' 2 - cosmic ray',header3
    ;sxaddhist,' 4 - saturated',header3
    ;sxaddhist,' 8 - unfixable',header3
    MWRFITS,mask,outfile,header3,/silent

    ; HDU #4 = Wavelength in units of A [DOUBLE]
    ;---------------------------------------------
    wave = dblarr(npix,3)
    for j=0,2 do wave[*,j]=double(frame.(j).wavelength[*,ifiber])
    MKHDR,header4,wave,/image
    sxaddpar,header4,'CTYPE1','Pixel'
    sxaddpar,header4,'CTYPE2','Chip'
    sxaddpar,header4,'BUNIT','Wavelength (Ang)'
    MWRFITS,wave,outfile,header4,/silent

    ; HDU #5 = Sky flux in units of 10^(-17) erg/s/cm^2/Ang [FLOAT]
    ;---------------------------------------------------------------
    sky = fltarr(npix,3)
    for j=0,2 do sky[*,j]=float(frame.(j).sky[*,ifiber])
    sky /= 1e-17
    MKHDR,header5,sky,/image
    sxaddpar,header5,'CTYPE1','Pixel'
    sxaddpar,header5,'CTYPE2','Chip'
    sxaddpar,header5,'BUNIT','Sky (10^-17 erg/s/cm^2/Ang)'
    MWRFITS,sky,outfile,header5,/silent

    ; HDU #6 = Sky error in units of 10^(-17) erg/s/cm^2/Ang [FLOAT]
    ;----------------------------------------------------------------
    skyerr = fltarr(npix,3)
    for j=0,2 do skyerr[*,j]=float(frame.(j).skyerr[*,ifiber])
    skyerr /= 1e-17
    MKHDR,header6,skyerr,/image
    sxaddpar,header6,'CTYPE1','Pixel'
    sxaddpar,header6,'CTYPE2','Chip'
    sxaddpar,header6,'BUNIT','Sky Error (10^-17 erg/s/cm^2/Ang)'
    MWRFITS,skyerr,outfile,header6,/silent

    ; HDU #7 = Telluric absorption flux in units of [FLOAT]
    ;--------------------------------------------------------
    telluric = fltarr(npix,3)
    for j=0,2 do telluric[*,j]=float(frame.(j).telluric[*,ifiber])
    MKHDR,header7,telluric,/image
    sxaddpar,header7,'CTYPE1','Pixel'
    sxaddpar,header7,'CTYPE2','Chip'
    sxaddpar,header7,'BUNIT','Telluric'
    MWRFITS,telluric,outfile,header7,/silent

    ; HDU #8 = Telluric error [FLOAT]
    ;-------------------------------------
    telerr = fltarr(npix,3)
    for j=0,2 do telerr[*,j]=float(frame.(j).telluricerr[*,ifiber])
    MKHDR,header8,telerr,/image
    sxaddpar,header8,'CTYPE1','Pixel'
    sxaddpar,header8,'CTYPE2','Chip'
    sxaddpar,header8,'BUNIT','Telluric Error'
    MWRFITS,telerr,outfile,header8,/silent

    ; HDU #9 = Wavelength solution coefficients [DOUBLE]
    ;-----------------------------------------------------
    nwcoef = n_elements(frame.(0).wcoef[0,*])
    wcoef = dblarr(nwcoef,3)
    for j=0,2 do wcoef[*,j]=double(frame.(j).wcoef[ifiber,*])
    MKHDR,header9,wcoef,/image
    sxaddpar,header9,'CTYPE1','Fiber'
    sxaddpar,header9,'CTYPE2','Parameters'
    sxaddpar,header9,'BUNIT','Wavelength Coefficients'
    sxaddhist,'Wavelength Coefficients to be used with PIX2WAVE.PRO:',header9
    sxaddhist,' 1 Global additive pixel offset',header9
    sxaddhist,' 4 Sine Parameters',header9
    sxaddhist,' 7 Polynomial parameters (first is a zero-point offset',header9
    sxaddhist,'                     in addition to the pixel offset)',header9
    MWRFITS,wcoef,outfile,header9,/silent

    ; HDU #10 = LSF coefficients [DOUBLE]
    ;-------------------------------------
    nlsfcoef = n_elements(frame.(0).lsfcoef[0,*])
    lsfcoef = dblarr(nlsfcoef,3)
    for j=0,2 do lsfcoef[*,j]=double(frame.(j).lsfcoef[ifiber,*])
    MKHDR,header10,lsfcoef,/image
    sxaddpar,header10,'CTYPE1','Fiber'
    sxaddpar,header10,'CTYPE2','Parameters'
    sxaddpar,header10,'BUNIT','LSF Coefficients'
    sxaddhist,'LSF Coefficients to be used with LSF_GH.PRO:',header10
    sxaddhist,'  binsize  The width of a pixel in X-units.  If this is non-zero',header10
    sxaddhist,'             then a "binned" Gauss-Hermite function is used.  If',header10
    sxaddhist,'             binsize=0 then a "normal, unbinned" Gauss-Hermite',header10
    sxaddhist,'             function is used.',header10
    sxaddhist,'  X0       An additive x-offset.  This is only used to',header10
    sxaddhist,'             evaluate the GH parameters that vary globally',header10
    sxaddhist,'             with X.',header10
    sxaddhist,'  Horder   The highest Hermite order, Horder=0 means',header10
    sxaddhist,'             only a constant term (i.e. only Gaussian).',header10
    sxaddhist,'             There are Horder Hermite coefficients (since we fix H0=1).',header10
    sxaddhist,'  Porder   This array gives the polynomial order for the',header10
    sxaddhist,'             global variation (in X) of each LSF parameter.',header10
    sxaddhist,'             That includes sigma and the Horder Hermite',header10
    sxaddhist,'             coefficients (starting with H1 because we fix H0=1)',header10
    sxaddhist,'             There will be Porder[i]+1 coefficients for',header10
    sxaddhist,'             parameter i.',header10
    sxaddhist,'  GHcoefs  The polynomial coefficients for sigma and the',header10
    sxaddhist,'             Horder Hermite parameters.  There are Porder[i]+1',header10
    sxaddhist,'             coefficients for parameter i.  The Hermite parameters',header10
    sxaddhist,'             start with H1 since we fix H0=1.',header10
    MWRFITS,lsfcoef,outfile,header10,/silent

    ; HDU #11 = Flux conversion factor
    ;-------------------------------------
    if tag_exist(frame,'FLUXCORR') and szflux[0] gt 1 then begin
      fluxcorr = frame.fluxcorr[*,*,ifiber]
      MKHDR,header11,fluxcorr,/image
      sxaddpar,header11,'CTYPE1','Pixel'
      sxaddpar,header11,'CTYPE2','Chip'
      sxaddpar,header11,'BUNIT','ADU to flux units conv factor (ergs/s/cm^2/A)'
      MWRFITS,fluxcorr,outfile,header11,/silent
    endif

    ;stop
    ;if keyword_set(single) then begin
    ;  ; save the first iteration visit spectrum if iterating for telluric
    ;  if iter eq 1 then begin
    ;   outfile1 = outdir+'apVisit-'+apred_vers+'-'+strtrim(mjd,2)+'-'+strtrim(tmass_name,2)+'_1.fits' 
    ;   spawn,'cp '+outfile+' '+outfile1
    ;  endif
    ;endif

  endif  ; object spectrum

  BOMB:

End

if keyword_set(stp) then stop

end
