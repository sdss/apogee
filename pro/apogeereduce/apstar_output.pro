pro apstar_output,starstr,visitstr,binstr,stars_dir,silent=silent,stp=stp,nolsf=nolsf,commiss=commiss,localdir=localdir,apstar_vers=apstar_vers,locationdir=locationdir,survey=survey

;+
;
; APSTAR_OUTPUT
;
; This outputs the final visit-combined object spectrum  and LSF array.
;
; INPUTS:
;  starstr     The star structure that contains the final
;               visit-combined spectrum.  This is created by
;               apvisitcomb.pro
;  visitstr    Structure with information for each visit
;  binstr      Binary check structure
;  stars_dir   The root output directory for the apStar and apLSF files.
;  /silent     Don't print anything to the screen.
;  /stp        Stop at the end of the program.
;  /commiss    Commissioning data, use apStarC instead of apStar
;
; OUTPUTS:
;  The apStar and apLSF files are written to the "star"
;  directory.
;
; USAGE:
;  IDL>apstar_output,starstr,visitstr,stars_dir
;
; By D.Nidever  July 2010
;-

cspeed = 2.99792458d5  ; speed of light in km/s

nstarstr = n_elements(starstr)
nvisitstr = n_elements(visitstr)
nstars_dir = n_elements(stars_dir)

; Not enough inputs
if nstarstr eq 0 or nvisitstr eq 0 or nstars_dir eq 0 then begin
  print,'Syntax - apstar_output,starstr,visitstr,stars_dir,silent=silent,stp=stp'
  return
endif

npix = n_elements(starstr.spec)
nvisits = starstr.nvisits

if ~keyword_set(survey) then survey='apogee'

;##############################################################
;  OUTPUT THE STAR FILE
;##############################################################

; apStar files contain the following:
;
;    * HDU #0 = Header only
;    * HDU #1 = Flux in units of 10^(-17) erg/s/cm^2/Ang [FLOAT]
;    * HDU #2 = Error in units of 10^(-17) erg/s/cm^2/Ang [FLOAT]
;    * HDU #3 = Flag mask [16-bit INT]
;    * HDU #4 = Average sky flux in units of 10^(-17) erg/s/cm^2/Ang [FLOAT]
;    * HDU #5 = Sky error in units of 10^(-17) erg/s/cm^2/Ang [FLOAT]
;    * HDU #6 = Telluric absorption [FLOAT]
;    * HDU #7 = Telluric error [FLOAT]
;    * HDU #8 = LSF coefficients [DOUBLE]


; Create the header
;--------------------
MKHDR,header,0

; Get some values from the 1st visit header
vhead1 = reform(starstr.header[0,*])
objid = sxpar(vhead1,'OBJID')
; some of the mags are 0, so get the max mag
jmag = max(visitstr.j)
hmag = max(visitstr.h)
kmag = max(visitstr.k)
jmagerr = max(visitstr.j_err)
hmagerr = max(visitstr.h_err)
kmagerr = max(visitstr.k_err)
locid = sxpar(vhead1,'LOCID')
if locid eq 0 then locid = sxpar(vhead1,'LOCATION')
plate = sxpar(vhead1,'PLATE')
targ1 = sxpar(vhead1,'TARG1')
targ2 = sxpar(vhead1,'TARG2')
targ3 = sxpar(vhead1,'TARG3')
ak_targ = sxpar(vhead1,'aktarg')
ak_targ_method = sxpar(vhead1,'akmethod')
ak_wise = sxpar(vhead1,'akwise')
sfd_ebv = sxpar(vhead1, 'sfd_ebv')
targ1=visitstr[0].apogee_target1
if nvisits gt 1 then for ivisit=1,nvisits-1 do targ1=targ1 OR visitstr[ivisit].apogee_target1
targ2=visitstr[0].apogee_target2
if nvisits gt 1 then for ivisit=1,nvisits-1 do targ2=targ2 OR visitstr[ivisit].apogee_target2
targ3=visitstr[0].apogee_target3
if nvisits gt 1 then for ivisit=1,nvisits-1 do targ3=targ3 OR visitstr[ivisit].apogee_target3

; Obj ID
sxaddpar,header,'OBJID',objid,' Object ID'
; magnitudes
sxaddpar,header,'J',jmag,' 2MASS J magnitude'
if jmagerr gt 0 then sxaddpar,header,'J_ERR',jmagerr,' 2MASS J magnitude uncertainty'
sxaddpar,header,'H',hmag,' 2MASS H magnitude'
if jmagerr gt 0 then sxaddpar,header,'H_ERR',hmagerr,' 2MASS H magnitude uncertainty'
sxaddpar,header,'K',kmag,' 2MASS Ks magnitude'
if jmagerr gt 0 then sxaddpar,header,'K_ERR',kmagerr,' 2MASS Ks magnitude uncertainty'

; ADDED THIS SO WE CAN GET TARGETTING INFORMATION INTO apStar file.
sxaddpar,header,'AKTARG',ak_targ,' Extinction used for targeting'
sxaddpar,header,'AKMETHOD',ak_targ_method,' Extinction method using for targeting'
sxaddpar,header,'AKWISE',ak_wise,' WISE all-sky extinction'
sxaddpar,header,'SFD_EBV',sfd_ebv,' SFD E(B-V)' 

sxaddpar,header,'TARG1',targ1,' First APOGEE targeting flag (bitwise, see docs)'
sxaddpar,header,'TARG2',targ2,' Second APOGEE targeting flag (bitwise, see docs)'
sxaddpar,header,'TARG3',targ3,' Third APOGEE targeting flag (bitwise, see docs)'
sxaddpar,header,'SURVEY',survey,' Survey name (for targeting flags)'
sxaddpar,header,'TELESCOP',visitstr[0].telescope,' Telescope'
sxaddpar,header,'LOCID',locid,' LocationID'
sxaddpar,header,'FIELD',apogee_field(locid,plate),' Field name'
; get catalogdb information
;allobj=getcat(objid,ra=visitstr[0].ra)
;if size(allobj,/type) eq 8 then begin
;  sxaddpar,header,'AKTARG',allobj.ak,' Extinction used for targeting'
;  sxaddpar,header,'AKMETHOD',allobj.method,' Extinction method using for targeting'
;  sxaddpar,header,'AKWISE',allobj.ak_wise,' WISE all-sky extinction'
;  sxaddpar,header,'SFD_EBV',allobj.sfd_ebv,' SFD E(B-V)' 
;  sxaddpar,header,'J',allobj.j,' 2MASS J magnitude'
;  sxaddpar,header,'J_ERR',allobj.j_err,' 2MASS J magnitude uncertainty'
;  sxaddpar,header,'H',allobj.h,' 2MASS H magnitude'
;  sxaddpar,header,'H_ERR',allobj.h_err,' 2MASS H magnitude uncertainty'
;  sxaddpar,header,'K',allobj.ks,' 2MASS Ks magnitude'
;  sxaddpar,header,'K_ERR',allobj.ks_err,' 2MASS Ks magnitude uncertainty'
;endif else begin
;  print,'Object not found in catalog file. You can continue, but no extinction in header'
;  sxaddpar,header,'AKMETHOD','NONE',' Extinction method'
;  ;stop
;endelse
; Coordinates
sxaddpar,header,'RA',visitstr[0].ra,' right ascension, deg, J2000'
sxaddpar,header,'DEC',visitstr[0].dec,' declination, deg, J2000'
sxaddpar,header,'GLON',visitstr[0].glon,' galactic longitude, deg, J2000'
sxaddpar,header,'GLAT',visitstr[0].glat,' galactic latitude, deg, J2000'

sxaddpar,header,'NVISITS',nvisits,' number of visit spectra combined'

sxaddpar,header,'COMBTYPE', starstr.combtype, ' velocities used in final combination. (0=VREL, 1=SYNTHVREL)'

; Weighting method
;wtind = where(stregex(vhead1,'APVISITCOMB: ',/boolean) eq 1,nwtind)
;if nwtind gt 0 then begin
;  line = vhead1[wtind[0]]
;  pos1 = strpos(line,'APVISITCOMB: ')+13
;  pos2 = strpos(line,'/')
;  len = pos2-pos1+1
;  if pos2 eq -1 then len=100
;  com = strmid(line,pos1,len)
;  sxaddhist,'Visit spectra combining '+com,header
;endif

; Final RVs
;-----------
;  take a weighted mean (with snr)
gdv = where(visitstr.vrelerr lt 1000,ngdv)
if ngdv gt 0 then begin
  vhelio = TOTAL(visitstr[gdv].vhelio*visitstr[gdv].snr)/TOTAL(visitstr[gdv].snr)
  verr = sqrt(TOTAL(visitstr[gdv].vrelerr^2*visitstr.snr^2)/TOTAL(visitstr[gdv].snr)^2)  ; weighted error
  verr_med = median(reform(visitstr[gdv].vrelerr),/even)
  if ngdv gt 1 then vscatter = stddev(visitstr[gdv].vhelio) else vscatter=0.0d0
  obsvhelio = TOTAL(visitstr[gdv].obsvhelio*visitstr[gdv].snr)/TOTAL(visitstr[gdv].snr)
  obsverr = sqrt(TOTAL(visitstr[gdv].obsvrelerr^2*visitstr.snr^2)/TOTAL(visitstr[gdv].snr)^2)  ; weighted error
  obsverr_med = median(reform(visitstr[gdv].obsvrelerr),/even)
  if ngdv gt 1 then obsvscatter = stddev(visitstr[gdv].obsvhelio) else obsvscatter=0.0d0
  synthvhelio = TOTAL(visitstr[gdv].synthvhelio*visitstr[gdv].snr)/TOTAL(visitstr[gdv].snr)
  synthverr = sqrt(TOTAL(visitstr[gdv].synthvrelerr^2*visitstr.snr^2)/TOTAL(visitstr[gdv].snr)^2)  ; weighted error
  synthverr_med = median(reform(visitstr[gdv].synthvrelerr),/even)
  if ngdv gt 1 then synthvscatter = stddev(visitstr[gdv].synthvhelio,/nan) else synthvscatter=0.0d0
  if ngdv gt 1 then synthscatter = stddev(visitstr[gdv].obsvhelio-visitstr[gdv].synthvhelio,/nan) else synthscatter=0.0d0
  VCONV,vhelio,visitstr[0].glon,visitstr[0].glat,vlsr,vgsr
endif else begin
  vhelio = 999999.0
  verr = 999999.0
  verr_med = 999999.0
  vscatter = -1.0
  obsvhelio = 999999.0
  obsverr = 999999.0
  obsverr_med = 999999.0
  obsvscatter = -1.0
  synthvhelio = 999999.0
  synthverr = 999999.0
  synthverr_med = 999999.0
  synthvscatter = -1.0
  synthscatter = -1.0
  vlsr = 999999.0
  vgsr = 999999.0
endelse

; STARFLAGS determine OR and AND flags, and add new bits
starflag=0L
andflag=visitstr[0].starflag
for i=0,nvisits-1 do begin
  starflag=starflag or visitstr[i].starflag
  andflag=andflag and visitstr[i].starflag
endfor
if ((synthscatter gt max([2.,2*verr_med])) OR (synthscatter LT 0)) then starflag = starflag or starflagval('SUSPECT_RV_COMBINATION') ;Recommend  star_warn set in aspcapflag if this is set
if ((synthscatter gt max([10.,10*verr_med])) OR (synthscatter LT 0)) then starflag = starflag or starflagval('BAD_RV_COMBINATION') ;Recomemnd star_bad set  in aspcapflag if this is set
if starstr.ccpfwhm gt 1.5*starstr.autofwhm then starflag = starflag or starflagval('SUSPECT_BROAD_LINES')

print,'Mean Vhelio=',stringize(vhelio,ndec=3),'+/-',stringize(verr,ndec=3),' km/s'
print,'Scatter = ',stringize(vscatter,ndec=3),' km/s'
sxaddhist,'MEAN VALUES:',header,/comment
sxaddpar,header,'VRAD',0.0,' Doppler shift (km/s) of this spectrum'  ; set to rest frame
sxaddpar,header,'VHELIO',vhelio,' mean barycentric velocity (km/s)'
sxaddpar,header,'VSCATTER',vscatter,' STDEV of VHELIO (km/s)'
sxaddpar,header,'VERR',verr,' weighted error in VHELIO (km/s)'
sxaddpar,header,'VERR_MED',verr_med,' median error in VHELIO (km/s)'
sxaddpar,header,'VLSR',vlsr,' mean LSR velocity (km/s)'
sxaddpar,header,'VGSR',vgsr,' mean GSR velocity (km/s)'
sxaddpar,header,'OVHELIO',obsvhelio,' mean barycentric velocity (km/s) from OBSVHELIO'
sxaddpar,header,'OVERR',obsverr,' weighted error in OBSVHELIO'
sxaddpar,header,'OVERR_ME',obsverr_med,' median error of OBSVHELIO velocity (km/s) from synth xcorr'
sxaddpar,header,'OVSCAT',obsvscatter,' STDEV of OBSVHELIO (km/s)'
sxaddpar,header,'SVHELIO',synthvhelio,' mean barycentric velocity (km/s) from SYNTHVHELIO'
sxaddpar,header,'SVERR',synthverr,' weighted error in SYNTHVHELIO'
sxaddpar,header,'SVERR_ME',synthverr_med,' median error of SYNTHVHELIO velocity (km/s) from synth xcorr'
sxaddpar,header,'SVSCAT',synthvscatter,' STDEV of SYNTHVHELIO (km/s)'
sxaddpar,header,'SYNTHSCA',synthscatter,' STDEV of OBSVHELIO-SYNTHVHELIO (km/s)'
sxaddpar,header,'SNR',median(starstr.spec[*,0]/starstr.err[*,0]),' final median Signal/Noise ratio'
sxaddpar,header,'STARFLAG',starflag,' bitwise OR of individual visit starflags'
sxaddpar,header,'ANDFLAG',andflag,' bitwise AND of individual visit starflags'
; mean fiber
if ngdv gt 0 then meanfib=TOTAL(visitstr[gdv].fiberid*visitstr[gdv].snr)/TOTAL(visitstr[gdv].snr) else meanfib=999999
if nvisits gt 1 and ngdv gt 1 then sigfib = stddev(visitstr[gdv].fiberid) else sigfib = 0.
sxaddpar,header,'MEANFIB',meanfib,' S/N weighted mean fiber number'
sxaddpar,header,'SIGFIB',sigfib,' standard deviation (unweighted) of fiber number'
; best stellar parameters from RV cross-correlations
if tag_exist(starstr,'rv_feh') then begin
  sxaddpar,header,'CHISQ',starstr.chisq, ' chi square from mini-grid cross correlation'
  sxaddpar,header,'RVFEH',starstr.rv_feh,' metallicity [Fe/H] from mini-grid cross correlation'
  sxaddpar,header,'RVTEFF',starstr.rv_teff,' effective temperature (K) from mini-grid cross correlation'
  sxaddpar,header,'RVLOGG',starstr.rv_logg,' surface gravity (dex) from mini-grid cross correlation'
  sxaddpar,header,'RVALPH',starstr.rv_alpha,' alpha abundance for mini-grid cross correlation'
  sxaddpar,header,'RVCARB',starstr.rv_carb,' carbon abundance for mini-grid cross correlation'
  sxaddpar,header,'CCPFWHM',starstr.ccpfwhm,' FWHM of RV CCF of star with template (km/s)'
  sxaddpar,header,'AUTOFWHM',starstr.autofwhm,' FWHM of RV CCF of template with template (km/s)'
endif else begin 
  sxaddpar,header,'RVFEH',visitstr[0].feh,' metallicity [Fe/H] from visit mini-grid cross correlation'
  sxaddpar,header,'RVTEFF',visitstr[0].teff,' effective temp (K) from visit mini-grid cross correlation'
  sxaddpar,header,'RVLOGG',visitstr[0].logg,' surface gravity from visit mini-grid cross correlation'
endelse

; Put in information for each visit
for i=0,nvisits-1 do begin
  num = strtrim(i+1,2)
  ; Get HJD
  hjd = helio_jd(visitstr[i].jd-2400000.0,visitstr[i].ra,visitstr[i].dec)
  sxaddhist,'VISIT '+num+' INFORMATION:',header,/comment
  sxaddpar,header,'SFILE'+num,file_basename(visitstr[i].file),' Visit #'+num+' spectrum file'
  sxaddpar,header,'DATE'+num,visitstr[i].dateobs,' DATE-OBS of visit '+num
  sxaddpar,header,'JD'+num,visitstr[i].jd,' Julian date of visit '+num
  sxaddpar,header,'HJD'+num,hjd,' Reduced Heliocentric JD of visit '+num
  sxaddpar,header,'FIBER'+num,visitstr[i].fiberid,' Fiber, visit '+num
  sxaddpar,header,'BC'+num,visitstr[i].bc,' Barycentric correction (km/s),visit '+num
  sxaddpar,header,'VTYPE'+num,visitstr[i].vtype,' RV type (1=chisq, 2=xcorr) from visit '+num
  sxaddpar,header,'VRAD'+num,visitstr[i].vrel,' Doppler shift (km/s) of visit '+num
  sxaddpar,header,'VERR'+num,visitstr[i].vrelerr,' error in VRAD (km/s)'
  sxaddpar,header,'VHELIO'+num,visitstr[i].vhelio,' Barycentric velocity (km/s), visit '+num
  sxaddpar,header,'CHISQ'+num,visitstr[i].chisq, ' chi square from visit mini-grid xcorr'
  sxaddpar,header,'RVTEFF'+num,visitstr[i].rv_teff,' effective temperature (K) from visit mini-grid xcorr'
  sxaddpar,header,'RVLOGG'+num,visitstr[i].rv_logg,' surface gravity (dex) from visit mini-grid xcorr'
  sxaddpar,header,'RVFEH'+num,visitstr[i].rv_feh,' metallicity [Fe/H] from visit mini-grid xcorr'
  sxaddpar,header,'RVALPH'+num,visitstr[i].rv_alpha,' alpha abundance from visit mini-grid xcorr'
  sxaddpar,header,'RVCARB'+num,visitstr[i].rv_carb,' carbon abundance from visit mini-grid xcorr'
  sxaddpar,header,'SNRVIS'+num,visitstr[i].snr,' Signal/Noise ratio, visit '+num
  sxaddpar,header,'FLAG'+num,visitstr[i].starflag,' STARFLAG for visit '+num
end

; Pixel limits of each chip
sxaddpar,header,'RMIN',starstr.pixlim[0,0],' Min pixel of red chip contrib, any frame'
sxaddpar,header,'RMAX',starstr.pixlim[1,0],' Maxmm pixel of red chip contrib, any frame'
sxaddpar,header,'GMIN',starstr.pixlim[0,1],' Min pixel of green chip contrib, any frame'
sxaddpar,header,'GMAX',starstr.pixlim[1,1],' Maxm pixel of green chip contrib, any frame'
sxaddpar,header,'BMIN',starstr.pixlim[0,2],' Min pixel of blue chip contrib, any frame'
sxaddpar,header,'BMAX',starstr.pixlim[1,2],' Maxm pixel of blue chip contrib, any frame'
sxaddpar,header,'ROVERMIN',starstr.pixlim_overlap[0,0],' Min pixel of red chip contrib, all frames'
sxaddpar,header,'ROVERMAX',starstr.pixlim_overlap[1,0],' Max pixel of red chip contrib, all frames'
sxaddpar,header,'GOVERMIN',starstr.pixlim_overlap[0,1],' Min pixel of green chip contrib, all frames'
sxaddpar,header,'GOVERMAX',starstr.pixlim_overlap[1,1],' Max pixel of green chip contrib, all frames'
sxaddpar,header,'BOVERMIN',starstr.pixlim_overlap[0,2],' Minimum pixel of blue chip contrib, all frames'
sxaddpar,header,'BOVERMAX',starstr.pixlim_overlap[1,2],' Maximum pixel of blue chip contrib, all frames'

; Wavelength coefficients
sxaddhist,'Wavelength polynomial coefficients (Ang):',header,/comment
nwcoef = n_elements(starstr.wcoef)
if nwcoef gt 2 then begin
 sxaddhist,' wave = wcoef0 + wcoef1*i + wcoef2*i^2 + ...',header,/comment
 sxaddhist,' where i is the pixel number starting with 1',header,/comment
 sxaddpar,header,'NWCOEF',nwcoef,' Number of wavelength coefficients'
 for i=0,nwcoef-1 do sxaddpar,header,'WCOEF'+strtrim(i,2),starstr.wcoef[i]
endif else begin
 sxaddpar,header,'NWAVE',n_elements(starstr.spec[*,0]),'Number of wavelengths in subsequent HDUs'
 sxaddpar,header,'CRVAL1',starstr.wcoef[0],'Start log10(wavelength) in subsequent HDUs'
 sxaddpar,header,'CDELT1',starstr.wcoef[1],'Dispersion in log10(wavelength) in subsequent HDUs'
 sxaddpar,header,'CRPIX1',1,'Pixel of starting wavelength in subsequent HDUs'
 sxaddpar,header,'CTYPE1','LOG-LINEAR','Logarithmic wavelength scale in subsequent HDUs'
 sxaddpar,header,'DC-FLAG',1,'Logarithmic wavelength scale in subsequent HDUs'
endelse

; Gap locations
;mask = total(starstr.comblsf,2) ge 0.1
;gapbeg = where(mask eq 0 and shift(mask,1) eq 1)
;gapend = where(mask eq 0 and shift(mask,-1) eq 1)
;if gapend[0] lt 1000 and n_elements(gapend) gt 1 then remove,0,gapend
;if n_elements(gapbeg) gt 1 and n_elements(gapend) gt 1 then begin
;  sxaddpar,header,'GAP1BEG',gapbeg[0]+1,' First pixel of gap 1 (first pix is 1)'
;  sxaddpar,header,'GAP1END',gapend[0]+1,' Last pixel of gap 1 (first pix is 1)'
;  sxaddpar,header,'GAP2BEG',gapbeg[1]+1,' First pixel of gap 2 (first pix is 1)'
;  sxaddpar,header,'GAP2END',gapend[1]+1,' Last pixel of gap 2 (first pix is 1)'
;endif else begin
;  sxaddhist,' NO GAPS FOUND',header
;endelse

; Information on the HDUs
leadstr = 'APSTAR: '
sxaddpar,head,'V_APSTAR',getvers()
sxaddhist,leadstr+systime(0),header
info = GET_LOGIN_INFO()
sxaddhist,leadstr+info.user_name+' on '+info.machine_name,header
sxaddhist,leadstr+'IDL '+!version.release+' '+!version.os+' '+!version.arch,header
sxaddhist,leadstr+' APOGEE Reduction Pipeline Version: '+getvers(),header
sxaddhist,leadstr+'The data are in separate extensions:',header
sxaddhist,leadstr+' HDU0 = Header only',header
sxaddhist,leadstr+' All image extensions have:',header
sxaddhist,leadstr+'   row 1: combined spectrum with individual pixel weighting',header
sxaddhist,leadstr+'   row 2: combined spectrum with global weighting',header
sxaddhist,leadstr+'   row 3-nvisits+2: individual resampled visit spectra ',header
sxaddhist,leadstr+'  unless nvisits=1, which only have a single row',header
sxaddhist,leadstr+' All spectra shifted to rest (vacuum) wavelength scale',header
sxaddhist,leadstr+' HDU1 - Flux (10^-17 ergs/s/cm^2/Ang)',header
sxaddhist,leadstr+' HDU2 - Error (10^-17 ergs/s/cm^2/Ang)',header
sxaddhist,leadstr+' HDU3 - Flag mask:',header
sxaddhist,leadstr+'   row 1: bitwise OR of all visits',header
sxaddhist,leadstr+'   row 2: bitwise AND of all visits',header
sxaddhist,leadstr+'   row 3-nvisits+2: individual visit masks ',header
sxaddhist,leadstr+' HDU4 - Sky (10^-17 ergs/s/cm^2/Ang)',header
sxaddhist,leadstr+' HDU5 - Sky Error (10^-17 ergs/s/cm^2/Ang)',header
sxaddhist,leadstr+' HDU6 - Telluric',header
sxaddhist,leadstr+' HDU7 - Telluric Error',header
sxaddhist,leadstr+' HDU8 - LSF coefficients',header
sxaddhist,leadstr+' HDU9 - RV and CCF structure',header

; Put files into subdirectories by Location ID
outfile=apogee_filename('Star',field=locationdir,obj=objid)
outdir=file_dirname(outfile)
file_mkdir,outdir

; Create filename
;-----------------
;  apStar-OBJID-MJD.fits
;  Need to append the last of the LAST visit
MJD = max(visitstr.jd - 2400000.5 )  ; last visit

if not keyword_set(silent) then $
  print,' Writing apStar frame to ',outfile

; HDU #0 = Header only
;----------------------
file_delete,outfile,/allow
FITS_WRITE,outfile,0,header
  
; HDU #1 = Flux in units of 10^(-17) erg/s/cm^2/Ang [FLOAT]
;-----------------------------------------------------------
flux = starstr.spec / 1e-17
bderr=where(starstr.err eq baderr(),nbd)
flux[bderr]=0.
MKHDR,header1,flux,/image
if nwcoef eq 2 then begin
 sxaddpar,header1,'CRVAL1',starstr.wcoef[0]
 sxaddpar,header1,'CDELT1',starstr.wcoef[1]
 sxaddpar,header1,'CRPIX1',1
 sxaddpar,header1,'CTYPE1','LOG-LINEAR'
 sxaddpar,header1,'DC-FLAG',1
endif else sxaddpar,header1,'CTYPE1','Pixel'
sxaddpar,header1,'BUNIT','Flux (10^-17 erg/s/cm^2/Ang)'
MWRFITS,flux,outfile,header1,/silent

; HDU #2 = Error in units of 10^(-17) erg/s/cm^2/Ang [FLOAT]
;------------------------------------------------------------
bderr=where(starstr.err eq baderr(),nbd)
err = starstr.err / 1e-17
if nbd gt 0 then err[bderr] = baderr()
MKHDR,header2,err,/image
if nwcoef eq 2 then begin
 sxaddpar,header2,'CRVAL1',starstr.wcoef[0]
 sxaddpar,header2,'CDELT1',starstr.wcoef[1]
 sxaddpar,header2,'CRPIX1',1
 sxaddpar,header2,'CTYPE1','LOG-LINEAR'
 sxaddpar,header2,'DC-FLAG',1
endif else sxaddpar,header2,'CTYPE1','Pixel'
sxaddpar,header2,'BUNIT','Error (10^-17 erg/s/cm^2/Ang)'
MWRFITS,errout(err),outfile,header2,/silent

; HDU #3 = Pixel mask [16-bit INT]
;---------------------------------
mask = fix(starstr.mask)
MKHDR,header3,mask,/image
if nwcoef eq 2 then begin
 sxaddpar,header3,'CRVAL1',starstr.wcoef[0]
 sxaddpar,header3,'CDELT1',starstr.wcoef[1]
 sxaddpar,header3,'CRPIX1',1
 sxaddpar,header3,'CTYPE1','LOG-LINEAR'
 sxaddpar,header3,'DC-FLAG',1
endif else sxaddpar,header3,'CTYPE1','Pixel'
sxaddpar,header3,'BUNIT','Flag Mask (bitwise)'
;sxaddhist,'Explanation of BITWISE flag mask (OR combined)',header3
;sxaddhist,' 1 - bad pixels',header3
;sxaddhist,' 2 - cosmic ray',header3
;sxaddhist,' 4 - saturated',header3
;sxaddhist,' 8 - unfixable',header3
MWRFITS,mask,outfile,header3,/silent

; HDU #4 = Wavelength in units of Ang [DOUBLE]
;------------------------------------------------------------
;wave = double(starstr.wave)
;MKHDR,header4,wave
;if nwcoef eq 2 then begin
; sxaddpar,header4,'CRVAL1',starstr.wcoef[0]
; sxaddpar,header4,'CDELT1',starstr.wcoef[1]
; sxaddpar,header4,'CRPIX1',1
; sxaddpar,header4,'CTYPE1','LOG-LINEAR'
; sxaddpar,header4,'DC-FLAG',1
;endif else sxaddpar,header4,'CTYPE1','Pixel'
;sxaddpar,header4,'BUNIT','Wavelength (Ang)'
;MWRFITS,wave,outfile,header4,/silent


; HDU #4 = Sky flux in units of 10^(-17) erg/s/cm^2/Ang [FLOAT]
;---------------------------------------------------------------
sky = starstr.sky / 1e-17
MKHDR,header4,sky,/image
if nwcoef eq 2 then begin
 sxaddpar,header4,'CRVAL1',starstr.wcoef[0]
 sxaddpar,header4,'CDELT1',starstr.wcoef[1]
 sxaddpar,header4,'CRPIX1',1
 sxaddpar,header4,'CTYPE1','LOG-LINEAR'
 sxaddpar,header4,'DC-FLAG',1
endif else sxaddpar,header4,'CTYPE1','Pixel'
sxaddpar,header4,'BUNIT','Sky (10^-17 erg/s/cm^2/Ang)'
MWRFITS,sky,outfile,header4,/silent

; HDU #5 = Sky error in units of 10^(-17) erg/s/cm^2/Ang [FLOAT]
;----------------------------------------------------------------
skyerr = starstr.skyerr / 1e-17
MKHDR,header5,skyerr,/image
if nwcoef eq 2 then begin
 sxaddpar,header5,'CRVAL1',starstr.wcoef[0]
 sxaddpar,header5,'CDELT1',starstr.wcoef[1]
 sxaddpar,header5,'CRPIX1',1
 sxaddpar,header5,'CTYPE1','LOG-LINEAR'
 sxaddpar,header5,'DC-FLAG',1
endif else sxaddpar,header5,'CTYPE1','Pixel'
sxaddpar,header5,'BUNIT','Sky error (10^-17 erg/s/cm^2/Ang)'
MWRFITS,skyerr,outfile,header5,/silent

; HDU #6 = Telluric absorption flux in units of [FLOAT]
;--------------------------------------------------------
tel = starstr.telluric
MKHDR,header6,tel,/image
if nwcoef eq 2 then begin
 sxaddpar,header6,'CRVAL1',starstr.wcoef[0]
 sxaddpar,header6,'CDELT1',starstr.wcoef[1]
 sxaddpar,header6,'CRPIX1',1
 sxaddpar,header6,'CTYPE1','LOG-LINEAR'
 sxaddpar,header6,'DC-FLAG',1
endif else sxaddpar,header6,'CTYPE1','Pixel'
sxaddpar,header6,'BUNIT','Telluric'
MWRFITS,tel,outfile,header6,/silent

; HDU #7 = Telluric error [FLOAT]
;---------------------------------
telerr = starstr.telerr
MKHDR,header7,tel,/image
if nwcoef eq 2 then begin
 sxaddpar,header7,'CRVAL1',starstr.wcoef[0]
 sxaddpar,header7,'CDELT1',starstr.wcoef[1]
 sxaddpar,header7,'CRPIX1',1
 sxaddpar,header7,'CTYPE1','LOG-LINEAR'
 sxaddpar,header7,'DC-FLAG',1
endif else sxaddpar,header7,'CTYPE1','Pixel'
sxaddpar,header7,'BUNIT','Telluric error'
MWRFITS,telerr,outfile,header7,/silent

; HDU #8 = LSF coefficients [DOUBLE]
;-------------------------------------
lsfcoef = starstr.lcoef
MKHDR,header8,lsfcoef,/image
sxaddpar,header8,'CTYPE1','Parameters'
sxaddpar,header8,'BUNIT','LSF Coefficients'
sxaddhist,'LSF Coefficients to be used with LSF_GH.PRO:',header8
sxaddhist,'  binsize  The width of a pixel in X-units.  If this is non-zero',header8
sxaddhist,'             then a "binned" Gauss-Hermite function is used.  If',header8
sxaddhist,'             binsize=0 then a "normal, unbinned" Gauss-Hermite',header8
sxaddhist,'             function is used.',header8
sxaddhist,'  X0       An additive x-offset.  This is only used to',header8
sxaddhist,'             evaluate the GH parameters that vary globally',header8
sxaddhist,'             with X.',header8
sxaddhist,'  Horder   The highest Hermite order, Horder=0 means',header8
sxaddhist,'             only a constant term (i.e. only Gaussian).',header8
sxaddhist,'             There are Horder Hermite coefficients (since we fix H0=1).',header8
sxaddhist,'  Porder   This array gives the polynomial order for the',header8
sxaddhist,'             global variation (in X) of each LSF parameter.',header8
sxaddhist,'             That includes sigma and the Horder Hermite',header8
sxaddhist,'             coefficients (starting with H1 because we fix H0=1)',header8
sxaddhist,'             There will be Porder[i]+1 coefficients for',header8
sxaddhist,'             parameter i.',header8
sxaddhist,'  GHcoefs  The polynomial coefficients for sigma and the',header8
sxaddhist,'             Horder Hermite parameters.  There are Porder[i]+1',header8
sxaddhist,'             coefficients for parameter i.  The Hermite parameters',header8
sxaddhist,'             start with H1 since we fix H0=1.',header8
MWRFITS,lsfcoef,outfile,header8,/silent

; HDU #9 = RV, CCF, and binary Structure
;-------------------------------
rvstr = {file:visitstr.file,fiber:visitstr.fiberid,plate:visitstr.plate,mjd:visitstr.mjd,$
         locid:visitstr.location_id,jd:visitstr.jd,$
         vrel:visitstr.vrel,vrelerr:visitstr.vrelerr,vhelio:visitstr.vhelio,bc:visitstr.bc, $
         binary: binstr.binary, stablerv_chi2: binstr.stablerv_chi2, $
         stablerv_rchi2: binstr.stablerv_rchi2, $
         chi2_threshold: binstr.chi2_threshold, stablerv_chi2_prob: binstr.stablerv_chi2_prob}
if tag_exist(starstr,'rv_teff') then begin
  rvstr = create_struct(rvstr,'teff',starstr.rv_teff,'logg',starstr.rv_logg,'feh',$
                           starstr.rv_feh,'alpha',starstr.rv_alpha,'carbon',starstr.rv_carb)
endif else begin
  rvstr = create_struct(rvstr,'teff',visitstr[0].teff,'logg',visitstr[0].logg,$
           'feh',visitstr[0].feh,'alpha',visitstr[0].alpha,'carbon',visitstr[0].carbon)
endelse
if tag_exist(visitstr,'CHISQ') then rvstr = create_struct(rvstr,'CHISQ',visitstr.chisq)
if tag_exist(visitstr,'SYNTHVREL') then rvstr=create_struct(rvstr,'SYNTHVREL',visitstr.synthvrel,$
                   'SYNTHVHELIO',visitstr.synthvhelio,'SYNTHVRELERR',visitstr.synthvrelerr)
if tag_exist(visitstr,'OBSVREL') then rvstr=create_struct(rvstr,'OBSVREL',visitstr.obsvrel,$
                   'OBSVHELIO',visitstr.obsvhelio,'OBSVRELERR',visitstr.obsvrelerr)
if tag_exist(starstr,'CCF') then begin
  rvstr = create_struct(rvstr,'ccf',starstr.ccf,'ccferr',starstr.ccferr,'ccflag',starstr.ccflag,'ccfw0',starstr.ccfw0,'ccfdw',starstr.ccfdw)
endif
if tag_exist(starstr,'AUTOCF') then rvstr=create_struct(rvstr,'autocf',starstr.autocf)
MWRFITS,rvstr,outfile,/silent

if not keyword_set(nolsf) then begin

;################################
;  OUTPUT THE COMBINED LSF FILE
;################################

; Create filename
;-----------------
;   apLSF-STAR8.fits
outlsffile=apogee_filename('StarLSF',field=locationdir,obj=objid)
file_delete,outlsffile,/allow

if not keyword_set(silent) then $
  print,' Writing combined LSF to ',outlsffile

; HDU #0 - Header only
;----------------------
MKHDR,lsfheader,0
sxaddpar,lsfheader,'OBJID',objid
sxaddpar,lsfheader,'SPECTRUM',outfile,' apSpec file this corresponds to'
sxaddhist,'HDU 1 = Combined LSF 2D array',lsfheader,/comment
sxaddhist,'HDU 2 = Fitted LSF coefficients',lsfheader,/comment
FITS_WRITE,outlsffile,0,lsfheader

; HDU #1 - Combined LSF 2D array
;--------------------------------
MKHDR,header1,starstr.comblsf,/image
sxaddpar,header1,'BUNIT','LSF'
sxaddhist,'Combined LSF 2D array',header1,/comment
MWRFITS,starstr.comblsf,outlsffile,header1,/silent

; HDU #2 - LSF coefficients [DOUBLE]
;-------------------------------------
lcoef = starstr.lcoef
MKHDR,header2,lcoef,/image
sxaddpar,header2,'CTYPE1','Parameters'
sxaddpar,header2,'BUNIT','LSF Coefficients'
sxaddhist,'These are the FITTED coefficients to the 2D combined',header2,/comment
sxaddhist,'LSF array in HDU #1',header2,/comment
sxaddhist,'LSF Coefficients to be used with LSF_GH.PRO:',header2
sxaddhist,'  binsize  The width of a pixel in X-units.  If this is non-zero',header2
sxaddhist,'             then a "binned" Gauss-Hermite function is used.  If',header2
sxaddhist,'             binsize=0 then a "normal, unbinned" Gauss-Hermite',header2
sxaddhist,'             function is used.',header2
sxaddhist,'  X0       An additive x-offset.  This is only used to',header2
sxaddhist,'             evaluate the GH parameters that vary globally',header2
sxaddhist,'             with X.',header2
sxaddhist,'  Horder   The highest Hermite order, Horder=0 means',header2
sxaddhist,'             only a constant term (i.e. only Gaussian).',header2
sxaddhist,'             There are Horder Hermite coefficients (since we fix H0=1).',header2
sxaddhist,'  Porder   This array gives the polynomial order for the',header2
sxaddhist,'             global variation (in X) of each LSF parameter.',header2
sxaddhist,'             That includes sigma and the Horder Hermite',header2
sxaddhist,'             coefficients (starting with H1 because we fix H0=1)',header2
sxaddhist,'             There will be Porder[i]+1 coefficients for',header2
sxaddhist,'             parameter i.',header2
sxaddhist,'  GHcoefs  The polynomial coefficients for sigma and the',header2
sxaddhist,'             Horder Hermite parameters.  There are Porder[i]+1',header2
sxaddhist,'             coefficients for parameter i.  The Hermite parameters',header2
sxaddhist,'             start with H1 since we fix H0=1.',header2
MWRFITS,lcoef,outlsffile,header2,/silent

endif

pl=1
if keyword_set(pl) then begin
  combspec=starstr.spec
  comberr=starstr.err
  combtel=starstr.telluric
  combsky=starstr.sky
  combmask=starstr.mask
  wave_final=starstr.wave/1.e4
  if keyword_set(localdir) then plotdir=localdir+'/'+locationdir+'/plots/' else plotdir=outdir+'/plots/'
  file_mkdir,plotdir
  set_plot,'PS'
  device,file=plotdir+file_basename(outfile,'.fits')+'.eps',/encap,xsize=48,ysize=12,/color
  smcolor
  ;!p.multi=[0,1,2]
  ;erase
  ;multiplot,[1,2],mxtitle='Wavelength'
  xs=1.5000 & xe=1.7000 & ys=0 & ye=2*median(combspec[*,0])
  
  plot,wave_final,combtel,color=3,xrange=[xs,xe],ystyle=4,charsize=2,thick=2,ytitle='Flux'
  ;!p.multi=[0,1,2]
  plot,wave_final,combspec[*,0],yrange=[0,2*median(combspec[*,0])],xrange=[xs,xe],charsize=2,/noerase,thick=2,ytitle='Flux'
  oplot,wave_final,combsky,color=3
  bd=where(combmask and maskval('PERSIST_LOW'),nmask)
  if nmask gt 0 then oplot,wave_final[bd],fltarr(nmask),psym=6,color=4
  bd=where(combmask and maskval('PERSIST_MED'),nmask)
  if nmask gt 0 then oplot,wave_final[bd],fltarr(nmask),psym=6,color=3
  bd=where(combmask and maskval('PERSIST_HIGH'),nmask)
  if nmask gt 0 then oplot,wave_final[bd],fltarr(nmask),psym=6,color=2
  bd=where(combmask and maskval('LITTROW_GHOST'),nmask)
  if nmask gt 0 then oplot,wave_final[bd],fltarr(nmask),psym=6,color=6
;  for i=2,n_elements(combspec[0,*])-1 do oplot,wave_final,combspec[*,i],color=(2+(i mod 6))
;  oplot,wave_final,combspec[*,1],thick=2,color=3
;  oplot,wave_final,combspec[*,0],thick=2
  ;!p.multi=[1,1,2]
  ;multiplot
  device,/close

  device,file=plotdir+file_basename(outfile,'.fits')+'SN.eps',/encap,xsize=48,ysize=12,/color
  plot,wave_final,combspec[*,0]/comberr[*,0],yrange=[0,2*median(combspec[*,0]/comberr[*,0])],charsize=2,thick=2,ytitle='S/N'
  for i=2,n_elements(comberr[0,*])-1 do oplot,wave_final,combspec[*,i]/comberr[*,i],color=(2+(i mod 6))
;  oplot,wave_final,combspec[*,1]/comberr[*,0],thick=2,color=3
;  oplot,wave_final,combspec[*,0]/comberr[*,0],thick=2
  ;multiplot,/reset
  device,/close
  set_plot,'X'
endif

if keyword_set(stp) then stop

end
