;+
;
; APQL_WAVESOL
;
; This program figures out the rough wavelength solution
; for the quicklook software
;
; INPUTS:
;  str      The quicklook structure.
;  allstr   The quicklook structure of all previous reads of this exposure.
;  /debug   Makes some diagnostic plots.
;  /silent  Don't print anything to the screen.
;
; OUTPUTS:
;  The wavelength solution values are updated in the STR structure.
;  =error   The error message if one occurred.
;
; USAGE:
;  IDL>apl_wavesol,str,allstr
;
; By D.Nidever  2010
;-
pro apql_wavesol,str,allstr,verbose=verbose,silent=silent,error=error

; Initialize values to BAD
str.wavefit_pars = 999999.
str.wavefit_rms = 999999
str.wavelength_status = 0

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
if n_elements(str) eq 0 then begin
  error = 'Not enough inputs'
  if not keyword_set(silent) then print,'Syntax - apql_wavesol,str,allstr,silent=silent,error=error'
  return
endif

; Default parameters
npix = 2048L

; make sure the system variable exists before we proceed
DEFSYSV,'!apql',exists=apql_exists
if not apql_exists then begin
  error = 'APQL_WAVESOL: system variable !apql is not defined'
  if not keyword_set(silent) then print,error
  return
endif

; Check that plugmap exists
if not PTR_VALID(!apql.plugmap.datastr) then begin
  error = 'APQL_WAVESOL: no plugmap data available'
  if not keyword_set(silent) then print,error
  return
endif

; Check that airglow lines are there
if not PTR_VALID(!apql.airglow.data) then begin
  error = 'APQL_WAVESOL: no airglow line data'
  if not keyword_set(silent) then print,error
  return
endif

; Check that the extracted spectra are there
if not PTR_VALID(str.frame) then begin
  error = 'APQL_WAVESOL: NO extracted spectra'
  if not keyword_set(silent) then print,error
  return
endif

skyind = where((*!apql.plugmap.datastr).fiberdata.objtype eq 'SKY' AND $
               (*!apql.plugmap.datastr).fiberdata.holetype eq 'OBJECT' AND $
               (*!apql.plugmap.datastr).fiberdata.SPECTROGRAPHID eq 2 AND $
               (*!apql.plugmap.datastr).fiberdata.fiberId ge 1,nskyind)
if nskyind eq 0 then begin
  error = 'APQL_WAVESOL: NO Sky fibers'
  if not keyword_set(silent) then print,error
  return
endif
skyfiberind = (*!apql.plugmap.datastr).fiberdata[skyind].fiberid-1

chipgap1 = 147
chipgap2 = 140

; Load the AIRGLOW linelist
airstr = (*!apql.airglow.data)

; Most of this code was taken from apwavecal_chip.pro

; Get initial wavelength guess
;The wavelength ranges for the three chips are:
;a: 1.5113 - 1.5778 microns
;b: 1.5824 - 1.6403 microns
;c: 1.6444 - 1.6803 microns

; Initialize wavelength range guesses
;wcoef0 = fltarr(3,2)
;wcoef0[0,*] = [15140.0d0, (15805.7d0-15140.0d0)/(2048L-1L)]
;wcoef0[1,*] = [15856.5d0, (16433.3d0-15856.5d0)/(2048L-1L)]
;wcoef0[2,*] = [16476.4d0, (16956.8d0-16476.4d0)/(2048L-1L)]

; VACUUM WAVELENGTHS
wcoef0arr = [ [ 16955.45703d0, -0.2128979266d0, -1.117692409d-05],$
              [ 16434.20508d0, -0.2613874376d0, -1.035568130d-05],$
              [ 15809.69238d0, -0.3065520823d0, -9.610030247d-06] ]
;coef0 = coef0arr[*,chipnum-1]

xpix = findgen(npix)

t0 = systime(1)

; Initialize the dispersion structure
npoly = 2
dispstr = REPLICATE({nlines:0L,nmatch:0L,coef:fltarr(npoly+1),perror:fltarr(npoly+1),rms:0.0,chisq:0.0,rchisq:0.0,xshift:0.0},3)

; Loop through the chips
For i=0,2 do begin

  coef0 = reform(wcoef0arr[*,i])

  skyframe = {flux:(*str.frame)[i].flux[*,skyfiberind],err:(*str.frame)[i].err[*,skyfiberind],$
              mask:(*str.frame)[i].mask[*,skyfiberind]}

  ; We need to be careful because there are slight shifts from
  ; fiber to fiber in the real data.  Might have to remove these

  ; Maybe just median sky spectrum
  medflux = median(skyframe.flux,dim=2)
  mederr = median(skyframe.err,dim=2)
  medframe = {flux:medflux,err:mederr,mask:lonarr(npix)}
  APPEAKFIT,medframe,linestr1,/nogauss,/silent
  nlines = n_elements(linestr1)
  ADD_TAG,linestr1,'CHIPNUM',i+1,linestr1
  ADD_TAG,linestr1,'WAVE_FIT',0.0d0,linestr1
  ADD_TAG,linestr1,'MODEL_WAVE',0.0d0,linestr1
  ADD_TAG,linestr1,'MODEL_FLUX',0.0,linestr1
  ADD_TAG,linestr1,'MODEL_NAME','',linestr1
  ADD_TAG,linestr1,'MODEL_MATCH',0,linestr1

  ; Wavelengths for this chip
  ;wcoef0 = wcoef0arr[*,i]
  w = poly(xpix,coef0)

  ; Wavelength range for this chip
  wr = poly([0.0d0,npix-1.0d0],coef0)
  ;gdlines = where(airstr.wave ge wr[1]-100 and airstr.wave le wr[0]+100,ngdlines)
  gdlines = where(airstr.wave ge wr[1]-50 and airstr.wave le wr[0]+50,ngdlines)
  linelist = airstr[gdlines]
  ADD_TAG,linelist,'X',0.0,linelist

  ; Get X-pixel positions for the lines
  x1 = findgen(npix+1601)-800
  w1 = poly(x1,coef0)
  si1 = sort(w1)
  si_lines = sort(linelist.wave)
  x_lines = spline(w1[si1],x1[si1],linelist[si_lines].wave)
  linelist[si_lines].x = x_lines

  ; Create the model spectrum
  mpar = dblarr(ngdlines*3)
  height_ratio = median([linestr1.height])/median([linelist.emission]>1) > 1
  mpar[indgen(ngdlines)*3] = linelist.emission*height_ratio > 1
  mpar[indgen(ngdlines)*3+1] = linelist.x
  mpar[indgen(ngdlines)*3+2] = 1.5
  mspec = GFUNC(xpix,mpar)

  ; Cross-correlation
  ;  positive xshift means that spec is shifted to the LEFT of refspec
  ;  need to ADD xshift to spec to get it to line up with refspec
  n = 150
  lag = findgen(n)-n/2
  xcorr = C_CORRELATE(medflux,mspec,lag)
  bestind = first_el(maxloc(xcorr))
  if bestind eq -1 then begin
     ; bad fit (medflux is probably 0 because of bad match of sky fibers
     continue
  endif
  xshift0 = lag[bestind]

  ; Fit with Gaussian
  yfit = MPFITPEAK(lag,xcorr,xpar,estimates=[max(xcorr),xshift0,3.0],/gaussian,/positive)
  xshift = xpar[1]
  ;print,'Xshift = ',strtrim(xshift,2),' pixels'


  ; Match the lines
  ;----------------
  ;For j=0,nlines-1 do begin
  ;  xdiff = (linestr1[j].gaussx + xshift) - linelist.x
  ;  closeind = where(abs(xdiff) lt 5.0,ncloseind)
  ;  if ncloseind gt 0 then begin
  ;    bestind = where(abs(xdiff) eq min(abs(xdiff)),nbestind)
  ;    linestr1[j].model_wave = linelist[bestind[0]].wave
  ;    linestr1[j].model_match = 1
  ;  end
  ;End
  ;match = where(linestr1.model_match eq 1,nmatch)
  ;if nmatch eq 0 then goto,BOMB

  SRCMATCH,linelist.x,linelist.x*0,linestr1.gaussx+xshift,linestr1.gaussx*0,5,ind1,ind2,count=nmatch
  if nmatch gt 0 then begin
    linestr1[ind2].model_wave = linelist[ind1].wave
    linestr1[ind2].model_match = 1
  endif
  ; Get better polynomial fit
  match = where(linestr1.model_match eq 1,nmatch)
  ; Make sure we have enough points
  if nmatch ge 3 then begin
    coef1 = ROBUST_POLY_FIT(linestr1[match].gaussx,linestr1[match].model_wave,2,/double)
  endif else begin
    print,'Not enough matched lines'
    goto,BOMB
  endelse
  ; Get better X-values for linelist
  x1 = findgen(npix+1601)-800
  w1 = poly(x1,coef1)
  si1 = sort(w1)
  si_lines = sort(linelist.wave)
  x_lines = spline(w1[si1],x1[si1],linelist[si_lines].wave)
  linelist[si_lines].x = x_lines
  ; Match again
  linestr1.model_wave = 0.0
  linestr1.model_match = 0
  SRCMATCH,linelist.x,linelist.x*0,linestr1.gaussx,linestr1.gaussx*0,5,ind1,ind2,count=nmatch
  if nmatch gt 0 then begin
    linestr1[ind2].model_wave = linelist[ind1].wave
    linestr1[ind2].model_match = 1
  endif else begin
    print,'No matches'
    goto,BOMB
  endelse
  ; Final polynomial coefficients
  match = where(linestr1.model_match eq 1,nmatch)
  if nmatch ge 3 then begin
    coef2 = ROBUST_POLY_FIT(linestr1[match].gaussx,linestr1[match].model_wave,2)
  endif else begin
    print,'Not enough matched lines'
    goto,BOMB
  endelse

  ;stop

  ;;; Fit by allowing the constant offset to vary
  ;;--------------------------------------------
  ;coef1 = ROBUST_POLY_FIT(linestr1[match].gaussx,linestr1[match].model_wave,2)
  ;
  ;parinfo1 = replicate({fixed:0,limited:[0,0],limits:[0.0d0,0.0d0]},n_elements(coef0))
  ;parinfo1[1:11].fixed = 1  ; only allow constant to vary
  ;; Put the X values on an "absolute" scale, need to add offsets
  ;model_wave = fiberlinestr[match].model_wave
  ;xobs = fiberlinestr[match].gaussx/xscale
  ;chip = fiberlinestr[match].chip

  dispstr[i].nlines = nlines
  dispstr[i].nmatch = nmatch
  dispstr[i].coef = coef2
  dispstr[i].xshift = xshift
  ;dispstr[i].rms = rms

  ; Add to ALL lines
  PUSH,linestr,linestr1

  BOMB:

  ;stop

ENDFOR

;medcoef = median(dispstr.coef,dim=2)
;if keyword_set(verbose) then print,'Final median parameters are = ',strtrim(medcoef,2)

; this takes ~0.45 sec so far
;print,systime(1)-t0

;stop

; Now fit a function to all of the lines
;---------------------------------------
; most of this taken from apwavecal.pro

if n_elements(linestr) eq 0 then begin
  error = 'APQL_WAVESOL: no airglow lines matchd'
  if not keyword_set(silent) then print,error
  return
endif

lines = linestr
; Only want lines with matches
gd = where(lines.model_match eq 1,ngd)
lines = lines[gd]
nlines = n_elements(lines)

; Add X
ADD_TAG,lines,'X',0.0d0,lines
lines.x = lines.gaussx


; From Hardware Technical Description document
; chip a: 15140.0 - 15805.7 A
; chip b: 15856.5 - 16433.3 A
; chip c: 16476.4 - 16956.8 A
;
; 2.9mm minimum spacing between the chips
; 18 micron pixels

; the chip gap between chips a+b is ~147.83 pixels
; the chip gap between chips b+c is ~139.50 pixels

; IN THE FUTURE WE SHOULD ADD IN CONSERVATIVE CONSTRAINTS ON THE PARAMETERS

; The parameters are:
;  4 sine parameters
;  2 chip gaps
;  7 poly parameters (first one is a zero-point offset)
;  The chip number must also be input (1, 2 or 3)
initsine = [ 9621.0591d0, -18738.972d0, -340.88099d0, 0.86038618d0 ]
initpoly = [ 0.0, -0.017533513d0, 0.55039528d0, -0.21868473d0, -0.83540194d0, 1.1044074d0,$
             0.37239986d0, -0.83149812d0 ]
initchipgap1 = 148.8
initchipgap2 = 150.7
;initpars = [initsine, initpoly, initxoffset, initchipgap1, initchipgap2]


; The parameters
; 1 xoffset
; 4 sine parameters
; 2 chip gaps
; 6 poly parameters
;initcoef1 = [foutstr.xoffset[ifiber], foutstr.sine, foutstr.chipgap1[ifiber],$
;             foutstr.chipgap2[ifiber], foutstr.poly]
initcoef1 = [median([dispstr.xshift]), initsine, initchipgap1, initchipgap2, initpoly]
;initcoef1 = lastpars
parinfo1 = replicate({fixed:0,limited:[0,0],limits:[0.0d0,0.0d0]},n_elements(initcoef1))
parinfo1[0].fixed = 1     ; fix xoffset for now
parinfo1[1:4].fixed = 0   ; allow sine params to float
;parinfo1[1].fixed = 1    ; fix sine Xoffset
parinfo1[5:6].fixed = 0   ; allow chip gap params to float
;parinfo1[[7,11]].fixed = 1  ; keep constant and cubic terms fixed to ZERO
;parinfo1[7:12].fixed = 1  ; fix poly coefficients to default values
parinfo1[7].fixed = 1     ; fix poly xoffset
;parinfo1[12].fixed = 1    ; fix 5th poly term
parinfo1.limited = 1
parinfo1.limits[0] = initcoef1 - 0.1*abs(initcoef1)
parinfo1.limits[1] = initcoef1 + 0.1*abs(initcoef1)
parinfo1[0].limits = [-1,1]*( abs(initcoef1[0]) > 30 )
;err = lines.gfit_perror[1]
err = lines.gaussx*0+1
bd = where(err eq 0,nbd,comp=gd)
if nbd gt 0 then err[bd] = median([err[gd]])
fa = {chipnum:lines.chipnum}
pars1 = mpfitfun('fit_pix2wave',lines.x,lines.model_wave,err,initcoef1,yfit=yfit1,parinfo=parinfo1,$
                 functargs=fa,status=status1,perror=perror1,bestnorm=chisq1,/quiet)
if status1 lt 1 then begin
  error = 'APQL_WAVESOL: Problem fitting the wavelength solution'
  if not keyword_set(silent) then print,error
  return
endif
yfit1 = fit_pix2wave(lines.x,pars1,chipnum=lines.chipnum,xb=xb1)

; Fit a second time
initcoef2 = pars1
parinfo2 = parinfo1
parinfo2[0].fixed = 0  ; allow xoffset to float
parinfo2.limits[0] = initcoef2 - 0.2*abs(initcoef2)
parinfo2.limits[1] = initcoef2 + 0.2*abs(initcoef2)
parinfo2[0].limits = [-20,20]
pars2 = mpfitfun('fit_pix2wave',lines.x,lines.model_wave,err,initcoef2,yfit=yfit2,parinfo=parinfo2,$
                 functargs=fa,status=status2,perror=perror2,bestnorm=chisq2,/quiet)
if status2 lt 1 then begin
  error = 'APQL_WAVESOL: Problem fitting the wavelength solution'
  if not keyword_set(silent) then print,error
  return
endif

; Remove outliers
diff = lines.model_wave-yfit2
sig = MAD(diff,/zero)
gd = where(abs(diff) lt 4*sig,ngd)
lines2 = lines[gd]

; Fit a final time
initcoef3 = pars2
parinfo3 = parinfo1
parinfo3[0].fixed = 0  ; allow xoffset to float
parinfo3.limits[0] = initcoef3 - 0.2*abs(initcoef3)
parinfo3.limits[1] = initcoef3 + 0.2*abs(initcoef3)
;err2 = lines2.gfit_perror[1]
err2 = lines2.gaussx*0+1
bd = where(err2 eq 0,nbd,comp=gd)
if nbd gt 0 then err2[bd] = median([err2[gd]])
fa2 = {chipnum:lines2.chipnum}
pars3 = mpfitfun('fit_pix2wave',lines2.x,lines2.model_wave,err2,initcoef3,yfit=yfit3,parinfo=parinfo3,$
                 functargs=fa2,status=status3,perror=perror3,bestnorm=chisq3,dof=dof3,/quiet)
if status3 lt 1 then begin
  error = 'APQL_WAVESOL: Problem fitting the wavelength solution'
  if not keyword_set(silent) then print,error
  return
endif
yfit3 = fit_pix2wave(lines2.x,pars3,chipnum=lines2.chipnum,xb=xb2)

; Final values
xb = xb2
yfit = yfit3
fpars = pars3
perror = perror3 * sqrt(chisq3/dof3)
rms = STDDEV(lines2.model_wave-yfit3)
sig = MAD(lines2.model_wave-yfit3,/zero)
rchisq = chisq3/dof3
nlines = n_elements(lines2)

; Sine ONLY part
fpars_sineonly = fpars
fpars_sineonly[7:*] = 0.0
yfit_sineonly = fit_pix2wave(lines2.x,fpars_sineonly,chipnum=lines2.chipnum)
lines2.wave_fit = yfit
ADD_TAG,lines2,'XB',0.0,lines2
lines2.xb = xb

;stop


; Final values
;fpars = pars2
;rms = STDDEV(lines.model_wave-yfit2)
;sig = MAD(lines.model_wave-yfit2)

str.wavefit_pars = fpars
str.wavefit_rms = rms

; Refit with Y errors
;ytop = fit_pix2wave(lines.y-lines.pixcenerr,pars2,chipnum=lines.chipnum)  ; convert X errors to Y errors
;yerr = abs(ytop-yfit2) > 0.001
;parinfo3 = parinfo1
;parinfo3.fixed = 0
;pars3 = mpfitfun('fit_pix2wave',lines.y,lines.model_wave,yerr,initcoef2,yfit=yfit3,parinfo=parinfo3,$
;                 functargs=fa,bestnorm=chisq3,dof=dof3,perror=perror3,status=status3,/quiet)
;perror = perror3*sqrt(chisq3/dof3)
;rchisq = chisq3/nlines
; Don't use the refitted parameters, the X-errors are sometimes
; WAY too small and will throw off the fit.

; Sine ONLY part
;pars2_sineonly = pars2
;pars2_sineonly[6:12] = 0.0
;yfit2_sineonly = fit_pix2wave(lines.y,pars2_sineonly,chipnum=lines.chipnum,yb=yb)
     

if keyword_set(verbose) then print,'Initial RMS = ',strtrim(rms,2)


; Plug values into the structure
;-------------------------------
coefstr = REPLICATE({fiber:0L,nlines:0L,coef:dblarr(15),coeferr:dblarr(15),rchisq:0.0,rms:0.0,niter:0L},1)
;coefstr.fiber = ifiber
coefstr.nlines = nlines
coefstr.coef = fpars
;coefstr.coeferr = perror
;coefstr.rchisq = rchisq
coefstr.rms = rms
;coefstr.niter = count


; Plug the values into the output arrays
;---------------------------------------
;  the parameters for pix2wave.pro are
;  [ Xoffset, 4 sine parameters, 7 poly parameters ]
; the only difference between the 3 chps is YOFFSET
;   chip1: xoffset = 0.
;   chip2: xoffset = 2048+chipgap1
;   chip3: xoffset = 4096+chipgap1+chipgap2
;chipgap1 = fpars[4]
;chipgap2 = fpars[5]
;coef_chip1 = [0.0, fpars[0:3], fpars[6:12] ]
;coef_chip2 = [2048.0+chipgap1, fpars[0:3], fpars[6:12] ]
;coef_chip3 = [4096.0+chipgap1+chipgap2, fpars[0:3], fpars[6:12] ]
xoffset = fpars[0]
chipgap1 = fpars[5]
chipgap2 = fpars[6]  
coef_chip1 = [-1023.5-2048-chipgap1+xoffset, fpars[1:4], fpars[7:*] ]
coef_chip2 = [-1023.5+xoffset, fpars[1:4], fpars[7:*] ]
coef_chip3 = [-1023.5+2048+chipgap2+xoffset, fpars[1:4], fpars[7:*] ]

;coefim1[i,*] = coef_chip1
;coefim2[i,*] = coef_chip2
;coefim3[i,*] = coef_chip3

; Making wavelengths for each chip
x1 = findgen(2048)
w1 = pix2wave(x1,coef_chip1)
w2 = pix2wave(x1,coef_chip2)
w3 = pix2wave(x1,coef_chip3)
;waveim1[i,*] = w1
;waveim2[i,*] = w2
;waveim3[i,*] = w3

;psfile = str.outdir+'apql_wavesol_'+strtrim(str.frameid,2)+'-'+str.readnum
;str.wavesol_figfile = psfile

; Put the wavelength ranges in the structure
;str.waverange_chipa = minmax(w1)  ; chip a
;str.waverange_chipb = minmax(w2)  ; chip b
;str.waverange_chipc = minmax(w3)  ; chip c


; Total fit
;yy = findgen(6432)
;ff = fit_pix2wave(;;,fpars,chipnum=;;*0+1)  ; all chip=1 since we are inputting YB
;  yr1 = [min([lines.model_wave,ff]),max([lines.model_wave,ff])]
;  yr1 = [yr1[0]-0.1*range(yr1),yr1[1]+0.1*range(yr1)]/1e4
xx = scale_vector(findgen(6432),-2198,4246)
;xx = scale_vector(findgen(6432),-3222,3222)
ff = fit_pix2wave(xx,fpars,chipnum=xx*0+2,xb=xx2)  ; all chip=2 since we are inputting XB

; Polynomial Component - Residuals 
ff_sineonly = fit_pix2wave(xx,fpars_sineonly,chipnum=xx*0+2) ; all chip=2 since we are inputting XB
ff_resid = ff-ff_sineonly
wdiff = lines2.model_wave - yfit_sineonly

;; Polynomial component
;ff_sineonly = fit_pix2wave(yy,pars2_sineonly,chipnum=yy*0+1) ; all chip=1 since we are inputting YB
;ff_resid = ff-ff_sineonly
;wdiff = lines.model_wave - yfit2_sineonly

;  yrange = max(wdiff)-min(wdiff)
;  yr2 = [min(wdiff)-0.1*yrange,max(wdiff)+0.1*yrange]


;
;; Wavelength ranges
;;angstrom = '!6!sA!r!u!9 %!6!n'   ; IDL Angstrom symbol
;angstrom = '&#8491;'             ; HTML Angstrom symbol
;PUSH,lines,'<table border=1>'
wr_meas = [w1[0],w1[npix-1], w2[0], w2[npix-1], w3[0], w3[npix-1]]
;;wr_exp = [ 15140.0d, 15805.7d,  15856.5d, 16433.3d,  16476.4d, 16956.8d ]
;wr_exp = [ 15140.060d, 15809.620d, 15855.114d, 16435.235d, 16471.889d, 16956.059d ]
wr_exp = [16955.779d, 16473.414d, 16434.587d, 15856.336d, 15809.828d, 15142.094d ]
dwr = wr_meas - wr_exp

str.waverange_exp[0,*] = wr_exp[0:1]
str.waverange_exp[1,*] = wr_exp[2:3]
str.waverange_exp[2,*] = wr_exp[4:5]
;str.waverange_exp[0,*] = [ 15140.060d, 15809.620d ]
;str.waverange_exp[1,*] = [ 15855.114d, 16435.235d ]
;str.waverange_exp[2,*] = [ 16471.889d, 16956.059d ]
;str.waverange_meas[0,*] = minmax(w1)
;str.waverange_meas[1,*] = minmax(w2)
;str.waverange_meas[2,*] = minmax(w3)
str.waverange_meas[0,*] = wr_meas[0:1]
str.waverange_meas[1,*] = wr_meas[2:3]
str.waverange_meas[2,*] = wr_meas[4:5]
str.waverange_diff = str.waverange_meas - str.waverange_exp


; WAVELENGTH STATUS
;--------------------
; 0-good, 1-bad
str.wavelength_status = 0   ; okay for now

dwr_thresh = 2.0  ; 1.0

; Chip A - Minimum Wavelength
if abs(dwr[0]) gt dwr_thresh then str.wavelength_status = 1
; Chip A - Maximum Wavelength
if abs(dwr[1]) gt dwr_thresh then str.wavelength_status = 1

; Chip B - Minimum Wavelength
if abs(dwr[2]) gt dwr_thresh then str.wavelength_status = 1
; Chip B - Maximum Wavelength
if abs(dwr[3]) gt dwr_thresh then str.wavelength_status = 1

; Chip C - Minimum Wavelength
if abs(dwr[4]) gt dwr_thresh then str.wavelength_status = 1
; Chip C - Maximum Wavelength
if abs(dwr[5]) gt dwr_thresh then str.wavelength_status = 1

; This takes about ~1 sec now.
;print,systime(1)-t0

print,'WAVEFIT_RMS = ',str.wavefit_rms
print,'WAVERANGE_DIFF = ',str.waverange_diff
print,'WAVELENGTH_STATUS = ',str.wavelength_status

;stop

end
