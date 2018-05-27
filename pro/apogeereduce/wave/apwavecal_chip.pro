function pix2wave_func,y,coef,linelist=linelist,gsigma=gsigma

npix = n_elements(y)
if n_elements(gsigma) eq 0 then gsigma = 10.0   ; Gaussian sigma
nlines = n_elements(linelist)

; Model
modarr = fltarr(npix)
for j=0,nlines-1 do begin
  ;yline = (linelist[j].wave-coef[0])/coef[1]
  yline = poly(linelist[j].wave,coef)
  lo = floor(yline-gsigma*5.0)>0
  hi = ceil(yline+gsigma*5.0)<(npix-1)
  if lo lt (npix-1) and hi gt 0 then $
    modarr[lo:hi] += gaussian(y[lo:hi],[1.0,yline,gsigma])
end

return,modarr
end


;#########################################

pro apwavecal_chip,frame,linestr,dispstr,verbose=verbose,silent=silent,$
                   pl=pl,error=error

;+
;
; APWAVECAL_CHIP
;
; Do the wavelength calibration on a ThAr frame.  For a single
; chip.
;
; INPUTS:
;  frame      The extracted 1D frame of a ThAr image (one chip)
;  /pl        Diagnostic plots.
;  /verbose   Print lots of information
;  /silent    Don't print anything to the screen.
;
; OUTPUTS:
;  linestr    Structure of all the lines.  One element
;               for each line on the chip.
;  dispstr    Dispersion function structure.  One element
;               for each fiber.
;  =error     The error message if one occurred.
;
; USAGE:
;  IDL>apwavecal_chip,lampfile,flatfile,linestr,dispstr 
;
; By D.Nidever  Feb. 2010
;-

;lampfile = '/net/stream/apogee/sp2d/ap2D-a-00000014.fits'
;flatfile = '/net/stream/apogee/sp2d/ap2D-a-00000003.fits'
;lampfile = '/net/stream/apogee/sp2d/ap2D-b-00000014.fits'
;flatfile = '/net/stream/apogee/sp2d/ap2D-b-00000003.fits'
;lampfile = '/net/stream/apogee/sp2d/ap2D-c-00000014.fits'
;flatfile = '/net/stream/apogee/sp2d/ap2D-c-00000003.fits'

apgundef,error
;setdisp,/silent

; Get APOGEE directories
dirs=getdir(apogee_dir,cal_dir,spectro_dir,apogee_vers,lib_dir)
if n_elements(direrr) gt 0 then return
linelist_dir = lib_dir+'arclines/'

;nlampfile = n_elements(lampfile)
;nflatfile = n_elements(flatfile)

; Not enough inputs
if n_elements(frame) eq 0 then begin
  print,'Syntax - apwavecal_chip,frame,linestr,dispstr,pl=pl'
  return
endif

npoly = 2  ; polynomial order to use
sz = size(frame.flux)
npix = sz[1]
nfibers = sz[2]
xpix = lindgen(npix)

; Lamp type
;lamptype = sxpar(frame1.header,'LAMPTYPE')
lampfile = file_basename(frame.filename,'.fits')
info = APFILEINFO(lampfile,/silent)

lampframeid = strmid(lampfile,7,8)
cmjd=getcmjd(long(lampframeid),mjd=mjd5)
if mjd5 ge 55982 and mjd5 le 55984 then lamptype = 'FFP' else $
lamptype = ''
if sxpar(frame.header,'LAMPUNE') eq 1 then lamptype='URANIUM'
if sxpar(frame.header,'LAMPTHAR') eq 1 then lamptype='THARNE'
lamptype = strupcase(lamptype)
print,'LAMPTYPE = ',lamptype

; Filename, chip number
filename = frame.filename
base = file_basename(filename,'.fits')    ; ap1D-a-00000014
dum = strsplit(base,'-',/extract)
chip = dum[1]

CASE chip of
'a': chipnum=1
'b': chipnum=2
'c': chipnum=3
else: stop,'CHIP '+chip+' not supported'
ENDCASE

; Load the linelist and reference spectrum
CASE lamptype of
'THARNE': begin
  lamplinelist = IMPORTASCII(linelist_dir+'tharne.lines.vac.apogee',/header,/silent)
  si = sort(lamplinelist.wave)
  lamplinelist = lamplinelist[si]

  ; Load the ThArNe median reference spectrum
  ;--------------------------------------------
  FITS_READ,linelist_dir+'ThArNe_medspectrum_'+dirs.instrument+'.fits',refspec3,refhead
  refspec = refspec3[*,chipnum-1]

end
'URANIUM': begin
  lamplinelist = IMPORTASCII(linelist_dir+'UNe.vac.apogee',/header,/silent)
  si = sort(lamplinelist.wave)
  lamplinelist = lamplinelist[si]

  ; Load the Uranium median reference spectrum
  ;--------------------------------------------
  FITS_READ,linelist_dir+'UNe_medspectrum_'+dirs.instrument+'.fits',refspec3,refhead
  refspec = refspec3[*,chipnum-1]

end
'FFP': begin
  ;lamplinelist=ffplines()
  lamplinelist = IMPORTASCII(linelist_dir+'ffp.dat',/header,/silent)
  FITS_READ,linelist_dir+'FFP_medspectrum_'+dirs.instrument+'.fits',refspec3,refhead
  refspec = refspec3[*,chipnum-1]

end
else: begin
  print,lamptype,' NOT SUPPORTED'
  return
end
ENDCASE

; Get initial wavelength guess
;The wavelength ranges for the three chips are:
;a: 1.5113 - 1.5778 microns
;b: 1.5824 - 1.6403 microns
;c: 1.6444 - 1.6803 microns
;
; /net/stream/apogee/data/20110117/bak/check_wavesol.pro
;  using ap1D-[abc]-00000049.fits (renamed 53) fiber 59 near the center
;coef0arr = [ [16972.776d0, -0.24268043d0 ],$
;             [16449.811d0, -0.28862290d0 ],$
;             [15815.436d0, -0.32527468d0 ] ]
; These were found using the X/WAVE info in the linelist itself
;   which were identified by Shetrone using frame 53.
; VACUUM WAVELENGTHS
if dirs.instrument eq 'apogee-n' then begin
  coef0arr = [ [ 16955.45703d0, -0.2128979266d0, -1.117692409d-05],$
             [ 16434.20508d0, -0.2613874376d0, -1.035568130d-05],$
             [ 15809.69238d0, -0.3065520823d0, -9.610030247d-06] ]
endif else begin
  coef0arr = [[  1.69577252e+04,  -2.14859462e-01,  -1.09959211e-05],$
              [  1.64324720e+04,  -2.63317139e-01,  -1.03074667e-05],$
              [  1.58023346e+04,  -3.08933509e-01,  -9.45618858e-06]]
endelse
coef0 = coef0arr[*,chipnum-1]

; Wavelength range for this chip
wr = poly([0.0d0,npix-1.0d0],coef0)
;gdlines = where(lamplinelist.wave ge min(wr)-100 and lamplinelist.wave le max(wr)+100,ngdlines)
gdlines = where(lamplinelist.wave ge min(wr)-100 and lamplinelist.wave le max(wr)+100 and $
                lamplinelist.flux gt 0,ngdlines)
linelist = lamplinelist[gdlines]
ADD_TAG,linelist,'X',0.0,linelist


;------------------------------------------
; Create median spectrum and fit the lines

; Make a straight median spectrum to get first guess of shift
;  in case there are large shifts
med0 = median(frame.flux,dim=2)
;mask0 = median(frame.mask,dim=2)
;mask0 = CONVOL(mask0,fltarr(41)+1)  ; grow the mask a bit
;mask0 = mask0/(mask0>1)
;mask0 = 1-mask0  ; want 1-good, 0-bad
;XCORLB,med0,refspec,100,xsh0,mask=mask0
;xsh0 = round(xsh0)

; Make "fake" spectra using identified lines and then x-correlate
;  this is more robust and won't get messed up by a few bright
;  lines.
x = findgen(npix)
temp1 = {flux:refspec,err:sqrt(refspec>1),mask:long(refspec)*0}
APPEAKFIT,temp1,reflinestr,/silent,count=nreflinestr
sigref = mad(refspec)
nsig=3
nsig=5
gdrefline = where(reflinestr.height gt nsig*sigref,ngdrefline)  ; only use lines well above the background
if ngdrefline lt 10 then begin   ; if there aren't many use all of the lines
  gdrefline = lindgen(nreflinestr)
  ngdrefline = nreflinestr
endif
refgauss = fltarr(ngdrefline*3)
refgauss[0:3*ngdrefline-3:3] = 1
refgauss[1:3*ngdrefline-2:3] = reflinestr[gdrefline].gaussx
refgauss[2:3*ngdrefline-1:3] = 2
refspec2 = gfunc(x,refgauss)

temp2 = {flux:med0,err:sqrt(med0>1),mask:long(med0*0)}
APPEAKFIT,temp2,medlinestr,/silent,count=nmedlinestr
sigmed = mad(med0)
gdmedline = where(medlinestr.height gt nsig*sigmed,ngdmedline)  ; only use lines well above the background
if ngdmedline lt 10 then begin   ; if there aren't many use all of the lines
  gdmedline = lindgen(nmedlinestr)
  ngdmedline = nmedlinestr
endif
medgauss = fltarr(ngdmedline*3)
medgauss[0:3*ngdmedline-3:3] = 1
medgauss[1:3*ngdmedline-2:3] = medlinestr[gdmedline].gaussx
medgauss[2:3*ngdmedline-1:3] = 2
medspec2 = gfunc(x,medgauss)

;XCORLB,medspec2,refspec2,100,xsh0
XCORLB,medspec2,refspec2,5,xsh0

; Cross-correlate all of the fibers against each other to get
; zero-point shifts and a "global" spectrum
;------------------------------------------
xshift = fltarr(nfibers)
;fiber_ref = frame.flux[*,nfibers/2]
;fiber_ref -= MEDFILT1D(fiber_ref,150,/edge)
fiber_arr = frame.flux*0
x = lindgen(npix)
for i=0,nfibers-1 do begin
  fiber1 = frame.flux[*,i]
  fiber1 -= MEDFILT1D(fiber1,150,/edge)
  err1 = frame.err[*,i]
  mask1 = frame.mask[*,i]

  ; Pre-shift fiber1
  fiber1b = shift(fiber1,-xsh0)
  err1b = shift(err1,-xsh0)
  mask1b = shift(mask1,-xsh0)

  ; grow the bad regions a bit
  mask1b = CONVOL(mask1b,fltarr(41)+1)
  mask1b = mask1b/(mask1b>1)

  mask1b = 1-mask1b  ; want 1-good, 0-bad

  ; Now measure the relative shift
  ;XCORLB,fiber1b,refspec,20,xsh,errspec=err1b
  ;XCORLB,fiber1b,refspec,20,xsh,mask=mask1b
  XCORLB,fiber1b,refspec,5,xsh,mask=mask1b

  xsh += xsh0 ; total shift

  xsh = -xsh  ; flip the sign to be consistent with previous version

  fiber2 = spline(x,fiber1,x-xsh)
  if xsh gt 0 then fiber2[0:ceil(xsh)]=0    ; fix the ends
  if xsh lt 0 then fiber2[npix-floor(abs(xsh))-1:npix-1] = 0
  fiber_arr[*,i] = fiber2
  xshift[i] = xsh

end

;stop

; Create MEDIAN spectrum
;------------------------
medfiber0 = MEDIAN(fiber_arr,dim=2)
; remove a smooth background component
cont_medfiber150 = MEDFILT1D(medfiber0,150,/edge)
cont_medfiber50 = MEDFILT1D(medfiber0,50,/edge)
cont_medfiber = cont_medfiber150
cont_medfiber[0:80] = cont_medfiber50[0:80]
cont_medfiber[npix-80:npix-1] = cont_medfiber50[npix-80:npix-1]
medfiber = medfiber0 - cont_medfiber
medfiber[0:4] = 0.0  ; reference pixels
medfiber[npix-4:npix-1] = 0.0


; Now find the "good" lines in this median spectrum
;--------------------------------------------------
medframe = {flux:medfiber,err:sqrt(medfiber>1),mask:bytarr(npix)}
APPEAKFIT,medframe,medlinestr0,/silent,/nocont
ADD_TAG,medlinestr0,'MODEL_WAVE',0.0d0,medlinestr0
ADD_TAG,medlinestr0,'MODEL_FLUX',0.0,medlinestr0
ADD_TAG,medlinestr0,'MODEL_NAME','',medlinestr0
ADD_TAG,medlinestr0,'MODEL_USEWAVE','',medlinestr0
ADD_TAG,medlinestr0,'MODEL_MATCH',0,medlinestr0
ADD_TAG,medlinestr0,'WAVE_FIT',0.0d0,medlinestr0
sig = sqrt( mad(medfiber)^2 + median(medfiber0) )  ; add in error from continuum
nsig=2
nsig=5
gdmedlines = where(medlinestr0.gpar[0] gt median(medfiber)+nsig*sig,ngdmedlines)
;gdmedlines = where(medlinestr0.gpar[0] gt median(medfiber)+5*sig,ngdmedlines)
medlinestr = medlinestr0[gdmedlines]
;medlinestr = medlinestr0

; Now cross-correlate with the reference spectrum
;--------------------------------------------------
;XCORLB,refspec,medfiber,20,medxsh
medxsh = 0         ; should be lined up now!!!!!

; Get linelist X-values
;w2p_coef = [-coef[0]/coef[1],1.0/coef[1]]   ; wave->pix
;for j=0,ngdlines-1 do linelist[j].x=poly(linelist[j].wave,w2p_coef)
; Get better X-positions for the linelist lines
xx = findgen(3000)-500
ww = poly(xx,coef0)
si = sort(ww)
for j=0,ngdlines-1 do linelist[j].x=spline(ww[si],xx[si],linelist[j].wave)
; Now shift the lines
linelist.x -= medxsh

; Now MATCH the lines
;---------------------
; Loop through the linelist lines
gsigma = 10.0   ; Gaussian sigma
for j=0,n_elements(linelist)-1 do begin
  ydiff = medlinestr.gaussx-linelist[j].x
  closeind = where(abs(ydiff) lt 2.0*gsigma,ncloseind)
  if ncloseind gt 0 then begin
    bestind = where(abs(ydiff) eq min(abs(ydiff)),nbestind)
    medlinestr[bestind[0]].model_wave = linelist[j].wave
    medlinestr[bestind[0]].model_flux = linelist[j].flux
    medlinestr[bestind[0]].model_name = linelist[j].name
    medlinestr[bestind[0]].model_usewave = linelist[j].usewave
    medlinestr[bestind[0]].model_match = 1
  end
end

; Fit the dispersion function
;----------------------------
;npoly = 2
gd = where(medlinestr.model_match eq 1,ngd)
med_fcoef1 = ROBUST_POLY_FIT(medlinestr[gd].gaussx,medlinestr[gd].model_wave,1)    ; linear
med_fcoef = ROBUST_POLY_FIT(medlinestr[gd].gaussx,medlinestr[gd].model_wave,npoly) ; quadratic
med_wpix = POLY(xpix,med_fcoef)
medlinestr.wave_fit = POLY(medlinestr.gaussx,med_fcoef)
med_wrms = MAD(medlinestr[gd].model_wave - medlinestr[gd].wave_fit)
if med_wrms eq 0.0 then med_wrms = ROBUST_SIGMA(medlinestr[gd].model_wave - medlinestr[gd].wave_fit)

; Reject outliers
gd = where(medlinestr.model_match eq 1 and abs(medlinestr.model_wave-medlinestr.wave_fit) lt 4*med_wrms,$
           ngd,complement=bd,ncomplement=nbd)
if nbd gt 0 then begin
  medlinestr[bd].model_match = 0
  medlinestr[bd].model_wave = 0.0
endif
med_wrms = MAD(medlinestr[gd].model_wave - medlinestr[gd].wave_fit)
if med_wrms eq 0.0 then med_wrms = ROBUST_SIGMA(medlinestr[gd].model_wave - medlinestr[gd].wave_fit)
medlinestr_match = medlinestr[gd]
nmatch = ngd

; Get parameter errors
xgerror1 = medlinestr[gd].gerror[1] > med_fcoef[1]/40.  ; lower error threshold
dum = POLY_FIT(medlinestr[gd].gaussx,medlinestr[gd].model_wave,npoly,measure_errors=xgerror1[gd],$
               sigma=perror1,chisq=chisq1,yerror=yerror1)
rchisq1 = chisq1/nmatch

; Plotting
;----------
;pl = 1 ;1
if keyword_set(pl) then begin
  !p.multi=[0,1,2]
  psym8
  plot,medlinestr[gd].gaussx,medlinestr[gd].model_wave/1e4,ps=8,xtit='X',ytit='Wavelength (Microns)',$
       tit='Median Spectrum - '+strtrim(nmatch,2)+' lines - '+strtrim(ngd,2)+' matches',$
       xr=[0,2047],yr=minmax(med_wpix)/1e4,xs=1,ys=1,charsize=1.5
  oplot,xpix,med_wpix/1e4,co=250
  plot,medlinestr[gd].gaussx,medlinestr[gd].model_wave-medlinestr[gd].wave_fit,ps=8,xtit='X',ytit='Residuals',$
       xr=[0,2047],xs=1,ys=1,charsize=1.5,tit='RMS = '+strtrim(med_wrms,2)
  oplot,[0,2047],[0,0],linestyle=2
  !p.multi=0
  wait,0.5
endif


if keyword_set(verbose) then begin
  print,'RMS(FIT) = ',strtrim(med_wrms,2),' Angstroms'
  print,'Dispersion = ',strtrim(med_fcoef1[1],2),' Ang/pix'
endif


;------------------------------------
; Fitting the lines in each fiber
;------------------------------------


if not keyword_set(silent) then begin
  ;print,''
  print,'Fitting lines in the Comparison Lamp Spectra'
  print,''
endif

; Starting LINE structure
dumstr = {fiber:0L,peakx:0.0,gaussx:0.0,pixcen:0.0,pixcenerr:0.0,gfit_pars:fltarr(4),gfit_perror:fltarr(4),gfit_height:0.0,gfit_center:0.0,$
                     gfit_chisq:0.0,gfit_rchisq:0.0,gfit_status:0L,lsffit_pars:fltarr(2),lsffit_perror:fltarr(2),lsffit_flux:0.0,$
                     lsffit_center:0.0,lsffit_chisq:0.0,lsffit_rchisq:0.0,lsffit_status:0L,dof:0L,rms:0.0,$
                     wave_fit:0.0,model_wave:0.0d0,model_flux:0.0,model_name:'',model_usewave:0,model_match:0,nsigdev:0.0,bad:0}
linestr = REPLICATE(dumstr,nfibers*nmatch) ;50*nfibers)
;cntlines = 0

; Find the peaks in the emission line spectra
For i=0L,nfibers-1 do begin

  fiber = reform(frame.flux[*,i])
  var = reform(frame.err[*,i])^2
  ;fiber = reform(wout[*,i,0])
  ;var = reform(wout[*,i,1])

  ; Remove the background
  med = MEDFILT1D(fiber,100,/edge)
  fiber -= med

  resid = fiber            ; lines will be removed

  ; Find peaks, must be greater than neighbors
  ;  set the threshold low
  ;smfiber = GSMOOTH(fiber,2)
  med = median(fiber)
  ;smsig = mad(smfiber)
  sig = mad(fiber)
  nsig = 2  ; 5
  nsig = 5  ; 5
  thresh = med+nsig*sig
  lft = [0.0,fiber[0:npix-2]]
  rgt = [fiber[1:npix-1],0.0]
  peak = where(fiber gt lft AND fiber gt rgt and fiber gt thresh and fiber gt nsig*sqrt(var) and $
               (xpix gt 3 and xpix lt npix-4),npeak)
;print,i,npeak
  ;lft = [0.0,smfiber[0:npix-2]]
  ;rgt = [smfiber[1:npix-1],0.0]
  ;peak = where(smfiber gt lft AND smfiber gt rgt and smfiber gt thresh and smfiber gt nsig*sqrt(var) and $
  ;             (xpix gt 3 and xpix lt npix-4),npeak)
  
  ;if keyword_set(verbose) then $
  ;  print,'Fiber ',strtrim(i+1,2),' - ',strtrim(npeak,2),' Lines found'

  ; Find the lines
  ;fiberframe = {flux:fiber,err:sqrt(var),mask:bytarr(npix)}
  ;APPEAKFIT,fiberframe,ilinestr,/silent,/nocont

  ; Loop through the peaks in the median spectrum
  For j=0L,nmatch-1 do begin

    ; Fit with Gaussians
    ;-------------------
    lineind = i*nmatch+j

    ; Get the shifted position from the median value
    xcen = medlinestr_match[j].gaussx - (xshift[i]-medxsh)
    xsig = medlinestr_match[j].gpar[2]

    ; Is there a peak here
    mindiff = min(abs(peak-xcen))
    ;if mindiff gt 5 then stop
;print,i,j,mindiff
    if mindiff gt 5 then goto,BOMB1
    ;if mindiff gt 10 then goto,BOMB1
    bestpeak = first_el(minloc(abs(peak-xcen)))
    xpeak = peak[bestpeak]

    ;print,i,j,xcen,xpeak

    ; X indices
    ;lo = ( peak[j]-5 ) > 0
    ;hi = ( peak[j]+5 ) < (npix-1)
    ;lo = ( floor(xcen-2.5*xsig) ) > 0
    ;hi = ( ceil(xcen+2.5*xsig) ) < (npix-1)
    lo = ( floor(xpeak-2.0*xsig) ) > 0
    hi = ( ceil(xpeak+2.0*xsig) ) < (npix-1)
    num = hi-lo+1
    xx = lindgen(num)+lo

    ; Fit Gaussian to MID peak
    ;  a = [Height, X, Sigma]
    apgundef,status,dof,chisq,rchisq,yfit,perror,yerror
    estimates1 = medlinestr_match[j].gpar
    estimates1[1] = xpeak ;xcen
    parinfo1 = replicate({fixed:0,limits:[0.0,0.0],limited:[1,1]},4)
    parinfo1[0].limits = [0.5,1.5]*estimates1[0]
    parinfo1[1].limits = [-1,1]+estimates1[1]
    parinfo1[2].limits = [0.2, 5.0]
    parinfo1[3].limits = [-2,2]*sig

    ; Only fit the center for now
    parinfo1[[0,2]].fixed = 1
    gpar1 = MPFITFUN('gaussbin',xx,resid[lo:hi],sqrt(var[lo:hi]),estimates1,parinfo=parinfo1,yfit=yfit1,$
                     perror=perror1,bestnorm=chisq1,dof=dof1,yerror=yerror1,status=status1,/quiet)
    if status1 lt 1 then goto,BOMB1
    pcerror1 = perror1 * sqrt(chisq1/dof1)

    ; Now fit height
    estimates2 = gpar1
    estimates2[0] = resid[xpeak]  ; start with the actual peak value in the spec
    parinfo2 = parinfo1
    parinfo2[0].fixed = 0
    parinfo2[0].limits = [0.5,1.5]*estimates2[0]
    parinfo2[1].fixed = 1
    parinfo2[1].fixed = 1
    gpar2 = MPFITFUN('gaussbin',xx,resid[lo:hi],sqrt(var[lo:hi]),estimates2,parinfo=parinfo2,yfit=yfit2,$
                     perror=perror2,bestnorm=chisq2,dof=dof2,yerror=yerror2,status=status2,/quiet)
    if status2 lt 1 then goto,BOMB1
    pcerror2 = perror2 * sqrt(chisq2/dof2)

    ; Allow all three parameters to vary, slightly
    estimates3 = gpar2
    parinfo3 = parinfo1
    parinfo3.fixed = 0
    parinfo3[0].limits = [0.7,1.3]*estimates3[0]
    parinfo3[1].limits = [-0.5,0.5]+estimates3[1]
    parinfo3[2].limits = [0.7, 1.3]*estimates3[2]
    gpar3 = MPFITFUN('gaussbin',xx,resid[lo:hi],sqrt(var[lo:hi]),estimates3,parinfo=parinfo3,yfit=yfit3,$
                     perror=perror3,bestnorm=chisq3,dof=dof3,yerror=yerror3,status=status3,/quiet)
    if status3 lt 1 then goto,BOMB1
    pcerror3 = perror3 * sqrt(chisq3/dof3)

    if status1 lt 1 or status2 lt 1 or status3 lt 1 then stop

    ; Final Gaussian parameters
    yfit = yfit3
    gpar = gpar3
    pcerror = pcerror3
    chisq = chisq3
    dof = dof3
    status = status3
    yerror = sqrt(mean( (resid[lo:hi]-yfit)^2 ))

    ; Fit with LSF
    ;--------------
    uselsf = tag_exist(frame,'LSFCOEF')
    if keyword_set(uselsf) then begin

      lsfcoef = reform(frame.lsfcoef[i,*])


      ; Convert to **dither combined pixels**
      x = xpix[lo:hi]*2
      xcenter = gpar[1]*2
      dw = 1.0  ;0.3    ; just a guess

      ; DO I NEED TO CONVERT TO DITHER COMBINED PIXELS??????
      stop

      ; Only use pixels well above the background
      ;usepix = where(fiber gt 3*sig,nusepix)
      ;if nusepix eq 0 then goto,BOMB1

      ; Fit the LSF to the line
      apgundef,status2,dof2,chisq2,rchisq2,yfit2,perror2,yerror2
      fa = {coef:lsfcoef}
      initpars = [ gpar[0]*gpar[2]*sqrt(2*!dpi)*2/dw, gpar[1]*2 ]
      parinfo2 = replicate({fixed:0,limits:[0.0,0.0],limited:[1,1]},2)
      parinfo2[0].limits = [ 0.8, 1.2 ]*initpars[0]
      parinfo2[1].limits = [-1.0,1.0]+initpars[1]
      pars4 = MPFITFUN('skyfit_lsf_gh_single',x,resid[lo:hi],sqrt(var[lo:hi]),initpars,functargs=fa,$
                       status=status4,dof=dof4,bestnorm=chisq4,parinfo=parinfo4,perror=perror4,$
                       yfit=yfit4,/quiet)      
      rchisq = chisq4/dof4
      yerror4 = sqrt(mean((resid[lo:hi]-yfit4)^2))
      pcerror4 = perror4 * sqrt(chisq4/dof4)


      ; Final fit and parameters
      yfit = yfit4    ; use LSF fit
      pars = pars4

      ; Remove line from spectrum
      resid_withline = resid  ; resid with this line
      resid[lo:hi] -= yfit

      ; Put in LINESTR
      linestr[lineind].fiber = i
      linestr[lineind].peakx = peak[j]
      linestr[lineind].gaussx = gpar[1]      ; Gaussian center
      linestr[lineind].pixcen = pars4[1]/2   ; use LSF center for FINAL CENTER
      linestr[lineind].pixcenerr = pcerror2[1]/2

      linestr[lineind].gfit_pars = gpar
      linestr[lineind].gfit_perror = pcerror
      linestr[lineind].gfit_height = gpar[0]
      linestr[lineind].gfit_center = gpar[1]
      linestr[lineind].gfit_chisq = chisq
      linestr[lineind].gfit_rchisq = chisq/dof
      linestr[lineind].gfit_status = status

      linestr[lineind].lsffit_pars = pars4      ; LSF fitting
      linestr[lineind].lsffit_perror = pcerror4  ; LSF fitting
      linestr[lineind].lsffit_flux = pars4[0]
      linestr[lineind].lsffit_center = pars4[1]
      linestr[lineind].lsffit_chisq = chisq4
      linestr[lineind].lsffit_rchisq = chisq4/dof4
      linestr[lineind].lsffit_status = status4

      linestr[lineind].dof = dof4
      linestr[lineind].rms = yerror4

    ; No LSF information, use Gaussian
    endif else begin

      ; Final fit and parameters
      ;yfit = yfit1
      ;pars = gpar

      ; Remove line from spectrum
      resid_withline = resid  ; resid with this line
      resid[lo:hi] -= yfit

      ; Put in LINESTR
      linestr[lineind].fiber = i
      linestr[lineind].peakx = xpeak  ;peak[j]
      linestr[lineind].gaussx = gpar[1]      ; Gaussian center
      linestr[lineind].pixcen = gpar[1]
      linestr[lineind].pixcenerr = pcerror[1]

      linestr[lineind].gfit_pars = gpar
      linestr[lineind].gfit_perror = pcerror
      linestr[lineind].gfit_height = gpar[0]
      linestr[lineind].gfit_center = gpar[1]
      linestr[lineind].gfit_chisq = chisq
      linestr[lineind].gfit_rchisq = chisq/dof
      linestr[lineind].gfit_status = status

      linestr[lineind].dof = dof
      linestr[lineind].rms = yerror

    endelse

    ; Add the linelist information
    linestr[lineind].model_wave = medlinestr_match[j].model_wave
    linestr[lineind].model_flux = medlinestr_match[j].model_flux
    linestr[lineind].model_name = medlinestr_match[j].model_name
    linestr[lineind].model_usewave = medlinestr_match[j].model_usewave
    linestr[lineind].model_match = 1


    ; Plotting
    ;pl = 0 ;1 
    lpl=0
    if keyword_set(lpl) and keyword_set(pl) then begin
      plot,xx,fiber[lo:hi],xs=1,ys=1,xtit='X',ytit='Counts',$
           tit='Fiber '+strtrim(i+1,2)+' Peak '+strtrim(j+1,2)+' Chisq='+stringize(chisq/dof,ndec=1)
      oplot,xx,sqrt(var[lo:hi]),co=150
      oplot,xx,resid_withline[lo:hi],co=200,thick=2
      oplot,xx,resid[lo:hi],co=100
      oplot,xx,yfit,co=180
      if keyword_set(uselsf) then begin
        oplot,xx,yfit2,co=250
        legend,['Spectrum','Error','Spectrum - Other lines Removed','Residual','Gaussian Fit','LSF Fit'],$
               textcolor=[255,150,200,100,180,250],/top,/left
      endif else begin
        legend,['Spectrum','Error','Spectrum - Other lines Removed','Residual','Gaussian Fit'],$
               textcolor=[255,150,200,100,180],/top,/left
      endelse
      ;wait,0.1
      ;stop
    endif

    ;stop

    BOMB1:

  End ; Loop through the lines

  ;cntlines += npeak

  ;gdlines = where(linestr.fiber eq i and linestr.gfit_status eq 1,ngdlines)
  ;print,ngdlines

  ;stop

End

; Trim extra lines
;linestr = linestr[0:cntlines-1]

; Remove crappy lines
;bdlines = where(linestr.gfit_pars[2] lt 0.4,nbdlines)
;if nbdlines gt 0 then REMOVE,bdlines,linestr
linestr0 = linestr
bdlines = where(linestr.gfit_status eq 0,nbdlines)
if nbdlines gt 0 then REMOVE,bdlines,linestr

; Remove lines with bad chisq values
;med_chisq = median(linestr.gfit_chisq)
;sig_chisq = mad(linestr.gfit_chisq) > 1
;bdlines = where(linestr.gfit_chisq gt med_chisq+10*sig_chisq,nbdlines)
;if nbdlines gt 0 then REMOVE,bdlines,linestr
; THE CHISQ VALUES ARE HIGH FOR BRIGHT LINES
;  BECAUSE GAUSSIAN IS NOT A GOOD FIT

; How many lines per fiber
nlines = histogram(linestr.fiber,bin=1,min=0,max=nfibers-1)

if not keyword_set(silent) then begin
  print,'Median NLINES per fiber = ',strtrim(long(median(nlines)),2)
  print,''
endif


;stop

;; Get a good estimate for the ZERO-POINT and DISPERSION
;;--------------------------------------------------------
;; using the middle fiber
;midfiber = nfibers/2
;gdobs = where(linestr.fiber eq midfiber,ngdobs)
;linestr1 = linestr[gdobs]
;
;; Observed
;gsigma = 10.0   ; Gaussian sigma
;obsarr = fltarr(npix)
;for j=0,ngdobs-1 do begin
;  lo = floor(linestr1[j].gfit_pars[1]-gsigma*5.0)>0
;  hi = ceil(linestr1[j].gfit_pars[1]+gsigma*5.0)<(npix-1)
;  obsarr[lo:hi] += gaussian(xpix[lo:hi],[1.0,linestr1[j].gfit_pars[1],gsigma])
;end
;
;; Try various dispersions
;ndisp = 31
;disparr = scale_vector(findgen(ndisp),0.7*coef0[1],1.3*coef0[1])
;xcorr_best = fltarr(ndisp)
;shift_best = fltarr(ndisp)
;for j=0,ndisp-1 do begin
;
;  disp1 = disparr[j]
;  tcoef = [coef0[0],disp1]
;  tcoef_inv = [-tcoef[0]/tcoef[1],1.0/tcoef[1]]
;  modarr = pix2wave_func(xpix,tcoef_inv,linelist=linelist,gsigma=gsigma)
;
;  ; Now cross-correlate
;  lag = lindgen(1001)-500
;  xcorr = C_CORRELATE(obsarr,modarr,lag)
;  bestind = where(xcorr eq max(xcorr),nbestind)
;  bestsh = lag[bestind[0]]
;
;  xcorr_best[j] = max(xcorr)
;  shift_best[j] = bestsh
;
;  ;plot,obsarr
;  ;oplot,shift(modarr,-bestsh),co=250
;  ;stop
;
;end
;
;bestind = where(max(xcorr_best) eq xcorr_best)
;dispersion = disparr[bestind[0]]
;bestshift = shift_best[bestind[0]]
;coef = [coef0[0]+dispersion*bestshift, dispersion]
;
;if not keyword_set(silent) then begin
;  print,'Initial wavelength parameters are:'
;  print,'Zero point = ',strtrim(coef[0],2),' Ang'
;  print,'Dispersion = ',strtrim(coef[1],2),' Ang/pix'
;  print,''
;endif
;
;;stop
;
;; AGAIN, Wavelength range for this chip
;wr = poly([0.0d0,npix-1.0d0],coef)
;;gdlines = where(str2.wave ge min(wr)-100 and str2.wave le max(wr)+100,ngdlines)
;gdlines = where(str2.wave ge min(wr)-100 and str2.wave le max(wr)+100 and str2.flux gt 150,ngdlines)
;linelist = str2[gdlines]
;ADD_TAG,linelist,'Y',0.0,linelist

;stop

; Get default linelist X-values
xx = findgen(3000)-500
ww = poly(xx,med_fcoef)
si = sort(ww)
for j=0,ngdlines-1 do linelist[j].x=spline(ww[si],xx[si],linelist[j].wave)
linelist.x += medxsh  ; subtract out the median xshift
linelist0 = linelist


if not keyword_set(silent) then begin
  print,'Matching lines and fitting dispersion functions'
  print,''
endif

; Initialize the dispersion structure
;npoly = 2
dispstr = REPLICATE({fiber:0L,nlines:0L,nmatch:0L,xshift:0.0,coef:fltarr(npoly+1),perror:fltarr(npoly+1),$
                     rms:0.0,chisq:0.0,rchisq:0.0},nfibers)
dispstr.xshift = xshift

; Match the lines and fit the dispersion function
;------------------------------------------------
For i=0,nfibers-1 do begin

  gdobs = where(linestr.fiber eq i,ngdobs)
  if ngdobs eq 0 then continue
  linestr1 = linestr[gdobs]

  ; ONLY SELECT LINES CLOSE TO THE ONES IN THE MEDIAN FIBER!!!!!!


  if keyword_set(verbose) then $
    print,'Fiber ',strtrim(i+1,2)

  ;; Observed
  ;gsigma = 10.0   ; Gaussian sigma
  ;obsarr = fltarr(npix)
  ;for j=0,ngdobs-1 do begin
  ;  lo = floor(linestr1[j].gfit_pars[1]-gsigma*5.0)>0
  ;  hi = ceil(linestr1[j].gfit_pars[1]+gsigma*5.0)<(npix-1)
  ;  obsarr[lo:hi] += gaussian(xpix[lo:hi],[1.0,linestr1[j].gfit_pars[1],gsigma])
  ;end
  ;
  ;; Model
  ;;modarr = fltarr(npix)
  ;;for j=0,ngdlines-1 do begin
  ;;  yline = (linelist[j].wave-coef[0])/coef[1]
  ;;  lo = floor(yline-gsigma*5.0)>0
  ;;  hi = ceil(yline+gsigma*5.0)<(npix-1)
  ;;  if lo lt (npix-1) and hi gt 0 then $
  ;;    modarr[lo:hi] += gaussian(xpix[lo:hi],[1.0,yline,gsigma])
  ;;end
  ;modarr = pix2wave_func(xpix,[-coef[0]/coef[1],1.0/coef[1]],linelist=linelist,gsigma=gsigma)
  ;
  ;; Now cross-correlate
  ;lag = lindgen(101)-50
  ;xcorr = C_CORRELATE(obsarr,modarr,lag)
  ;bestind = where(xcorr eq max(xcorr),nbestind)
  ;bestsh = lag[bestind[0]]
  ;
  ;;xcorr_best[j] = max(xcorr)
  ;;shift_best[j] = bestsh
  ;
  ;; Apply the constant offset
  ;dw = coef[1]
  ;coef1 = coef
  ;coef1[0] += dw*bestsh
  ;
  ;; Now figure out the linear coefficients
  ;func = 'pix2wave_func'
  ;initpars1 = [-coef1[0]/coef1[1],1.0/coef1[1]]
  ;;pars0 = coef1
  ;fa = {linelist:linelist,gsigma:gsigma}
  ;parinfo1 = replicate({fixed:0,limits:[0.0,0.0],limited:[0,0]},2)
  ;parinfo1[0].fixed = 0
  ;parinfo1[0].limits = [-5.0,5.0]+initpars1[0]
  ;parinfo1[0].limited = [1,1]
  ;parinfo1[1].fixed = 0
  ;parinfo1[1].limits = [initpars1[1]/2.0, initpars1[1]*2.0]
  ;parinfo1[1].limited = [1,1]
  ;err = obsarr*0.0 + 1.0
  ;pars1 = MPFITFUN(func, xpix, obsarr, err, initpars1, functargs=fa, bestnorm=bestnorm1, parinfo=parinfo1,$
  ;                dof=dof1,perror=perror1,status=status1,yfit=yfit1,/quiet)
  ;
  ;; Whether this works or not depends a lot on the initial COEF value.
  ;; It might be more robust to loop through ~10 DISPERSION values
  ;; to make sure that we've got the right one.
  ;
  ;; Now fit the quadratic term
  ;initpars2 = [ pars1, 1d-8]
  ;pars2 = MPFITFUN(func, xpix, obsarr, err, initpars2, functargs=fa, bestnorm=bestnorm2,$
  ;                dof=dof2,perror=perror2,status=status2,yfit=yfit2,/quiet)
  ;
  ;; Now fit the cubic term
  ;;initpars3 = [ pars2, 1d-12]
  ;;pars3 = MPFITFUN(func, xpix, obsarr, err, initpars3, functargs=fa, bestnorm=bestnorm3,$
  ;;                dof=dof3,perror=perror3,status=status3,yfit=yfit3)

  ; Use the coefficients from the median spectrum and cross-correlation
  ; to get a first estimate of the coefficients
  ;--------------------------------------------------------------------  
  ;coef1 = med_fcoef1
  ;coef1[0] += xshift[i]*med_fcoef1[1]
  ;pars2 = [-coef1[0]/coef1[1],1.0/coef1[1]]
  ;
  ;; Get X for model lines
  ;linelist.x = 0.0  ; initialize to zero
  ;for j=0,ngdlines-1 do linelist[j].x=poly(linelist[j].wave,pars2)

  ; Get linelist X-values
  linelist = linelist0    ; use the linelist with zero-shift X-values
  ; Now shift the lines
  linelist.x -= xshift[i]

  ; Now MATCH the lines
  ;------------------------
  ;; Loop through the linelist lines
  ;for j=0,n_elements(linelist)-1 do begin
  ;  ydiff = linestr1.gaussx-linelist[j].x
  ;  closeind = where(abs(ydiff) lt 2.0*gsigma,ncloseind)
  ;  if ncloseind gt 0 then begin
  ;    bestind = where(abs(ydiff) eq min(abs(ydiff)),nbestind)
  ;    linestr1[bestind[0]].model_wave = linelist[j].wave
  ;    linestr1[bestind[0]].model_flux = linelist[j].flux
  ;    linestr1[bestind[0]].model_name = linelist[j].name
  ;    linestr1[bestind[0]].model_match = 1
  ;  end
  ;end

  ; Fit the dispersion function
  ;----------------------------
  gd = where(linestr1.model_match eq 1,ngd)
  if ngd lt npoly+1 then begin
    print,'Not enough lines to fit polynomial solution, ifiber: ',i,' ngd: ',ngd
    goto,BOMB
  endif
  fcoef1 = ROBUST_POLY_FIT(linestr1[gd].pixcen,linestr1[gd].model_wave,npoly)
  wpix = POLY(xpix,fcoef1)
  linestr1.wave_fit = POLY(linestr1.pixcen,fcoef1)
  wrms1 = MAD(linestr1[gd].model_wave - linestr1[gd].wave_fit)
  if wrms1 eq 0.0 then wrms1 = ROBUST_SIGMA(linestr1[gd].model_wave - linestr1[gd].wave_fit)

  ; Reject outliers
  gd = where(linestr1.model_match eq 1 and abs(linestr1.model_wave-linestr1.wave_fit) lt 4*wrms1,$
             ngd,complement=bd,ncomplement=nbd)
  if nbd gt 0 then begin
    linestr1[bd].model_match = 0
    linestr1[bd].model_wave = 0.0
  endif
  if ngd eq 0 then begin
    print,'No good lines left'
    continue
  endif
  wrms1 = MAD(linestr1[gd].model_wave - linestr1[gd].wave_fit)
  if wrms1 eq 0.0 then wrms1 = ROBUST_SIGMA(linestr1[gd].model_wave - linestr1[gd].wave_fit)
  nmatch = ngd

  ; Get parameter errors
  ;dum = where(linestr1.model_match eq 1,nmatch)
  ;ygerror1 = linestr1[gd].gerror[1] > fcoef1[1]/40.  ; lower error threshold
  ygerror1 = linestr1[gd].pixcenerr > abs(fcoef1[1])/40.  ; lower error threshold
  ;dum = POLY_FIT(linestr1[gd].gaussx,linestr1[gd].model_wave,npoly,measure_errors=ygerror1[gd],$
  dum = POLY_FIT(linestr1[gd].pixcen,linestr1[gd].model_wave,npoly,measure_errors=ygerror1[gd],$
                 sigma=perror1,chisq=chisq1,yerror=yerror1)
  rchisq1 = chisq1/nmatch


  if keyword_set(verbose) then begin
    print,'RMS(FIT) = ',strtrim(wrms1,2),' Angstroms'
    print,'Dispersion = ',strtrim(fcoef1[1],2),' Ang/pix'
  endif

  ; Plotting
  ;----------
  ;pl=0 ;1
  if keyword_set(pl) then begin
    !p.multi=[0,1,2]
    psym8
    plot,linestr1[gd].pixcen,linestr1[gd].model_wave/1e4,ps=8,xtit='X',ytit='Wavelength (Microns)',$
         tit='Fiber '+strtrim(i+1,2)+' - '+strtrim(ngdobs,2)+' lines - '+strtrim(ngd,2)+' matches',xr=[0,2047],$
         yr=minmax(wpix)/1e4,xs=1,ys=1,charsize=1.5
    oplot,xpix,wpix/1e4,co=250
    plot,linestr1[gd].pixcen,linestr1[gd].model_wave-linestr1[gd].wave_fit,ps=8,xtit='X',ytit='Residuals',$
         xr=[0,2047],xs=1,ys=1,charsize=1.5,tit='RMS = '+strtrim(wrms1,2)+' Ang'
    oplot,[0,2047],[0,0],linestyle=2
    !p.multi=0
  endif

  linestr1.nsigdev = abs(linestr1.model_wave-linestr1.wave_fit)/wrms1
  bad = where( abs(linestr1.model_wave-linestr1.wave_fit) gt 4.0*wrms1,nbad)
  if nbad gt 0 then begin
    linestr1[bad].model_match = 0
    linestr1[bad].bad = 1
  endif
  missing = where( linestr1.model_match eq 0,nmissing)
  if keyword_set(verbose) then begin
    print,strtrim(nbad,2),' outlier(s)'
    print,strtrim(nmissing,2),' missing'
  endif

  ;stop

  ; SECOND ITERATION
  ;------------------
  iter = 0
  if (nbad gt 0 or nmissing gt 0) and keyword_set(iter) then begin

    if keyword_set(verbose) then $
      print,'REMATCHING lines'

    ; Get a new wave->pix coefficients
    pars3 = robust_poly_fit(linestr1[gd].model_wave,linestr1[gd].gaussx,1)
    linelist.x = 0.0  ; initialize to zero
    for j=0,ngdlines-1 do linelist[j].x=poly(linelist[j].wave,pars3)

    ; Loop through the linelist lines
    for j=0,n_elements(linelist)-1 do begin
      wdiff = linestr1.wave_fit-linelist[j].wave
      thresh = 5.0*wrms1 > linestr1.gfit_pars[2]*4.0
      closeind = where(abs(wdiff) lt thresh,ncloseind)
      if ncloseind gt 0 then begin
        bestind = where(abs(ydiff) eq min(abs(ydiff)),nbestind)
        linestr1[bestind[0]].model_wave = linelist[j].wave
        linestr1[bestind[0]].model_flux = linelist[j].flux
        linestr1[bestind[0]].model_name = linelist[j].name
        linestr1[bestind[0]].model_usewave = linelist[j].usewave
        linestr1[bestind[0]].model_match = 1
      end
    end

    ; HOW COME THE CRs ARE COMING THROUGH?????
    ;stop,'how come the CRs are coming through??? the threshold should limit this'
    ; In chip b there are two lines that are very close to each other and
    ; are sometimes blended.

    ; Refit
    gd2 = where(linestr1.model_match eq 1,ngd)
    fcoef2 = ROBUST_POLY_FIT(linestr1[gd2].pixcen,linestr1[gd2].model_wave,npoly)
    wpix = POLY(xpix,fcoef2)
    linestr1.wave_fit = POLY(linestr1.pixcen,fcoef2)
    wrms2 = MAD(linestr1[gd2].model_wave - linestr1[gd2].wave_fit)
    if wrms2 eq 0.0 then wrms2 = ROBUST_SIGMA(linestr1[gd2].model_wave - linestr1[gd2].wave_fit)

    ; You can still have outliers here if a CR is close enough to
    ; a real line to find a "match".
    ; SELECT only "good" lines for a final fit

    ; Get parameter errors
    dum = where(linestr1.model_match eq 1,nmatch)
    ygerror2 = linestr1[gd2].pixcenerr > abs(fcoef2[1])/40.  ; lower error threshold
    dum = POLY_FIT(linestr1[gd2].pixcen,linestr1[gd2].model_wave,npoly,measure_errors=ygerror2[gd2],$
                   sigma=perror2,chisq=chisq2,yerror=yerror2)
    rchisq2 = chisq2/nmatch

    if keyword_set(verbose) then begin
      print,'RMS(FIT) = ',strtrim(wrms2,2),' Angstroms'
      print,'Dispersion = ',strtrim(fcoef2[1],2),' Ang/pix'
    endif

    ; Plotting
    ;----------
    if keyword_set(pl) then begin
      !p.multi=[0,1,2]
      psym8
      plot,linestr1[gd2].pixcen,linestr1[gd2].model_wave/1e4,ps=8,xtit='X',ytit='Wavelength (Microns)',$
           tit='Fiber '+strtrim(i+1,2)+' - '+strtrim(ngdobs,2)+' lines - '+strtrim(ngd2,2)+' matches',$
           xr=[0,2047],yr=minmax(wpix)/1e4,xs=1,ys=1,charsize=1.5
      oplot,xpix,wpix/1e4,co=250
      plot,linestr1[gd2].pixcen,linestr1[gd2].model_wave-linestr1[gd2].wave_fit,ps=8,xtit='X',ytit='Residuals',$
           xr=[0,2047],xs=1,ys=1,charsize=1.5,tit='RMS = '+strtrim(wrms2,2)+' Ang'
      oplot,[0,2047],[0,0],linestyle=2
      !p.multi=0
      wait,0.5
    endif

    linestr1.nsigdev = abs(linestr1.model_wave-linestr1.wave_fit)/wrms2
    bad = where( abs(linestr1.model_wave-linestr1.wave_fit) gt 4.0*wrms2,nbad)
    if nbad gt 0 then begin
      linestr1[bad].model_match = 0
      linestr1[bad].bad = 1
    endif
    missing = where( linestr1.model_match eq 0,nmissing)
    if keyword_set(verbose) then begin
      print,strtrim(nbad,2),' outlier(s)'
      print,strtrim(nmissing,2),' missing'
    endif

    ; Plug it into DISPSTR
    dispstr[i].fiber = i
    dispstr[i].nlines = ngdobs
    dispstr[i].nmatch = nmatch
    dispstr[i].coef = fcoef2
    dispstr[i].perror = perror2
    dispstr[i].rms = wrms2
    dispstr[i].chisq = chisq2
    dispstr[i].rchisq = rchisq2


  ; No iteration
  endif else begin

    ; Plug it into DISPSTR
    dispstr[i].fiber = i
    dispstr[i].nlines = ngdobs
    dispstr[i].nmatch = nmatch
    dispstr[i].coef = fcoef1
    dispstr[i].perror = perror1
    dispstr[i].rms = wrms1
    dispstr[i].chisq = chisq1
    dispstr[i].rchisq = rchisq1

  endelse

  BOMB:


  if ngdobs eq 0 then stop

  ; Put line information back in to LINESTR
  if ngdobs gt 0 then linestr[gdobs] = linestr1

  if keyword_set(verbose) then begin
    print,'Final parameters: ',strtrim(dispstr[i].coef,2)
    print,''
  endif

  if keyword_set(pl) then wait,0.3

  ;if i eq 225 then stop

  ;stop

End

; Get matches for bad fibers
bdfibers = where(dispstr.nmatch eq 0,nbdfibers)
if nbdfibers gt 0 then begin

  print,'Fixing ',strtrim(nbdfibers,2),' bad fibers'

  medcoef = dispstr.coef
  medcoef[*,bdfibers] = !values.f_nan
  medcoef = MEDFILT2D(medcoef,15,dim=2,/edge)
  
  ; Loop through the bad fibers
  For i=0,nbdfibers-1 do begin

    iline = bdfibers[i]

    ; Find the lines
    temp = {flux:frame.flux[*,iline],err:frame.err[*,iline],mask:frame.mask[*,iline]}
    temp.flux -= MEDFILT1D(temp.flux,100,/edge)
    APPEAKFIT,temp,linestr0,count=nlinestr0,/silent
    if nlinestr0 eq 0 then continue

    linestr1 = REPLICATE(dumstr,nlinestr0)
    STRUCT_ASSIGN,linestr0,linestr1
    linestr1.gfit_pars = linestr0.gpar
    linestr1.gfit_perror = linestr0.gerror
    linestr1.gfit_height = linestr0.gpar[0]
    linestr1.gfit_center = linestr0.gpar[1]
    linestr1.gfit_chisq = linestr0.chisq
    linestr1.gfit_rchisq = linestr0.rchisq
    linestr1.gfit_status = linestr0.status
    linestr1.pixcen = linestr0.gpar[1]

    ; Get wavelength from medium coef solution
    wave = poly(linestr1.gaussx,reform(medcoef[*,iline]))

    ; Match to linelist
    SRCOR2,wave,wave*0.,linelist.wave,linelist.wave*0,1.0,ind1,ind2,opt=1,/silent
    dum = where(ind1 ne -1,nmatch)
    if nmatch eq 0 then continue

    linestr1[ind1].model_wave = linelist[ind2].wave
    linestr1[ind1].model_flux = linelist[ind2].flux
    linestr1[ind1].model_name = linelist[ind2].name
    linestr1[ind1].model_usewave = linelist[ind2].usewave
    linestr1[ind1].model_match = 1
    wrms1 = MAD(wave[ind1]-linelist[ind2].wave)

    ; Only keep the matches
    linestr1 = linestr1[ind1]

    ; Fill in the information
    dispstr[iline].fiber = iline
    dispstr[iline].nlines = nlinestr0
    dispstr[iline].nmatch = nmatch
    dispstr[iline].coef = reform(medcoef[*,iline])
    dispstr[iline].xshift = xshift[iline]
    dispstr[iline].rms = wrms1

    ; If we have enough lines then get new dispersion solution
    if nmatch ge npoly+1 then begin
      coef = robust_poly_fit(linestr1.gaussx,linestr1.model_wave,npoly)
      wave2 = poly(linestr1.gaussx,coef)
      linestr1.wave_fit = wave2

      wrms = MAD(linestr1.wave_fit-linestr1.model_wave)
      linestr1.nsigdev = abs(linestr1.model_wave-linestr1.wave_fit)/wrms
      bad = where( abs(linestr1.model_wave-linestr1.wave_fit) gt 4.0*wrms,nbad)
      if nbad gt 0 then begin
        linestr1[bad].model_match = 0
        linestr1[bad].bad = 1
      endif

      dispstr[iline].coef = coef
      dispstr[iline].rms = wrms
    endif

    ; Add the new lines
    PUSH,linestr,linestr1    

  Endfor  ; bad fiber loop

endif ; some bad fibers


medcoef = median(dispstr.coef,dim=2)
if not keyword_set(silent) then $
  print,'Final median parameters are = ',strtrim(medcoef,2)

; We could fit a function to the cefficients as a function of Fiber #
; to improve things.

;stop

end
