pro ap1dwavecal_refit,frame,outframe,plugmap=plugmap,silent=silent,verbose=verbose,plot=plot,stp=stp,$
                noshift=noshift,pfile=pfile,chipshift=chipshift

;+
;
; AP1DWAVECAL
;
; This wavelength calibrates a single APOGEE frame by using the
; airglow lines.  A constant offset is found and applied.
;
; INPUTS:
;  frame     A structure with the header/data information for an
;                undersampled frame.  This must have the LSF and
;                WAVELENGTH information already appended to the structure.
;  /chipshift  Fit shifts for each chip separately
;  /noshift  Don't shift, just use the wavelength coefficients as is.
;  /silent   Don't print anything to the screen.
;  /verbose  Print lots of information to the screen
;  /plot     Make some plots
;  /stp      Stop at the end of the program.
;
; OUTPUTS:
;  outframe  The same frame but now with wavelength information added.
;
; USAGE:
;  IDL>ap1dwavecal,frame,outframe
;
; By D. Nidever  March 2010
;-

apgundef,outframe

; Not enough inputs
if n_elements(frame) eq 0 then begin
  print,'Syntax - ap1dwavecal,frame,outframe,noshift=noshift,silent=silent,verbose=verbose,pl=pl,stp=stp'
  return
endif

; Get APOGEE directories
dirs=getdir(apogee_dir,cal_dir,spectro_dir,apogee_vers,lib_dir)
linelist_dir = lib_dir+'skylines/'
if FILE_TEST(linelist_dir,/directory) eq 0 then begin
  print,'LINELISTS Directory ',linelist_dir,' NOT FOUND'
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
needtags2 = ['HEADER','FLUX','ERR','MASK','WCOEF']
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
nfibers = sz[2]
; Is this a dither-combined spectrum?
xscale = 1    ; assume original non-dither combined spectrum
if npix eq 4096 then xscale = 2

; Get the wavelength calibration frame for this dither pair
;---------------------------------------------------------------
wcoef1 = frame.chipa.wcoef
wcoef2 = frame.chipb.wcoef
wcoef3 = frame.chipc.wcoef
ncoefpar = n_elements(wcoef1[0,*])

; Wavelength parameters
;wr = MINMAX( [wim1,wim2,wim3] )
;dw = MEDIAN( [wim1-shift(wim1,0,1), wim2-shift(wim2,0,1), wim3-shift(wim3,0,1)] ) / xscale
w0all = fltarr(nfibers)
for i=0,nfibers-1 do w0all[i]=pix2wave(0.0,reform(wcoef1[i,*]))
w1all = fltarr(nfibers)
for i=0,nfibers-1 do w1all[i]=pix2wave(2047.0,reform(wcoef3[i,*]))
;wr = [ MIN(w0all), MAX(w1all) ]
wr = [ MAX(w0all), MIN(w1all) ]
midw = [ pix2wave(findgen(2048),reform(wcoef1[nfibers/2,*])),$
       pix2wave(findgen(2048),reform(wcoef2[nfibers/2,*])),$
       pix2wave(findgen(2048),reform(wcoef3[nfibers/2,*])) ]
dw = abs(MEDIAN(slope(midw)))
nw = ceil((range(wr)/dw)) + 1
;psig = MEDIAN([linestr1.gpar[2],linestr2.par[2],linestr3.gpar[2]])
;wsig = psig*dw
wsig = 3.0*dw
wave = findgen(nw)*dw+min(wr)


chiptag = ['a','b','c']

; Initialize outframe, add WAVELENGTH
For i=0,2 do begin
  chstr0 = frame.(i)
  tags = tag_names(chstr0)

  ; Make the new chip structure
  for j=0,n_elements(tags)-1 do begin
    if j eq 0 then begin
      chstr = CREATE_STRUCT(tags[j],chstr0.(j))
    endif else begin
      chstr = CREATE_STRUCT(chstr,tags[j],chstr0.(j))
    endelse

    ; Add WAVELENGTH after MASK
    if tags[j] eq 'MASK' then $
      chstr = CREATE_STRUCT(chstr,'WAVELENGTH',dblarr(npix,nfibers))
  end

  ; Add to the final OUTFRAME
  if i eq 0 then begin
    outframe = CREATE_STRUCT('chip'+chiptag[i],chstr)
  endif else begin
    outframe = CREATE_STRUCT(outframe,'chip'+chiptag[i],chstr)
  endelse

end


; No SHIFT, just use wavelength coefficients as is
if keyword_set(noshift) then begin

  if not keyword_set(silent) then print,'NOT Shifting Wavelengths'

  x = dindgen(npix) / xscale

  ; Add in the new wavelength arrays
  For i=0,nfibers-1 do begin

    ; Loop through the chips
    For j=0,2 do begin
      ; Calculate new wavelength array and add to the output structure
      wcoef1 = outframe.(j).wcoef[i,*]
      w = PIX2WAVE(x,wcoef1)
      outframe.(j).wavelength[*,i] = w
    Endfor ; chip loop

  Endfor  ; fiber loop

  return

endif  ; /noshift



;; Remove the continuum from the spectra and set lower error boundary
;subframe = frame
;for i=0,2 do begin
;  med = MEDFILT2D(subframe.(i).flux,101,dim=1,/edge)
;  subframe.(i).flux -= med
;  ;sig = MAD(subframe.(i).flux,/zero)
;  ;mederr = median(subframe.(i).err)
;  ;adderr = sqrt(sig^2 - mederr^2)  ; want total median err to be SIG
;  ;subframe.(i).err = sqrt( subframe.(i).err^2 + adderr^2 )
;end

; Find the airglow lines
;---------------------------
;  We know basically where they are, so we could just go look for
;  them.

if keyword_set(plugmap) then begin
  skyplugind = where(plugmap.fiberdata.spectrographid eq 2 and $
                   plugmap.fiberdata.holetype eq 'OBJECT' and $
                   plugmap.fiberdata.objtype eq 'SKY',nsky)
  if nsky eq 0 then begin
    error = 'halt: NO SKY fibers.  CANNOT do sky wavelength calibration.'
    if not keyword_set(silent) then print,error
    stop
  endif else print,'NSKY: ', nsky
  skyfiberid = plugmap.fiberdata[skyplugind].fiberid
  skyindex = 300-skyfiberid
endif else skyindex=indgen(nfibers)
nsky=n_elements(skyindex)

; Check if there is a previously saved linelist to use
base = file_basename(frame.(0).filename,'.fits.fz')
base = file_basename(base,'.fits')
framenum = first_el(strsplit(base,'-',/extract),/last)
;savefile = frame.(0).wave_dir+dirs.prefix+'1dwavecal_'+framenum+'.dat'
;refitlines = 1
;if file_test(savefile) and not keyword_set(refitlines) then begin
;
;  print,'Restoring previously saved linelist file ',savefile
;  restore,savefile
;
;  ind1 = where(linestr.chip eq 1) & linestr1=linestr[ind1]
;  ind2 = where(linestr.chip eq 2) & linestr2=linestr[ind2]
;  ind3 = where(linestr.chip eq 3) & linestr3=linestr[ind3]
;
;endif else begin

;fibers = lindgen(11)*30
;fibers[0] = 1
;fibers[10] = 299

  print,'Finding airglow lines'
  APPEAKFIT,frame.chipa,linestr1,nsigthresh=7,fibers=skyindex
  APPEAKFIT,frame.chipb,linestr2,nsigthresh=7,fibers=skyindex
  APPEAKFIT,frame.chipc,linestr3,nsigthresh=7,fibers=skyindex
  ADD_TAG,linestr1,'CHIP',1,linestr1
  ADD_TAG,linestr2,'CHIP',2,linestr2
  ADD_TAG,linestr3,'CHIP',3,linestr3
  ADD_TAG,linestr1,'WCEN',0.0d0,linestr1
  ADD_TAG,linestr2,'WCEN',0.0d0,linestr2
  ADD_TAG,linestr3,'WCEN',0.0d0,linestr3
  linestr = [linestr1,linestr2,linestr3]
;endelse

; At some point we should change this so that it only looks for peaks
; where the airglow lines should be

; NEED TO FIT THE LINES WITH THE PROPER LSF!!!

; Load the AIRGLOW linelist
airstr = IMPORTASCII(linelist_dir+'airglow_oct18a.txt',/header,/silent)
;airstr = IMPORTASCII(linelist_dir+'airglow.txt',/header,/silent)
;gd = where(airstr.doublet eq 0 and airstr.emission gt 20,ngd)
;airstr[gd].usewave = 1
nairstr = n_elements(airstr)

; Setting up fake spectrum
nbright = nairstr
bright = indgen(nbright)  ; use all for now
;bright = where(airstr.bright eq 1,nbright)
;;nbright = nlinewave
;;bright = indgen(nbright)
refpar = dblarr(nbright*3)
refpar[indgen(nbright)*3] = 1
refpar[indgen(nbright)*3+1] = airstr[bright].wave
refpar[indgen(nbright)*3+2] = wsig
refspec = GFUNC(wave,refpar)

chiptag = ['a','b','c']

; Initialize fiber shift structure
if not keyword_set(chipshift) then begin
  shiftstr = REPLICATE({fiber:0,xscale:xscale,noriglines:0L,nlines:0L,nmatch:0L,xshift0:999999L,$
                        xshift:999999.0,pixshift:999999.0,sig:999999.0,err:999999.0,$
                        interpshift:0.0,useshift:0.0,medoff:0.0,pars1:fltarr(11)},nfibers)
endif else begin
  shiftstr = REPLICATE({fiber:0,xscale:xscale,noriglines:0L,nlines:0L,nmatch:lonarr(3),xshift0:fltarr(3)+999999L,$
                        xshift:fltarr(3)+999999.0,pixshift:fltarr(3)+999999.0,sig:fltarr(3)+999999.0,err:fltarr(3)+999999.0,$
                        interpshift:fltarr(3),useshift:fltarr(3),medoff:fltarr(3),pars1:fltarr(11,3)},nfibers)
endelse

if keyword_set(verbose) then begin
  print,'---------------------'
  print,'FIBER NMATCH PIXSHIFT'
  print,'---------------------'
endif

x = findgen(npix)

apgundef,allfiberlinestr

; Now match the airglow lines 
For i=0,nsky-1 do begin

  ifiber = skyindex[i]

  ;print,'Fiber ',strtrim(ifiber+1,2),'/',strtrim(nsky,2)

  apgundef,fiberlinestr

  shiftstr[i].fiber = ifiber

  spec = fltarr(nw)  ; initializing fake spectrum

  ; Get all of the lines for this fiber
  for j=0,2 do begin
    chiplinestr = (SCOPE_VARFETCH('linestr'+strtrim(j+1,2)))
    wcoef = reform( (SCOPE_VARFETCH('wcoef'+strtrim(j+1,2)))[ifiber,*] )
    ind = where(chiplinestr.fiber eq ifiber,nind)
    if nind eq 0 then goto,BOMB1
    ilinestr = chiplinestr[ind]

    ; Convert pixel positions to wavelengths
    pcen = ilinestr.gpar[1] / xscale    ; convert to undersampled pixel values
    lwcen = PIX2WAVE(pcen,wcoef)
    ilinestr.wcen = lwcen

    ; Add to fiber linelist
    PUSH,fiberlinestr,ilinestr

    BOMB1:

  endfor

  if n_elements(fiberlinestr) eq 0 then goto,BOMB   ; no lines found

  ; Number of original lines
  shiftstr[i].noriglines = n_elements(fiberlinestr)

  ; Remove lines that are too close together
  fiberlinestr0 = fiberlinestr    ; backup
  nlines = n_elements(fiberlinestr)
  bad = lonarr(nlines)
  for j=0,nlines-2 do begin
    wdist = abs(fiberlinestr[j].wcen-fiberlinestr[j+1:nlines-1].wcen)/dw
    bd = where(wdist lt 2,nbd)
    if nbd gt 0 then bad[bd+j+1] = 1
  end
  badind = where(bad eq 1,nbadind)
  if nbadind gt 0 then REMOVE,badind,fiberlinestr

  ; remove lines from superpersistence region
  ;; Add MJD upper limit for when the blue chip was replaced
  ;if ifiber gt 200 then begin
  ;  badind=where(fiberlinestr.chip eq 3,nbadind)
  ;  if nbadind gt 0 then REMOVE,badind,fiberlinestr
  ;endif

  nlines = n_elements(fiberlinestr)

  ; Number of lines
  shiftstr[i].nlines = nlines

  ; Make the "fake" spectrum in wavelength space
  ipar = dblarr(nlines*3)
  ipar[indgen(nlines)*3] = 1.0
  ipar[indgen(nlines)*3+1] = fiberlinestr.wcen
  ipar[indgen(nlines)*3+2] = wsig
  spec = GFUNC(wave,ipar)

  ; Cross-correlation
  ;  positive xshift means that spec is shifted to the LEFT of refspec
  ;  need to ADD xshift to spec to get it to line up with refspec
  n = 50
  lag = findgen(n)-n/2
  xcorr = C_CORRELATE(spec,refspec,lag)
  bestind = first_el(maxloc(xcorr))
  xshift0 = lag[bestind]

  ; Fit with Gaussian
  yfit = MPFITPEAK(lag,xcorr,xpar,estimates=[max(xcorr),xshift0,3.0],/gaussian,/positive,nterms=5)
  xshift = xpar[1]


  ; Match the lines
  ;----------------
  ADD_TAG,fiberlinestr,'MODEL_XPIX',0.0d0,fiberlinestr
  ADD_TAG,fiberlinestr,'MODEL_WAVE',0.0d0,fiberlinestr
  ADD_TAG,fiberlinestr,'MODEL_ID',0,fiberlinestr
  ADD_TAG,fiberlinestr,'MODEL_USEWAVE',0,fiberlinestr
  ADD_TAG,fiberlinestr,'MODEL_MATCH',0,fiberlinestr
  For j=0,nlines-1 do begin
    wdiff = (fiberlinestr[j].wcen + xshift*dw) - airstr.wave
    wthresh = 2.0*dw*fiberlinestr[j].gpar[2] < 1.0
    closeind = where(abs(wdiff) lt wthresh,ncloseind)
    if ncloseind gt 0 then begin
      bestind = where(abs(wdiff) eq min(abs(wdiff)),nbestind)
      fiberlinestr[j].model_xpix = airstr[bestind[0]].xpix
      fiberlinestr[j].model_wave = airstr[bestind[0]].wave
      fiberlinestr[j].model_id = airstr[bestind[0]].id
      fiberlinestr[j].model_usewave = airstr[bestind[0]].usewave
      fiberlinestr[j].model_match = 1
    end
  End
  ;match = where(fiberlinestr.model_match eq 1,nmatch)
  match = where(fiberlinestr.model_match eq 1 and fiberlinestr.model_usewave eq 1,nmatch)
  if nmatch lt 2 then goto,BOMB

  ;; Chip loop
  if keyword_set(chipshift) then nchipshift=3 else nchipshift=1
  for j=0,nchipshift-1 do begin

    if keyword_set(chipshift) then begin
       gdchiplines = where(fiberlinestr.chip eq j+1,ngdchiplines)
       fiberlinestr1 = fiberlinestr[gdchiplines]
    endif else fiberlinestr1=fiberlinestr

    match = where(fiberlinestr1.model_match eq 1 and fiberlinestr1.model_usewave eq 1,nmatch)
    if nmatch lt 2 then goto,BOMB2

    ; Fit by allowing the constant offset to vary
    ;--------------------------------------------
    parinfo1 = replicate({fixed:0,limited:[0,0],limits:[0.0d0,0.0d0]},n_elements(wcoef1[0,*]))
    parinfo1.fixed = 1
    parinfo1[0].fixed = 0  ; only allow constant X-shift to vary

    ; Put the X values on an "absolute" scale, need to add offsets
    model_wave = fiberlinestr1[match].model_wave
    xobs = fiberlinestr1[match].gaussx/xscale
    chip = fiberlinestr1[match].chip

    ; Construct XB, chipnum is 1, 2, or 3
    ; Fix the zero-point at the center of the middle chip
    ;XB = X + (chipnum eq 1)*(-1023.5-2048-chipgap1) + $
    ;         (chipnum eq 2)*(-1023.5) + $
    ;         (chipnum eq 3)*(-1023.5+2048+chipgap2)
    ind1 = where(chip eq 1,nind1)
    if nind1 gt 0 then xobs[ind1] = xobs[ind1]+wcoef1[ifiber,0]
    ind2 = where(chip eq 2,nind2)
    if nind2 gt 0 then xobs[ind2] = xobs[ind2]+wcoef2[ifiber,0]
    ind3 = where(chip eq 3,nind3)
    if nind3 gt 0 then xobs[ind3] = xobs[ind3]+wcoef3[ifiber,0]
    ADD_TAG,fiberlinestr1,'XB',0.0,fiberlinestr1
    fiberlinestr1[match].xb = xobs

    initcoef1 = reform(wcoef2[ifiber,*])
    initcoef1[0] = 0.0
    initcoef1[0] += xshift/xscale  ; add the pixel shift


    ; First attempt  
    ;err = xobs*0.+1.0
    ;err = fiberlinestr1[match].gerror[1] > 0.1
    err = fiberlinestr1[match].gerror[1]*0.26
    bderr = where(err lt 0.001,nbderr)
    if nbderr gt 0 then err[bderr]=0.05
    pars0 = mpfitfun('pix2wave',xobs,model_wave,err,initcoef1,yfit=yfit0,$
                     parinfo=parinfo1,status=status0,perror=perror0,/quiet)

    ; Second attempt. Remove outliers and refit
    ww0 = pix2wave(xobs,pars0)
    diff0 = model_wave - ww0
    med0 = median(diff0)
    sig0 = mad(diff0)
    gdpts0 = where(abs(diff0-med0) lt 3.5*sig0,ngdpts0)
    pars1 = mpfitfun('pix2wave',xobs[gdpts0],model_wave[gdpts0],err[gdpts0],initcoef1,yfit=yfit1,$
                     parinfo=parinfo1,status=status1,perror=perror1,/quiet)

    ; Third attempt. Remove outliers and refit
    ww1 = pix2wave(xobs,pars1)
    diff1 = model_wave - ww1
    med1 = median(diff1)
    sig1 = mad(diff1)
    gdpts1 = where(abs(diff1-med1) lt 3.5*sig1,ngdpts)
    pars2 = mpfitfun('pix2wave',xobs[gdpts1],model_wave[gdpts1],err[gdpts1],initcoef1,yfit=yfit2,$
                     parinfo=parinfo1,status=status2,perror=perror2,bestnorm=chisq2,dof=dof2,/quiet)
    pcerror = perror2*sqrt(chisq2/dof2)  ; scaled uncertainties

    pixshift = pars2[0]
    fww = pix2wave(xobs,pars2)
    fsig = MAD(model_wave[gdpts1] - fww[gdpts1],/zero)

    ; Plotting
    if keyword_set(plot) and not keyword_set(pfile) then begin
      plot,xobs,model_wave-fww,ps=8
      oplot,xobs[gdpts1],model_wave[gdpts1]-yfit2,ps=4,co=250
      oplot,[-10000,10000],[0,0],linestyle=2
    endif

    ; Print out the information
    if keyword_set(chipshift) then begin
      fmt = '(I4,I6,I5,3F9.3)'
      if keyword_set(verbose) then print,ifiber,ngdchiplines,j+1,pixshift,fsig,pcerror[0],format=fmt
    endif else begin
      fmt = '(I4,I6,3F9.3)'
      if keyword_set(verbose) then print,ifiber,nmatch,pixshift,fsig,pcerror[0],format=fmt
    endelse

    ; Put values in the shift structure
    shiftstr[ifiber].nmatch[j] = nmatch
    shiftstr[ifiber].xshift0[j] = xshift0
    shiftstr[ifiber].xshift[j] = xshift
    shiftstr[ifiber].pixshift[j] = pixshift
    shiftstr[ifiber].sig[j] = fsig
    shiftstr[ifiber].err[j] = pcerror[0]

    ; Add residuals
    ADD_TAG,fiberlinestr1,'WDIFF',999999.0,fiberlinestr1
    fiberlinestr1[match].wdiff = fiberlinestr1[match].model_wave-fww

    PUSH,allfiberlinestr,fiberlinestr1
    BOMB2:
   ;stop
  endfor ; chip loop
  BOMB:
Endfor

if keyword_set(verbose) then print,'---------------------'

;; Pick out individual good lines
;; red  X    ID   WAVE
;;   238.8   363  16903.682
;;   526.8   350  16840.495
;;  1732.08  303  16553.822
;;  1934.9   291  16502.367
;; green X   ID   WAVE
;;    82.5   ??   ??
;;   180.9   285  16388.502
;;   895.7   239  16194.626
;;  1298.6   205  16079.769
;;  1467.2   191  16030.848
;;  1666.1   181  15972.610
;; blue  X   ID   WAVE
;;   350.4   118  15702.538
;;   575.5   94   15631.522
;;   682.6   86   15597.638
;;   930.2   63   15517.799
;;   955.7   62   15509.791
;;   982.1   61   15500.873
;;  1303.8   42   15395.335
;;  1492.6   24   15332.407
;;  1624.6   17   15287.794
;;  1763.7    8   15240.961
;;  1922.1    6   15187.145
goodid = [6,8,17,24,42,61,62,63,86,94,118,181,191,205,239,285,291,303,350,363]
undefine,ind
for i=0,n_elements(goodid)-1 do begin
  ind1 = where(allfiberlinestr.model_id eq goodid[i],nind1)
  if nind1 gt 0 then push,ind,ind1
endfor
;glines = allfiberlinestr[ind]
glines = allfiberlinestr

;dw = [0.21,0.26,0.31]
dw = dblarr(2048,3)
for i=0,2 do dw[0:2046,i] = abs(slope(frame.(i).wave[*,150]))
dw[2047,*] = dw[2046,*]
skyshift = fltarr(nsky,3)
set_plot,'X'
!p.multi=[0,1,3]
erase
for i=0,nsky-1 do begin
  ind1 = where(glines.fiber eq skyindex[i] and glines.chip eq 1 and glines.model_wave gt 1,nind1)
  if nind1 gt 1 then begin
    dw1 = dw[round(glines[ind1].gaussx),0]
    skyshift[i,0] = median((glines[ind1].wcen-glines[ind1].model_wave)/dw1)
    !p.multi=[3,1,3]
    plot,glines[ind1].model_wave,glines[ind1].wcen-glines[ind1].model_wave,yr=[-0.2,0.2]
  endif
  ind2 = where(glines.fiber eq skyindex[i] and glines.chip eq 2 and glines.model_wave gt 1,nind2)
  if nind2 gt 1 then begin
    dw2 = dw[round(glines[ind2].gaussx),1]
    skyshift[i,1] = median((glines[ind2].wcen-glines[ind2].model_wave)/dw2)
    !p.multi=[2,1,3]
    plot,glines[ind2].model_wave,glines[ind2].wcen-glines[ind2].model_wave,yr=[-0.2,0.2]
  endif
  ind3 = where(glines.fiber eq skyindex[i] and glines.chip eq 3 and glines.model_wave gt 1,nind3)
  if nind3 gt 3 then begin
    dw3 = dw[round(glines[ind3].gaussx),2]
    skyshift[i,2] = median((glines[ind3].wcen-glines[ind3].model_wave)/dw3)
    !p.multi=[1,1,3]
    plot,glines[ind3].model_wave,glines[ind3].wcen-glines[ind3].model_wave,yr=[-0.2,0.2]
  endif
endfor
ind=where(glines.model_wave gt 1)
openu,1,'sky.dat',/append
for i=0,n_elements(ind)-1 do begin
  j=ind[i]
  disp=dw[round(glines[j].gaussx),glines[j].chip-1]
  fid=strsplit(file_basename(sxpar(frame.chipa.header,'FILENAME')),'-',/ext)
  printf,1,glines[j].chip,glines[j].fiber,glines[j].model_wave,glines[j].gaussx,(glines[j].wcen-glines[j].model_wave)/disp,fid[1],sxpar(frame.chipa.header,'PLATEID'),format='(i5,i5,f12.3,f12.3,f12.3,a12,i6)'
endfor
close,1


set_plot,'PS'
!p.multi=[0,0,0]
device,decomp=0
loadct,39,/silent
psym8
!p.font = 0
file = 'plots/ap1dwavecal_'+framenum+'_pixshift_fiber'
file = pfile
ps_open,file,/color,thick=4,/encap
yr = [min(skyshift),max(skyshift)]
yr = [yr[0]-0.05*range(yr),yr[1]+0.05*range(yr)]
plot,skyindex,skyshift[*,0],ps=1,xr=[0,300],yr=yr,xs=1,ys=1,xtit='Fiber',ytit='Pixel Shift (pix)',tit=framenum,charsize=1.2
oplot,skyindex,skyshift[*,0],ps=1,co=250
coef1 = robust_poly_fitq(skyindex,skyshift[*,0],1)
oplot,skyindex,poly(skyindex,coef1),co=250
oplot,skyindex,skyshift[*,1],ps=1,co=150
coef2 = robust_poly_fitq(skyindex,skyshift[*,1],1)
oplot,skyindex,poly(skyindex,coef2),co=150
oplot,skyindex,skyshift[*,2],ps=1,co=80
coef3 = robust_poly_fitq(skyindex,skyshift[*,2],1)
oplot,skyindex,poly(skyindex,coef3),co=80
ps_close
ps2gif,file+'.eps',/delete,/eps,chmod='664'o
;spawn,['epstopdf',file+'.eps'],/noshell

return

stop

plotc,glines.fiber,glines.wcen-glines.model_wave,glines.xb,ps=1,xtit='Fiber Number',ytit='Wave-Model Wave (A)',tit='Color-coded by XB',charsize=1.2
;plotc,glines.xb,glines.wcen-glines.model_wave,glines.fiber,ps=1,xtit='XB',ytit='Wave-Model Wave (A)',tit='Color-coded by Fiber',charsize=1.2

stop

;; INTERPOLATION
if keyword_set(chipshift) then nchipshift=3 else nchipshift=1
for i=0,nchipshift-1 do begin

  if keyword_set(chipshift) then print,'chip = ',strtrim(i+1,2)

  ; Fit the shifts
  ;x = findgen(nfibers)
  x = shiftstr.fiber
  gd = where(shiftstr.nmatch[i] gt 0,ngd)
  ;shcoef = ROBUST_POLY_FIT(x[gd],shiftstr[gd].pixshift,2)
  shcoef = AP_ROBUST_POLY_FIT(x[gd],shiftstr[gd].pixshift[i],1)
  dum = POLY_FIT(x[gd],shiftstr[gd].pixshift[i],1,measure_errors=shiftstr[gd].err[i],sigma=shcoeferr)
  print,'Fit coefficients: ', shcoef
  print,'Fit coeff errors: ',shcoeferr

  ; Get the interpolated shift at all fibers
  x = findgen(nfibers)
  interpshift = POLY(x,shcoef)
  shiftstr.interpshift[i] = interpshift
  sig = MAD(shiftstr[gd].pixshift[i]-interpshift[gd],/zero)

  ; Which shift to use
  ;   use interpolated shift for fibers with no measured lines
  ;   also use interpolated shifts for fibers with measured shifts
  ;     that are inconsistent with the poly fit (i.e. outliers)
  shiftstr.useshift[i] = shiftstr.pixshift[i]  ; measure pixel shift
  bd = where(shiftstr.nmatch[i] eq 0 OR $
             abs(shiftstr.pixshift[i]-shiftstr.interpshift[i]) gt 3.0*sig,nbd)
  if nbd gt 0 then shiftstr[bd].useshift[i] = shiftstr[bd].interpshift[i]
  if nsky lt nfibers then shiftstr.useshift[i] = shiftstr.interpshift[i]

  ; ALWAYS USE THE POLYNOMIAL FIT!!! REQUIRED IF WE ARE ONLY USING SKY FIBERS
  shiftstr.useshift[i] = shiftstr.interpshift[i]

  ; Plot the shifts
  ;pl=1
  if keyword_set(plot) then begin
    yr = minmax(shiftstr[gd].pixshift[i])
    yr = [yr[0]-0.1*range(yr),yr[1]+0.1*range(yr)]
    if keyword_set(pfile) then begin
      set_plot,'PS'
      file_mkdir,file_dirname(pfile)
      device,file=pfile+'.eps',/encap,/color,xsize=16,ysize=16
      smcolor,/ps
    endif else begin
      smcolor
    endelse
    x = shiftstr.fiber
    plot,x,shiftstr.pixshift[i],ps=1,xr=[0,nfibers],yr=[-0.55,0.25],xs=1,ys=1,$
         xtit='Fiber Number',ytit='Pixel Shift',tit='Measured Pixel Shifts'
    x = findgen(nfibers)
    oplot,x,interpshift,co=3
    if nbd gt 0 then oplot,[x[bd]],[shiftstr[bd].pixshift[i]],ps=4,co=2
    al_legend,['Fit','Bad values','Zero '+string(format='(f8.3)',shcoef[0]),'Slope '+string(format='(e10.2)',shcoef[1])],textcolor=[3,2,1,1],/top,/left
    if keyword_set(pfile) then begin
      device,/close
      ps2gif,pfile+'.eps',/delete,/eps,chmod='664'o
    endif
  endif
stop
endfor ; chip loop

stop

; Add median shift to the header
apaddpar,outframe,'AP1DWAVECAL: Wavelength calibration',/history
apaddpar,outframe,'MEDWSH',median(shiftstr.useshift),' Median wave zero shift (pixels)'
apaddpar,outframe,'AP1DWAVECAL: median SIG of airglow fits = '+stringize(median(shiftstr.sig),ndec=5)+' pix',/history
apaddpar,outframe,'AP1DWAVECAL: linear fit to pixel offsets vs. fiber #',/history
apaddpar,outframe,'AP1DWAVECAL:   COEF = [ '+stringize(shcoef[0],ndec=5)+', '+stringize(shcoef[1],ndec=5)+' ]',/history
apaddpar,outframe,'SHCOEF0',shcoef[0],' pixel shift fit coef 0'
apaddpar,outframe,'SHCOEF1',shcoef[1],' pixel shift fit coef 1'
apaddpar,outframe,'AP1DWAVECAL: SIG of pixel shifts = '+stringize(MAD(shiftstr.pixshift),ndec=5)+' pix',/history

; Add in the new wavelength arrays
if keyword_set(verbose) then print,'Computing new Wavelength arrays'
For i=0,nfibers-1 do begin

  pixshift = shiftstr[i].useshift

  ; Calculate new wavelength array and add to the output structure
  x = dindgen(npix) / xscale
  coef1 = wcoef1[i,*]
  coef1[0] = coef1[0]+pixshift  ; add pixel shift
  warr1 = PIX2WAVE(x,coef1)
  ;outframe.chipa.data[*,i,1] = warr1  ; wave in second plane
  outframe.chipa.wavelength[*,i] = warr1  ; wave in second plane
  outframe.chipa.wcoef[i,*] = coef1

  coef2 = wcoef2[i,*]
  coef2[0] = coef2[0]+pixshift  ; add pixel shift
  warr2 = PIX2WAVE(x,coef2)
  ;outframe.chipb.data[*,i,1] = warr2  ; wave in second plane
  outframe.chipb.wavelength[*,i] = warr2  ; wave in second plane
  outframe.chipb.wcoef[i,*] = coef2

  coef3 = wcoef3[i,*]
  coef3[0] = coef3[0]+pixshift  ; add pixel shift
  warr3 = PIX2WAVE(x,coef3)
  ;outframe.chipc.data[*,i,1] = warr3  ; wave in second plane
  outframe.chipc.wavelength[*,i] = warr3  ; wave in second plane
  outframe.chipc.wcoef[i,*] = coef3

  ;stop

End

; Add the structures to FRAME
frame = create_struct(temporary(frame),'LINESTR',allfiberlinestr,'WSHIFTSTR',shiftstr)

;; Save the fitting information
;if not keyword_set(refitlines) then begin
;  base = file_basename(frame.(0).filename,'.fits.fz')
;  base = file_basename(base,'.fits')
;  framenum = first_el(strsplit(base,'-',/extract),/last)
;  savefile = frame.(0).wave_dir+dirs.prefix+'1Dwavecal-'+framenum+'.dat'
;  save,shiftstr,linestr,allfiberlinestr,file=savefile
;  ; save the shift information in a table
;  mwrfits,shiftstr,frame.(0).wave_dir+dirs.prefix+'1Dwavecal-'+framenum+'.fits',/create
;endif


stop

if keyword_set(stp) then stop

end
