;+
;
; APWAVECAL_MEAN
;
; Take the mean of many apWave calibration files
;
; INPUTS:
;  ids        Array of 8-digit IDs for wavelength calibration files.
;  outname    Name of output file.
;  /shift     Shift each wavelength solution in X to align them before
;                averaging.
;  /verbose   Give verbose output.
;
; OUTPUTS:
;  outfile  The averaged wavelength calibration file is output.
;
; USAGE:
;  IDL>apwavecal_mean,ids,outfile,verbose=verbose
;
; Written by D. Nidever  Aug 2018
;-

pro apwavecal_mean,ids,outname,verbose=verbose,shift=shift

;; Inputs
nfiles = n_elements(ids)
if nfiles eq 0 then begin
  print,'Syntax - apwavecal_mean,ids,outfile,shift=shift,verbose=verbose'
  return
endif  

chiptag = ['a','b','c']

print,'Averaging wavelength calibration for files: ',ids

;; Get the filenames
dirs = getdir(apodir,caldir,spectrodir,vers)
wave_dir = dirs.caldir+'/wave/'  
files = wave_dir+'apWave-a-'+ids+'.fits'
test = file_test(files)
bd = where(test eq 0,nbd)
if nbd gt 0 then begin
  print,'Files not found: ',strjoin(files[bd],' ')
  return
endif

;; Load the coefficients
wstr = replicate({id:'',file:'',mjd:0L,coef:dblarr(300,14,3)},nfiles)
for i=0,nfiles-1 do begin
  cmjd = getcmjd(long(ids[i]),mjd=mjd)
  wfiles = wave_dir+'/apWave-'+chiptag+'-'+ids[i]+'.fits'
  FITS_READ,wfiles[0],coef1,head1,exten=1
  FITS_READ,wfiles[1],coef2,head2,exten=1
  FITS_READ,wfiles[2],coef3,head3,exten=1
  wstr[i].id = ids[i]
  wstr[i].file = wfiles[0]
  wstr[i].mjd = mjd
  wstr[i].coef[*,*,0] = coef1
  wstr[i].coef[*,*,1] = coef2
  wstr[i].coef[*,*,2] = coef3
endfor

;; Initialize the final arrays/headers
  

;; Fiber loop
npix = 2048L
nfibers = 300L
x1 = dindgen(npix)
sharr = dblarr(nfibers,nfiles)
npoly = 8
coefstr = replicate({fiber:0,multi_rms:0.0,chipgap1:0.0,chipgap2:0.0,coef:fltarr(npoly),perror:fltarr(npoly+2),sig:0.0,rms:0.0,rchisq:0.0},nfibers)
For i=0,nfibers-1 do begin

  ;; Construct mean wavelength array
  allwave = dblarr(npix,3,nfiles)
  For j=0,nfiles-1 do begin
    coef1 = reform(wstr[j].coef[i,*,0])
    coef2 = reform(wstr[j].coef[i,*,1])
    coef3 = reform(wstr[j].coef[i,*,2])
    w1 = pix2wave(x1,coef1)
    w2 = pix2wave(x1,coef2)
    w3 = pix2wave(x1,coef3)
    allwave[*,0,j] = w1
    allwave[*,1,j] = w2
    allwave[*,2,j] = w3
  Endfor
  medwave = median(allwave,dim=3,/even)
  ;; Check for outliers
  ;rms = mad(reform(allwave[1000,1,*]))
  meddiff = fltarr(nfiles)
  for j=0,nfiles-1 do meddiff[j]=median(allwave[*,*,j]-medwave)
  rms = mad(meddiff)
  coefstr[i].multi_rms = rms
  ;gdfiles = where(abs(reform(allwave[1000,1,*])-medwave[1000,1]) lt (5*rms>0.1),ngdfiles,comp=bdfiles,ncomp=nbdfiles)
  gdfiles = where(abs(meddiff) lt (5*rms>0.1),ngdfiles,comp=bdfiles,ncomp=nbdfiles)
  if nbdfiles gt 0 then print,strtrim(nbdfiles,2),' outliers found: ',ids[bdfiles]
  ;; Make the mean wavelength solution using "good" files
  mnwave = total(allwave[*,*,gdfiles],3)/ngdfiles

  ;; Measure the X shift for each file
  ;;  respect to the mean wavelengths
  if keyword_set(shift) then begin
    xind1 = lindgen(100)*20+10   ; use every 20th pixel to compare
    xind = [xind1, xind1+2048, xind1+4096]
    shwave = dblarr(npix,3,ngdfiles)
    For j=0,ngdfiles-1 do begin
      coef1 = reform(wstr[gdfiles[j]].coef[i,*,0])
      coef2 = reform(wstr[gdfiles[j]].coef[i,*,1])
      coef3 = reform(wstr[gdfiles[j]].coef[i,*,2])
      w1 = pix2wave(x1,coef1)
      w2 = pix2wave(x1,coef2)
      w3 = pix2wave(x1,coef3)
      w = [w1,w2,w3]

      ;; Use the derivative to the find the shift
      ;; Take the simplest case of a line that is shifted
      ;; y = mx+b, y2 = m(x+dx)+b = mx+b + mdx = y + mdx
      ;; The derivataive is just y' = m, so we have
      ;; y2-y1 = dx*y' or dx = (y2-y1)/y' 
      xsh = median( (w[xind]-mnwave[xind])/(mnwave[xind+1]-mnwave[xind]) )
      sharr[i,j] = xsh
      if keyword_set(verbose) then print,'Fiber=',strtrim(i+1,2),' ID=',ids[j],' Xsh=',strtrim(xsh,2)
      ;; Shifted wavelength arrays
      shwave[*,0,j] = pix2wave(x1-xsh,coef1)
      shwave[*,1,j] = pix2wave(x1-xsh,coef2)
      shwave[*,2,j] = pix2wave(x1-xsh,coef3)
    Endfor
    ;; Find the mean wavelength solution
    mnwave = total(shwave,3)/ngdfiles    
  Endif  ; shift

  ;; Re-fit the coefficients
  xx = [x1,x1,x1]
  yy = [mnwave[*,0],mnwave[*,1],mnwave[*,2]]
  err = yy*0+rms
  chipnum = [lonarr(npix)+1,lonarr(npix)+2,lonarr(npix)+3]
  chipgap1 = -1023.5d0-2048-coef1[0]
  chipgap2 = coef3[0]-2048+1023.5d0
  polypars = coef1[6:*]
  ; X values get divided by 3000 in pix2wave.pro but not func_chipgap_poly
  ; divide polypars by 3000^power
  pow = indgen(n_elements(polypars)-1)+1
  polypars[1:*] /= 3000.^pow
  initpars = [chipgap1, chipgap2, polypars]
  fa = {chipnum:chipnum}
  pars = MPFITFUN('func_chipgap_poly',xx,yy,err,initpars,status=status,dof=dof,$
                  functargs=fa,bestnorm=chisq,perror=perror,yfit=yfit,/quiet)  
  yfit = func_chipgap_poly(xx,pars,chipnum=chipnum,xb=xb)
  rchisq = chisq/dof
  diff = yy-yfit
  sig = mad(diff)
  std = stddev(diff)  

  coefstr[i].fiber = i
  coefstr[i].chipgap1 = pars[0]
  coefstr[i].chipgap2 = pars[1]
  coefstr[i].coef = pars[2:*]
  coefstr[i].perror = perror
  coefstr[i].rchisq = rchisq
  coefstr[i].sig = sig
  coefstr[i].rms = std

  if keyword_set(verbose) then print,i+1,rms,std,format='(I5,2G10.4)'
  
  ;if std gt 1e-5 then stop

  ;stop

Endfor


; Get the wavelength coeffcients and arrays
;------------------------------------------
wavestr = REPLICATE({coef:dblarr(nfibers,6+npoly),wave:dblarr(npix,nfibers)},3)
x1 = findgen(npix)
for i=0,nfibers-1 do begin
  polypars = coefstr[i].coef
  sinepars = [0.0d,0.0,1.0,0.0]
  
  ; X values get divided by 3000 in pix2wave.pro
  ; multiply polypars by 3000^power
  pow = indgen(n_elements(polypars)-1)+1
  polypars[1:*] *= 3000.^pow
  
  ; Plug the values into the output arrays
  ;---------------------------------------
  ;  the parameters for pix2wave.pro are
  ;  [ Yoffset, 4 sine parameters, 7 poly parameters ]
  ; the only difference between the 3 chps is YOFFSET
  ;   chip1: yoffset = 0.
  ;   chip2: yoffset = 2048+chipgap1
  ;   chip3: yoffset = 4096+chipgap1+chipgap2
  ;chipgap1 = fpars[4]
  ;chipgap2 = fpars[5]
  ;coef_chip1 = [0.0, fpars[0:3], fpars[6:12] ]
  ;coef_chip2 = [2048.0+chipgap1, fpars[0:3], fpars[6:12] ]
  ;coef_chip3 = [4096.0+chipgap1+chipgap2, fpars[0:3], fpars[6:12] ]
  xoffset = 0.0
  chipgap1 = coefstr[i].chipgap1
  chipgap2 = coefstr[i].chipgap2
  coef_chip1 = [-1023.5d -2048-chipgap1+xoffset, sinepars, 0.0, polypars ]
  coef_chip2 = [-1023.5d +xoffset, sinepars, 0.0, polypars ]
  coef_chip3 = [-1023.5d +2048+chipgap2+xoffset, sinepars, 0.0, polypars ]
  
  ;XB = X + (chipnum eq 1)*(-1023.5-2048-chipgap1) + $
  ;         (chipnum eq 2)*(-1023.5) + $
  ;         (chipnum eq 3)*(-1023.5+2048+chipgap2)
  
  ; Making the chip-specific coefficients
  wavestr[0].coef[i,*] = coef_chip1
  wavestr[1].coef[i,*] = coef_chip2
  wavestr[2].coef[i,*] = coef_chip3
  
  ; Making wavelengths for each chip
  w1 = pix2wave(x1,coef_chip1)
  w2 = pix2wave(x1,coef_chip2)
  w3 = pix2wave(x1,coef_chip3)
  wavestr[0].wave[*,i] = w1
  wavestr[1].wave[*,i] = w2
  wavestr[2].wave[*,i] = w3
endfor
  
  
; Output the coefficients array and wavelengths per fiber for each chip
;----------------------------------------------------------------------
if not keyword_set(silent) then $
  print,'Writing to = '+wave_dir+dirs.prefix+'Wave-[abc]-'+outname+'.fits'
  
;; Start with the header of the first good file
wfiles0 = wave_dir+'/apWave-'+chiptag+'-'+ids[gdfiles[0]]+'.fits'
  
; Loop through the chips
for i=0,2 do begin
  outfile = wave_dir+dirs.prefix+'Wave-'+chiptag[i]+'-'+outname+'.fits'
  head0 = headfits(wfiles0[i],exten=0)
  head1 = headfits(wfiles0[i],exten=1)
  head2 = headfits(wfiles0[i],exten=2)
  head3 = headfits(wfiles0[i],exten=3)
  ;; Update HDU0
  leadstr = 'APWAVECAL_MEAN'
  sxaddhist,'APWAVECAL_MEAN:  Combined '+strtrim(ngdfiles,2)+' wave cal files',head0
  sxaddhist,leadstr+systime(0),head0
  info = GET_LOGIN_INFO()
  sxaddhist,leadstr+info.user_name+' on '+info.machine_name,head0
  sxaddhist,leadstr+'IDL '+!version.release+' '+!version.os+' '+!version.arch,head0
  ; Put all IDS into header
  for j=0,ngdfiles-1 do sxaddhist,leadstr+'WAVEID'+strtrim(j+1,2)+'='+ids[gdfiles[j]],head0
  for j=0,ngdfiles-1 do sxaddpar,head0,'WAVEID'+strtrim(j+1,2),ids[gdfiles[j]]

  medrms = median(coefstr.rms)
  sxaddhist,leadstr+' Median RMS='+stringize(medrms,ndec=4),head0

  ; HDU0 - header only
  FITS_WRITE,outfile,0,head0  ; write the coeffient array
  ; HDU1 - chip-specific coefficients
  MKHDR,head1,wavestr[i].coef,/image
  MWRFITS,wavestr[i].coef,outfile,head1,/silent  ; write the coeffient array
  ; HDU2 - wavelength arrays
  MKHDR,head2,wavestr[i].wave,/image
  MWRFITS,wavestr[i].wave,outfile,head2,/silent   ; add the wavelength array to the first extension
  ; HDU3 - full coefficients
  MKHDR,head3,transpose(coefstr.coef),/image
  MWRFITS,transpose(coefstr.coef),outfile,head3,/silent  ; the actual coefficients
  
endfor ; chip loop

if keyword_set(stp) then stop

end
