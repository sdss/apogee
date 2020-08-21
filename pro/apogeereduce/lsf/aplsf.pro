pro aplsf,lampframes,waveid,lsfim,fitmethod=fitmethod,nsigfit=nsigfit,gauss=gauss,$
          psfid=psfid,Porder=Porder,WPorder=WPorder,Wproftype=Wproftype,pl=pl,$
          verbose=verbose,silent=silent,outdir=outdir,clobber=clobber,stp=stp,fibers=fibers

;+
;
; APLSF
;
; This program computes the 1D LSF for a dither combined pair of ThAr
; observations.
;
; INPUTS:
;  lampframes   The directory and ID8 numbers (concatenated) for a
;                dither pair of ThAr observations (a two-element
;                array).  Their 1D extracted spectra must already
;                have been created.
;  waveid      The directory and MJD5 number (concatenated) for the
;                wavelength calibration file to use.
;  =fitmethod  The method to fit the LSF (per fiber):
;                (1) Fit each line separately
;                (2) Fit all lines in a chip together (the default)
;                (3) Fit all lines together.
;  =nsigfit    The number of standard deviations (+/-) around a line
;                to fit.  The default is nsigfit=7.
;  =Porder     The polynomial order for how the Gauss-Hermite
;                coefficients are allowed to globally vary with Y.
;                The default is Porder=[2,1,1,1,1,0]
;  =WPorder    The polynomial order for the two wing parameter
;                coefficients are allowed to globally vary with Y.
;                The default is WPorder=[0,0]
;  =Wproftype  The wing profile type.  The supported types are:
;                 1-Gaussian (default)
;                 2-Lorentzian
;                 3-Exponential
;                 4-1/r^2
;                 5-1/r^3
;  /singlefit  Fit each line separately.  The default is to fit all
;                lines in a fiber (all3 chips at the same time).
;  /clobber    Overwrite any existing files.
;  =outdir     The output directory.
;  /pl         Plot the fits
;  /verbose    Verbose output to the screen.
;  /silent     Don't print anything to the screen
;  /stp        Stop at the end of the program
;
; OUTPUTS:
;  lsfim       A 2D array with the LSF parameters for each fiber
;  A file is also written to the lsf/ directory called apLSF-MJD5.fits
;
; USAGE:
;  IDL>aplsf,lampframes,waveid,lsfim
;
; By D.Nidever  March 2010
;-

t0_0 = systime(1)

; Get APOGEE directories
dirs=getdir(apogee_dir,cal_dir,spectro_dir,apogee_vers,lib_dir)
lsf_dir = cal_dir+'/lsf/'
linelist_dir = lib_dir+'skylines/'

apgundef,lsfstr

nlampframes = n_elements(lampframes)
nwaveid = n_elements(waveid)

; Not enough inputs
if nlampframes eq 0 then begin
  print,'Syntax - aplsf,lampframes,waveid,lsfim,psfid=psfid,fitmethod=fitmethod,pl=pl,'
  print,'               outdir=outdir,porder=porder,wporder=wporder,wproftype=wproftype,'
  print,'               clobber=clobber,verbose=verbose,silent=silent,stp=stp'
  return
endif

chiptag = ['a','b','c']

; Settings
if n_elements(nsigfit) eq 0 then nsigfit=8  ;8
if n_elements(verbose) eq 0 then verbose=1  ; verbose by default
if n_elements(forcepositive) eq 0 then forcepositive = 1  ; LSF must be positive

; if nsigfit is too small then the GH wings aren't constrained enough.

; Fitting METHOD
if not keyword_set(pl) then pl=0
if n_elements(fitmethod) eq 0 then fitmethod=2  ; /chipfit by default
if fitmethod lt 1 or fitmethod gt 3 then begin
  print,'FITMETHOD must be 1-3'
  return
endif
case fitmethod of
1: fitmethod_str='Fit each line separately'
2: fitmethod_str='Fit all lines in a chip together'
3: fitmethod_str='Fit all lines together'
endcase

; Allow the GH coefficients to vary as polynomials
;  this is the ORDER, so 0 means constant
;if n_elements(Porder) eq 0 then Porder = [1, 0, 1, 1, 1]
;if n_elements(Porder) eq 0 then Porder = [2, 1, 1, 1, 1,0]
if n_elements(Porder) eq 0 then Porder = [1, 1, 1, 1, 1,0]
nPorder = n_elements(Porder)
Horder = nPorder-1    ; The highest Hermite order
nGHcoefs = TOTAL(Porder+1)
if n_elements(WPorder) eq 0 then WPorder = [0,0]
if n_elements(Wproftype) eq 0 then Wproftype=1
nWpar = n_elements(WPorder)
nWcoefs = TOTAL(WPorder+1)

; Construct chip names
lampframeid1 = string(long(file_basename(lampframes[0])),format='(I08)')
files1 = file_dirname(lampframes[0])+'/'+dirs.prefix+'1D-'+['a','b','c']+'-'+lampframeid1+'.fits'
files = files1
if nlampframes gt 1 then begin
  lampframeid2 = string(long(file_basename(lampframes[1])),format='(I08)')
  files2 = file_dirname(lampframes[1])+'/'+dirs.prefix+'1D-'+['a','b','c']+'-'+lampframeid2+'.fits'
  PUSH,files,files2
endif
;files = [files1,files2]

; Get file info
info = APFILEINFO(files,/silent)

; Check the files
okay = (info.exists AND info.sp1dfmt AND info.allchips AND ((info.naxis eq 3) OR (info.exten eq 1)))
if total(okay) ne nlampframes*3 then begin
  print,'Files not okay'
  stop
  return
endif

; Check if the output files already exist
if n_elements(outdir) eq 0 then outdir=lsf_dir
outfile = outdir+dirs.prefix+'LSF-a-'+lampframeid1+'.fits'
if file_test(outfile) eq 1 and not keyword_set(clobber) then begin
  print,'Output files already exist and CLOBBER not set'
  return
endif


; Load the DATA
APLOADFRAME,files1,frame1,/exthead  ; loading frame 1

; Deal with NANs and ERR=0
for k=0,2 do begin
  bdnan = where(finite(frame1.(k).flux) eq 0 or finite(frame1.(k).err) eq 0,nbdnan)
  if nbdnan gt 0 then begin
    frame1.(k).flux[bdnan] = 0.0
    frame1.(k).err[bdnan] = baderr()
    frame1.(k).mask[bdnan] = 1
  endif

  bdzero = where(frame1.(k).err le 0,nbdzero)
  if nbdzero gt 0 then begin
    frame1.(k).flux[bdzero] = 0.0
    frame1.(k).err[bdzero] = baderr()
    frame1.(k).mask[bdzero] = 1
  endif
endfor

; Lamp type
lamptype = ''
if sxpar(frame1.(0).header,'LAMPUNE') eq 1 then lamptype='URANIUM'
if sxpar(frame1.(0).header,'LAMPTHAR') eq 1 then lamptype='THARNE'
if lamptype eq '' and info[0].exptype eq 'OBJECT' then lamptype='SKY'
lamptype = strupcase(lamptype)
if lamptype eq '' then begin
  print,'No LAMPTYPE for this exposure'
  return
endif
lamptype = strupcase(lamptype)


; Wavelength information
;-----------------------
;if tag_exist(frame1.(0),'WAVE') eq 0 then begin
;
;
;  ; Wavelength calibration file not input
;  if n_elements(waveid) eq 0 then begin
;    print,'Need a wavelength calibration file.  Wavelengths not in the ap1D file'
;    return
;  endif
;
;  ; Get the Wave file
;  ;wavefiles = wave_dir+dirs.prefix+'Wave-'+['a','b','c']+'-'+wavemjd5[0]+'.fits'
;  wavefiles = file_dirname(waveid)+'/'+dirs.prefix+'Wave-'+['a','b','c']+'-'+file_basename(waveid)+'.fits'
;  waveinfo = APFILEINFO(wavefiles,/silent)
;  okay = (waveinfo.exists)
;  if total(okay) ne 3 then begin
;    print,'Wave File not okay'
;    return
;  endif
;
;  ; Get the wavelength information
;  ;APLOADFRAME,wavefiles,waveframe
;  fits_read,wavefiles[0],wcoef1,head1,exten=1
;  fits_read,wavefiles[0],wim1,exten=2
;  fits_read,wavefiles[1],wcoef2,head2,exten=1
;  fits_read,wavefiles[1],wim2,exten=2
;  fits_read,wavefiles[2],wcoef3,head3,exten=1
;  fits_read,wavefiles[2],wim3,exten=2
;  waveframe = {chipa:{header:head1,coef:wcoef1,wave:wim1, wave2:0.*wim1},$
;               chipb:{header:head2,coef:wcoef2,wave:wim2, wave2:0.*wim2},$
;               chipc:{header:head3,coef:wcoef3,wave:wim3, wave2:0.*wim3}}
;  ; It looks like we only need these for the chip offsets/gaps
;
;; Use wavelengths from the ap1D file
;endif else begin
  waveframe = {chipa:{coef:frame1.(0).wcoef,wave:frame1.(0).wavelength, wave2:0.*frame1.(0).wavelength},$
               chipb:{coef:frame1.(1).wcoef,wave:frame1.(1).wavelength, wave2:0.*frame1.(1).wavelength},$
               chipc:{coef:frame1.(2).wcoef,wave:frame1.(2).wavelength, wave2:0.*frame1.(2).wavelength}}
;endelse



; Get the PSF information
if n_elements(psfid) gt 0 then begin

  ; Construct chip names
  psfframeid = string(long(file_basename(psfid)),format='(I08)')
  psffiles = file_dirname(psfid)+'/'+dirs.prefix+'PSF-'+chiptag+'-'+psfframeid+'.fits'
  tstr1 = MRDFITS(psffiles[0],1,/silent)
  tstr2 = MRDFITS(psffiles[1],1,/silent)
  tstr3 = MRDFITS(psffiles[2],1,/silent)
  tracestr = {tstr1:tstr1, tstr2:tstr2, tstr3:tstr3}

; PSFID not input, check the header
endif else begin

  lampfiles0 = file_dirname(lampframes[0])+'/'+dirs.prefix+'1D-'+chiptag+'-'+file_basename(lampframes[0])+'.fits'

  info = apfileinfo(lampfiles0,/silent)
  if min(info.exists) eq 0 then begin
    print,lampframes[0],' NOT FOUND'
    return
  endif
  head0 = headfits(lampfiles0[0])

  ; check the header for trace filename
  psffile = sxpar(head0,'PSFFILE',count=npsffile)
  if npsffile eq 0 then begin
    print,'CANNOT get PSFID from header.  Please input'
    return
  endif
  psfbase = file_basename(psffile,'.fits')
  dum = strsplit(psfbase,'-',/extract)
  psfid = file_dirname(psffile)+'/'+strtrim(dum[2],2)

  ; Construct chip names
  psfframeid = string(long(file_basename(psfid)),format='(I08)')
  psffiles = file_dirname(psfid)+'/'+dirs.prefix+'PSF-'+chiptag+'-'+psfframeid+'.fits'
  tstr1 = MRDFITS(psffiles[0],1,/silent)
  tstr2 = MRDFITS(psffiles[1],1,/silent)
  tstr3 = MRDFITS(psffiles[2],1,/silent)
  tracestr = {tstr1:tstr1, tstr2:tstr2, tstr3:tstr3}

endelse


print,'#############'
print,'Running APLSF'
print,'#############'
print,'Lamp Frame ID1 = ',lampframes[0],' LAMPTYPE = ',lamptype

print,'Porder = ',strtrim(porder,2)
print,'WPorder = ',strtrim(wporder,2)
print,'Wproftype = ',strtrim(wproftype,2)
print,'Fitmethod = ',strtrim(fitmethod,2),'  ',fitmethod_str

lsfid=file_basename(lampframes[0])

; DITHER COMBINING
;----------------------
if nlampframes gt 1 then begin

  ; Load the second frame
  print,'Lamp Frame ID2 = ',lampframes[1]
  APLOADFRAME,files2,frame2,/exthead  ; loading frame 1

  if abs(sxpar(frame1.chipa.header,'DITHPIX')-sxpar(frame2.chipa.header,'DITHPIX')) lt 0.01 then begin
    combframe = frame1
    for j=0,2 do begin
      combframe.(j).flux = frame1.(j).flux + frame2.(j).flux
      combframe.(j).err = sqrt(frame1.(j).err^2 + frame2.(j).err^2)
    endfor
    binsize = 1

  endif else begin
    stop,'dither combination for LSF not yet tested'

  ; Deal with NANs and ERR=0
  for k=0,2 do begin
    bdnan = where(finite(frame2.(k).flux) eq 0 or finite(frame2.(k).err) eq 0,nbdnan)
    if nbdnan gt 0 then begin
      frame2.(k).flux[bdnan] = 0.0
      frame2.(k).err[bdnan] = baderr()
      frame2.(k).mask[bdnan] = 1
    endif

    bdzero = where(frame2.(k).err le 0,nbdzero)
    if nbdzero gt 0 then begin
      frame2.(k).flux[bdzero] = 0.0
      frame2.(k).err[bdzero] = baderr()
      frame2.(k).mask[bdzero] = 1
    endif
  endfor

  ;------------------------------
  ; STEP I:  Measure Dither Shift
  ;------------------------------
  shiftstr = REPLICATE({index:-1L,framenum:'',shift:999999.0,shifterr:999999.0},2)
  ;print,'STEP 1: Measuring the DITHER SHIFT with APDITHERSHIFT'
  print,'Measuring the DITHER SHIFT with APDITHERSHIFT'
  APDITHERSHIFT,frame1,frame2,shift,shifterr,/xcorr
  ;shift = 0.50
  ;shifterr = 0.002
  shiftstr.index = [0,1]
  shiftstr.framenum = file_basename(lampframes)
  shiftstr.shift = [0.0, shift]
  shiftstr.shifterr = [0.0, shifterr]
  ;print,'UNCOMMENT THIS LINE!!!'

  ;--------------------------
  ; STEP II:  Dither Combine
  ;--------------------------
  ;print,'STEP 2: Combining DITHER FRAMES with APDITHERCOMB'
  print,'Combining DITHER FRAMES with APDITHERCOMB'
  plugmap = {fiberdata:replicate({objid:-1L,holetype:'OBJECT',objtype:'SKY',fiberid:0L},300)} ; make fake plugmap file
  ; fiberid=1 is at the top of the detector or index=299
  ; index = 300-fiberid
  ;plugmap.fiberdata.fiberid = indgen(300)+1
  plugmap.fiberdata.fiberid = 300-indgen(300)
  APDITHERCOMB,[frame1,frame2],shiftstr,plugmap,combframe,error=error
  if n_elements(error) ne 0 then begin
    print,'Error in combining the frames'
    return
  endif

  ;xmulti = 2       ; dither combined X units
  binsize = 2

; NO DITHER COMBINING
;--------------------
  endelse

Endif else begin
  combframe = frame1
  ;xmulti = 1        ; single frame X units
  binsize = 1
Endelse

sz = size(combframe.chipa.flux)
npix = sz[1]
nfibers = sz[2]

; Remove the continuum from the spectra and set lower error boundary
combframe0 = combframe
for i=0,2 do begin
  med = MEDFILT2D(combframe.(i).flux,101,dim=1,/edge)
  combframe.(i).flux -= med
  ;sig = MAD(combframe.(i).flux,/zero)
  ;mederr = median(combframe.(i).err)
  ;adderr = sqrt(sig^2 - mederr^2)  ; want total median err to be SIG
  ;combframe.(i).err = sqrt( combframe.(i).err^2 + adderr^2 )
end




;---------------------------
; STEP III:  Find the Lines
;---------------------------
;if keyword_set(gauss) or lamptype ne 'SKY' or fitmethod ne 2 then begin
if keyword_set(gauss) then begin

  print,'Find the Lines with APPEAKFIT'
  if keyword_set(gauss) then begin
    APPEAKFIT,combframe.chipa,linestr1,fibers=[10,150,290],nsigthresh=10
    APPEAKFIT,combframe.chipb,linestr2,fibers=[10,150,290],nsigthresh=10
    APPEAKFIT,combframe.chipc,linestr3,fibers=[10,150,290],nsigthresh=10
  endif else begin
    APPEAKFIT,combframe.chipa,linestr1
    APPEAKFIT,combframe.chipb,linestr2
    APPEAKFIT,combframe.chipc,linestr3
  endelse
  ;restore,'aplsf_linestr.dat'
  print,'UNCOMMENT THESE LINES!!!'
  ADD_TAG,linestr1,'CHIP',1,linestr1
  ADD_TAG,linestr2,'CHIP',2,linestr2
  ADD_TAG,linestr3,'CHIP',3,linestr3
  linestr = [linestr1,linestr2,linestr3]


  ; Add wavelengths for the lines
  ADD_TAG,linestr,'WAVE',0.0d0,linestr
  ADD_TAG,linestr,'X',0.0,linestr
  ADD_TAG,linestr,'DISP',0.0,linestr

  ind1 = where(linestr.chip eq 1,nind1)
  x1 = linestr[ind1].gaussx-1023.5-2048-150
  w1 = dblarr(nind1)
  disp1 = dblarr(nind1)
  for i=0,nind1-1 do begin
    w1[i] = pix2wave(linestr[ind1[i]].gaussx,waveframe.(0).coef[linestr[ind1[i]].fiber,*])
    x1[i] = linestr[ind1[i]].gaussx + waveframe.(0).coef[linestr[ind1[i]].fiber,0]
    disp1[i] = pix2wave(linestr[ind1[i]].gaussx+0.5,waveframe.(0).coef[linestr[ind1[i]].fiber,*])-$
               pix2wave(linestr[ind1[i]].gaussx-0.5,waveframe.(0).coef[linestr[ind1[i]].fiber,*])
  end
  linestr[ind1].wave = w1
  linestr[ind1].x = x1
  linestr[ind1].disp = disp1

  ind2 = where(linestr.chip eq 2,nind2)
  x2 = linestr[ind2].gaussx-1023.5
  w2 = dblarr(nind2)
  disp2 = dblarr(nind2)
  for i=0,nind2-1 do begin
    w2[i] = pix2wave(linestr[ind2[i]].gaussx,waveframe.(1).coef[linestr[ind2[i]].fiber,*])
    x2[i] = linestr[ind2[i]].gaussx + waveframe.(1).coef[linestr[ind2[i]].fiber,0]
    disp2[i] = pix2wave(linestr[ind2[i]].gaussx+0.5,waveframe.(1).coef[linestr[ind2[i]].fiber,*])-$
               pix2wave(linestr[ind2[i]].gaussx-0.5,waveframe.(1).coef[linestr[ind2[i]].fiber,*])
  end
  linestr[ind2].wave = w2
  linestr[ind2].x = x2
  linestr[ind2].disp = disp2

  ind3 = where(linestr.chip eq 3,nind3)
  x3 = linestr[ind3].gaussx-1023.5+2048+150
  w3 = dblarr(nind3)
  disp3 = dblarr(nind3)
  for i=0,nind3-1 do begin
    w3[i] = pix2wave(linestr[ind3[i]].gaussx,waveframe.(2).coef[linestr[ind3[i]].fiber,*])
    x3[i] = linestr[ind3[i]].gaussx + waveframe.(2).coef[linestr[ind3[i]].fiber,0]
    disp3[i] = pix2wave(linestr[ind3[i]].gaussx+0.5,waveframe.(2).coef[linestr[ind3[i]].fiber,*])-$
               pix2wave(linestr[ind3[i]].gaussx-0.5,waveframe.(2).coef[linestr[ind3[i]].fiber,*])
  end
  linestr[ind3].wave = w3
  linestr[ind3].x = x3
  linestr[ind3].disp = disp3

  ; For Uranium remove these lines that are blended
  if lamptype eq 'URANIUM' then begin
    badwave = [15295.5, 15369.9, 15590.15, 15190.9, 15687.6d0, 16045.6, 16061.0, 16540.6, 16667.9]
    ; with the height thresholds we only lose about 2 lines in the blue
    for i=0,n_elements(badwave)-1 do begin
      bd = where(abs(linestr.wave-badwave[i]) lt 1.0,nbd)
      if nbd gt 0 then REMOVE,bd,linestr
    end
  endif
  ; Remove BAD lines using a polynomial fit to the SIGMA values
  if keyword_set(pl) then begin
    set_plot,'PS'
    device,file='aplsf.eps',/encap,/color
    plot,linestr.x,linestr.gpar[2],ps=1,xtit='X',ytit='Gaussian Sigma (pixels)',tit='Removing Bad Lines'
  endif

  if keyword_set(gauss) then begin
    if keyword_set(pl) then smcolor
    outfile = lsf_dir+dirs.prefix+'LSF-'+lampframeid1+'.dat'
    openw,lun,/get_lun,outfile
    ; initialize splines for wave2pix
    for ichip=0,2 do $
      for ifiber=0,299 do waveframe.(ichip).wave2[*,ifiber] = spl_init(reverse(waveframe.(ichip).wave[*,ifiber]),reverse(indgen(n_elements(waveframe.(ichip).wave[*,ifiber]))))

    if lamptype eq 'SKY' then begin
      ;readcol,linelist_dir+'airglow.new',iair,airw,airf,airs,airsep,airname,airuse,format='(i,f,f,i,f,a,x,i)'
      readcol,linelist_dir+'airglow.txt',iair,airw,airf,airs,airsep,airname,airuse,format='(i,x,x,f,f,i,f,a,i)'
      for i=0,n_elements(airw)-1 do begin
        if (airuse[i]) then begin
          ;good=where(abs(linestr.wave-airw[i]) lt 1 and linestr.height gt 2000)
          good=where(abs(linestr.wave-airw[i]) lt 1)
          if good[0] ge 0 then begin
            if keyword_set(pl) then oplot,[linestr[good].x],[linestr[good].gpar[2]],color=3,ps=1
            g150=where(linestr[good].fiber eq 150 or linestr[good].fiber eq 10 or linestr[good].fiber eq 290)
            for j=0,n_elements(g150)-1 do begin
              k=good[g150[j]]
              fwhm=2.354*linestr[k].gpar[2]
              ichip=linestr[k].chip-1
              ifiber=linestr[k].fiber
              xpred=spl_interp(reverse(waveframe.(ichip).wave[*,ifiber]),reverse(indgen(n_elements(waveframe.(ichip).wave[*,ifiber]))),waveframe.(ichip).wave2[*,ifiber],airw[i])
              printf,lun,getmjd5(combframe.chipa.header),linestr[k].wave,linestr[k].x,fwhm,fwhm*linestr[k].disp,$
                     linestr[k].gaussx,airw[i],xpred,linestr[k].chip,linestr[k].fiber,linestr[k].height,airw[i],lamptype,$
                     format='(i8,f12.2,6f10.2,i6,i6,f12.1,f12.2,2x,a)'
            endfor
          endif
        endif
      endfor
    endif else begin
      good=where(linestr.height gt 2000)
      if good[0] ge 0 then begin
        if keyword_set(pl) then oplot,[linestr[good].x],[linestr[good].gpar[2]],color=3,ps=1
        g150=where(linestr[good].fiber eq 150 or linestr[good].fiber eq 10 or linestr[good].fiber eq 290)
        for j=0,n_elements(g150)-1 do begin
          k=good[g150[j]]
          fwhm=2.354*linestr[k].gpar[2]
          printf,lun,getmjd5(combframe.chipa.header),linestr[k].wave,linestr[k].x,fwhm,fwhm*linestr[k].disp,$
                 linestr[k].gaussx,linestr[k].chip,linestr[k].fiber,linestr[k].height,lamptype,format='(i8,f12.2,4f10.2,i6,i6,f12.1,2x,a)'
        endfor
      endif
    endelse
    free_lun,lun
    return
  endif

  si = sort(linestr.x)
  sigcoef = robust_poly_fit(linestr[si].x,linestr[si].gpar[2],5)
  xx = scale_vector(findgen(1000),-3222,3222)
  if keyword_set(pl) then oplot,xx,poly(xx,sigcoef),co=250
  std = MAD(linestr.gpar[2]-poly(linestr.x,sigcoef),/zero)
  gd = where(abs(linestr.gpar[2]-poly(linestr.x,sigcoef)) lt 4*std,ngd)
  if keyword_set(pl) then oplot,linestr[gd].x,linestr[gd].gpar[2],ps=1,co=150
  linestr0 = linestr
  linestr = linestr[gd]

endif

; Load the airglow lines
if lamptype eq 'SKY' then begin
  airstr = IMPORTASCII(linelist_dir+'airglow.txt',/header,/silent)
  nairstr = n_elements(airstr)
endif
if lamptype eq 'THARNE' then begin
  airstr = IMPORTASCII(linelist_dir+'tharne.lines.vac.apogee',/header,/silent)
  nairstr = n_elements(lampstr)
  pl=1
endif


;------------------------
; STEP IV:  Fit the LSF
;------------------------
print,'Fit the LSF'
;sz = size(combframe.chipa.flux)
;npix = sz[1]
;nfibers = sz[2]
nlsfpar = 3+nPorder+TOTAL(Porder+1)
if nWpar gt 0 then nlsfpar += 2+nWpar+TOTAL(WPorder+1)
lsfim1 = fltarr(nfibers,nlsfpar)
lsfim2 = fltarr(nfibers,nlsfpar)
lsfim3 = fltarr(nfibers,nlsfpar)


apgundef,allfitstr


if keyword_set(verbose) then begin
  print,'--------------------------------------------------------------------'
  case fitmethod of
  1:  print,'Fiber Nlines Chisq  Sigma     H0      H1      H2      H3      H4 '
  2:  print,'Fiber Chip  Nlines Chisq  Sigma     H0      H1      H2      H3      H4 '
  3:  print,'Fiber Nlines Chisq  Sigma     H0      H1      H2      H3      H4 '
  else:
  endcase
  print,'--------------------------------------------------------------------'
endif

if fitmethod eq 2 then begin
  ntotcoef = total(Porder+1)
  if nWpar gt 0 then ntotcoef+=total(WPorder+1)
  ;nlsfpars = 3+n_elements(Porder)+total(Porder+1)
  ;if nWpar gt 0 then nlsfpars += 2+n_elements(WPorder)+total(WPorder+1)
  fitstr = REPLICATE({fiber:0L,chip:0L,xoffset:0.0,par:dblarr(ntotcoef),perror:dblarr(ntotcoef),ghcoefs:fltarr(total(Porder+1)),$
                      wcoefs:fltarr(total(WPorder+1)),lsfpars:fltarr(nlsfpar),status:0L,chisq:1.e10,rchisq:1.e10,npts:0L,rms:0.0},nfibers*3)
  apgundef,linestr
endif


if ~keyword_set(fibers) then ifibers=indgen(nfibers) else ifibers=fibers

; Loop through the Fibers
;For i=0,nfibers-1 do begin
for ii=0,n_elements(ifibers)-1 do begin
  i=ifibers[ii]

  if not keyword_set(verbose) then print,'Fiber ',strtrim(i+1,2)

  apgundef,fiberlinestr

  ; FITTING METHOD
  ;-----------------
  CASE FITMETHOD OF

  ;========================================================================
  ; ALL LINES FIT SEPARATELY
  ;========================================================================
  ;  This does NOT find a solution can be used to make the apLSF files
  ;  This is only useful for seeing what the GH parameters are for
  ;   individual lines
  1: BEGIN


    spec = fltarr(npix*3)
    errspec = fltarr(npix*3)

    ; The measured chip gaps
    ;  want them in dither combined units, multiply by 2
    ;xoff1 = waveframe.(1).data[i,0] * xmulti  ; Xoffset for chip b
    ;xoff2 = waveframe.(2).data[i,0] * xmulti  ; Xoffset for chip c
    xoff1 = waveframe.(0).coef[i,0] * binsize  ; Xoffset for chip a
    xoff2 = waveframe.(1).coef[i,0] * binsize  ; Xoffset for chip b
    xoff3 = waveframe.(2).coef[i,0] * binsize  ; Xoffset for chip c
    xoff = [xoff1, xoff2, xoff3]

    ; Construct the Y array
    x = fltarr(npix*3)
    x[0:npix-1] = findgen(npix)+xoff1
    x[npix:2*npix-1] = findgen(npix)+xoff2
    x[2*npix:3*npix-1] = findgen(npix)+xoff3

    ; Loop through the Chips
    For j=0,2 do begin

      ;fiber = combframe.(j).data[*,i,0]
      ;var = combframe.(j).data[*,i,1]
      fiber = combframe.(j).flux[*,i]
      err = combframe.(j).err[*,i]

      ; Remove the continuum
      ;cont = MEDFILT1D(fiber,101,/edge)
      ;fiber -= cont

      chiplinestr = (SCOPE_VARFETCH('linestr'+strtrim(j+1,2)))
      ind = where(chiplinestr.fiber eq i,nind)
      if nind eq 0 then goto,BOMB1
      ilinestr = chiplinestr[ind]

      ilinestr.peakx += xoff[j]
      ilinestr.gaussx += xoff[j]
      ilinestr.gpar[1] += xoff[j]
      ilinestr.gpar0[1] += xoff[j]

      spec[j*npix:(j+1)*npix-1] = fiber
      errspec[j*npix:(j+1)*npix-1] = err
      ;spec[yoff:yoff+npix-1] = fiber
      ;errspec[yoff:yoff+npix-1] = sqrt(var)

      ; Add to fiber linelist
      PUSH,fiberlinestr,ilinestr

      BOMB1:
    End
    nlines = n_elements(fiberlinestr)
    nspec = n_elements(spec)
    ;y = findgen(nspec)

    ; gausshermitebin.pro works for ONE line
    ; need a function like GHM_all_1d_decom  that runs it on all lines

    ; Only fit pixels that are close to a line, within +/-10 sigma or so.

    resid = spec  ; we will subtract the fitted lines from this as we go

    fitstr = REPLICATE({fiber:0L,chip:0L,line:0L,par:fltarr(8),perror:fltarr(8),flux:0.0,center:0.0,ghpar:fltarr(nPorder),$
                        status:0L,chisq:0.0,rchisq:0.0,npts:0L,rms:0.0},nlines)
    For j=0,nlines-1 do begin

      ; Loading info into the structure
      fitstr[j].fiber = i
      fitstr[j].chip = fiberlinestr[j].chip
      fitstr[j].line = j

      ; Select pixels near this line
      gd = where(abs(x-fiberlinestr[j].gaussx) lt nsigfit*fiberlinestr[j].gpar[2],ngd)

      ; Check for contamination in the wings
      ;  if so, reduce to +/-2.5 sigma.
      inner = where(abs(x-fiberlinestr[j].gaussx) le 2.5*fiberlinestr[j].gpar[2],ninner)
      outer = where(abs(x-fiberlinestr[j].gaussx) lt nsigfit*fiberlinestr[j].gpar[2] AND $
                    abs(x-fiberlinestr[j].gaussx) gt 2.5*fiberlinestr[j].gpar[2],nouter)
      totinner = TOTAL(resid[inner])
      totouter = TOTAL(resid[outer])
      ;if totouter/totinner*100 gt 5 or totouter gt 3*TOTAL(errspec[outer]) then begin
      if  totouter/totinner*100 gt 5 and totouter gt 5*TOTAL(errspec[outer]) then begin
        gd = where(abs(x-fiberlinestr[j].gaussx) le 2.5*fiberlinestr[j].gpar[2],ngd)
      endif

      ;;  6 Gauss-Hermitze LSF parameters, the height and center aren't neded
      ;initpar = [fiberlinestr[j].gpar[0:2], 1.0, fltarr(4)]
      ;npars = n_elements(initpar)
      ;parinfo = REPLICATE({limited:[0,0],limits:[0.0,0.0],fixed:0},npars)
      ;binsize = 2
      ;fa = {binsize:binsize}
      ;initpar[0] = initpar[0]/binsize

      xin = x[gd]
      ;specin = spec[gd]
      specin = resid[gd]
      errspecin = errspec[gd]

      Porder_single = Porder*0  ; all constant
      nPorder = n_elements(Porder_single)
      ;binsize = xmulti  ;2
      ;Xglobalcenter = 6430
      Xglobalcenter = 0
      WPorder_single = WPorder*0  ; all constant

      initpar = dblarr(2+1+TOTAL(Porder_single+1)+TOTAL(WPorder_single+1))
      ; We need to multiply the heights by sig*sqrt(2*pi) because the LSFs
      ;  are normalized but the Gaussians that APPEAKFIT.PRO have an
      ;  area under the curve of ht*sig*sqrt(2*pi).
      initpar[0] = fiberlinestr[j].gpar[0] * fiberlinestr[j].gpar[2] * sqrt(2*!dpi)   ; flux
      initpar[1] = fiberlinestr[j].gpar[1] ; -Xoffset        ; center
      initpar[2] = -Xglobalcenter   ;Xoffset
      initpar[3] = fiberlinestr[j].gpar[2]
      npars = n_elements(initpar)

      parinfo = REPLICATE({limited:[0,0],limits:[0.0,0.0],fixed:0},npars)
      parinfo[0].limited=1 & parinfo[0].limits=[0.9,1.1]*initpar[0]
      parinfo[1].limited=1 & parinfo[1].limits=[-1,1]+initpar[1]
      parinfo[2].fixed = 1  ; Keep X-center FIXED
      parinfo[3].limited=1 & parinfo[3].limits=[0.8,1.2]*initpar[3]
      parinfo[4].limited=1 & parinfo[4].limits=[-5,5]
      parinfo[5].limited=1 & parinfo[5].limits=[-5,5]
      parinfo[6].limited=1 & parinfo[6].limits=[-5,5]
      parinfo[7].limited=1 & parinfo[7].limits=[-5,5]

      ; kludge, keeping sigma=1 fixed
      ;initpar[3] = 1.0
      ;parinfo[3].fixed = 1
      ;parinfo[3].limited=0

      loarr = 0
      hiarr = ngd-1
      fa = {binsize:binsize,nlines:1,loarr:loarr,hiarr:hiarr,porder:porder_single,$
            wproftype:wproftype,wporder:wporder_single}

      ; Initial guess fit
      ;yfit0 = fit_lsf_gh(yin,initpar,binsize=binsize,nlines=1,loarr=loarr,hiarr=hiarr,porder=porder_single)

      ; Maybe I should constrain the center/width a bit more.

      ; Fit the parameters
      par = MPFITFUN('fit_lsf_gh',xin,specin,errspecin,initpar,$
                     parinfo=parinfo,yfit=yfit,dof=dof,status=status,bestnorm=chisq,$
                     perror=perror,functargs=fa,niter=niter,/quiet)
      if status lt 1 then begin
        if keyword_set(verbose) then print,'Fitting problems'
        fitstr[j].status = status
        fitstr[j].rchisq = 999999.
        goto,BOMB
      endif

      pcerror = perror * sqrt(chisq/dof)  ; scaled uncertainties
      rms = sqrt(mean((specin-yfit)^2))
      rchisq = chisq/dof
      if finite(rchisq) eq 0 then rchisq=999999.

      ; Remove line from spectrum
      resid_withline = resid  ; resid with this line
      resid[gd] -= yfit

      ; Stuff the information in the structure
      fitstr[j].par = par
      fitstr[j].perror = pcerror
      fitstr[j].flux = par[0]
      fitstr[j].center = par[1]
      fitstr[j].ghpar = par[3:*]
      fitstr[j].status = status
      fitstr[j].chisq = chisq
      fitstr[j].rchisq = rchisq
      fitstr[j].npts = ngd
      fitstr[j].rms = rms

      ; Print the parameters
      if keyword_set(verbose) then begin
         fmt = '(I4,F9.4,F10.1,F10.3,'+strtrim(nPorder,2)+'G10.3)'
         print,format=fmt,i+1,rchisq,par[0],par[1],par[3:*]
      endif

      ; Plotting
      if keyword_set(pl) then begin
        xr = minmax(xin)
        yr = minmax([specin,yfit])
        yr = [yr[0]-range(yr)*0.15,yr[1]+range(yr)*0.1]
        plot,xin,specin,tit='Fiber '+strtrim(i+1,2)+' Line '+strtrim(j+1,2),xr=xr,yr=yr,xs=1,ys=1
        oplot,xr,[0,0],linestyle=2
        oplot,xin,yfit,co=250
        off = range([specin,yfit])*0.075
        oplot,xin,specin-yfit-off,co=150
        oplot,xr,[0,0]-off,linestyle=2
        legend,['Data','Model','Residuals'],textcolor=[255,250,150],/top,/left,charsize=1.2
        xyouts,mean(xr),yr[1]-0.05*range(yr),'Chisq = '+strtrim(rchisq,2),align=0.5,charsize=1.2
        wait,0.3
      endif

      ;if ngd lt 15 then stop
      stop

      BOMB:

    End  ; line loop

    ; Add to the final fit structure for all lines
    PUSH,allfitstr,fitstr

    if i eq 50 then stop

    ;stop


  END  ; all lines fit separately



  ;========================================================================
  ; FIT ALL LINES PER CHIP TOGETHER
  ;========================================================================
  2: BEGIN

    ; THIS NOW USES A MODEL OF THE LINES AND FITS IT TO THE SPECTRUM
    ; ITSELF.  IT'S NOT NECESSARY TO FIND THE LINES BEFOREHAND.
    ; THIS CURRENTLY ONLY WORKS FOR SKY/AIRGLOW LINES/EXPOSURES.

    if lamptype ne 'SKY' then begin
      print,'This only works for SKY at the moment, continue at your peril...'
      ;stop
      ;return
    endif

    set_plot,'ps'
    file_mkdir,outdir+'/html'
    file_mkdir,outdir+'/plots'
    ; Loop through the Chips
    For j=0,2 do begin

      spec = combframe.(j).flux[*,i]
      errspec = combframe.(j).err[*,i]
      maskspec = combframe.(j).mask[*,i]
      wave = waveframe.(j).wave[*,i]
      wcoef = waveframe.(j).coef[i,*]
      modelflux = spec*0.0

      ; set these here in case we bomb for this fiber and need to use interpolated parameters
      fitstr[3*i+j].fiber = i
      fitstr[3*i+j].chip = j+1
      fitstr[3*i+j].lsfpars[0] = binsize
      fitstr[3*i+j].lsfpars[1] = npix/2

      ; Fix NANs
      bd = where(finite(spec) eq 0,nbd)
      if nbd gt 0 then begin
        spec[bd] = 0.0
        errspec[bd] = baderr()
      endif

      ;; Fix very low errspec values
      ;mederr = median(errspec)
      ;sigerr = mad(errspec)
      ;bderr = where(errspec lt mederr-2*sigerr,nbderr)
      ;if nbderr gt 0 then errspec[bderr]=mederr

      ; Increase the noise in the error spectrum
      ;   NOT SURE THIS IS NEEDED ANYMORE!!!
      diff = shift(spec,1)-spec
      sig = MAD(diff,/zero)
      ;sig = MAD(diff,/zero)/sqrt(2.)  ; correct for added noise from sub
      mederr = median(errspec)
      adderr = sqrt((sig^2 - mederr^2) > 0)  ; want total median err to be SIG
      if adderr gt 5 then $
        errspec = sqrt( errspec^2 + adderr^2 )
      ;print,'ADDERR=',adderr  

      ; Do a better job of removing the continuum
      med0 = MEDFILT1D(spec,101,/edge)
      temp = spec
      sig1 = MAD(temp)
      mask = long(abs(temp-med0) gt 2*sig1)
      mask = convol(mask,indgen(5))
      mask = mask/(mask>1)
      bd1 = where(mask eq 1,nbd1)
      if nbd1 gt 0 then temp[bd1]=!values.f_nan
      med1 = MEDFILT1D(temp,101,/edge)

      temp2 = spec
      mask2 = long(abs(temp2-med1) gt 2*sig1)
      mask2 = convol(mask2,indgen(5))
      mask2 = mask2/(mask2>1)
      bd2 = where(mask2 eq 1,nbd2)
      if nbd2 gt 0 then temp2[bd2]=!values.f_nan
      med2 = MEDFILT1D(temp2,101,/edge)

      spec -= med2  ; remove the continuum

      ; Increase the error spectrum if necessary
      sig = MAD(spec,/zero)
      mederr = median(errspec)
      adderr = sqrt(sig^2 - mederr^2)  ; want total median err to be SIG
      ;if adderr gt 0.5*mederr then errspec=sqrt(errspec^2 + adderr^2)
      if mederr gt sig*3 then begin
        print,'The error spectrum does not match the NOISE'
        goto,BOMB2
        stop
      endif

      ;; chip-dependent Porder
      ;case j of
      ;0: Porder = [2,1,1,1,1]
      ;1: Porder = [1,1,1,1,1]
      ;2: Porder = [1,1,1,1,1]
      ;endcase
      ;nPorder = n_elements(Porder)
      ;Horder = nPorder-1 


      ; Get all AIRGLOW lines for this chip
      gdlines = where(airstr.wave ge min(wave) and airstr.wave le max(wave) and airstr.uselsf eq 1,ngdlines)

      
      ; Start the CHIPLINESTR structure
      chiplinestr = REPLICATE({fiber:0L,chip:0L,x:0.0,model_id:0L,model_wave:0.0d0,model_emission:0.0,$
                               model_doublet:0L,model_dbl_wsep:0.0d0,$
                               model_dbl_xsep:0.0d0,model_type:'',lsffit_pars:dblarr(3),$
                               lsffit_perror:dblarr(3),lsffit_wave:0.0d0,lsffit_flux:0.0,lsffit_chisq:0.0,$
                               lsffit_status:0L,lsffit_ghcoefs:dblarr(nghcoefs+nwcoefs)},ngdlines)
      chiplinestr.fiber = i
      chiplinestr.chip = j
      chiplinestr.model_id = airstr[gdlines].id
      chiplinestr.model_wave = airstr[gdlines].wave
      chiplinestr.model_type = airstr[gdlines].name
      if tag_exist(airstr,'EMISSION') then chiplinestr.model_emission=airstr[gdlines].emission
      if tag_exist(airstr,'DOUBLET') then chiplinestr.model_doublet=airstr[gdlines].doublet
      if tag_exist(airstr,'DBL_WSEP') then chiplinestr.model_dbl_wsep=airstr[gdlines].dbl_wsep

      ; Get X-values for the airglow lines
      ;ADD_TAG,chiplinestr,'X',0.0,chiplinestr
      si = sort(wave)
      x = findgen(npix)
      for k=0,ngdlines-1 do chiplinestr[k].x = spline(wave[si],x[si],chiplinestr[k].model_wave)

      ; Convert DOUBLET WSEP to pixels
      dwall = abs(wave[1:*]-wave[0:npix-2])
      dblind = where(chiplinestr.model_doublet eq 1,ndblind)
      for k=0,ndblind-1 do begin
        ind1 = dblind[k]
        dw1 = dwall[round(chiplinestr[ind1].x)]
        chiplinestr[ind1].model_dbl_xsep = chiplinestr[ind1].model_dbl_wsep / dw1
      end

      ; Only keep lines with decent flux
     ; maxflux = fltarr(ngdlines)
     ; for k=0,ngdlines-1 do maxflux[k] = max(spec[ (round(chiplinestr[k].x)-1)>0:(round(chiplinestr[k].x)+1)<(npix-1)])
     ; case j of
     ;  0: thresh=200
     ;  1: thresh=500
     ;  2: thresh=300
     ; endcase
     ; gdlines = where(maxflux gt thresh,ngdlines)
     ; if ngdlines lt 5 then gdlines = where(maxflux gt 2*thresh,ngdlines)
     ; if ngdlines lt 5 then gdlines = where(maxflux gt 4*thresh,ngdlines)
     ; chiplinestr = chiplinestr[gdlines]

      ; Need to fit all lines together
      ;-------------------------------
      nsigfit = 8
      nsigfit = 6

      ; Get all pixels close to a line
      mask = intarr(npix)
      xloarr = lonarr(ngdlines)
      xhiarr = lonarr(ngdlines)
      for k=0,ngdlines-1 do begin
        ;gd = where(abs(x-chiplinestr[k].x) lt nsigfit*1.0,ngd)
        if chiplinestr[k].model_doublet eq 0 then begin
          gd = where(abs(x-chiplinestr[k].x) lt nsigfit*2.0,ngd)
        endif else begin
          gd = where( x ge (chiplinestr[k].x-0.5*chiplinestr[k].model_dbl_xsep-nsigfit*2.0) and $
                      x le (chiplinestr[k].x+0.5*chiplinestr[k].model_dbl_xsep+nsigfit*2.0),ngd)
        endelse
        xlo = min(gd)
        xhi = max(gd)
        mask[xlo:xhi] = 1
        xloarr[k] = xlo
        xhiarr[k] = xhi

        ; uncomment to see individual lines
        if keyword_set(pl) and pl gt 1 then begin
          plot,x[gd],spec[gd]
          peak=max(spec[gd],imax)
          initpar = [peak, x[gd[imax]], 1.30]
          yfit = mpfitpeak(x[gd],spec[gd],par,nterms=3,estimates=initpar,/gaussian,/positive)
          fwhm = 2.35482 * par[2]  ; FWHM in pixels
          print,(wave[xhi]+wave[xlo])/2.,(wave[xhi]+wave[xlo])/2./(fwhm*(wave[xhi]-wave[xlo])/(n_elements(gd)-1.))
          oplot,x[gd],yfit
          legend,[string(format='(f6.0)',(wave[xhi]+wave[xlo])/2.),'R = '+string(format='(f6.0)',(wave[xhi]+wave[xlo])/2./(fwhm*(wave[xlo]-wave[xhi])/(n_elements(gd)-1.)))]
        stop
        endif
      end
      useind = where(mask eq 1,nuseind)

      ; Input values
      xin = x[useind]
      specin = spec[useind]
      errspecin = errspec[useind]

      bd=where((maskspec[useind] and badmask()) gt 0,nbd)
      if float(nbd)/n_elements(useind) gt 0.5 then goto,BOMB2

      ; "Downweight" pixel with low S/N
      ;snr_thresh = 2 ;5 ; 1
      ;bdpix = where(specin/(errspecin>1) lt snr_thresh,nbdpix)
      ;if nbdpix gt 0 then errspecin[bdpix] = ( errspecin[bdpix]*5 < 1e30 )  ; 10


      ; lo/hi indices for the "input" arrays
      npixarr = xhiarr-xloarr+1
      loarr = lonarr(ngdlines)
      hiarr = lonarr(ngdlines)
      ;hiarr[0] = npixarr[0]-1
      for k=0,ngdlines-1 do begin
        loarr[k] = where(useind eq xloarr[k])
        hiarr[k] = where(useind eq xhiarr[k])
      end

      ; Allow the GH coefficients to vary as polynomials
      ;  this is the ORDER, so 0 means constant
      ;porder = [2, 2, 2, 2, 2, 2]
      ;Porder = [0, 0, 0, 0, 0]
      nPorder = n_elements(Porder)

      ; Now fit the LSF
      ;  the flux of each line is allowed to vary
      ; The parameters are:
      ;  Nline flux and Nline center parameters (these are interleaved
      ;     e.g. flux1, cen1, flux1, cen2, etc.)
      ;  6 Gauss-Hermitze LSF parameters, the height and center aren't needed

      ; The parameters are:
      ;  Nlines   flux/center for the Nlines
      ;  X0       Xoffset
      ;  GHcoefs  The Gauss-Hermite polynomial variation coefficients
      ;binsize = 2
      ;binsize = xmulti

      ; The "global" center
      ; I want the coefficients that vary globally with Y to be fit around
      ; the center of the array, so that the constant term is the value
      ; at the center of the array, near 6430 from the left edge.
      ;Xglobalcenter = 6430
      Xglobalcenter = npix/2

      ; ****GHcoefs means GH and Wing coefficients here****

      ; Initial parameters
      initpar = dblarr(3*ngdlines+1+nGHcoefs+nWcoefs)
      ; We need to multiply the heights by sig*sqrt(2*pi) because the LSFs
      ;  are normalized but the Gaussians that APPEAKFIT.PRO have an
      ;  area under the curve of ht*sig*sqrt(2*pi).
      initpar[0:3*ngdlines-3:3] = (spec[chiplinestr.x]>1) * 1.0 * sqrt(2*!dpi)   ; flux
      initpar[1:3*ngdlines-2:3] = chiplinestr.x
      initpar[2:3*ngdlines-1:3] = 0.0      ;            constant Y-offset
      initpar[3*ngdlines] = Xglobalcenter
      ;initpar[3*ngdlines+1:*] = GHcoefs
      ;initpar[3*ngdlines+1:3*ngdlines+nGHcoefs] = GHcoefs

      ;initpar[3*nlines+1] = median(chiplinestr.gpar[2])  ; sigma

      ; Start with GHcoefs from previous fibers
      prevind = where(fitstr.chip eq (j+1) and fitstr.fiber lt i,nprevind)
      if nprevind gt 0 then begin
        init_GHcoefs = fitstr[0].par*0
        case 1 of
        ; Only one, just use it
        (nprevind eq 1): init_par=fitstr[prevind].par
        ; Two, take mean
        (nprevind eq 2): init_par=TOTAL(fitstr[prevind].par,2)/float(nprevind)
        ; Three or more, robust linear fit
        (nprevind ge 3): begin
          for k=0,n_elements(init_par)-1 do begin
            med_init_par1 = MEDIAN(fitstr[prevind].par[k],/even)
            init_par[k] = med_init_par1
          end
        end
        else: stop,'This should never happen'
        endcase
        initpar[3*ngdlines+1:*] = init_par

      ; First fiber
      endif else begin

        initpar[3*ngdlines+1] = 1.5  ; sigma

        ;Wcoefs = [0.03, 3.5]
        Wcoefs = fltarr(nWcoefs)
        Wcoefs[0] = 0.03
        Wcoefs[WPorder[0]+1] = 3.5
        initpar[3*ngdlines+nGHcoefs+1:3*ngdlines+nGHcoefs+nWcoefs] = Wcoefs

      endelse

      ;if nWcoefs gt 0 then $
      ;  initpar[3*ngdlines+nGHcoefs+1:3*ngdlines+nGHcoefs+nWcoefs] = Wcoefs

      npars = n_elements(initpar)
      parinfo = REPLICATE({limited:[0,0],limits:[0.0,0.0],fixed:0},npars)
      parinfo[3*ngdlines:*].fixed = 1  ; Keep all LSF params fixed
      parinfo[0:3*ngdlines-3:3].limited = [1,0] ; flux must be positive
      parinfo[0:3*ngdlines-3:3].limits[0] = 0.0
      parinfo[1:3*ngdlines-2:3].fixed = 1       ; keep x-values fixed
      parinfo[2:3*ngdlines-1:3].fixed = 1       ; keep y-offsets fixed
      ;parinfo[2:3*ngdlines-1:3].limited = 1
      ;parinfo[2:3*ngdlines-1:3].limits[0] = -1.0*mad(fiberflux)   ; can't move too much
      ;parinfo[2:3*ngdlines-1:3].limits[1] =  1.0*mad(fiberflux)   ; can't move too much
      fa = {x:xin,y:specin,err:errspecin,binsize:binsize,nlines:ngdlines,loarr:loarr,$
            hiarr:hiarr,porder:porder,wproftype:wproftype,wporder:wporder,$
            doublet:chiplinestr.model_doublet,dbl_sep:chiplinestr.model_dbl_xsep}

      ; Initial guess fit
      yfit0 = SKYFIT_LSF_GH(xin,initpar,_extra=fa)

      ; First fit, only heights allowed to vary
      ;----------------------------------------
      initpar1 = initpar
      parinfo1 = parinfo
      t0 = systime(1)
      par1 = MPFIT('skyfit_lsf_gh_dev',initpar1,$
                    parinfo=parinfo1,dof=dof,status=status1,bestnorm=chisq1,$
                    perror=perror1,functargs=fa,niter=niter1,autoderivative=0,maxiter=3,/quiet)
      yfit1 = skyfit_lsf_gh(xin,par1,_extra=fa)
      t1 = systime(1)
      ;print,t1-t0


      ; Vary heights and sigma zero-point
      ;----------------------------------
      initpar2 = par1
      parinfo2 = parinfo
      initpar2[0:3*ngdlines-3:3] = initpar2[0:3*ngdlines-3:3] > 0.001  ; make sure heights are not 0.0
      parinfo2[0:3*ngdlines-3:3].fixed = 0
      if chisq1/dof lt 5 then begin  ; only constrain if the fit was good
        parinfo2[0:3*ngdlines-3:3].limited = [1,1]
        parinfo2[0:3*ngdlines-3:3].limits[0] = 0.5*initpar2[0:3*ngdlines-3:3] > 0.0001  ; constrain heights
        parinfo2[0:3*ngdlines-3:3].limits[1] = 1.5*initpar2[0:3*ngdlines-3:3]
      endif
      parinfo2[1:3*ngdlines-2:3].fixed = 1     ; centers fixed!
      parinfo2[2:3*ngdlines-1:3].fixed = 1     ; keep y-offsets fixed
      parinfo2[3*ngdlines:*].fixed = 1
      parinfo2[3*ngdlines].fixed = 1  ; keep Xglobalcenter fixed
      parinfo2[3*ngdlines+1].fixed = 0 ; allow sigma zero-point to vary
      parinfo2[3*ngdlines+1].limited[0] = 1
      parinfo2[3*ngdlines+1].limits[0] = 0.01  ; must be positive
      ;parinfo2[npars-2].limited=[1,0] & parinfo2[npars-2].limits=0.001  ; wing area must be positive
      ;parinfo2[npars-1].limited=[1,0] & parinfo2[npars-1].limits=1.5    ; minimum wing sigma

      t0 = systime(1)
      par2 = MPFIT('skyfit_lsf_gh_dev',initpar2,$
                       parinfo=parinfo2,dof=dof2,status=status2,bestnorm=chisq2,$
                       perror=perror2,functargs=fa,niter=niter2,autoderivative=0,$
                       nfev=nfev2,ftol=1e-10,gtol=1e-10,xtol=1e-10,/quiet,maxiter=50) ;,maxiter=15)
      t1 = systime(1)
      yfit2 = skyfit_lsf_gh(xin,par2,_extra=fa)

      ; If the fit is bad try letting the y-offsets float
      if chisq2/dof2 gt 5 then begin
        parinfo2[2:3*ngdlines-1:3].fixed = 0
        parinfo2[2:3*ngdlines-1:3].limited = [1,1]
        parinfo2[2:3*ngdlines-1:3].limits[0] = -mad(spec)
        parinfo2[2:3*ngdlines-1:3].limits[1] = mad(spec)

        t0 = systime(1)
        par2 = MPFIT('skyfit_lsf_gh_dev',initpar2,$
                       parinfo=parinfo2,dof=dof2,status=status2,bestnorm=chisq2,$
                       perror=perror2,functargs=fa,niter=niter2,autoderivative=0,$
                       nfev=nfev2,ftol=1e-10,gtol=1e-10,xtol=1e-10,/quiet,maxiter=50) ;,maxiter=15)
        t1 = systime(1)
        yfit2 = skyfit_lsf_gh(xin,par2,_extra=fa)
      endif

      ; Let LSF float more
      ;-------------------
      initpar3 = par2
      parinfo3 = parinfo2
      parinfo3[0:3*ngdlines-1].fixed = 0
      parinfo3[0:3*ngdlines-3:3].limits[0] = 0.8*initpar3[0:3*ngdlines-3:3]  ; constrain heights
      parinfo3[0:3*ngdlines-3:3].limits[1] = 1.2*initpar3[0:3*ngdlines-3:3]
      parinfo3[0:3*ngdlines-3:3].limited = [1,0]
      parinfo3[0:3*ngdlines-3:3].limits[0] = 0   ; flux must be positive
      parinfo3[1:3*ngdlines-2:3].fixed = 1       ; centers fixed!
      parinfo3[2:3*ngdlines-1:3].fixed = 1       ; keep y-offsets fixed
      parinfo3[3*ngdlines:*].fixed = 0
      parinfo3[3*ngdlines:*].limited = 0
      parinfo3[3*ngdlines].fixed = 1  ; keep Xglobalcenter fixed
      parinfo3[npars-2].limited=[1,0] & parinfo3[npars-2].limits=0.001  ; wing area must be positive
      parinfo3[npars-1].limited=[1,0] & parinfo3[npars-1].limits=1.5    ; minimum wing sigma

      t0 = systime(1)
      par3 = MPFIT('skyfit_lsf_gh_dev',initpar3,$
                       parinfo=parinfo3,dof=dof3,status=status3,bestnorm=chisq3,$
                       perror=perror3,functargs=fa,niter=niter3,/quiet,autoderivative=0,$
                       nfev=nfev3,ftol=1e-7,gtol=1e-10,xtol=1e-10,maxiter=50) ;maxiter=15)
      t1 = systime(1)
      yfit3 = skyfit_lsf_gh(xin,par3,_extra=fa)

      ; Final parameters
      fpar = par3
      fperror = perror3 * sqrt(chisq3/dof3) 
      fchisq = chisq3/dof3
      status = status3
      fyfit = yfit3

      ; Plotting
      if keyword_set(pl) then begin
        plotfile=outdir+'/plots/'+dirs.prefix+string(format='("LSF-",a,"_",i1,"_",i3.3)',lsfid,j,i) 
        device,file=plotfile+'.eps',/color,/encap
        loadct,39
        !p.multi=[0,1,3]
        xr = [0,n_elements(xin)-1]
        yr = [-3000,max([specin,yfit3])*1.2]/1e4
        plot,specin/1.e4,xr=xr,yr=yr,xs=1,ys=1,ytit='Counts / 1e4',tit='Sky Fiber='+strtrim(i,2)+' Chip='+strtrim(j,2)+$
             ' Chisq='+stringize(fchisq,ndec=3)+' Nlines='+strtrim(ngdlines,2),charsize=1.3
        oplot,fyfit/1.e4,co=250,linestyle=2
        oplot,(specin-fyfit-2000)/1.e4
        for k=0,ngdlines-1 do begin
          tmp= min(abs(xin-chiplinestr[k].x),ix)
          oplot,[ix,ix],[0.,1e5]
        endfor
        legend,['Data','Final Fit: FLUX, CENTER, and LSF params fit'],$
               textcolor=[200,250],/top,/left,charsize=1.2
        xi=indgen(n_elements(xin))
        plot,xr,xr*0,xr=xr,yr=[-5,5],xs=1,ys=1,tit='fractional residual',linestyle=1,ytit='sigma resid'
        oplot,xr,xr*0,color=150,linestyle=1
        ;iplot=where(fyfit gt 80)
        ;oplot,xi[iplot],(specin[iplot]-fyfit[iplot])/fyfit[iplot]
        oplot,xi,(specin-fyfit)/errspecin
        for k=0,ngdlines-1 do begin
          tmp= min(abs(xin-chiplinestr[k].x),ix)
          oplot,[ix,ix],[-1e5,1e5]/1e4
        endfor

        xr = [0,2048]
        yr=[-500,1000]
        plot,xin,specin,xr=xr,yr=yr,xs=1,ys=1,ytit='Counts',tit='Sky Fiber='+strtrim(i,2)+' Chip='+strtrim(j,2)+$
             ' Chisq='+stringize(fchisq,ndec=3)+' Nlines='+strtrim(ngdlines,2),charsize=1.3
        oplot,spec,color=150
        oplot,xin,specin
        oplot,xin,fyfit,co=250,linestyle=2
        oplot,xin,specin-fyfit-200
        device,/close
        ps2jpg,plotfile+'.eps',/eps,chmod='664'o,/delete
        !p.multi=[0,0,0]
      endif
      combframe.(j).err[useind,i] = combframe.(j).flux[useind,i]
      combframe.(j).err[useind,i] -= fyfit

      ; Print the parameters
      fGHcoefs = fpar[3*ngdlines+1:*]  ; Gauss-Hermite parameters
      ;if keyword_set(verbose) then $
      print,i,j,ngdlines,fchisq,fGHcoefs,format='(I5,I5,I5,F9.2,'+strtrim(n_elements(fGHcoefs),2)+'g9.2)'

      ; Save the fitted parameters
      chiplinestr.lsffit_flux = fpar[0:3*ngdlines-3:3]
      chiplinestr.lsffit_pars[0] = fpar[0:3*ngdlines-3:3]
      chiplinestr.lsffit_pars[1] = fpar[1:3*ngdlines-2:3]
      chiplinestr.lsffit_pars[2] = fpar[2:3*ngdlines-1:3]
      chiplinestr.lsffit_perror[0] = fperror[0:3*ngdlines-3:3]
      chiplinestr.lsffit_perror[1] = fperror[1:3*ngdlines-2:3]
      chiplinestr.lsffit_perror[2] = fperror[2:3*ngdlines-1:3]
      chiplinestr.lsffit_chisq = fchisq
      for k=0,ngdlines-1 do chiplinestr[k].lsffit_wave = PIX2WAVE(chiplinestr[k].lsffit_pars[1],wcoef)
      chiplinestr.lsffit_status = status
      chiplinestr.lsffit_ghcoefs = fGHcoefs

      ; Put the parameters in the output file
      ;lsfpars = [binsize, fpar[3*ngdlines], Horder, Porder, fGHcoefs]
      GHcoefs = fGHcoefs[0:nGHcoefs-1]
      Wcoefs = fGHcoefs[nGHcoefs:*]
      lsfpars = [binsize, fpar[3*ngdlines], Horder, Porder, GHcoefs, Wproftype, nWpar, WPorder, Wcoefs]
      ;lsfpars[1] += xoff  ; add offset to Xoffset
      case j of
      0: lsfim1[i,*] = lsfpars
      1: lsfim2[i,*] = lsfpars
      2: lsfim3[i,*] = lsfpars
      endcase
      fitstr[3*i+j].lsfpars = lsfpars

      ; fiber and chip loaded above in case fit bombs before here (e.g. missing traces)
      ;fitstr[3*i+j].fiber = i
      ;fitstr[3*i+j].chip = j+1
      ;fitstr[3*i+j].xoffset = xoff
      ; Put the fitting information in the structure      
      fitstr[3*i+j].par = fGHcoefs
      fitstr[3*i+j].perror = fperror[3*ngdlines+1:*]
      fitstr[3*i+j].ghcoefs = GHcoefs
      fitstr[3*i+j].wcoefs = Wcoefs
      fitstr[3*i+j].status = status3
      fitstr[3*i+j].chisq = fchisq
      fitstr[3*i+j].npts = n_elements(specin)
      fitstr[3*i+j].rms = sqrt(mean( (specin-yfit3)^2))

      modelflux[useind] = yfit3

      BOMB2:

      ; Stuff the information back into the large structure
      ;linestr[fibchind] = chiplinestr
      PUSH,linestr,chiplinestr

      ; Plot the chip spectrum and the model
      ;pl = 0
      ;if keyword_set(pl) then begin
      ;  plot,spec
      ;  oplot,modelflux,co=250
      ;  stop
      ;endif

      ; Check the "centered-ness" using a Gaussian fit
      xx = findgen(71)*0.2-35*0.2
      lsf = lsf_gh(xx+1024,1024,lsfpars)
      lsf_yfit = mpfitpeak(xx,lsf,pp,/gaussian,/positive,nterms=3)
      print,'LSF off-centerness = ',stringize(pp[1],ndec=4),' pixels'
      ;window,2
      ;plot,xx,lsf
      ;oplot,xx,lsf_yfit,co=250
      ;wset,0

      ;stop

    End ; chip loop

  END ; all lines in a chip fit together


  ;========================================================================
  ; FIT ALL LINES TOGETHER
  ;========================================================================
  3: BEGIN


    spec = fltarr(npix*3)
    errspec = fltarr(npix*3)

    ; The measured chip gaps
    ;  want them in dither combined units, multiply by 2
    ;xoff1 = waveframe.(1).data[i,0] * xmulti  ; Xoffset for chip b
    ;xoff2 = waveframe.(2).data[i,0] * xmulti  ; Xoffset for chip c
    xoff1 = waveframe.(0).coef[i,0] * binsize  ; Xoffset for chip a
    xoff2 = waveframe.(1).coef[i,0] * binsize  ; Xoffset for chip b
    xoff3 = waveframe.(2).coef[i,0] * binsize  ; Xoffset for chip c
    xoff = [xoff1, xoff2, xoff3]

    ; Construct the X array
    x = fltarr(npix*3)
    x[0:npix-1] = findgen(npix)+xoff1
    x[npix:2*npix-1] = findgen(npix)+xoff2
    x[2*npix:3*npix-1] = findgen(npix)+xoff3

    ; Loop through the Chips
    For j=0,2 do begin

      ;fiber = combframe.(j).data[*,i,0]
      ;var = combframe.(j).data[*,i,1]
      fiber = combframe.(j).flux[*,i]
      err = combframe.(j).err[*,i]

      ; Remove the continuum
      ;cont = MEDFILT1D(fiber,101,/edge)
      ;fiber -= cont

      ;chiplinestr = (SCOPE_VARFETCH('linestr'+strtrim(j+1,2)))
      min_heightthresh = 10*median(errspec)
      max_heightthresh = 5e4
      ind = where(linestr.chip eq (j+1) and linestr.fiber eq i and linestr.height gt min_heightthresh and $
                  linestr.height lt max_heightthresh,nind)
      ;ind = where(chiplinestr.fiber eq i,nind)
      if nind eq 0 then goto,BOMB3
      ;ilinestr = chiplinestr[ind]
      ilinestr = linestr[ind]

      ; Convert pixel positions to wavelengths
      ;pcen = ilinestr.gpar[1] / xscale    ; convert to undersampled pixel values
      ;lwcen = PIX2WAVE(pcen,wcoef)
      ;ilinestr.wcen = lwcen

      ilinestr.peakx += xoff[j]
      ilinestr.gaussx += xoff[j]
      ilinestr.gpar[1] += xoff[j]
      ilinestr.gpar0[1] += xoff[j]

      spec[j*npix:(j+1)*npix-1] = fiber
      errspec[j*npix:(j+1)*npix-1] = err

      ; Add to fiber linelist
      PUSH,fiberlinestr,ilinestr

      BOMB3:
    End
    nlines = n_elements(fiberlinestr)
    nspec = n_elements(spec)
    ;y = findgen(nspec)

    ; gausshermitebin.pro works for ONE line
    ; need a function like GHM_all_1d_decom  that runs it on all lines

    ; Only fit pixels that are close to a line, within +/-10 sigma or so.

    ; KLUDGE
    ;  make sure the errors on the continuum are high enough
    ;sig0 = MAD(spec,/zero)
    ;errspec = sqrt(errspec^2 + sig0^2)

    ; Get all pixels close to a line
    mask = intarr(nspec)
    xloarr = lonarr(nlines)
    xhiarr = lonarr(nlines)
    for j=0,nlines-1 do begin
      ;ylo = round(fiberlinestr[j].gaussx-nsigfit.0*fiberlinestr[j].gpar[2]) > 0
      ;yhi = round(fiberlinestr[j].gaussx+nsigfitfiberlinestr[j].gpar[2]) < (nspec-1)
      gd = where(abs(x-fiberlinestr[j].gaussx) lt nsigfit*fiberlinestr[j].gpar[2],ngd)
      xlo = min(gd)
      xhi = max(gd)
      mask[xlo:xhi] = 1
      ;gd = where(abs(y-fiberlinestr[j].gaussx) lt 10.*fiberlinestr[i].gpar[2],ngd)
      ;if ngd gt 0 then mask[gd] = 1
      xloarr[j] = xlo
      xhiarr[j] = xhi    
    end
    useind = where(mask eq 1,nuseind)

    ; Input values
    ;Xoffset = 6430
    xin = x[useind] ;-Xoffset
    specin = spec[useind]
    errspecin = errspec[useind]

    ; lo/hi indices for the "input" arrays
    npixarr = xhiarr-xloarr+1
    loarr = lonarr(nlines)
    hiarr = lonarr(nlines)
    ;hiarr[0] = npixarr[0]-1
    for j=0,nlines-1 do begin
      ;loarr[j] = hiarr[j-1]+1
      ;hiarr[j] = loarr[j]+npixarr[j]-1
      loarr[j] = where(useind eq xloarr[j])
      hiarr[j] = where(useind eq xhiarr[j])
    end

    ;stop

    ; Allow the GH coefficients to vary as polynomials
    ;  this is the ORDER, so 0 means constant
    ;porder = [2, 2, 2, 2, 2, 2]
    ;Porder = [0, 0, 0, 0, 0]
    nPorder = n_elements(Porder)

    ; Now fit the LSF
    ;  the flux of each line is allowed to vary
    ; The parameters are:
    ;  Nline flux and Nline center parameters (these are interleaved
    ;     e.g. flux1, cen1, flux1, cen2, etc.)
    ;  6 Gauss-Hermitze LSF parameters, the height and center aren't needed

    ; The parameters are:
    ;  Nlines   flux/center for the Nlines
    ;  X0       Xoffset
    ;  GHcoefs  The Gauss-Hermite polynomial variation coefficients
    ;binsize = 2
    ;binsize = xmulti

    ; The "global" center
    ; I want the coefficients that vary globally with Y to be fit around
    ; the center of the array, so that the constant term is the value
    ; at the center of the array, near 6430 from the left edge.
    ;Xglobalcenter = 6430
    Xglobalcenter = 0

    ; The Xoffset:  this is used to set the global center and is used
    ; to apply the offsets required to chips b and c pixels to be on 
    ; the "global" pixel scale (X=0 at the left edge of chip a).

    initpar = dblarr(3*nlines+1+TOTAL(porder+1))
    ; We need to multiply the heights by sig*sqrt(2*pi) because the LSFs
    ;  are normalized but the Gaussians that APPEAKFIT.PRO have an
    ;  area under the curve of ht*sig*sqrt(2*pi).
    initpar[0:3*nlines-3:3] = fiberlinestr.gpar[0] * fiberlinestr.gpar[2] * sqrt(2*!dpi)   ; flux
    initpar[1:3*nlines-2:3] = fiberlinestr.gpar[1] ; -Xoffset        ; center
    initpar[2:3*nlines-1:3] = 0.0          ; Y offset
    initpar[3*nlines] = -Xglobalcenter   ;Xoffset
    initpar[3*nlines+1] = median(fiberlinestr.gpar[2])
    if i gt 0 then initpar[3*nlines+1:*]=GHcoefs  ; Use fit from previous fiber
    ; leave the GHcoefs at ZERO for first fiber
    npars = n_elements(initpar)

    ;if i eq 0 then ghpars = [median(fiberlinestr.gpar[2]),1.0,fltarr(4)]
    ;; use i-1 ghpars as initial guesses
    ;initpar = [(fiberlinestr.gpar[0:1])(*), ghpars]
    ;npars = nlines*2+6
    ;initpar[2*indgen(nlines)] = fiberlinestr.gpar[0]/binsize
    ;   ; appeakfit used binsize=1
    parinfo = REPLICATE({limited:[0,0],limits:[0.0,0.0],fixed:0},npars)
    parinfo[3*nlines].fixed = 1  ; Keep X-center FIXED
    parinfo[0:3*nlines-3:3].limited = [1,0] ; flux must be positive
    parinfo[0:3*nlines-3:3].limits[0] = 0.0
    parinfo[2:3*nlines-1:3].limited = [1,1]
    parinfo[2:3*nlines-1:3].limits = [-1,1]*mad(spec)   ; can't move too much
    parinfo[3*nlines+1].limited = [1,0]
    parinfo[3*nlines+1].limits[0] = 0.001     ; sigma must be positive
    fa = {binsize:binsize,nlines:nlines,loarr:loarr,hiarr:hiarr,porder:porder}

    ; Kludge, fix sigma=1
    ;initpar[2*nlines+1] = 1.0
    ;parinfo[2*nlines+1].fixed = 1
    ;parinfo[2*nlines+1].limited = 0


    ; Initial guess fit
    yfit0 = fit_lsf_gh(xin,initpar,binsize=binsize,nlines=nlines,loarr=loarr,hiarr=hiarr,porder=porder)

    ; ITERATE
    endflag = 0
    count = 0
    WHILE (endflag eq 0) do begin

      if count eq 0 then initpar1=initpar else initpar1=lastpar

      ; First fit the GH parameters
      parinfo1 = parinfo
      parinfo1[0:3*nlines-1].fixed = 1
      par1 = MPFITFUN('fit_lsf_gh',xin,specin,errspecin,initpar1,$
                     parinfo=parinfo1,yfit=yfit1,dof=dof,status=status1,bestnorm=chisq1,$
                     perror=perror1,functargs=fa,niter=niter1,/quiet)

      ;stop

      ; Now fit the height/centers, keep GH pars fixed
      parinfo2 = parinfo
      parinfo2[3*nlines+1:*].fixed = 1
      par2 = MPFITFUN('fit_lsf_gh',xin,specin,errspecin,par1,$
                     parinfo=parinfo2,yfit=yfit2,dof=dof,status=status2,bestnorm=chisq2,$
                     perror=perror2,functargs=fa,niter=niter2,/quiet) ;,autoderivative=0) ;,/quiet
      ; the autoderivative isn't working properly in lsf_ghb.pro

      ; Increase the errors for outlier points
      ;  don't remove the points or it will mess up the index arrays
      ;  should probably reduce DOF by nbd
      if count gt 3 then begin
        diff = (specin-yfit2)/(errspecin>1)
        bd = where(abs(diff) gt 7,nbd)
        if nbd gt 0 then errspecin[bd] *= 100
      endif

      ;stop

      ; Check if we should quit
      if count gt 0 then begin
        ;if count gt 5 then endflag = 1   ; max 5 interations
        if count gt 10 then endflag = 1   ; max 5 interations
        diffpar = abs(par2-lastpar)/abs(par2)
        if max(diffpar)*100 lt 3 then endflag = 1
      endif
      lastpar = par2

      count++

      ;stop

    ENDWHILE

    ; Final parameters
    fpar = par2
    yfit = yfit2
    fchisq = chisq2/dof
    pcerror = perror2 * sqrt(chisq2/dof)    ; scaled uncertainties
    ;print,'Chisq = ',strtrim(fchisq,2)

    ;; Now allow all to float
    ;parinfo3 = parinfo
    ;par3 = MPFITFUN('fit_lsf_gh',yin,specin,errspecin,par2,$
    ;               parinfo=parinfo3,yfit=yfit3,dof=dof,status=status3,bestnorm=chisq3,$
    ;               perror=perror3,functargs=fa,/quiet)

    ;dum = lsf_gh(yin,initpar,binsize=binsize,nlines=nlines,loarr=loarr,hiarr=hiarr)

    ; Plotting
    if keyword_set(pl) then begin
      plot,specin,xs=1,ytit='Counts',tit='Fiber '+strtrim(i+1,2)+' Chisq = '+stringize(fchisq,ndec=3)+' Nlines='+strtrim(nlines,2)
      oplot,yfit,co=250,linestyle=2
      ;xyouts,n_elements(specin)*0.5,0.9*max(specin),'Chisq = '+strtrim(fchisq,2),align=0.5,charsize=1.3
    endif

    ; I can see a slight difference between the fit on the blue/red ends 

    ; Print the parameters
    GHcoefs = fpar[3*nlines+1:*]  ; Gauss-Hermite parameters
    nGHcoefs = n_elements(GHcoefs)
    if keyword_set(verbose) then begin
      ;fmt = '(I4,I6,6F8.4)'
      ;print,format=fmt,i+1,nlines,ghpars
       ;fmt = '(I4,I6,F8.4,'+strtrim(nGHcoefs,2)+'G10.3)'
       fmt = '(I4,I6,F9.4,'+strtrim(nGHcoefs,2)+'G10.3)'
       print,format=fmt,i+1,nlines,fchisq,GHcoefs
    endif

    ; Put the parameters in the output file
    lsfpars = [binsize, fpar[3*nlines], Horder, Porder, GHcoefs]
    lsfim1[i,*] = lsfpars
    lsfim1[i,1] += xoff[0]  ; add offset to Xoffset
    lsfim2[i,*] = lsfpars
    lsfim2[i,1] += xoff[1]  ; add offset to Xoffset
    lsfim3[i,*] = lsfpars
    lsfim3[i,1] += xoff[2]  ; add offset to Xoffset

    ;stop

  END ; fitting all lines together

  else: stop,'This fitting method not supported'

  ENDCASE  ; different fitting methods

  ;stop

End ; fiber loop


mwrfits,combframe.chipa.flux,outdir+dirs.prefix+string(format='("LSF-",a)',lsfid)+'_flux.fits',/create
mwrfits,combframe.chipb.flux,outdir+dirs.prefix+string(format='("LSF-",a)',lsfid)+'_flux.fits'
mwrfits,combframe.chipc.flux,outdir+dirs.prefix+string(format='("LSF-",a)',lsfid)+'_flux.fits'
mwrfits,combframe.chipa.err,outdir+dirs.prefix+string(format='("LSF-",a)',lsfid)+'_res.fits',/create
mwrfits,combframe.chipb.err,outdir+dirs.prefix+string(format='("LSF-",a)',lsfid)+'_res.fits'
mwrfits,combframe.chipc.err,outdir+dirs.prefix+string(format='("LSF-",a)',lsfid)+'_res.fits'
htmlfile=outdir+'/html/'+dirs.prefix+string(format='("LSF-",a)',lsfid)
openw,html,htmlfile+'.html',/get_lun
printf,html,'<HTML><BODY><TABLE BORDER=2>'
for ii=0,n_elements(ifibers)-1 do begin
  i=ifibers[ii]
  printf,html,'<TR>'
  for j=0,2 do begin
      plotfile='../plots/'+dirs.prefix+string(format='("LSF-",a,"_",i1,"_",i3.3)',lsfid,j,i) 
      printf,html,'<TD><IMG SRC='+plotfile+'.jpg>'
  endfor
endfor
  free_lun,html


if keyword_set(verbose) then $
  print,'------------------------------------------------------------'

;stop

;save,lsfim,file='aplsfim.dat'
;stop
OUTPUT:
;restore,'aplsfim.dat'
;lsfim[*,1] = -6430   ; there was a bug
;lsfim1 = lsfim
;lsfim2 = lsfim
;lsfim2[*,1] += 4391.5182
;lsfim3 = lsfim
;lsfim3[*,1] += 8766.7898
;restore,'aplsfim2.dat'

skiptohere:

; You do see variations of the GH coefficients with fiber #


; Fit polynomials to the different coefficients with YPOS
; use these at least for outliers.

; Fit the parameters as a function of YPOS
;------------------------------------------
print,'Fitting parameters as functions of YPOS'
;ypos = fltarr(nfibers)
;for i=0,nfibers-1 do ypos[i]=poly(npix/2,tstr2[i].coef)

; since we may not have all fibers in the PSF, just use row number to interpolate 
;  we could instead get these from a PSF that has all fibers, i.e., the fiberid frame
ypos=indgen(nfibers)

;fitind = [0,1,2,3,4,5,6,7,8,9,10,11,12,13]
npar = n_elements(fitstr[0].par)
fitind = indgen(npar)
fitorder = 3 ;2
fitcoef = fltarr(n_elements(fitind),3,fitorder+1)
fitstr1 = fitstr
fitstr2 = fitstr

for chip=0,2 do begin
  ind_chip = where(fitstr1.chip eq chip+1,nind)
  ; get fibers with good chisq values
  med_chisq = median(fitstr1[ind_chip].chisq)
  sig_chisq = mad(fitstr1[ind_chip].chisq)
  gd = where( fitstr1[ind_chip].chisq lt 3.5*sig_chisq+med_chisq,ngd)
  ind = ind_chip[gd]
  for i=0,n_elements(fitind)-1 do begin
    fitcoef[i,chip,*] = ROBUST_POLY_FIT(ypos[fitstr1[ind].fiber],fitstr1[ind].par[fitind[i]],fitorder)
    fitstr2[ind_chip].par[fitind[i]] = POLY(ypos[fitstr[ind_chip].fiber],fitcoef[i,chip,*])
  endfor
endfor

; Median smoothing
fitstr3 = fitstr
for chip=0,2 do begin
  ind = where(fitstr.chip eq chip+1,nind)
  for i=0,n_elements(fitind)-1 do begin
    smlen = 31 < nind
    smpar = MEDFILT1D(fitstr[ind].par[fitind[i]],smlen,/edge)
    fitstr3[ind].par[fitind[i]] = smpar
  endfor
endfor

; Use Median smoothed parameters
;print,'USING MEDIAN SMOOTHED PARAMETERS'
;fitstr = fitstr3

; Using fitted parameters only for fibers with BAD fits
print,'Using fitted parameters ONLY for fibers with BAD Chisq'
for i=0,2 do begin
  ind_chip = where(fitstr.chip eq i+1,nind)

  ; Look at chisq
  med_chisq = median(fitstr[ind_chip].chisq)
  sig_chisq = mad(fitstr[ind_chip].chisq)
  print,'median chisq: ',med_chisq,' sig_chisq: ', sig_chisq

  ; Look at difference in parameters compared to MEDIAN smoothed
  diff = fitstr[ind_chip].par-fitstr3[ind_chip].par
  sigpar = MAD(diff,dim=2)
  sigdiff = diff/(sigpar#replicate(1,nind))
  totsigdiff = sqrt(total(sigdiff^2,1))
  maxsigdiff = max(sigdiff,dim=1)


  bd = where( (fitstr[ind_chip].chisq gt 3.0*sig_chisq+med_chisq),nbd,comp=gd)
  ;bd = where( (fitstr[ind_chip].chisq gt 3.0*sig_chisq+med_chisq) OR $
  ;            (maxsigdiff gt median(maxsigdiff)+3.5*mad(maxsigdiff)),nbd,comp=gd)
  if nbd gt 0 then begin
    print,'replacing all pars:',fitstr[ind_chip[bd]].fiber,fitstr[ind_chip[bd]].chip
    ; replace with nearest fiber with good fit
    for ibd=0,nbd-1 do begin
      junk=min(abs(fitstr[ind_chip[gd]].fiber-fitstr[ind_chip[bd[ibd]]].fiber),imin)
      fitstr[ind_chip[bd[ibd]]].par = fitstr[ind_chip[gd[imin]]].par
      print,'replace: ', fitstr[ind_chip[bd[ibd]]].fiber,fitstr[ind_chip[bd[ibd]]].chip,fitstr[ind_chip[bd[ibd]]].chisq
      print,' with: ', fitstr[ind_chip[gd[imin]]].fiber,fitstr[ind_chip[gd[imin]]].chip
    endfor
    ; replace with smoothed parameters
    ;fitstr[ind_chip[bd]].par = fitstr3[ind_chip[bd]].par
    ;stop
  endif

  ; Holtz., I think the following is a bad idea, since parameters can be correlated

  ; Change fitted values for parameters that are FAR from the median values
  ;for j=0,n_elements(fitind)-1 do begin
  ;  par = fitstr[ind_chip].par[fitind[j]]
  ;  medpar = fitstr3[ind_chip].par[fitind[j]]
  ;  sig = MAD(par-medpar)
  ;  bdpar = where(abs(par-medpar) gt 4.*sig,nbdpar)
  ;  if nbdpar gt 0 then begin
  ;    print,'replacing indiv pars: ',fitstr[ind_chip[bdpar]].fiber,fitstr[ind_chip[bdpar]].chip
  ;    fitstr[ind_chip[bdpar]].par[fitind[j]]=fitstr3[ind_chip[bdpar]].par[fitind[j]]
  ;  endif
  ;endfor

endfor



;lsfpars = [binsize, fpar[3*ngdlines], Horder, Porder, GHcoefs, Wproftype, nWpar, WPorder, Wcoefs]
;nlsfpars = 3+n_elements(Porder)+total(Porder+1)+2+n_elements(WPorder)+total(WPorder+1)
;add_tag,fitstr,'LSFPARS',fltarr(nlsfpars),fitstr
nghcoefs = total(Porder+1)
nwcoefs = total(WPorder+1)
for i=0,n_elements(fitstr)-1 do begin
  ghcoefs = fitstr[i].par[0:nghcoefs-1]
  wcoefs = fitstr[i].par[nghcoefs:*]
  binsize = fitstr[i].lsfpars[0]
  Xglobalcenter = fitstr[i].lsfpars[1]
  lsfpars = [binsize, Xglobalcenter, Horder, Porder, ghcoefs, Wproftype, nWpar, WPorder, wcoefs]
  fitstr[i].lsfpars = lsfpars
endfor


; Maybe make separate plots for each chip

; Plot the fits
if keyword_set(pl) then begin
  psym8
  !p.multi=[0,4,4]
  ;erase

  coarr = [250,150,80]

  for i=0,n_elements(fitind)-1 do begin
    y = fitstr.par[fitind[i]]
    med = median(y)
    sig = mad(y)
    yr = [-5,5]*sig+med
    yr[0] = yr[0] > min(y)
    yr[1] = yr[1] < max(y)

    plot,fitstr.par[fitind[i]],xr=[0,2048],yr=yr,xs=1,ys=1,$
      tit='Parameter='+strtrim(fitind[i],2),charsize=1.8,xtit='YPOS',ytit='Value',/nodata
    for j=0,2 do begin
      ind = where(fitstr.chip eq j+1,nind)
      oplot,ypos[fitstr1[ind].fiber],fitstr1[ind].par[fitind[i]],co=coarr[j],ps=8,sym=0.7
      oplot,ypos[fitstr[ind].fiber],fitstr[ind].par[fitind[i]],co=coarr[j],thick=1,linestyle=2
      ;oplot,ypos[fitstr[ind].fiber],fitstr3[ind].par[fitind[i]],co=coarr[j],thick=1 ;,linestyle=2
    endfor
  endfor
  if pl gt 1 then stop
  !p.multi=0
endif

; SAVE DIAGNOSTIC PLOTS LIKE IN APWAVECAL.PRO

;stop

;;------------------------------------
;; Fit a second time, but constrained
;;-----------------------------------
;if FITMETHOD eq 2 then begin
;
;  ; Use fitted parameters as initial guess
;  fitstr4 = fitstr
;  add_tag,fitstr4,'INITPAR',fitstr[0].par*0,fitstr4
;  fitstr4.initpar = fitstr2.par
;
;  ;lsfpars = [binsize, fpar[3*ngdlines], Horder, Porder, GHcoefs, Wproftype, nWpar, WPorder, Wcoefs]
;  nlsfpars = 3+n_elements(Porder)+total(Porder+1)+2+n_elements(WPorder)+total(WPorder+1)
;  add_tag,fitstr4,'LSFPARS',fltarr(nlsfpars),fitstr4
;
;  x = findgen(npix)
;
;  ; Loop through the Fibers
;  For i=0,nfibers-1 do begin
;
;    if not keyword_set(verbose) then print,'Fiber ',strtrim(i+1,2)
;
;    apgundef,fiberlinestr
;
;    ; Loop through the Chips
;    For j=0,2 do begin
;
;      spec = combframe.(j).flux[*,i]
;      errspec = combframe.(j).err[*,i]
;      wave = waveframe.(j).wave[*,i]
;      wcoef = waveframe.(j).coef[i,*]
;      modelflux = spec*0.0
;
;      ; Do a better job of removing the continuum
;      med0 = MEDFILT1D(spec,101,/edge)
;      temp = spec
;      sig1 = MAD(temp)
;      mask = long(abs(temp-med0) gt 2*sig1)
;      mask = convol(mask,indgen(5))
;      mask = mask/(mask>1)
;      bd1 = where(mask eq 1,nbd1)
;      if nbd1 gt 0 then temp[bd1]=!values.f_nan
;      med1 = MEDFILT1D(temp,101,/edge)
;
;      temp2 = spec
;      mask2 = long(abs(temp2-med1) gt 2*sig1)
;      mask2 = convol(mask2,indgen(5))
;      mask2 = mask2/(mask2>1)
;      bd2 = where(mask2 eq 1,nbd2)
;      if nbd2 gt 0 then temp2[bd2]=!values.f_nan
;      med2 = MEDFILT1D(temp2,101,/edge)
;
;      spec -= med2  ; remove the continuum
;
;
;      ; Get all lines for this fiber/chip from linestr
;      gdlines = where(linestr.fiber eq i and linestr.chip eq j and linestr.lsffit_pars[0] gt 0.0,ngdlines)
;      chiplinestr = linestr[gdlines]
;
;
;      ; Need to fit all lines together
;      ;-------------------------------
;      nsigfit = 8
;
;
;      ; Get all pixels close to a line
;      mask = intarr(npix)
;      xloarr = lonarr(ngdlines)
;      xhiarr = lonarr(ngdlines)
;      for k=0,ngdlines-1 do begin
;        ;gd = where(abs(x-chiplinestr[k].x) lt nsigfit*1.0,ngd)
;        if chiplinestr[k].model_doublet eq 0 then begin
;          gd = where(abs(x-chiplinestr[k].x) lt nsigfit*2.0,ngd)
;        endif else begin
;          gd = where( x ge (chiplinestr[k].x-0.5*chiplinestr[k].model_dbl_xsep-nsigfit*2.0) and $
;                      x le (chiplinestr[k].x+0.5*chiplinestr[k].model_dbl_xsep+nsigfit*2.0),ngd)
;        endelse
;        xlo = min(gd)
;        xhi = max(gd)
;        mask[xlo:xhi] = 1
;        xloarr[k] = xlo
;        xhiarr[k] = xhi
;      end
;      useind = where(mask eq 1,nuseind)
;
;      ; Input values
;      xin = x[useind]
;      specin = spec[useind]
;      errspecin = errspec[useind]
;
;      ; lo/hi indices for the "input" arrays
;      npixarr = xhiarr-xloarr+1
;      loarr = lonarr(ngdlines)
;      hiarr = lonarr(ngdlines)
;      ;hiarr[0] = npixarr[0]-1
;      for k=0,ngdlines-1 do begin
;        loarr[k] = where(useind eq xloarr[k])
;        hiarr[k] = where(useind eq xhiarr[k])
;      end
;
;
;     ; Allow the GH coefficients to vary as polynomials
;     ;  this is the ORDER, so 0 means constant
;     ;porder = [2, 2, 2, 2, 2, 2]
;     ;Porder = [0, 0, 0, 0, 0]
;     nPorder = n_elements(Porder)
;
;     ; Now fit the LSF
;     ;  the flux of each line is allowed to vary
;     ; The parameters are:
;     ;  Nline flux and Nline center parameters (these are interleaved
;     ;     e.g. flux1, cen1, flux1, cen2, etc.)
;     ;  6 Gauss-Hermitze LSF parameters, the height and center aren't needed
;
;     ; The parameters are:
;     ;  Nlines   flux/center for the Nlines
;     ;  X0       Xoffset
;     ;  GHcoefs  The Gauss-Hermite polynomial variation coefficients
;     ;binsize = 2
;     ;binsize = xmulti
;
;     ; The "global" center
;     ; I want the coefficients that vary globally with Y to be fit around
;     ; the center of the array, so that the constant term is the value
;     ; at the center of the array, near 6430 from the left edge.
;     ;Xglobalcenter = 6430
;     Xglobalcenter = npix/2
;
;     ; ****GHcoefs means GH and Wing coefficients here****
;
;     ind = where(fitstr4.fiber eq i and fitstr4.chip eq j+1,nind)
;     init_ghcoefs = fitstr4[ind].initpar
;
;     ; Initial parameters
;     initpar = dblarr(3*ngdlines+1+nGHcoefs+nWcoefs)
;     ; We need to multiply the heights by sig*sqrt(2*pi) because the LSFs
;     ;  are normalized but the Gaussians that APPEAKFIT.PRO have an
;     ;  area under the curve of ht*sig*sqrt(2*pi).
;     initpar[0:3*ngdlines-3:3] = chiplinestr.lsffit_pars[0]
;     initpar[1:3*ngdlines-2:3] = chiplinestr.lsffit_pars[1]
;     initpar[2:3*ngdlines-1:3] = 0.0      ;            constant Y-offset
;     initpar[3*ngdlines] = Xglobalcenter
;     initpar[3*ngdlines+1:*] = init_ghcoefs
;
;     fa = {x:xin,y:specin,err:errspecin,binsize:binsize,nlines:ngdlines,loarr:loarr,$
;           hiarr:hiarr,porder:porder,wproftype:wproftype,wporder:wporder,$
;           doublet:chiplinestr.model_doublet,dbl_sep:chiplinestr.model_dbl_xsep}
;
;     ; Constraints
;     npars = n_elements(initpar)
;     parinfo = REPLICATE({limited:[0,0],limits:[0.0,0.0],fixed:0},npars)
;     ;parinfo[1:3*ngdlines-2:3].fixed = 1       ; keep x-values fixed
;     parinfo[0:3*ngdlines-1].fixed = 0
;     parinfo[0:3*ngdlines-1].limited = [1,0]
;     ;parinfo[0:3*ngdlines-1].limits[0] = 0   ; flux must be positive
;     parinfo[0:3*ngdlines-3:3].limited = 1
;     parinfo[0:3*ngdlines-3:3].limits[0] = initpar[0:3*ngdlines-1:3]-0.2*initpar[0:3*ngdlines-1:3]
;     parinfo[0:3*ngdlines-3:3].limits[1] = initpar[0:3*ngdlines-1:3]+0.2*initpar[0:3*ngdlines-1:3]
;     parinfo[1:3*ngdlines-2:3].fixed = 1
;     ;parinfo[1:3*ngdlines-2:3].limited = 1
;     ;parinfo[1:3*ngdlines-2:3].limits[0] = -0.3 + initpar[1:3*ngdlines-2:3]
;     ;parinfo[1:3*ngdlines-2:3].limits[1] = 0.3 + initpar[1:3*ngdlines-2:3]
;     parinfo[2:3*ngdlines-1:3].fixed = 1       ; keep y-offsets fixed
;
;     ; Let the heights AND centers vary a little bit
;     initpar1 = initpar
;     parinfo1 = parinfo
;     parinfo1[0:3*ngdlines-3:3].fixed = 0
;     parinfo1[0:3*ngdlines-3:3].limited = 1
;     parinfo1[0:3*ngdlines-3:3].limits[0] = initpar1[0:3*ngdlines-1:3]-0.2*initpar1[0:3*ngdlines-1:3]
;     parinfo1[0:3*ngdlines-3:3].limits[1] = initpar1[0:3*ngdlines-1:3]+0.2*initpar1[0:3*ngdlines-1:3]
;     parinfo1[1:3*ngdlines-2:3].fixed = 0
;     parinfo1[1:3*ngdlines-2:3].limited = 1
;     parinfo1[1:3*ngdlines-2:3].limits[0] = -3 + initpar1[1:3*ngdlines-2:3]
;     parinfo1[1:3*ngdlines-2:3].limits[1] = 3 + initpar1[1:3*ngdlines-2:3]
;     parinfo1[2:3*ngdlines-1:3].fixed = 1       ; keep y-offsets fixed
;     parinfo1[3*ngdlines:*].fixed = 1   ; the rest are fixed
;     t0 = systime(1)
;     par1 = MPFIT('skyfit_lsf_gh_dev',initpar1,$
;                      parinfo=parinfo1,dof=dof,status=status1,bestnorm=chisq1,$
;                      perror=perror1,functargs=fa,niter=niter1,/quiet,autoderivative=0,$
;                      nfev=nfev1,ftol=1e-5,gtol=1e-5,xtol=1e-5,maxiter=3) ;,/fastnorm)
;
;
;     ; Second fit, let all float
;     parinfo2 = parinfo1
;     initpar2 = par1
;     parinfo2[0:3*ngdlines-3:3].limits[0] = initpar2[0:3*ngdlines-1:3]-0.1*initpar2[0:3*ngdlines-1:3]
;     parinfo2[0:3*ngdlines-3:3].limits[1] = initpar2[0:3*ngdlines-1:3]+0.1*initpar2[0:3*ngdlines-1:3]
;     parinfo2[1:3*ngdlines-2:3].limits[0] = -0.3 + initpar2[1:3*ngdlines-2:3]
;     parinfo2[1:3*ngdlines-2:3].limits[1] = 0.3 + initpar2[1:3*ngdlines-2:3]
;
;     ; Keep all higher order terms or constant parameters fixed
;     ;  keep wing parameters fixed
;     startind = 3*ngdlines+1
;     parinfo[startind:*].fixed = 1  ; all fixed for now
;     hi_ind = total(Porder+1,/cum)-1
;     lo_ind = [0,hi_ind[0:Horder-1]+1]
;
;     for k=0,Horder do begin
;       if Porder[k] gt 0 then begin
;         indx = lo_ind[k]+startind
;         parinfo[indx].fixed = 0
;         parinfo[indx].limited = 1
;         parinfo[indx].limits[0] = initpar[indx]-0.1*abs(initpar[indx])
;         parinfo[indx].limits[1] = initpar[indx]+0.1*abs(initpar[indx])
;       endif        
;     end
;
;
;     t0 = systime(1)
;     par2 = MPFIT('skyfit_lsf_gh_dev',initpar2,$
;                      parinfo=parinfo2,dof=dof2,status=status2,bestnorm=chisq2,$
;                      perror=perror2,functargs=fa,niter=niter2,autoderivative=0,$
;                      nfev=nfev2,ftol=1e-7,gtol=1e-7,xtol=1e-7,maxiter=10,/quiet)
;     yfit2 = skyfit_lsf_gh(xin,par2,_extra=fa)
;
;     ; Final parameters
;     fpar = par2
;     fperror = perror2 * sqrt(chisq2/dof2) 
;     fchisq = chisq2/dof2
;     status = status2
;
;     ; Plotting
;     if keyword_set(pl) then begin
;       xr = [0,n_elements(xin)-1]
;       plot,specin,xs=1,xr=xr,ytit='Counts',tit='Sky Fiber='+strtrim(i,2)+' Chip='+strtrim(j,2)+$
;            ' Chisq='+stringize(fchisq,ndec=3)+' Nlines='+strtrim(ngdlines,2),charsize=1.3
;       oplot,yfit2,co=250,linestyle=2
;       legend,['Data','Final fit'],textcolor=[200,250],/top,/left,charsize=1.2
;     endif
;
;     ; Print the parameters
;     fGHcoefs = fpar[3*ngdlines+1:*]  ; Gauss-Hermite parameters
;     ;if keyword_set(verbose) then $
;     print,i,j,ngdlines,fchisq,fGHcoefs,format='(I5,I5,I5,F9.2,'+strtrim(n_elements(fGHcoefs),2)+'g9.2)'
;
;     ; Save the fitted parameters
;     chiplinestr.lsffit_flux = fpar[0:3*ngdlines-3:3]
;     chiplinestr.lsffit_pars[0] = fpar[0:3*ngdlines-3:3]
;     chiplinestr.lsffit_pars[1] = fpar[1:3*ngdlines-2:3]
;     chiplinestr.lsffit_pars[2] = fpar[2:3*ngdlines-1:3]
;     chiplinestr.lsffit_perror[0] = fperror[0:3*ngdlines-3:3]
;     chiplinestr.lsffit_perror[1] = fperror[1:3*ngdlines-2:3]
;     chiplinestr.lsffit_perror[2] = fperror[2:3*ngdlines-1:3]
;     chiplinestr.lsffit_chisq = fchisq
;     for k=0,ngdlines-1 do chiplinestr[k].lsffit_wave = PIX2WAVE(chiplinestr[k].lsffit_pars[1],wcoef)
;     chiplinestr.lsffit_status = status
;     chiplinestr.lsffit_ghcoefs = fGHcoefs
;
;     ; Put the parameters in the output file
;     GHcoefs = fGHcoefs[0:nGHcoefs-1]
;     Wcoefs = fGHcoefs[nGHcoefs:*]
;     lsfpars = [binsize, fpar[3*ngdlines], Horder, Porder, GHcoefs, Wproftype, nWpar, WPorder, Wcoefs]
;     ;lsfpars[1] += xoff  ; add offset to Xoffset
;     ;case j of
;     ;0: lsfim1[i,*] = lsfpars
;     ;1: lsfim2[i,*] = lsfpars
;     ;2: lsfim3[i,*] = lsfpars
;     ;endcase
;     fitstr4[3*i+j].lsfpars = lsfpars
;
;     ; Put the fitting information in the structure      
;     fitstr4[3*i+j].fiber = i
;     fitstr4[3*i+j].chip = j+1
;     ;fitstr4[3*i+j].xoffset = xoff
;     fitstr4[3*i+j].par = fGHcoefs
;     fitstr4[3*i+j].perror = fperror[3*ngdlines+1:*]
;     fitstr4[3*i+j].ghcoefs = GHcoefs
;     fitstr4[3*i+j].wcoefs = Wcoefs
;     fitstr4[3*i+j].status = status2
;     fitstr4[3*i+j].chisq = fchisq
;     fitstr4[3*i+j].npts = n_elements(specin)
;     fitstr4[3*i+j].rms = sqrt(mean( (specin-yfit2)^2))
;
;     modelflux[useind] = yfit2
;
;     ; Stuff the information back into the large structure
;     ;linestr[fibchind] = chiplinestr
;     PUSH,linestr2,chiplinestr
;
;     ; Plot the chip spectrum and the model
;     ;pl = 0
;     ;if keyword_set(pl) then begin
;     ;  plot,spec
;     ;  oplot,modelflux,co=250
;     ;  stop
;     ;endif
;
;     ;wait,0.5
;     if fchisq gt 40 then stop
;
;     ;stop
;
;   End ; chip loop
;
; Endfor ; fiber loop
;
;
; ; Plot the fits
; if keyword_set(pl) then begin
;   psym8
;   !p.multi=[0,4,4]
;   ;erase
;
;   coarr = [250,150,80]
;
;   for i=0,n_elements(fitind)-1 do begin
;     y = fitstr4.par[fitind[i]]
;     med = median(y)
;     sig = mad(y)
;     yr = [-4,4]*sig+med
;     yr[0] = yr[0] > min(y)
;     yr[1] = yr[1] < max(y)
;
;     plot,fitstr4.par[fitind[i]],xr=[0,2048],yr=yr,xs=1,ys=1,$
;       tit='Parameter='+strtrim(fitind[i],2),charsize=1.8,xtit='YPOS',ytit='Value',/nodata
;     for j=0,2 do begin
;       ind = where(fitstr.chip eq j+1,nind)
;       oplot,ypos[fitstr4[ind].fiber],fitstr4[ind].par[fitind[i]],co=coarr[j],ps=8,sym=0.7
;       ;oplot,ypos[fitstr[ind].fiber],fitstr2[ind].par[fitind[i]],co=coarr[j],thick=1,linestyle=2
;       ;oplot,ypos[fitstr[ind].fiber],fitstr3[ind].par[fitind[i]],co=coarr[j],thick=1 ;,linestyle=2
;     endfor
;     ;stop
;   endfor
;   !p.multi=0
; endif
;
;
;endif ; fitmethod=2




;*************
;another option is to fit each line separately and then
;fit a 2D poly function for each hermite coefficient.
;at least it would be nice to see what we're fitting.
;*************

; Create the output LSF parameter array
for i=0,nfibers-1 do begin
  for j=0,2 do begin

    ind = where(fitstr.fiber eq i and fitstr.chip eq j+1,nind)

    if nind gt 0 then begin
      case j of
      0: lsfim1[fitstr[ind].fiber,*] = fitstr[ind].lsfpars
      1: lsfim2[fitstr[ind].fiber,*] = fitstr[ind].lsfpars
      2: lsfim3[fitstr[ind].fiber,*] = fitstr[ind].lsfpars
      endcase
    endif

  endfor
endfor



; How many LSF pixels to use for the LSF
sigind = Horder+4
sigarr = lsfim1[*,sigind]
maxsigma = MAX(sigarr)
nlsfpix = 2*ceil(maxsigma*5)+1  ; odd pixels

ylsf = findgen(nlsfpix)-nlsfpix/2
if n_elements(verbose) eq 0 then verbose = 1


;------------------------------
; Step 5: Make the LSF array
;------------------------------
;print,'Step 5: Create the LSF array'
print,'Creating the LSF array'
; Loop through the Chips
For i=0,2 do begin

  print,'Chip ',chiptag[i]

  ; The LSF parameter array for this chip
  ilsfim = (SCOPE_VARFETCH('lsfim'+strtrim(i+1,2)))

  ; Initialize the lsf array
  lsfarr = fltarr(npix,nfibers,nlsfpix)   ; FLOAT for now

  ; Loop through the fibers
  ;For j=0,nfibers-1 do begin
  for jj=0,n_elements(ifibers)-1 do begin
    j=ifibers[jj]

    if (j+1) mod 50 eq 0 then print,'Fiber ',strtrim(j+1,2),'/',strtrim(long(nfibers),2)

    ; Get the LSF for this chip/fiber
    lsfpars = reform(ilsfim[j,*])

    ; Make 2D LSF array
    xlsf = REPLICATE(1.0d0,npix)#(dindgen(nLSFpix)-nLSFpix/2)
    xlsf += dindgen(npix)#REPLICATE(1.0d0,nLSFpix)
    xcenter = dindgen(npix)
    lsf2d = LSF_GH(xlsf,xcenter,lsfpars)


    ; Stuff it in the array
    lsfarr[*,j,*] = lsf2d

    ; Loop through the pixels
    ;For k=0,npix-1 do begin
    ;
    ;   ; Get the LSF for this chip/fiber/pixel
    ;   lsfpars = reform(ilsfim[j,*])
    ;   lsf = LSF_GH(ylsf+k,k,lsfpars)
    ;
    ;   lsfarr[k,j,*] = lsf    ; stuff it in the array
    ;
    ;End ; pixel loop

    ; if this LSF is discrepant with previous one, replace it
    if j gt 0 then begin
      maxdiff=max(abs(lsfarr[1000,j,*]-lsfarr[1000,j-1,*]))
      if maxdiff gt .03 then begin
        print,'not halted: maxdiff exceeds 0.03 for fiber: ', j,' chip: ', i
        ilsfim[j,*]=ilsfim[j-1,*]
        lsfarr[*,j,*] = lsf_gh(xlsf,xcenter,reform(ilsfim[j,*]))
      endif
    endif

  End ; fiber loop


  ; Creating header for the main data unit, the LSF parameters
  head1 = combframe.(0).header
  sxaddpar,head1,'LAMPTYPE',lamptype
  sxaddpar,head1,'LSFMETHD',fitmethod,' Fitting method: 1-single, 2-chips, 3-all'
  sxaddhist,' ',head1
  sxaddhist,'LSF parameters from APLSF:',head1
  sxaddhist,'The Y-units are in dither combined pixels (starting with 0)',head1
  sxaddhist,'There are LSF parameters for each fiber [Nfibers,Npars]',head1
  sxaddhist,'The parameters are:',head1
  sxaddhist,' 1: Binsize - The width of one pixel',head1
  sxaddhist,' 2: Yoffset - The offset to apply to Y',head1
  sxaddhist,' 3: Horder - The highest Hermite order.  We fix H0=1',head1
  sxaddhist,' 4->Horder+4: Porder - The polynomial order for the global variation',head1
  sxaddhist,'                     with Y for Sigma and the Horder Hermite parameters',head1
  sxaddhist,'                     There are Porder[i]+1 coefficients for parameter i',head1
  sxaddhist,' Horder+5->Horder+4+SUM(Porder+1): GHcoefs - The polynomial coefficients',head1
  sxaddhist,'                       for Sigma and the Horder Hermite paramters',head1
  sxaddhist,'The LSF parameters are to be used with the LSF_GH.PRO function',head1
  sxaddhist,'For example, if you want the LSF at position Y=250 for fiber 10,',head1
  sxaddhist," IDL>fits_read,'apLsf-a-5529.fits',lsfim  ; load the parameters",head1
  sxaddhist,' IDL>pars = lsfim[9,*]                    ; parameters for fiber 10',head1
  sxaddhist,' IDL>y = findgen(21)-10 + 250             ; create pixel array',head1
  sxaddhist,' IDL>lsf = lsf_gh(y,250,pars)             ; get the LSF',head1
  sxaddhist,'',head1
  sxaddhist,'LAMPTYPE = '+lamptype,head1
  sxaddhist,'LSF Fitting Method = '+strtrim(fitmethod,2),head1
  sxaddhist,'1-every line fit separately',head1
  sxaddhist,'2-all lines per chip/fiber fit together',head1
  sxaddhist,'3-all lines per fiber (all three chips) fir together',head1

  ; Header for the first extension with the LSF array
  ;head2 = combframe.(0).header
  MKHDR,head2,lsfarr
  sxaddpar,head2,'CTYPE1','Pixels'
  sxaddpar,head2,'CRPIX1',1
  sxaddpar,head2,'CRVAL1',1
  sxaddpar,head2,'CDELT1',1
  sxaddpar,head2,'CTYPE2','Fibers'
  sxaddpar,head2,'CRPIX2',1
  sxaddpar,head2,'CRVAL2',1
  sxaddpar,head2,'CDELT2',1
  sxaddpar,head2,'CTYPE3','LSF pixels'
  sxaddpar,head2,'CRPIX3',1
  sxaddpar,head2,'CRVAL3',-nlsfpix/2
  sxaddpar,head2,'CDELT3',1
  sxaddhist,' ',head2
  sxaddhist,'LSFs from APLSF:',head2
  sxaddhist,'The LSFs are normalized and are by default to be used with dither',head2
  sxaddhist,'  combined pixels. There is an LSF for each fiber and pixel.',head2
  sxaddhist,'The array is [Npix, Nfibers, Nlsfpix].  There are pixels to cover',head2
  sxaddhist,'approximately +/-5 sigma => '+strtrim(nlsfpix,2)+' pixels.',head2

  ; Output the file
  ;----------------
  ;outfile = lsf_dir+'apLSF-'+chiptag[i]+'-'+strtrim(info[0].mjd5,2)+'.fits'
  ;outfile = lsf_dir+'apLSF-'+chiptag[i]+'-'+strtrim(info[0].mjd5,2)+'-'+file_basename(lampframes[0])+'.fits'
  ;outfile = lsf_dir+'apLSF-'+chiptag[i]+'-'+strtrim(info[0].mjd5,2)+'-'+lampframeid1+'.fits'
  if n_elements(outdir) eq 0 then outdir=lsf_dir
  if file_test(outdir,/directory) eq 0 then file_mkdir,outdir
  if keyword_set(fibers) then  $
  outfile = outdir+dirs.prefix+'LSF-'+chiptag[i]+'-'+lampframeid1+'fibers.fits'  else $
  outfile = outdir+dirs.prefix+'LSF-'+chiptag[i]+'-'+lampframeid1+'.fits'
  print,'Writing LSF file to ',outfile
  FITS_WRITE,outfile,ilsfim,head1                ; Write the parameter array
  MWRFITS,lsfarr,outfile,head2,/silent                ; add the LSF array

  ; with 19 pixels and using float these files are ~92MB

  print,''

  ;stop

End ; chip loop

; Save the final structure
if keyword_set(fibers) then $
savefile = outdir+dirs.prefix+'LSF-'+lampframeid1+'fibers.sav' else $
savefile = outdir+dirs.prefix+'LSF-'+lampframeid1+'.sav'
print,'Saving final structure to ',savefile
save,linestr,fitstr1,fitstr,file=savefile


print,'APLSF FINISHED'
dt = systime(1)-t0_0
print,'Time elapsed = ',strtrim(dt/60.,2),' min.'


;stop

if keyword_set(stp) then stop

end
