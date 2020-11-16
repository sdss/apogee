;+
;
; APWAVECAL
;
; Do the wavelength calibration on a ThAr frame.  Fit all
; three chips together.
;
; INPUTS:
;  lampid       The directory and ID8 string for the ThAr 1D frame concatenated
;  lsfid        The direcotry and ID8 string for the LSF calibration file.
;  =psfid       The directory and ID8 string for the PSF calibration file.
;  =fitmethod   The fitting method to use: (1) separately poly fit to
;                 each fiber (default), (2) 2D polynomial fit to all fibers,
;                 or (3) sine+poly fit with constraints across all
;                 fibers (original method).
;  /refitlines  Refit the lines if the SAVE file already exists,
;                 otherwise the lines in the SAVE file will be used.
;  /pl          Show plots.
;  /save        Save diagnostic plots.
;  /verbose     Verbose output to the screen.
; 
; OUTPUTS:
;  linestr      Structure of all the lines.  One element
;                 for each line on the chip.
;  dispstr      Dispersion function structure.  One element
;                 for each fiber.
;
; USAGE:
;  IDL>apwavecal,lampid8,lsfid,psfid=psfid
;
; By D.Nidever  Feb. 2010
;-

;function pix2wave_sine,x,par
;radeg = 180.0d0/!dpi
;return,par[0]*(sin((X/par[1]+par[2])/radeg)+par[3])
;end
;
;;########################
;
;function pix2wave_poly,x,par
;; NO cubic terms
;return,poly(x,[par[0], par[1], par[2], 0.0,  par[3], par[4]])
;end
;
;;########################

pro apwavecal,lampid,lsfid,coefstr,psfid=psfid,verbose=verbose,pl=pl,refitlines=refitlines,fitmethod=fitmethod,$
              save=save,stp=stp,name=name,medfilt=medfilt


; Error Handling
;------------------
; Establish error handler. When errors occur, the index of the  
; error is returned in the variable Error_status:  
CATCH, Error_status 

;This statement begins the error handler:  
if (Error_status ne 0) then begin 
   error = !ERROR_STATE.MSG  
   if not keyword_set(silent) then print,error
   CATCH, /CANCEL 
   return
endif

;setdisp,/silent

; Get APOGEE directories
dirs=getdir(apogee_dir,cal_dir,spectro_dir,apogee_vers,lib_dir)
psf_dir = spectro_dir+'cal/psf/'
wave_dir = spectro_dir+'cal/wave/'
lsf_dir = spectro_dir+'cal/lsf/'
linelist_dir = lib_dir+'arclines/'

;print,'TEMPORARY KLUDGE WITH THE WAVE_DIR DIRECTORY!!!!'
;wave_dir = '/net/apogee/spectro/t0.91/cal/wave/'

nlampid = n_elements(lampid)

; Not enough inputs
if nlampid eq 0 then begin
  print,'Syntax - apwavecal,lampid,lsfid,psfid=psfid,refitlines=refitlines'
  return
endif

chiptag = ['a','b','c']

; Defaults
if n_elements(pl) eq 0 then pl=0                ; no plotting by default
if n_elements(save) eq 0 then save=1            ; save by default
if n_elements(fitmethod) eq 0 then fitmethod=1  ; separately poly fit to each fiber

if keyword_set(pl) or keyword_set(save) then psym8

print,''
print,'-----------------'
print,'Running APWAVECAL'
print,'-----------------'
print,strtrim(nlampid,2),' Lamp Frames input'
print,''

; Lamp Frame IDs
lampframeid = string(long(file_basename(lampid)),format='(I08)')

; Create the output name
allinfo = APFILEINFO(file_dirname(lampid)+'/'+dirs.prefix+'1D-a-'+lampframeid+'.fits',/silent)
mjd5 = strtrim(long(allinfo[0].mjd5),2)
if keyword_set(name) then outname = string(name,format='(i8.8)') $
  else outname = strjoin(lampframeid,'-')

; Starting the output plot file
if keyword_set(save) then begin
  psfile = ''
  plots_dir = wave_dir+'plots/'
  if FILE_TEST(plots_dir,/directory) eq 0 then FILE_MKDIR,plots_dir
  !p.font = 0
  ;psfile = plots_dir+'apWave_'+outname
  ;file_delete,[psfile+'.ps',psfile+'.psf'],/allow
  ;orig_p=!p & orig_x=!x & orig_y=!y & orig_z=!z
  ;ps_open,psfile,/color,thick=4
  ;device,/inches,xsize=11,ysize=8.5
  ;ps_p=!p & ps_x=!x & ps_y=!y & ps_z=!z
  ;set_plot,'X' & !p=orig_p & !x=orig_x & !y=orig_y & !z=orig_z
endif

; Get the PSF information
if n_elements(psfid) gt 0 then begin

  ; Construct chip names
  psfframeid = string(long(file_basename(psfid)),format='(I08)')
  psffiles = apogee_filename('PSF',chip=chiptag,num=psfframeid)
  tstr1 = MRDFITS(psffiles[0],1,/silent)
  tstr2 = MRDFITS(psffiles[1],1,/silent)
  tstr3 = MRDFITS(psffiles[2],1,/silent)
  tracestr = {tstr1:tstr1, tstr2:tstr2, tstr3:tstr3}

; PSFID not input, check the header
endif else begin

  lampfiles0= apogee_filename('1D',chip=chiptag,num=lampframeid[0])

  info = apfileinfo(lampfiles0,/silent)
  if min(info.exists) eq 0 then begin
    print,lampid[0],' NOT FOUND'
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
  psffiles = apogee_filename('PSF',chip=chiptag,num=psfframeid)
  tstr1 = MRDFITS(psffiles[0],1,/silent)
  tstr2 = MRDFITS(psffiles[1],1,/silent)
  tstr3 = MRDFITS(psffiles[2],1,/silent)
  tracestr = {tstr1:tstr1, tstr2:tstr2, tstr3:tstr3}

endelse

;; Fix PSF structure if some fibers are missing
if n_elements(tstr1) lt 300 then begin
  print,'Less than 300 fibers in PSF file.  Temporary Kludge!'
  ;; Use a PSF file that has 300 fibers/traces
  psffiles_fid = file_dirname(psffiles[0])+'/apPSF-'+['a','b','c']+'-02830055.fits'
  print,'Using fiducial file apPSF-[abc]-02830055.fits'
  ;; CHIP A
  tstr1_fid = MRDFITS(psffiles_fid[0],1,/silent)
  ;; Assume the first few are okay
  mnoff0 = median(tstr1_fid[0:4].gaussy-tstr1[0:4].gaussy)  ; get offset
  ;; Now match them up
  srcor,tstr1_fid.gaussy,tstr1_fid.gaussy*0,tstr1.gaussy+mnoff0,tstr2.gaussy*0,2.0,ind1,ind2,opt=1,/silent,count=nmatch
  if nmatch lt n_elements(tstr1) then print,'Not enough matches'
  ;; Create new trace structure using fiducial as template
  mnoff = median(tstr1_fid[ind1].gaussy-tstr1[ind2].gaussy)
  mnoff_coef0 = median(tstr1_fid[ind1].coef[0]-tstr1[ind2].coef[0])
  tstr1_orig = tstr1   ; the original
  tstr1 = tstr1_fid    ; initialize the new array
  tstr1.peaky -= round(mnoff)   ; offset the quantities from the fiducial
  tstr1.gaussy -= mnoff
  tstr1.gcoef[1] -= mnoff
  tstr1.coef[0] -= mnoff_coef0
  tstr1.modelcoef[0] -= mnoff_coef0
  for i=0,n_tags(tstr1)-1 do tstr1[ind1].(i)=tstr1_orig[ind2].(i)  ; copy over original
  ;; -- CHIP B --
  tstr2_fid = MRDFITS(psffiles_fid[1],1,/silent)
  ;; Assume the first few are okay
  mnoff0 = median(tstr2_fid[0:4].gaussy-tstr2[0:4].gaussy)
  srcor,tstr2_fid.gaussy,tstr2_fid.gaussy*0,tstr2.gaussy+mnoff0,tstr2.gaussy*0,2.0,ind1,ind2,opt=1,/silent,count=nmatch
  if nmatch lt n_elements(tstr2) then print,'Not enough matches'
  ;; Create new trace structure using fiducial as template
  mnoff = median(tstr2_fid[ind1].gaussy-tstr2[ind2].gaussy)
  mnoff_coef0 = median(tstr2_fid[ind1].coef[0]-tstr2[ind2].coef[0])
  tstr2_orig = tstr2   ; the original
  tstr2 = tstr2_fid    ; initialize the new array
  tstr2.peaky -= round(mnoff)   ; offset the quantities from the fiducial
  tstr2.gaussy -= mnoff
  tstr2.gcoef[1] -= mnoff
  tstr2.coef[0] -= mnoff_coef0
  tstr2.modelcoef[0] -= mnoff_coef0
  for i=0,n_tags(tstr2)-1 do tstr2[ind1].(i)=tstr2_orig[ind2].(i)  ; copy over original
  ;; -- CHIP C --
  tstr3_fid = MRDFITS(psffiles_fid[2],1,/silent)
  ;; Assume the first few are okay
  mnoff0 = median(tstr3_fid[0:4].gaussy-tstr3[0:4].gaussy)
  srcor,tstr3_fid.gaussy,tstr3_fid.gaussy*0,tstr3.gaussy+mnoff0,tstr2.gaussy*0,2.0,ind1,ind2,opt=1,/silent,count=nmatch
  if nmatch lt n_elements(tstr3) then print,'Not enough matches'
  ;; Create new trace structure using fiducial as template
  mnoff = median(tstr3_fid[ind1].gaussy-tstr3[ind2].gaussy)
  mnoff_coef0 = median(tstr3_fid[ind1].coef[0]-tstr3[ind2].coef[0])
  tstr3_orig = tstr3   ; the original
  tstr3 = tstr3_fid    ; initialize the new array
  tstr3.peaky -= round(mnoff)   ; offset the quantities from the fiducial
  tstr3.gaussy -= mnoff
  tstr3.gcoef[1] -= mnoff
  tstr3.coef[0] -= mnoff_coef0
  tstr3.modelcoef[0] -= mnoff_coef0
  for i=0,n_tags(tstr3)-1 do tstr3[ind1].(i)=tstr3_orig[ind2].(i)  ; copy over original
  ;; Concatenate the structures into TRACESTR
  tracestr = {tstr1:tstr1, tstr2:tstr2, tstr3:tstr3}
endif

; Check if the savefile exists
savefile = wave_dir+dirs.prefix+'Wave-'+outname+'.dat'
test = file_test(savefile)

; Fitting the lines
if file_test(savefile) eq 0 or keyword_set(refitlines) then begin

  ; Loop through the lamp frames
  ;-------------------------------
  ; This only makes sense if they are different spectra (i.e. Uranium
  ;   and ThArNe)
  apgundef,linestr,dispstr,allinfo
  For i=0,nlampid-1 do begin

    apgundef,lampfiles,info
    apgundef,linestr1,linestr2,linestr3
    apgundef,dispstr1,dispstr2,dispstr3

    ; Construct chip names
    lampfiles = apogee_filename('1D',chip=chiptag,num=lampframeid[i])

    ; Get file info
    info = APFILEINFO(lampfiles,/silent)

    ; Check the files
    ;okay = (info.exists AND info.sp2dfmt AND info.allchips AND info.sz2048 AND (info.naxis eq 3))
    okay = (info.exists AND info.sp1dfmt AND info.allchips AND ((info.naxis eq 3) OR (info.exten eq 1)))
    if total(okay) ne 3 then begin
      print,'Files not okay'
      return
    endif

    ; Load the frame
    APLOADFRAME,lampfiles,frame0

    ; Deal with NANs and ERR=0
    for k=0,2 do begin
      bdnan = where(finite(frame0.(k).flux) eq 0 or finite(frame0.(k).err) eq 0,nbdnan)
      if nbdnan gt 0 then begin
        frame0.(k).flux[bdnan] = 0.0
        frame0.(k).err[bdnan] = baderr()
        frame0.(k).mask[bdnan] = 1
      endif

      bdzero = where(frame0.(k).err le 0,nbdzero)
      if nbdzero gt 0 then begin
        frame0.(k).flux[bdzero] = 0.0
        frame0.(k).err[bdzero] = baderr()
        frame0.(k).mask[bdzero] = 1
      endif
      ;bdmask = where(frame0.(k).mask gt 0,nbdmask)
      bdmask = where(frame0.(k).mask and badmask(),nbdmask)
      frame0.(k).flux[bdmask] = 0.0
      frame0.(k).err[bdmask] = baderr()

      ; median filtering across fibers will really mess up the individual fiber positions!
      ; perhaps this was introduced to deal with missing fibers?
      ;medfilt=5
      if n_elements(medfilt) gt 0 then begin
        f=medfilt2D(frame0.(k).flux,medfilt,dim=2)
        frame0.(k).flux = f
      endif
    endfor

    ; Lamp type
    if info[i].mjd5 ge 55982 and info[i].mjd5 le 55984 then lamptype = 'FFP' else $
    lamptype = ''
    if sxpar(frame0.(0).header,'LAMPUNE') eq 1 then lamptype='URANIUM'
    if sxpar(frame0.(0).header,'LAMPTHAR') eq 1 then lamptype='THARNE'
    lamptype = strupcase(lamptype)

    if lamptype eq '' then begin
      print,'No arc lamp was on for this exposure'
      return
    endif

    ; Add LSF information to the frame structure
    ;--------------------------------------------
    if n_elements(lsfid) gt 0 then begin

      lsfframeid = string(long(file_basename(lsfid[0])),format='(I08)')
      lsffiles = apogee_filename('LSF',chip=chiptag,num=lsfframeid)

      ; Loop through the chips
      for k=0,2 do begin
        ; Get the LSF calibration data
        FITS_READ,lsffiles[k],lsfcoef,lhead
        ; Add to the chip structure
        chstr = frame0.(k)
        chstr = CREATE_STRUCT(temporary(chstr),'LSFFILE',lsffiles[k],'LSFCOEF',lsfcoef)
        ; Now add this to the final FRAME structure
        if k eq 0 then begin
          frame = CREATE_STRUCT('chip'+chiptag[k],chstr)
        endif else begin
          frame = CREATE_STRUCT(frame,'chip'+chiptag[k],chstr)
        endelse
      end

    ; No LSF information
    endif else begin
      frame = frame0
    endelse
    apgundef,frame0   ; free up memory

    print,'Lamp Frame ID = ',lampid[i]
    print,'LAMPTYPE = ',lamptype
    ;print,'Flat Frame ID = ',flatid8

    ; Get the chip data
    print,''
    print,'-------------------------------------------------'
    print,'Getting WAVELENGTH DATA and SOLUTIONS for CHIP a'
    print,'-------------------------------------------------'
    print,''
    APWAVECAL_CHIP,frame.(0),linestr1,dispstr1,pl=pl

    print,''
    print,'-------------------------------------------------'
    print,'Getting WAVELENGTH DATA and SOLUTIONS for CHIP b'
    print,'-------------------------------------------------'
    print,''
    APWAVECAL_CHIP,frame.(1),linestr2,dispstr2,pl=pl

    print,''
    print,'-------------------------------------------------'
    print,'Getting WAVELENGTH DATA and SOLUTIONS for CHIP c'
    print,'-------------------------------------------------'
    print,''
    APWAVECAL_CHIP,frame.(2),linestr3,dispstr3,pl=pl

    ; KLUDGE for missing fiber in the red
    ; just copy values for fiber 299 to fiber 300
    if lampid[i] eq '/net/stream/apogee/data/20110325/00000899' then begin

      ;save,file='/net/stream/apogee/data/20110325/apwavecal_899.dat'
      restore,'/net/stream/apogee/data/20110325/apwavecal_899.dat'

      ind = where(linestr1.fiber eq 298)
      temp300 = linestr1[ind]
      temp300.fiber = 299
      linestr1 = [linestr1,temp300]

      ind = where(dispstr1.fiber eq 298)
      temp300 = dispstr1[ind]
      temp300.fiber = 299
      dispstr1 = [dispstr1,temp300]

      ; Fix the trace
      temp = tracestr
      tstr1 = temp.tstr1
      temp300 = tstr1[298]
      temp300.peaky += 7
      temp300.gaussy += 6.66
      temp300.modelcoef[0] += 7.0288
      temp300.coef[0] += 6.721
      tstr1 = [tstr1,temp300]
      tracestr = {tstr1:tstr1,tstr2:temp.tstr2,tstr3:temp.tstr3}

    endif


    nfibers1 = n_elements(dispstr1)
    nfibers2 = n_elements(dispstr2)
    nfibers3 = n_elements(dispstr3)

    if nfibers2 ne nfibers1 or nfibers3 ne nfibers1 then begin
      print,'Number of fibers in three chips are NOT the same'
      return
    endif

    ;if nfibers1 ne 300 then begin
    ;  print,'ONLY ',strtrim(nfibers1,2),' fibers in chip a.  Need 300'
    ;  return
    ;endif
    ;if nfibers2 ne 300 then begin
    ;  print,'ONLY ',strtrim(nfibers2,2),' fibers in chip b.  Need 300'
    ;  return
    ;endif
    ;if nfibers3 ne 300 then begin
    ;  print,'ONLY ',strtrim(nfibers3,2),' fibers in chip c.  Need 300'
    ;  return
    ;endif
    ;nfibers = 300

    nfibers = nfibers1
    sz = size(frame.(0).flux)
    npix = sz[1]


    ; Add extra tags to the line structures
    ADD_TAG,linestr1,'CHIPNUM',1L,linestr1
    ADD_TAG,linestr2,'CHIPNUM',2L,linestr2
    ADD_TAG,linestr3,'CHIPNUM',3L,linestr3
    frame_linestr = [linestr1,linestr2,linestr3]
    ADD_TAG,frame_linestr,'X',0.0,frame_linestr
    frame_linestr.x = frame_linestr.pixcen
    ADD_TAG,frame_linestr,'LAMPTYPE',lamptype,frame_linestr
    ADD_TAG,frame_linestr,'FRAME',lampframeid[i],frame_linestr

    ; Combine with previous linelist
    PUSH,linestr,frame_linestr

    ; Add extra tags to the dispersion structures
    ADD_TAG,dispstr1,'CHIPNUM',1L,dispstr1
    ADD_TAG,dispstr2,'CHIPNUM',2L,dispstr2
    ADD_TAG,dispstr3,'CHIPNUM',3L,dispstr3
    frame_dispstr = [dispstr1,dispstr2,dispstr3]
    ADD_TAG,frame_dispstr,'LAMPTYPE',lamptype,frame_dispstr
    ADD_TAG,frame_dispstr,'FRAME',lampframeid[i],frame_dispstr

    ; Combine with previous dispersion structures
    PUSH,dispstr,frame_dispstr

    ; Combine the frame info structures
    PUSH,allinfo,info

    ;stop

  Endfor


  ; Get Y-positions from the trace in the middle of the GREEN chip
  ;  ADD TO LINESTR
  ypos1 = fltarr(nfibers)
  for j=0,nfibers-1 do ypos1[j] = poly(0,tstr2[j].coef)
  ADD_TAG,linestr,'FIBER_YPOS1',0.0,linestr
  linestr.fiber_ypos1 = ypos1[linestr.fiber]
  ypos2 = fltarr(nfibers)
  for j=0,nfibers-1 do ypos2[j] = poly(2047,tstr2[j].coef)
  ADD_TAG,linestr,'FIBER_YPOS2',0.0,linestr
  linestr.fiber_ypos2 = ypos2[linestr.fiber]
  ymid = fltarr(nfibers)
  for j=0,nfibers-1 do ymid[j] = poly(1023.5,tstr2[j].coef)
  ADD_TAG,linestr,'FIBER_YMID',0.0,linestr
  linestr.fiber_ymid = ymid[linestr.fiber]

  ADD_TAG,linestr,'YPOS',0.0,linestr
  for j=0L,n_elements(linestr)-1 do $
      linestr[j].ypos = poly(linestr[j].x,(tracestr.(linestr[j].chipnum-1))[linestr[j].fiber].coef)

  indmatch = where(linestr.model_match eq 1 and linestr.model_usewave eq 1,nindmatch)
  mlinestr = linestr[indmatch]  ; lines with matches

; Use the lines in the SAVE file
Endif else begin

  print,'Using lines in SAVE file ',savefile
  restore,savefile

  ; Load the last frame
  lastlampid = first_el(lampid,/last)
  lastlampid1d = file_dirname(lastlampid)+'/'+dirs.prefix+'1D-'+file_basename(lastlampid)
  APLOADFRAME,lastlampid1d,frame
  sz = size(frame.(0).flux)
  npix = sz[1]
  nfibers = sz[2]

  ;apgundef,flinestr,coefstr

Endelse




;============================
; DIFFERENT FITTING METHODS
;============================
CASE fitmethod of

;################################################
;# 1.) Polynomial fits to each Fiber separately
;################################################

1: Begin  
  
  print,''
  print,'Using Fitting Method 1: Polynomial fits to each Fiber separately'
  print,'---------------------------------------------------------------'
  print,''

  if tag_exist(mlinestr,'XB') eq 0 then add_tag,mlinestr,'XB',0.0,mlinestr
  
  ; Step 1 - Get the chipgaps by fitting a polynomial to each fiber
  ;------------------------------------------------------------------
  resid = mlinestr.model_wave*0
  npoly = 8 ;6 ; 4
  parstr = replicate({i:0,pars:fltarr(npoly+2),perror:fltarr(npoly+2),ypos:0.0,sig:0.0},nfibers)
  for i=0,nfibers-1 do begin
    ind = where(mlinestr.fiber eq i,nind)
    if nind gt 0 then begin
      xx = mlinestr[ind].x
      yy = mlinestr[ind].model_wave
      ;err = yy*0.0+1.0
      err = mlinestr[ind].gfit_perror[1] > 0.1 
      chipnum = mlinestr[ind].chipnum
      initpars = [140.0d, 150.,fltarr(npoly)]
      fa = {chipnum:chipnum}
      pars1 = MPFITFUN('func_chipgap_poly',xx,yy,err,initpars,status=status,dof=dof,$
                         functargs=fa,bestnorm=chisq,perror=perror,yfit=yfit,/quiet)  
      yfit1 = func_chipgap_poly(xx,pars1,chipnum=chipnum,xb=xb)
  
      ; Remove outliers and refit
      diff = yy-yfit1
      sig = mad(diff)
      gd = where(abs(diff) lt 2.5*sig,nbd)
  
      fa2 = {chipnum:chipnum[gd]}
      pars2 = MPFITFUN('func_chipgap_poly',xx[gd],yy[gd],err[gd],pars1,status=status2,dof=dof2,$
                       functargs=fa2,bestnorm=chisq2,perror=perror2,yfit=yfit2,/quiet)      
      yfit2 = func_chipgap_poly(xx,pars2,chipnum=chipnum)
      rchisq2 = chisq2/dof2
      diff2 = mlinestr[ind].model_wave-yfit2
      sig2 = mad(diff2)
  
      resid[ind] = yy-yfit2
      parstr[i].i = i
      parstr[i].pars = pars2
      parstr[i].perror = perror2
      parstr[i].ypos = median(mlinestr[ind].ypos)
      parstr[i].sig = mad(yy-yfit2)
    endif else begin
    ;; no lines for this fiber
      parstr[i].i = i
      parstr[i].pars = 999999.
      parstr[i].perror = 999999.
      ;parstr[i].ypos = median(mlinestr[ind].ypos)
      parstr[i].sig = 999999.
    endelse
  endfor
  
  ; fit chipgaps with YPOS
  gdfib = where(parstr.sig lt 1000,ngdfib)
  chipgap1_coef = ap_robust_poly_fit(parstr[gdfib].ypos,parstr[gdfib].pars[0],1)
  chipgap1 = poly(parstr.ypos,chipgap1_coef)
  chipgap2_coef = ap_robust_poly_fit(parstr[gdfib].ypos,parstr[gdfib].pars[1],1)
  chipgap2 = poly(parstr.ypos,chipgap2_coef)
  
  ;plot,parstr.ypos,parstr.pars[0],ps=8,/ysty
  ;oplot,parstr.ypos,chipgap1,co=250
  ;plot,parstr.ypos,parstr.pars[1],ps=8,/ysty
  ;oplot,parstr.ypos,chipgap2,co=250
  
  ; Get XB using the linear fits to the chipgaps
  for i=0,nfibers-1 do begin
    ind = where(mlinestr.fiber eq i,nind)
    xx = mlinestr[ind].x
    chipnum = mlinestr[ind].chipnum
    pars = [chipgap1[i],chipgap2[i],fltarr(5)]
    yfit = func_chipgap_poly(xx,pars,chipnum=chipnum,xb=xb)
    mlinestr[ind].xb = xb
  endfor
  
  ; These residuals are already at 0.033A
  ;  even if we use npoly=4
  print,'Initial poly fits to each fiber. Sig = ',stringize(mad(resid),ndec=4),' A'
  
  
  ; Step 2 - Fit polynomial for each fiber holding chipgaps constant
  ;------------------------------------------------------------------
  resid = mlinestr.model_wave*0
  coefstr = replicate({i:0,chipgap1:0.0,chipgap2:0.0,coef:fltarr(npoly),perror:fltarr(npoly),ypos:0.0,sig:0.0,rms:0.0},nfibers)
  flinestr = mlinestr
  for i=0,nfibers-1 do begin
    ind = where(mlinestr.fiber eq i,nind)
    if nind gt 0 then begin
      xx = mlinestr[ind].x
      yy = mlinestr[ind].model_wave
      ;err = yy*0.0+1.0
      err = mlinestr[ind].gfit_perror[1] > 0.1 
      chipnum = mlinestr[ind].chipnum
      initpars = [chipgap1[i],chipgap2[i],parstr[i].pars[2:*]]
      fa = {chipnum:chipnum}
      parinfo = replicate({fixed:0},n_elements(initpars))
      parinfo[0:1].fixed = 1
      pars1 = MPFITFUN('func_chipgap_poly',xx,yy,err,initpars,status=status,dof=dof,$
                       functargs=fa,bestnorm=chisq,parinfo=parinfo,perror=perror,yfit=yfit,/quiet)  
      yfit1 = func_chipgap_poly(xx,pars1,chipnum=chipnum,xb=xb)
  
      ; Remove outliers and refit
      diff = yy-yfit1
      sig = mad(diff)
      gd = where(abs(diff) lt 2.5*sig,nbd)
  
      fa2 = {chipnum:chipnum[gd]}
      pars2 = MPFITFUN('func_chipgap_poly',xx[gd],yy[gd],err[gd],pars1,status=status2,dof=dof2,$
                       functargs=fa2,bestnorm=chisq2,parinfo=parinfo,perror=perror2,yfit=yfit2,/quiet)      
      yfit2 = func_chipgap_poly(xx,pars2,chipnum=chipnum)
      rchisq2 = chisq2/dof2
      diff2 = mlinestr[ind].model_wave-yfit2
      sig2 = mad(diff2)
  
      ;plot,xb,diff2,ps=8,xs=1,ys=1,tit='Fiber = '+strtrim(i,2)
  
      flinestr[ind].wave_fit = yfit2
  
      resid[ind] = yy-yfit2
      coefstr[i].i = i
      coefstr[i].chipgap1 = chipgap1[i]
      coefstr[i].chipgap2 = chipgap2[i]
      coefstr[i].coef = pars2[2:*]
      coefstr[i].perror = perror2[2:*]
      coefstr[i].ypos = median(mlinestr[ind].ypos)
      coefstr[i].sig = sig2
      coefstr[i].rms = stddev(diff2)
      ;stop
   endif else begin
   ;; No lines for this fiber
      coefstr[i].i = i
      coefstr[i].chipgap1 = chipgap1[i]
      coefstr[i].chipgap2 = chipgap2[i]
      coefstr[i].coef = 999999.
      coefstr[i].perror = 999999.
      ;coefstr[i].ypos = median(mlinestr[ind].ypos)
      coefstr[i].sig = 999999.
      coefstr[i].rms = 999999.
   endelse  

  endfor
  
  print,'Poly fits to each fiber. Sig = ',stringize(mad(resid),ndec=4),' A'
  
  
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
  
  ;stop
  
  ; Output the coefficients array and wavelengths per fiber for each chip
  ;----------------------------------------------------------------------
  mjd5 = strtrim(long(allinfo[0].mjd5),2)
  ;outname = mjd5+'-'+strjoin(lampframeid,'-')
  if keyword_set(name) then outname = name $
    else outname = strjoin(lampframeid,'-')
  if not keyword_set(silent) then $
    print,'Writing output files to = '+wave_dir+dirs.prefix+'Wave-[abc]-'+outname+'.fits'
  
  lampfiles0 = file_dirname(lampid[0])+'/'+dirs.prefix+'1D-'+chiptag+'-'+lampframeid[0]+'.fits'
  
  ; Loop through the chips
  for i=0,2 do begin
  
    outfile = wave_dir+dirs.prefix+'Wave-'+chiptag[i]+'-'+outname+'.fits'
    head0 = headfits(lampfiles0[i])
    ;sxaddpar,head0,'LAMPTYPE',lamptype
  
    ; Put LAMP filenames into header
    for j=0,n_elements(lampid)-1 do $
      sxaddpar,head0,'LAMPFIL'+strtrim(j+1,2),lampid[j]
  
  
    ; Initialize the Output header
    ;------------------------------
    leadstr = 'APWAVECAL: '
    sxaddhist,leadstr+systime(0),head0
    info = GET_LOGIN_INFO()
    sxaddhist,leadstr+info.user_name+' on '+info.machine_name,head0
    sxaddhist,leadstr+'IDL '+!version.release+' '+!version.os+' '+!version.arch,head0
    sxaddhist,leadstr+'Output File:',head0
    sxaddhist,leadstr+' HDU0 - Header only',head0
    sxaddhist,leadstr+' HDU1 - Chip-specific wavelength coefficients [Nfiber,Ncoef]',head0
    sxaddhist,leadstr+' HDU2 - Wavelength arrays [Npix,Nfiber]',head0
    sxaddhist,leadstr+' HDU3 - Full coefficients [Nfiber,Ncoef]',head0
    sxaddhist,leadstr+'',head0
  
    ; Put LAMP filenames into header
    for j=0,n_elements(lampid)-1 do begin
      lampfil = file_dirname(lampid[j])+'/'+dirs.prefix+'1D-a-'+lampframeid[j]+'.fits'
      hd = headfits(lampfil)
      lamptype = ''
      if sxpar(hd,'LAMPUNE') eq 1 then lamptype='URANIUM'
      if sxpar(hd,'LAMPTHAR') eq 1 then lamptype='THARNE'
      lamptype = strupcase(lamptype)
  
      sxaddhist,leadstr+' LAMPFILE'+strtrim(j+1,2)+'='+lampid[j],head0
      sxaddhist,leadstr+'  LAMPTYPE='+lamptype,head0
    endfor
  
    medrms = median(coefstr.rms)
    fsig = mad(flinestr.model_wave-flinestr.wave_fit,/zero)
    sxaddhist,leadstr+' Median RMS='+stringize(medrms,ndec=4),head0
    sxaddhist,leadstr+' Final SIG='+stringize(fsig,ndec=4),head0
  
  
    ; HDU0 - header only
    FITS_WRITE,outfile,0,head0  ; write the coeffient array
    ; HDU1 - chip-specific coefficients
    MKHDR,head1,wavestr[i].coef,/image
    sxaddpar,head1,'CTYPE1','Fiber'
    sxaddpar,head1,'CTYPE2','Coefficient'
    ;sxaddpar,head1,'BUNIT',''
    MWRFITS,wavestr[i].coef,outfile,head1,/silent  ; write the coeffient array
    ; HDU2 - wavelength arrays
    MKHDR,head2,wavestr[i].wave,/image
    sxaddpar,head2,'CTYPE1','Pixel'
    sxaddpar,head2,'CTYPE2','Fiber'
    sxaddpar,head2,'BUNIT','Wavelength (Angstroms)'
    MWRFITS,wavestr[i].wave,outfile,head2,/silent   ; add the wavelength array to the first extension
    ; HDU3 - full coefficients
    MKHDR,head3,transpose(coefstr.coef),/image
    sxaddpar,head3,'CTYPE1','Fiber'
    sxaddpar,head3,'CTYPE2','Coefficient'
    ;sxaddpar,head3,'BUNIT',''
    MWRFITS,transpose(coefstr.coef),outfile,head3,/silent  ; the actual coefficients
  
    ;stop
  
  endfor ; chip loop
  
  
  ; Save the structures
  savefile = wave_dir+dirs.prefix+'Wave-'+outname+'.dat'
  print,'Saving fitting information to ',savefile
  SAVE,linestr,mlinestr,flinestr,dispstr,parstr,coefstr,wavestr,fil=savefile
  
  medrms = median(coefstr.rms)
  print,'Median RMS = ',strtrim(medrms,2)
  
End  ; poly fits to each fiber sepatately




;#######################################################
;# 2.) 2D poly fit to all fibers
;#######################################################


2: Begin

  print,''
  print,'Using Fitting Method 2: 2D poly fit to all fibers'
  print,'-------------------------------------------------'
  print,''

  if tag_exist(mlinestr,'XB') eq 0 then add_tag,mlinestr,'XB',0.0,mlinestr
  
  ; Step 1 - Get the chipgaps by fitting a polynomial to each fiber
  ;------------------------------------------------------------------
  resid = mlinestr.model_wave*0
  npoly = 5 ; 4
  parstr = replicate({i:0,pars:fltarr(npoly+2),ypos:0.0,sig:0.0},nfibers)
  for i=0,nfibers-1 do begin
    ind = where(mlinestr.fiber eq i,nind)
    xx = mlinestr[ind].x
    yy = mlinestr[ind].model_wave
    ;err = yy*0.0+1.0
    err = mlinestr[ind].gfit_perror[1] > 0.1 
    chipnum = mlinestr[ind].chipnum
    initpars = [140.0,150.,fltarr(npoly)]
    fa = {chipnum:chipnum}
    pars = MPFITFUN('func_chipgap_poly',xx,yy,err,initpars,status=status,dof=dof,$
                       functargs=fa,bestnorm=chisq,parinfo=parinfo,perror=perror,yfit=yfit,/quiet)  
    yfit = func_chipgap_poly(xx,pars,chipnum=chipnum,xb=xb)
    resid[ind] = yy-yfit
    parstr[i].i = i
    parstr[i].pars = pars
    parstr[i].ypos = median(mlinestr[ind].ypos)
    parstr[i].sig = mad(yy-yfit)
  endfor
  
  ; fit chipgaps with YPOS
  chipgap1_coef = ap_robust_poly_fit(parstr.ypos,parstr.pars[0],1)
  chipgap1 = poly(parstr.ypos,chipgap1_coef)
  chipgap2_coef = ap_robust_poly_fit(parstr.ypos,parstr.pars[1],1)
  chipgap2 = poly(parstr.ypos,chipgap2_coef)
  
  ;plot,parstr.ypos,parstr.pars[0],ps=8,/ysty
  ;oplot,parstr.ypos,chipgap1,co=250
  ;plot,parstr.ypos,parstr.pars[1],ps=8,/ysty
  ;oplot,parstr.ypos,chipgap2,co=250
  
  ; Get XB using the linear fits to the chipgaps
  for i=0,nfibers-1 do begin
    ind = where(mlinestr.fiber eq i,nind)
    xx = mlinestr[ind].x
    chipnum = mlinestr[ind].chipnum
    pars = [chipgap1[i],chipgap2[i],fltarr(5)]
    yfit = func_chipgap_poly(xx,pars,chipnum=chipnum,xb=xb)
    mlinestr[ind].xb = xb
  endfor
  
  ; These residuals are already at 0.033A
  ;  even if we use npoly=4
  print,'Initial poly fits to each fiber. Sig = ',stringize(mad(resid),ndec=4),' A'
  
  ; Step 2 - Do a 2D global fit
  ;------------------------------
  npars = 21 ;15 ;11 ;8 ;6 ;8 ;11
  
  xx = mlinestr.xb
  yy = mlinestr.ypos
  zz = mlinestr.model_wave
  ;err = zz*0+1
  err = mlinestr.gfit_perror[1] > 0.1  ; this is important!!
  ;bd = where(err eq 0,nbd,comp=gd)
  ;if nbd gt 0 then err[bd] = median(err[gd])
  ;toogood = where(err lt 0.001,ntoogood)
  ;if ntoogood gt 0 then err[toogood]=0.001
  
  ; Use the same YPOS for all lines in a fiber
  for i=0,nfibers-1 do begin
    ind = where(mlinestr.fiber eq i,nind)
    ;yy[ind] = median(mlinestr[ind].ypos)
    yy[ind] = parstr[i].ypos
  endfor
  ; This makes no difference!! FOR NOW
  
  ; Iterate and remove outliers
  count = 0
  endflag = 0
  oldrms = 999999.
  WHILE (endflag eq 0) do begin
  
    apgundef,status,dof,chisq,rchisq,yfit
    initpars = dblarr(npars)
    initpars[0] = 1.0
    parinfo = REPLICATE({limited:[0,0],limits:[0.0,0.0],fixed:0},npars)
    pars1 = MPFIT2DFUN('func_poly2d',xx,yy,zz,err,initpars,status=status1,dof=dof1,$
                      bestnorm=chisq1,parinfo=parinfo,perror=perror1,yfit=yfit1,/quiet)  
    rchisq1 = chisq1/dof1
  
    ; Remove outliers and refit
    diff = zz-yfit1
    sig = mad(diff)
    gd = where(abs(diff) lt 2.5*sig,nbd)
  
    pars2 = MPFIT2DFUN('func_poly2d',xx[gd],yy[gd],zz[gd],err[gd],pars1,status=status2,dof=dof2,$
                       bestnorm=chisq2,parinfo=parinfo,perror=perror2,yfit=yfit2,/quiet)      
    rchisq2 = chisq2/dof2
    yfit2 = func_poly2d(xx,yy,pars2)
    diff2 = mlinestr.model_wave-yfit2
    sig2 = mad(diff2)
    std2 = stddev(diff2)
  
    ;plotc,xx,yy,zz-yfit2,ps=3,min=-0.3,max=0.3,xs=1,ys=1
  
    print,'reduced Chi2 = ',stringize(rchisq2,ndec=5)
    print,'RMS = ',stringize(sig2,ndec=4),' A'
    ;print,'STD = ',stringize(std2,ndec=4),' A'
  
  
    ; Loop through each fiber and find the median residual
    fiber_waveoff = fltarr(nfibers)
    resid = diff2
    zz = mlinestr.model_wave
    for i=0,nfibers-1 do begin
      indall = where(mlinestr.fiber eq i,nindall)
      ind = where(mlinestr.fiber eq i and abs(diff2) lt 3.5*sig2,nind)
      fiber_waveoff[i] = median(mlinestr[ind].model_wave-yfit2[ind])
      resid[indall] -= fiber_waveoff[i]
      zz[indall] -= fiber_waveoff[i]
    endfor
  
    ;plotc,xx,yy,resid,ps=3,min=-0.10,max=0.10,xs=1,ys=1
    rms = mad(resid)
    print,'RMS = ',stringize(rms,ndec=4)
  
    drms = rms-oldrms
    if abs(drms/rms) lt 0.01 or count gt 5 then endflag=1
    oldrms = rms
    count++
  
    ;stop
  
  ENDWHILE
  fpars = pars2
  
  
  ;stop
  
  ; Step 3 - Convert 2D coefficients to 1D coeffs for each fiber
  ;--------------------------------------------------------------
  npoly = 6
  flinestr = mlinestr
  coefstr = REPLICATE({fiber:-1L,nlines:0L,chipgap1:0.0,chipgap2:0.0,coef:dblarr(npoly),$
                       coeferr:dblarr(5),rchisq:99.99,rms:99.99,sig:99.99,niter:0L},nfibers)
  for i=0,nfibers-1 do begin
  
    ;15:  a = p[0] + p[1]*x + p[2]*x^2 + p[3]*x^3 + p[4]*x^4 + p[5]*y + p[6]*x*y + $
    ;       p[7]*(x^2)*y + p[8]*(x^3)*y + p[9]*y^2 + p[10]*x*y^2 + p[11]*(x^2)*y^2 + $
    ;       p[12]*y^3 + p[13]*x*y^3 + p[14]*y^4
  
    ; There are four types of terms: constant, X, Y, cross terms
    ; Y terms become constant terms, cross terms become Y terms
    ; constant terms: p[0], p[5]*y, p[9]*y^2, p[12]*y^3, p[14]*y^4
    ; X^1 terms: p[1]*x, p[6]*x*y, p[10]*x*y^2, p[13]*x*y^3
    ; X^2 terms: p[2]*x^2, p[7]*(x^2)*y, p[11]*(x^2)*y^2
    ; X^3 terms: p[3]*x^3, p[8]*(x^3)*y
    ; X^4 terms: p[4]*x^4
  
    ;ypos = parstr[i].ypos
    ;coef = fltarr(5)
    ;; Add fiber-to-fiber constant wavelength offset
    ;coef[0] = fiber_waveoff[i] + fpars[0] + fpars[5]*ypos + fpars[9]*ypos^2 + fpars[12]*ypos^3 + fpars[14]*ypos^4
    ;coef[1] = fpars[1] + fpars[6]*ypos + fpars[10]*ypos^2 + fpars[13]*ypos^3
    ;coef[2] = fpars[2] + fpars[7]*ypos + fpars[11]*ypos^2
    ;coef[3] = fpars[3] + fpars[8]*ypos
    ;coef[4] = fpars[4]
  
  
    ;21:  a = p[0] + p[1]*x + p[2]*x^2 + p[3]*x^3 + p[4]*x^4 + p[5]*x^5 + p[6]*y + p[7]*x*y + $
    ;         p[8]*(x^2)*y + p[9]*(x^3)*y + p[10]*(x^4)*y + p[11]*y^2 + p[12]*x*y^2 + $
    ;         p[13]*(x^2)*y^2 + p[14]*(x^3)*y^2 + p[15]*y^3 + p[16]*x*y^3 + p[17]*(x^2)*y^3 + $ 
    ;         p[18]*y^4 + p[19]*x*y^4 + p[20]*y^5
  
    ; There are four types of terms: constant, X, Y, cross terms
    ; Y terms become constant terms, cross terms become Y terms
    ; constant terms: p[0], p[6]*y, p[11]*y^2, p[15]*y^3, p[18]*y^4, p[20]*y^5
    ; X^1 terms: p[1]*x,   p[7]*x*y,     p[12]*x*y^2,     p[16]*x*y^3, p[19]*x*y^4
    ; X^2 terms: p[2]*x^2, p[8]*(x^2)*y, p[13]*(x^2)*y^2, p[17]*(x^2)*y^3
    ; X^3 terms: p[3]*x^3, p[9]*(x^3)*y, p[14]*(x^3)*y^2
    ; X^4 terms: p[4]*x^4, p[10]*(x^4)*y
    ; X^5 terms: p[5]*x^5
  
    ypos = parstr[i].ypos
    coef = fltarr(6)
    ; Add fiber-to-fiber constant wavelength offset
    coef[0] = fiber_waveoff[i] + fpars[0] + fpars[6]*ypos + fpars[11]*ypos^2 + fpars[15]*ypos^3 + $
              fpars[18]*ypos^4 + fpars[20]*ypos^5
    coef[1] = fpars[1] + fpars[7]*ypos + fpars[12]*ypos^2 + fpars[16]*ypos^3 + fpars[19]*ypos^4
    coef[2] = fpars[2] + fpars[8]*ypos + fpars[13]*ypos^2 + fpars[17]*ypos^3
    coef[3] = fpars[3] + fpars[9]*ypos + fpars[14]*ypos^2
    coef[4] = fpars[4] + fpars[10]*ypos
    coef[5] = fpars[5]
  
    ind = where(mlinestr.fiber eq i,nind)
  
    chipgap1 = parstr[i].pars[0]
    chipgap2 = parstr[i].pars[1]
    yfit = func_chipgap_poly(mlinestr[ind].x,[chipgap1,chipgap2,coef],chipnum=mlinestr[ind].chipnum,xb=xb)
    ;plot,xb,mlinestr[ind].model_wave-yfit,ps=8,tit='Fiber = '+strtrim(i,2),yr=[-0.2,0.2],xs=1,ys=1
    ;print,mad(mlinestr[ind].model_wave-yfit)
    diff = mlinestr[ind].model_wave-yfit
    rms = stddev(diff)
    sig = mad(diff)

    flinestr[ind].wave_fit = yfit

    coefstr[i].fiber = i
    coefstr[i].nlines = nind
    coefstr[i].chipgap1 = chipgap1
    coefstr[i].chipgap2 = chipgap2
    coefstr[i].coef = coef
    ;coefstr[i].coeferr = coeferr
    ;coefstr[i].rchisq = rchisq
    coefstr[i].rms = rms
    coefstr[i].sig = sig
  
  endfor

  
  ; Get the wavelength coeffcients and arrays
  ;------------------------------------------
  ;wavestr = REPLICATE({coef:dblarr(nfibers,npar-2),wave:dblarr(npix,nfibers)},3)
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
  mjd5 = strtrim(long(allinfo[0].mjd5),2)
  ;outname = mjd5+'-'+strjoin(lampframeid,'-')
  if keyword_set(name) then outname = name $
    else outname = strjoin(lampframeid,'-')
  if not keyword_set(silent) then $
    print,'Writing output files to = '+wave_dir+dirs.prefix+'Wave-[abc]-'+outname+'.fits'
  
  lampfiles0 = file_dirname(lampid[0])+'/'+dirs.prefix+'1D-'+chiptag+'-'+lampframeid[0]+'.fits'
  
  ; Loop through the chips
  for i=0,2 do begin
  
    outfile = wave_dir+dirs.prefix+'Wave-'+chiptag[i]+'-'+outname+'.fits'
    head0 = headfits(lampfiles0[i])
    ;sxaddpar,head0,'LAMPTYPE',lamptype
  
    ; Put LAMP filenames into header
    for j=0,n_elements(lampid)-1 do $
      sxaddpar,head0,'LAMPFIL'+strtrim(j+1,2),lampid[j]
  
  
    ; Initialize the Output header
    ;------------------------------
    leadstr = 'APWAVECAL: '
    sxaddhist,leadstr+systime(0),head0
    info = GET_LOGIN_INFO()
    sxaddhist,leadstr+info.user_name+' on '+info.machine_name,head0
    sxaddhist,leadstr+'IDL '+!version.release+' '+!version.os+' '+!version.arch,head0
    sxaddhist,leadstr+'Output File:',head0
    sxaddhist,leadstr+' HDU0 - Header only',head0
    sxaddhist,leadstr+' HDU1 - Chip-specific wavelength coefficients [Nfiber,Ncoef]',head0
    sxaddhist,leadstr+' HDU2 - Wavelength arrays [Npix,Nfiber]',head0
    sxaddhist,leadstr+' HDU3 - Full coefficients [Nfiber,Ncoef]',head0
    sxaddhist,leadstr+'',head0
    sxaddhist,leadstr+' Fitting Method 2: 2D polynomial fit',head0
    sxaddpar,head0,'FITMETH',fitmethod,' Wavelength fitting method'  

    ; Put LAMP filenames into header
    for j=0,n_elements(lampid)-1 do begin
      lampfil = file_dirname(lampid[j])+'/'+dirs.prefix+'1D-a-'+lampframeid[j]+'.fits'
      hd = headfits(lampfil)
      lamptype = ''
      if sxpar(hd,'LAMPUNE') eq 1 then lamptype='URANIUM'
      if sxpar(hd,'LAMPTHAR') eq 1 then lamptype='THARNE'
      lamptype = strupcase(lamptype)
  
      sxaddhist,leadstr+' LAMPFILE'+strtrim(j+1,2)+'='+lampid[j],head0
      sxaddhist,leadstr+'  LAMPTYPE='+lamptype,head0
    endfor
  
    medrms = median(coefstr.rms)
    fsig = mad(flinestr.model_wave-flinestr.wave_fit,/zero)
    sxaddhist,leadstr+' Median RMS='+stringize(medrms,ndec=4),head0
    sxaddhist,leadstr+' Final SIG='+stringize(fsig,ndec=4),head0
  
  
    ; HDU0 - header only
    FITS_WRITE,outfile,0,head0  ; write the coeffient array
    ; HDU1 - chip-specific coefficients
    MKHDR,head1,wavestr[i].coef,/image
    sxaddpar,head1,'CTYPE1','Fiber'
    sxaddpar,head1,'CTYPE2','Coefficient'
    ;sxaddpar,head1,'BUNIT',''
    MWRFITS,wavestr[i].coef,outfile,head1,/silent  ; write the coeffient array
    ; HDU2 - wavelength arrays
    MKHDR,head2,wavestr[i].wave,/image
    sxaddpar,head2,'CTYPE1','Pixel'
    sxaddpar,head2,'CTYPE2','Fiber'
    sxaddpar,head2,'BUNIT','Wavelength (Angstroms)'
    MWRFITS,wavestr[i].wave,outfile,head2,/silent   ; add the wavelength array to the first extension
    ; HDU3 - full coefficients
    MKHDR,head3,transpose(coefstr.coef),/image
    sxaddpar,head3,'CTYPE1','Fiber'
    sxaddpar,head3,'CTYPE2','Coefficient'
    ;sxaddpar,head3,'BUNIT',''
    MWRFITS,transpose(coefstr.coef),outfile,head3,/silent  ; the actual coefficients
  
  endfor ; chip loop

  
    ; Save the structures
  savefile = wave_dir+dirs.prefix+'Wave-'+outname+'.dat'
  print,'Saving fitting information to ',savefile
  SAVE,linestr,mlinestr,flinestr,dispstr,coefstr,wavestr,fil=savefile

  medrms = median(coefstr.rms)
  print,'Median RMS = ',strtrim(medrms,2)

  ;stop

End ; 2D poly fitting method



;########################################################
;# 3.) Sine+poly fits with constraints across all fibers
;########################################################
; This was the original method

3: Begin

  print,''
  print,'Using Fitting Method 3: Sine+poly fits with constraints across all fibers'
  print,'--------------------------------------------------------------------------'
  print,''

  ; Measuring the Fiber X offsets
  ;-------------------------------
  print,'Measuring the Fiber X offsets'
  
  ; Cross-correlate all of the fibers against each other to get
  ; zero-point shifts
  ;------------------------------------------
  mxshift0 = fltarr(nfibers)
  ;fiber_ref = frame.(1).flux[*,nfibers/2]
  lopix = 500
  hipix = 1800
  fiber_ref = median(frame.(1).flux[lopix:hipix,*],dim=2)
  fiber_ref = MEDFILT1D(fiber_ref,3,/edge)     ; smooth to get rid of bad pixels
  fiber_ref -= MEDFILT1D(fiber_ref,150,/edge) ; get rid of continuum
  sig = mad(fiber_ref)
  bdpix = where(fiber_ref lt 5*sig,nbdpix)
  if nbdpix gt 0 then fiber_ref[bdpix]=0
  ;mask = long(fiber_ref gt 5*sig)
  ;mask = convol(mask,fltarr(21)+1)
  ;mask = mask/(mask>1)
  x = lindgen(npix)
  for i=0,nfibers-1 do begin
    fiber1 = frame.(1).flux[lopix:hipix,i]
    fiber1 = MEDFILT1D(fiber1,3,/edge)    ; smooth to get rid of bad pixels
    fiber1 -= MEDFILT1D(fiber1,150,/edge) ; get rid of continuum
    sig = mad(fiber1)
    bdpix = where(fiber1 lt 5*sig,nbdpix)
    if nbdpix gt 0 then fiber1[bdpix]=0
    ; Now measure the relative shift
    XCORLB,fiber_ref,fiber1,20,xsh  ;,mask=mask
    mxshift0[i] = xsh
  end
  mxshift0 -= median(mxshift0)
  xshift = mxshift0
  
  ; Use the lines to find more accurate Xoffsets
  ;gdlines = where(mlinestr.fiber eq nfibers/2 and mlinestr.chipnum eq 1,ngdlines)
  gdlines = where(mlinestr.fiber eq nfibers/2,ngdlines)
  midlines = mlinestr[gdlines]
  mxshift = fltarr(nfibers)
  for i=0,nfibers-1 do begin
    ;gdlines1 = where(mlinestr.fiber eq i and mlinestr.chipnum eq 1,ngdlines1)
    gdlines1 = where(mlinestr.fiber eq i,ngdlines1)
    x1 = midlines.gaussx - mxshift0[i] + mxshift0[nfibers/2]
    x2 = mlinestr[gdlines1].gaussx
    SRCOR2,x1,x1*0,x2,x2*0.,4,ind1,ind2,opt=1,/silent
    xsh = median([midlines[ind1].gaussx-x2[ind2]])
    mxshift[i] = xsh
  end
  mxshift -= median(mxshift)
  xshift = mxshift

  if keyword_set(pl) then begin
    plot,mxshift0,ps=-1,xtit='Fiber',ytit='X Offset',tit='Fiber offset'
    oplot,mxshift,ps=-4,co=250
  endif
  if keyword_set(save) then begin
    ;set_plot,'PS' & !p=ps_p & !x=ps_x & !y=ps_y & !z=ps_z
    psfile1 = plots_dir+dirs.prefix+'Wave_'+outname+'_xoffset'
    PUSH,psfiles,psfile1
    ps_open,psfile1,thick=4,/color,/encap
    loadct,39,/silent
    plot,mxshift0,ps=-1,xtit='Fiber',ytit='X Offset',tit='Fiber offset',charsize=1.3
    oplot,mxshift,ps=-4,co=250
    ps_close
    ;set_plot,'X' & !p=orig_p & !x=orig_x & !y=orig_y & !z=orig_z
  endif
  
  print,''
  print,'Fitting single solution to all lines and all fibers'
  print,'----------------------------------------------------'
  print,''
  
  
  ; No fitting, get initial structure
  polychipgaps = 2 ;1
  FIT_WAVESOL_ALLFIBERS,mlinestr,outstr0,/allfixed,polychipgaps=polychipgaps,/silent
  outstr0.xoffset = xshift  ; use the middle-chip x-correlation values for the xoffsets
  outstr0.poly = 0           ; start with no poly terms
  
  ; Could use the difference in XSHIFT between the 1st/2nd chips and
  ; 2nd/3rd chips to get first estimate of chipgap variation with Y
  
  
  ; Fit a robust polynomial to all of the lines and remove BAD outliers
  coef0 = AP_ROBUST_POLY_FIT(outstr0.xb,mlinestr.model_wave,5) 
  resid = mlinestr.model_wave - poly(outstr0.xb,coef0)
  if keyword_set(pl) then begin
    plot,outstr0.xb,resid,ps=3,/ysty,xtit='X',ytit='Resid',tit='Poly fit to all lines'
  endif
  sig = MAD(resid) > 0.5
  bad = where(abs(resid) gt 7*sig,nbad)
  if nbad gt 0 then begin
    if keyword_set(pl) then oplot,[outstr0.xb[bad]],[resid[bad]],ps=1,co=250
    print,'Removing ',strtrim(nbad,2),' outliers'
    REMOVE,bad,mlinestr
  endif
  
  
  ; Allow ONLY sine to vary
  print,'Getting initial parameters - Allowing ONLY sine to vary'
  FIT_WAVESOL_ALLFIBERS,mlinestr,outstr1,initstr=outstr0,/fixoffset,/fixpoly,/fixchipgap,polychipgaps=polychipgaps,/silent,pl=pl
  
  ; Allow ONLY chipgaps to vary
  print,'Getting initial parameters - Allowing ONLY chipgaps to vary'
  FIT_WAVESOL_ALLFIBERS,mlinestr,outstr2,initstr=outstr1,/fixoffset,/fixpoly,/fixsine,polychipgaps=polychipgaps,/silent,pl=pl
  
  
  ; Outlier rejection loop
  print,'OUTLIER REJECTION LOOP'
  count = 0
  endflag = 0
  last_outstr = outstr2
  apgundef,outstr1,outstr2
  WHILE (endflag eq 0) do begin
  
    print,'Iteration=',strtrim(count+1,2)
  
    polychipgaps = 2 ;1
  
    ; 1.) First fit
    ; allow sine and chipgaps to vary
    FIT_WAVESOL_ALLFIBERS,mlinestr,outstr1,initstr=last_outstr,/fixoffset,polychipgaps=polychipgaps,/fixpoly,/silent,pl=pl
  
    ; allow poly and chipgaps to vary
    FIT_WAVESOL_ALLFIBERS,mlinestr,outstr2,initstr=outstr1,/fixoffset,polychipgaps=polychipgaps,/fixsine,pl=pl
  
    ; Remove outliers
    ;mlinestr.wave_fit = outstr2.yfit
    diff = mlinestr.model_wave-outstr2.yfit
    gd = where(abs(diff) lt 5*outstr2.sig,ngd)
    mlinestr0 = mlinestr
    mlinestr = mlinestr[gd]
    print,'SIG=',stringize(outstr2.sig,ndec=3),' NGOOD=',strtrim(ngd,2)
  
    nlast = n_elements(mlinestr0)  ; last linelist
    percdiff = abs(float(nlast-ngd)/float(ngd))*100.
    if (nlast eq ngd) or (percdiff lt 0.1) or (count ge 5) then endflag=1
  
    last_outstr = outstr2
  
    count++
  
    ;stop
  
  ENDWHILE
  
  
  ; Final output structure
  foutstr = last_outstr
  
  if keyword_set(save) then begin
    ;set_plot,'PS' & !p=ps_p & !x=ps_x & !y=ps_y & !z=ps_z
    psfile1 = plots_dir+dirs.prefix+'Wave_'+outname+'_allfibers'
    PUSH,psfiles,psfile1
    ps_open,psfile1,thick=4,/color,/encap
    loadct,39,/silent
    plotc,foutstr.xb,mlinestr0.model_wave-foutstr.yfit,foutstr.ypos,ps=1,/ysty,xr=[-3222,3222],xs=1,$
          xtit='X',ytit='Residuals (Ang)',tit='SIG='+stringize(foutstr.sig,ndec=4)+' Ang',charsize=1.3
    oplot,[0,0]-1023.5-75.,[-100,100],linestyle=2
    oplot,[0,0]+1023.5+75.,[-100,100],linestyle=2
    oplot,[-4000,4000],[0,0],linestyle=2
    ps_close
    ;spawn,['convert',psfile1+'.eps',psfile1+'.jpg'],/noshell
    ;spawn,['convert',psfile1+'.eps',psfile1+'.pdf'],/noshell
    ;set_plot,'X' & !p=orig_p & !x=orig_x & !y=orig_y & !z=orig_z
  endif
  
  ;stop
  
  ; Remove lines that are consistently bad
  smodel_wave = strtrim(mlinestr.model_wave,2)
  ui = uniq(smodel_wave,sort(smodel_wave))
  model_linestr = replicate({swave:'',wave:0.0d0,medoff:0.0,sig:0.0,bad:0},n_elements(ui))
  model_linestr.swave = smodel_wave[ui]
  model_linestr.wave = mlinestr[ui].model_wave
  for i=0,n_elements(ui)-1 do begin
    ind = where(smodel_wave eq model_linestr[i].swave,nind)
  
    medoff = median([mlinestr[ind].model_wave-mlinestr[ind].wave_fit])
    sig = mad([mlinestr[ind].model_wave-mlinestr[ind].wave_fit])
  
    model_linestr[i].medoff = medoff
    model_linestr[i].sig = sig
  
  end
  
  sig_medoff = mad(model_linestr.medoff)
  sig_sig = mad(model_linestr.sig)
  med_sig = median(model_linestr.sig)
  
  bdlines = where( (abs(model_linestr.medoff) gt 3.5*sig_medoff) or ( abs(model_linestr.sig-med_sig) gt 3.5*sig_sig),nbdlines)
  if nbdlines gt 0 then model_linestr.bad=1
  print,'Removing ',strtrim(nbdlines,2),' bad lines'
  for i=0,nbdlines-1 do begin
    bdind = where( strtrim(model_linestr.wave,2) eq model_linestr[bdlines[i]].swave,nbdind)
    mlinestr[bdind].model_match = 0
  end
  
  ;stop
  
  
  ; MAKE A LARGE FIGURE AT THE END THAT HAS LOTS OF USEFUL
  ; DIAGNOSTIC PLOTS
  
  ; When plotting up residuals vs. XB color-coded with Ypos
  ; -if there are large differences in the residuals as a function YPOS
  ;   (color):
  ;  red/blue  - rotation
  ;  green     - xoffsets
  
  
  ; FITTING THE CHIPGAPS AND POLY SEPARATELY IS PROBLEMATIC.  fitting
  ; chipgaps only will try to keep the mean zero, but then it doesn't
  ; have the shape for the poly terms to fit out.
  ;
  ;
  
  skiptohere1:
  
  npix = 2048L
  npar = n_elements(foutstr.sine)+n_elements(foutstr.poly)+1+2
  coefstr0 = REPLICATE({fiber:-1L,nlines:0L,coef:dblarr(npar),coeferr:dblarr(npar),rchisq:99.99,rms:99.99,sig:99.99,niter:0L},nfibers)
  print,'Fitting sine components and chipgaps for all fibers individually'
  
  ; Just fit sine and chipgaps, no poly
  for i=0,nfibers-1 do begin
  
    ifiber = i
    if (i+1) mod 50 eq 0 then print,strtrim(i+1,2),'/',strtrim(nfibers,2)
  
    ; Get the lines for this fiber
    gd1 = where(mlinestr.fiber eq ifiber and mlinestr.model_match eq 1,ngd1)
    if ngd1 eq 0 then goto,BOMB
    lines = mlinestr[gd1]
  
    ; The parameters
    ; 1 xoffset
    ; 4 sine parameters
    ; 2 chip gaps
    ; 6 poly parameters
    ;initcoef1 = [-foutstr.xoffset[ifiber], foutstr.sine, foutstr.chipgap1[ifiber],$
    initcoef1 = [foutstr.xoffset[ifiber], foutstr.sine, foutstr.chipgap1[ifiber],$
                 foutstr.chipgap2[ifiber], foutstr.poly]
    ;initcoef1 = lastpars
    parinfo1 = replicate({fixed:0,limited:[0,0],limits:[0.0d0,0.0d0]},n_elements(initcoef1))
    parinfo1[0].fixed = 1     ; fix xoffset for now
    parinfo1[1:4].fixed = 0   ; allow sine params to float
    parinfo1[5:6].fixed = 0   ; allow chip gap params to float
    ;parinfo1[[7,11]].fixed = 1  ; keep constant and cubic terms fixed to ZERO
    ;parinfo1[7:12].fixed = 1  ; fix poly coefficients to default values
    parinfo1[7].fixed = 1     ; fix poly xoffset
    ;parinfo1[12].fixed = 1    ; fix 5th poly term
    parinfo1.limited = 1
    parinfo1.limits[0] = initcoef1 - 0.5*abs(initcoef1)
    parinfo1.limits[1] = initcoef1 + 0.5*abs(initcoef1)
    parinfo1[8:*].fixed = 1
    err = lines.gfit_perror[1]
    bd = where(err eq 0,nbd,comp=gd)
    if nbd gt 0 then err[bd] = median(err[gd])
    toogood = where(err lt 0.001,ntoogood)
    if ntoogood gt 0 then err[toogood]=0.001
    fa = {chipnum:lines.chipnum}
  
    ; Fix sine Xoffset
    ; Fix sine Yoffset to zero
    ; sinpars[0]*( SIN( (XB+sinpars[1])/sinpars[2]/radeg ) + sinpars[3])
    parinfo1[2].fixed = 1
    ;parinfo1[4].fixed = 1
    ;initcoef1[4] = 0.0
  
    pars1 = mpfitfun('fit_pix2wave',lines.x,lines.model_wave,err,initcoef1,yfit=yfit1,parinfo=parinfo1,$
                     functargs=fa,status=status1,perror=perror1,bestnorm=chisq1,dof=dof1,/quiet)
    yfit1 = fit_pix2wave(lines.x,pars1,chipnum=lines.chipnum,xb=xb1)
  
    ;plot,xb1,lines.model_wave-yfit1,ps=8,/ysty
    ;oplot,[-4000,4000],[0,0],linestyle=2
  
    ; Final values
    yfit = yfit1
    perror = perror1 * sqrt(chisq1/dof1)
    rms = STDDEV(lines.model_wave-yfit)
    sig = MAD(lines.model_wave-yfit,/zero)
    rchisq = chisq1/dof1
    nlines = n_elements(lines)
  
      ; Plug values into the structure
    ;-------------------------------
    coefstr0[i].fiber = ifiber
    coefstr0[i].nlines = nlines
    coefstr0[i].coef = pars1
    coefstr0[i].coeferr = perror
    coefstr0[i].rchisq = chisq1/dof1
    coefstr0[i].rms = rms
    coefstr0[i].sig = sig
  
    BOMB:
  
    ;stop
  
  endfor
  
  
  ; Fit the parameters as a function of YPOS
  ;------------------------------------------
  print,'Fitting SINE and CHIPGAP parameters as functions of YPOS'
  ypos = fltarr(nfibers)
  for i=0,nfibers-1 do ypos[i]=poly(npix/2,tstr2[i].coef)
  
  ; Don't fit 0-xoffset. 2-sine xoffset, or poly coefficients
  fitind = [1,3,4,5,6]
  fitorder = 2
  fitcoef = fltarr(n_elements(fitind),fitorder+1)
  for i=0,n_elements(fitind)-1 do $
    fitcoef[i,*] = ROBUST_POLY_FIT(ypos,coefstr0.coef[fitind[i]],fitorder)
  
  ; Plot the fits
  if keyword_set(pl) then begin
    !p.multi=[0,3,2]
    erase
    for i=0,n_elements(fitind)-1 do begin
      plot,ypos,coefstr0.coef[fitind[i]],ps=8,sym=0.7,/xsty,/ysty,$
        tit='Parameter='+strtrim(fitind[i],2),charsize=1.8,xtit='YPOS',ytit='Value'
      oplot,ypos,poly(ypos,fitcoef[i,*]),co=250,thick=2
    endfor
    !p.multi=0
  endif
  ; Save the plot
  if keyword_set(save) then begin
    ;set_plot,'PS' & !p=ps_p & !x=ps_x & !y=ps_y & !z=ps_z
    psfile1 = plots_dir+dirs.prefix+'Wave_'+outname+'_sinegapfit'
    PUSH,psfiles,psfile1
    ps_open,psfile1,thick=4,/color,/encap
    loadct,39,/silent
    !p.multi=[0,3,2]
    ;erase
    for i=0,n_elements(fitind)-1 do begin
      plot,ypos,coefstr0.coef[fitind[i]],ps=8,sym=0.7,/xsty,/ysty,$
        tit='Parameter='+strtrim(fitind[i],2),charsize=1.8,xtit='YPOS',ytit='Value'
      oplot,ypos,poly(ypos,fitcoef[i,*]),co=250,thick=2
    endfor
    !p.multi=0
    ps_close
    ;set_plot,'X' & !p=orig_p & !x=orig_x & !y=orig_y & !z=orig_z
  endif
  
  ; Use fitted values
  coefstr1 = REPLICATE({fiber:0L,nlines:0L,coef:dblarr(npar),coeferr:dblarr(npar),initcoef:dblarr(npar),rchisq:0.0,rms:0.0,sig:0.0,niter:0L},nfibers)
  coefstr1.initcoef[0] = coefstr0.coef[0]
  coefstr1.initcoef[2] = coefstr0.coef[2]
  coefstr1.initcoef[7] = coefstr0.coef[7]
  for i=0,n_elements(fitind)-1 do $
    coefstr1.initcoef[fitind[i]] = poly(ypos,fitcoef[i,*])
  
  ;wait,0.5
  ;stop
  
  
  ; Now keep sine components fixed and only fit polynomial components
  print,'Fitting poly components and fixing sine/chipgaps to poly fitted values'
  
  ; Just fit poly
  for i=0,nfibers-1 do begin
  
    ifiber = i
    if (i+1) mod 50 eq 0 then print,strtrim(i+1,2),'/',strtrim(nfibers,2)
  
    ; Get the lines for this fiber
    gd1 = where(mlinestr.fiber eq ifiber and mlinestr.model_match eq 1,ngd1)
    lines = mlinestr[gd1]
  
    ; The parameters
    ; 1 xoffset
    ; 4 sine parameters
    ; 2 chip gaps
    ; 6 poly parameters
    initcoef1 = coefstr1[i].initcoef
    initcoef1[7:*] = foutstr.poly
  
    parinfo1 = replicate({fixed:0,limited:[0,0],limits:[0.0d0,0.0d0]},n_elements(initcoef1))
    parinfo1[0].fixed = 1     ; fix xoffset
    parinfo1[1:4].fixed = 1   ; fix sine params
    parinfo1[5:6].fixed = 1   ; fix chip gap params
    parinfo1[7].fixed = 1     ; fix poly xoffset
    parinfo1.limited = 1
    parinfo1.limits[0] = initcoef1 - 0.2*abs(initcoef1)
    parinfo1.limits[1] = initcoef1 + 0.2*abs(initcoef1)
    parinfo1[8:*].fixed = 0   ; allow poly components to float
    err = lines.gfit_perror[1]
    bd = where(err eq 0,nbd,comp=gd)
    if nbd gt 0 then err[bd] = median(err[gd])
    toogood = where(err lt 0.001,ntoogood)
    if ntoogood gt 0 then err[toogood]=0.001
    fa = {chipnum:lines.chipnum}
  
    ; Remove the sine component from the Y-values
    initcoef1_sineonly = initcoef1
    initcoef1_sineonly[7:*] = 0.0
    yfit1_sineonly = fit_pix2wave(lines.x,initcoef1_sineonly,chipnum=lines.chipnum,xb=xb1)
    wave_nosine = lines.model_wave - yfit1_sineonly
    initcoef1_nosine = initcoef1
    initcoef1_nosine[1:4] = 0.0
    initcoef1_nosine[3] = 1   ; must be 1 or it will be an error
    parinfo1[1:4].limits[0] = -1
    parinfo1[1:4].limits[1] = 2
  
    ; Use robust polyfit to get initial fit
    npoly = 6
    xb = xb1/3000.   ; need to divide by 3000, that's how pix2wave does it
    pars0 = AP_ROBUST_POLY_FIT(xb,wave_nosine,npoly)
    yfit0 = poly(xb,pars0)
  
    ; Use MPFIT to refine solution
    parinfo1[8:*].limited = 1
    parinfo1[8:*].limits[0] = pars0 - 0.05*abs(pars0)
    parinfo1[8:*].limits[1] = pars0 + 0.05*abs(pars0)
    initcoef1_nosine[8:*] = pars0
    pars1 = mpfitfun('fit_pix2wave',lines.x,wave_nosine,err,initcoef1_nosine,yfit=yfit1,parinfo=parinfo1,$
                     functargs=fa,status=status1,perror=perror1,bestnorm=chisq1,dof=dof1,/quiet)
    yfit1 = fit_pix2wave(lines.x,pars1,chipnum=lines.chipnum,xb=xb1)
  
    ; Reject outliers and refit
    diff = wave_nosine-yfit1
    sig = MAD(diff,/zero)
    gd = where(abs(diff) lt 3*sig,ngd)
    fa2 = {chipnum:lines[gd].chipnum}
    pars2 = mpfitfun('fit_pix2wave',lines[gd].x,wave_nosine[gd],err[gd],pars1,yfit=yfit2,parinfo=parinfo1,$
                     functargs=fa2,status=status2,perror=perror2,bestnorm=chisq2,dof=dof2,/quiet)
    yfit2 = fit_pix2wave(lines.x,pars2,chipnum=lines[gd].chipnum,xb=xb2)
  
    ; Add the sine terms back in
    fpars = pars2  ; pars1
    fpars[1:4] = initcoef1_sineonly[1:4]
  
    ; Final values
    yfit = fit_pix2wave(lines.x,fpars,chipnum=lines.chipnum,xb=xb)
    perror = perror2 * sqrt(chisq2/dof2)
    rms = STDDEV(lines.model_wave-yfit)
    sig = MAD(lines.model_wave-yfit,/zero)
    rchisq = chisq2/dof2
    nlines = n_elements(lines)
  
    ;plot,xb1,wave_nosine,ps=8,/ysty
    ;oplot,xb1,yfit0,ps=1,co=250
    ;oplot,xb1,yfit1,ps=1,co=150
    ;oplot,[-4000,4000],[0,0],linestyle=2
  
      ; Plug values into the structure
    ;-------------------------------
    coefstr1[i].fiber = ifiber
    coefstr1[i].nlines = nlines
    coefstr1[i].coef = fpars
    coefstr1[i].coeferr = perror
    coefstr1[i].rchisq = rchisq
    coefstr1[i].rms = rms
    coefstr1[i].sig = sig
  
    ;stop
  
  endfor
  
  
  ; Fit the parameters as a function of YPOS
  ;------------------------------------------
  print,'Fitting POLYNOMIAL parameters as functions of YPOS'
  ypos = fltarr(nfibers)
  for i=0,nfibers-1 do ypos[i]=poly(npix/2,tstr2[i].coef)
  
  ; Fit poly parameters
  ;fitind2 = [8,9,10,11,12]
  ;fitind2 = [10,11,12,13,14]
  ;fitind2 = [10,11,12,13,14]
  ;fitorder = 4 ;3 ;2
  ;fitcoef2 = fltarr(n_elements(fitind2),fitorder+1)
  ;for i=0,n_elements(fitind2)-1 do $
  ;  fitcoef2[i,*] = ROBUST_POLY_FIT(ypos,coefstr1.coef[fitind2[i]],fitorder)
  
  ; Use fitted values
  coefstr2 = REPLICATE({fiber:0L,nlines:0L,coef:dblarr(npar),coeferr:dblarr(npar),initcoef:dblarr(npar),rchisq:0.0,rms:0.0,sig:0.0,niter:0L},nfibers)
  coefstr2.initcoef[0:7] = coefstr1.coef[0:7]
  for i=0,n_elements(fitind2)-1 do $
    coefstr2.initcoef[fitind2[i]] = poly(ypos,fitcoef2[i,*])
  
  ; Parameter 8 and 9 (constant and linear terms) show higher-order
  ; structure.  Use Median smoothing for them
  par8 = coefstr1.coef[8]
  sm_par8 = MEDFILT1D(par8,31,/edge)
  par9 = coefstr1.coef[9]
  sm_par9 = MEDFILT1D(par9,11,/edge)
  coefstr2.initcoef[8] = sm_par8
  coefstr2.initcoef[9] = sm_par9
  ; 10-14 as well
  par10 = coefstr1.coef[10]
  sm_par10 = MEDFILT1D(par10,31,/edge)
  coefstr2.initcoef[10] = sm_par10
  par11 = coefstr1.coef[11]
  sm_par11 = MEDFILT1D(par11,31,/edge)
  coefstr2.initcoef[11] = sm_par11
  par12 = coefstr1.coef[12]
  sm_par12 = MEDFILT1D(par12,31,/edge)
  coefstr2.initcoef[12] = sm_par12
  par13 = coefstr1.coef[13]
  sm_par13 = MEDFILT1D(par13,31,/edge)
  coefstr2.initcoef[13] = sm_par13
  par14 = coefstr1.coef[14]
  sm_par14 = MEDFILT1D(par14,31,/edge)
  coefstr2.initcoef[14] = sm_par14
  
  
  ; Plot the fits
  if keyword_set(pl) then begin
    !p.multi=[0,4,2]
    erase
    plot,ypos,coefstr1.coef[8],ps=8,sym=0.7,/xsty,/ysty,$
        tit='Parameter=8 MED SMOOTHED',charsize=1.8,xtit='YPOS',ytit='Value'
    oplot,ypos,coefstr2.initcoef[8],co=250,thick=2
    plot,ypos,coefstr1.coef[9],ps=8,sym=0.7,/xsty,/ysty,$
        tit='Parameter=9 MED SMOOTHED',charsize=1.8,xtit='YPOS',ytit='Value'
    oplot,ypos,coefstr2.initcoef[9],co=250,thick=2
    plot,ypos,coefstr1.coef[10],ps=8,sym=0.7,/xsty,/ysty,$
        tit='Parameter=10 MED SMOOTHED',charsize=1.8,xtit='YPOS',ytit='Value'
    oplot,ypos,coefstr2.initcoef[10],co=250,thick=2
    plot,ypos,coefstr1.coef[11],ps=8,sym=0.7,/xsty,/ysty,$
        tit='Parameter=11 MED SMOOTHED',charsize=1.8,xtit='YPOS',ytit='Value'
    oplot,ypos,coefstr2.initcoef[11],co=250,thick=2
    plot,ypos,coefstr1.coef[12],ps=8,sym=0.7,/xsty,/ysty,$
        tit='Parameter=12 MED SMOOTHED',charsize=1.8,xtit='YPOS',ytit='Value'
    oplot,ypos,coefstr2.initcoef[12],co=250,thick=2
    plot,ypos,coefstr1.coef[13],ps=8,sym=0.7,/xsty,/ysty,$
        tit='Parameter=13 MED SMOOTHED',charsize=1.8,xtit='YPOS',ytit='Value'
    oplot,ypos,coefstr2.initcoef[13],co=250,thick=2
    plot,ypos,coefstr1.coef[14],ps=8,sym=0.7,/xsty,/ysty,$
        tit='Parameter=14 MED SMOOTHED',charsize=1.8,xtit='YPOS',ytit='Value'
    oplot,ypos,coefstr2.initcoef[14],co=250,thick=2
  
    for i=0,n_elements(fitind2)-1 do begin
      plot,ypos,coefstr1.coef[fitind2[i]],ps=8,sym=0.7,/xsty,/ysty,$
        tit='Parameter='+strtrim(fitind2[i],2),charsize=1.8,xtit='YPOS',ytit='Value'
      oplot,ypos,coefstr2.initcoef[fitind2[i]],co=250,thick=2
      ;oplot,ypos,poly(ypos,fitcoef2[i,*]),co=250,thick=2
    endfor
    !p.multi=0
  endif
  ; Save the plots
  if keyword_set(save) then begin
    ;set_plot,'PS' & !p=ps_p & !x=ps_x & !y=ps_y & !z=ps_z
    psfile1 = plots_dir+dirs.prefix+'Wave_'+outname+'_polyfit'
    PUSH,psfiles,psfile1
    ps_open,psfile1,thick=4,/color,/encap
    loadct,39,/silent
    device,/inches,xsize=10,ysize=10
    !p.multi=[0,3,3]
    ;erase
    plot,ypos,coefstr1.coef[8],ps=8,sym=0.7,/xsty,/ysty,$
        tit='Parameter=8 MEDIAN SMOOTHED',charsize=1.8,xtit='YPOS',ytit='Value'
    oplot,ypos,coefstr2.initcoef[8],co=250,thick=2
    plot,ypos,coefstr1.coef[9],ps=8,sym=0.7,/xsty,/ysty,$
        tit='Parameter=9 MEDIAN SMOOTHED',charsize=1.8,xtit='YPOS',ytit='Value'
    oplot,ypos,coefstr2.initcoef[9],co=250,thick=2
    plot,ypos,coefstr1.coef[10],ps=8,sym=0.7,/xsty,/ysty,$
        tit='Parameter=10 MEDIAN SMOOTHED',charsize=1.8,xtit='YPOS',ytit='Value'
    oplot,ypos,coefstr2.initcoef[10],co=250,thick=2
    plot,ypos,coefstr1.coef[11],ps=8,sym=0.7,/xsty,/ysty,$
        tit='Parameter=11 MEDIAN SMOOTHED',charsize=1.8,xtit='YPOS',ytit='Value'
    oplot,ypos,coefstr2.initcoef[11],co=250,thick=2
    plot,ypos,coefstr1.coef[12],ps=8,sym=0.7,/xsty,/ysty,$
        tit='Parameter=12 MEDIAN SMOOTHED',charsize=1.8,xtit='YPOS',ytit='Value'
    oplot,ypos,coefstr2.initcoef[12],co=250,thick=2
    plot,ypos,coefstr1.coef[13],ps=8,sym=0.7,/xsty,/ysty,$
        tit='Parameter=13 MEDIAN SMOOTHED',charsize=1.8,xtit='YPOS',ytit='Value'
    oplot,ypos,coefstr2.initcoef[13],co=250,thick=2
    plot,ypos,coefstr1.coef[14],ps=8,sym=0.7,/xsty,/ysty,$
        tit='Parameter=14 MEDIAN SMOOTHED',charsize=1.8,xtit='YPOS',ytit='Value'
    oplot,ypos,coefstr2.initcoef[14],co=250,thick=2
  
    for i=0,n_elements(fitind2)-1 do begin
      plot,ypos,coefstr1.coef[fitind2[i]],ps=8,sym=0.7,/xsty,/ysty,$
        tit='Parameter='+strtrim(fitind2[i],2),charsize=1.8,xtit='YPOS',ytit='Value'
      oplot,ypos,coefstr2.initcoef[fitind2[i]],co=250,thick=2
      ;oplot,ypos,poly(ypos,fitcoef2[i,*]),co=250,thick=2
    endfor
    !p.multi=0
    ps_close
    ;set_plot,'X' & !p=orig_p & !x=orig_x & !y=orig_y & !z=orig_z
  endif
  
  ;wait,0.5
  ;stop
  
  ; There are signs of the fiber-specific Xoffset in the linear
  ; polynomial term.  Why is that happening?
  
  
  ; SHOULD WE TRY TO FIT THE POLYNOMIAL COEFFICIENTS AS A 2D FUNCTION???
  ; OR FIT 1D FUNCTION TO ALL RESIDUALS???
  
  ;for i=0,nfibers-1 do begin
  ;
  ;  ifiber = i
  ;  if (i+1) mod 50 eq 0 then print,strtrim(i+1,2),'/',strtrim(nfibers,2)
  ;
  ;  ; Get the lines for this fiber
  ;  gd1 = where(mlinestr.fiber eq ifiber and mlinestr.model_match eq 1,ngd1)
  ;  lines = mlinestr[gd1]
  ;
  ;  pars = coefstr2[i].initcoef
  ;  pars_sineonly = pars
  ;  pars_sineonly[7:*] = 0.0
  ;  yfit1 = fit_pix2wave(lines.x,pars_sineonly,chipnum=lines.chipnum,xb=xb1)
  ;
  ;  resid = lines.model_wave-yfit1
  ;  PUSH,xball,xb1
  ;  PUSH,residall,resid
  ;  PUSH,fiberall,resid*0+i
  ;
  ;end
  
  
  print,''
  print,'Final tweaking of individual fibers'
  print,'------------------------------------'
  print,''
  
  apgundef,flinestr
  
  ; Loop through each fiber
  FOR i=0,nfibers-1 do begin
  
    ifiber = i
    if (i+1) mod 50 eq 0 then print,strtrim(i+1,2),'/',strtrim(nfibers,2)
  
    ; Get the lines for this fiber
    gd1 = where(mlinestr.fiber eq ifiber and mlinestr.model_match eq 1,ngd1)
    lines = mlinestr[gd1]
  
  
    ; From Hardware Technical Description document
    ; chip a: 15140.0 - 15805.7 A
    ; chip b: 15856.5 - 16433.3 A
    ; chip c: 16476.4 - 16956.8 A
    ;
    ; 2.9mm minimum spacing between the chips
    ; 18 micron pixels
  
  
    ; The parameters are:
    ;  4 sine parameters
    ;  2 chip gaps
    ;  7 poly parameters (first one is a zero-point offset)
    ;;  The chip number must also be input (1, 2 or 3)
    ;;initsine = [ 10824.783d0, 17110.9d0, 365.15342d0,  0.65602507d0]
    ;initsine = [ 10824.783d0, -17110.9d0-6800., -365.15342d0,  0.65602507d0]
    ;initchipgaps = [ 147.83, 139.50 ]
    ;;initpoly = [-3088.4787d0, -0.046553049d0, 1.6631944d-07, 3.7208707d-08, 0.0, -3.3738407d-15, -1.9423719d-19]
    ;initpoly = [-3088.4787d0, -0.046553049d0, -1.6631944d-07, 3.7208707d-08, 0.0, -3.3738407d-15, 1.9423719d-19]
    ;initpars = [ initsine, initchipgaps, initpoly ]
  
  
  
    ; The parameters
    ; 1 xoffset
    ; 4 sine parameters
    ; 2 chip gaps
    ; 6 poly parameters
    ;initcoef1 = [foutstr.xoffset[ifiber], foutstr.sine, foutstr.chipgap1[ifiber],$
    ;             foutstr.chipgap2[ifiber], foutstr.poly]
    ;initcoef1 = lastpars
    initcoef1 = coefstr2[i].initcoef
    parinfo1 = replicate({fixed:0,limited:[0,0],limits:[0.0d0,0.0d0]},n_elements(initcoef1))
    parinfo1[0].fixed = 1     ; fix xoffset for now
    parinfo1[1:4].fixed = 0   ; allow sine params to float
    parinfo1[5:6].fixed = 0   ; allow chip gap params to float
    ;parinfo1[[7,11]].fixed = 1  ; keep constant and cubic terms fixed to ZERO
    ;parinfo1[7:12].fixed = 1  ; fix poly coefficients to default values
    parinfo1[2].fixed = 1     ; fix sine xoffset
    parinfo1[7].fixed = 1     ; fix poly xoffset
    ;parinfo1[12].fixed = 1    ; fix 5th poly term
    parinfo1[8:*].fixed = 0    ; allow poly terms to float
    parinfo1.limited = 1
    parinfo1.limits[0] = initcoef1 - 0.05*abs(initcoef1)
    parinfo1.limits[1] = initcoef1 + 0.05*abs(initcoef1)
    err = lines.gfit_perror[1]
    bd = where(err eq 0,nbd,comp=gd)
    if nbd gt 0 then err[bd] = median(err[gd])
    toogood = where(err lt 0.001,ntoogood)
    if ntoogood gt 0 then err[toogood]=0.001
    fa = {chipnum:lines.chipnum}
  
    parinfo1[0].fixed = 1     ; fix xoffset
    parinfo1[2].fixed = 1
    parinfo1[5:6].fixed = 1   ; fix chip gap params
    parinfo1[7].fixed = 1     ; fix poly xoffset
  
    parinfo1[1:4].fixed = 1  ; keep sine components fixed
    parinfo1[11:*].fixed = 1 ; keep higher poly terms fixed
  
  
    ;parinfo1[0].fixed = 0     ; allow xoffset to vary
  
    ; testing
    parinfo1[0].fixed = 0     ; fix xoffset for now
    parinfo1[1:4].fixed = 0   ; allow sine params to float
    parinfo1[5:6].fixed = 1   ; fix chip gap params
    parinfo1[2].fixed = 1     ; fix sine xoffset
    parinfo1[7].fixed = 1     ; fix poly xoffset
    parinfo1[8:*].fixed = 0    ; allow poly terms to float
    parinfo1.limited = 1
    parinfo1.limits[0] = initcoef1 - 0.2*abs(initcoef1)
    parinfo1.limits[1] = initcoef1 + 0.2*abs(initcoef1)
  
  
    ; ONLY ALLOW LOW-ORDER POLY TERMS TO FLOAT!!!!
  
  
    ; First fit
    pars1 = mpfitfun('fit_pix2wave',lines.x,lines.model_wave,err,initcoef1,yfit=yfit1,parinfo=parinfo1,$
                     functargs=fa,status=status1,perror=perror1,bestnorm=chisq1,/quiet)
    yfit1 = fit_pix2wave(lines.x,pars1,chipnum=lines.chipnum,xb=xb1)
  
    ;plot,xb1,lines.model_wave-yfit1,ps=8
    ;oplot,[-4000,4000],[0,0],linestyle=2
  
    ;stop
  
  
    ; Remove outliers
    diff = lines.model_wave-yfit1
    sig = MAD(diff,/zero)
    gd = where(abs(diff) lt 3*sig,ngd)
    lines2 = lines[gd]
  
    ; Second Fit
    initcoef2 = pars1
    parinfo2 = parinfo1
    parinfo2.limits[0] = initcoef2 - 0.1*abs(initcoef2)
    parinfo2.limits[1] = initcoef2 + 0.1*abs(initcoef2)
  
    parinfo2[0].fixed = 1
    parinfo2[2].fixed = 1
    parinfo2[5:6].fixed = 1
    parinfo2[7].fixed = 1
  
    err2 = lines2.gfit_perror[1]
    bd = where(err2 eq 0,nbd,comp=gd)
    if nbd gt 0 then err2[bd] = median(err2[gd])
    toogood = where(err2 lt 0.001,ntoogood)
    if ntoogood gt 0 then err2[toogood]=0.001
    fa2 = {chipnum:lines2.chipnum}
    pars2 = mpfitfun('fit_pix2wave',lines2.x,lines2.model_wave,err2,initcoef2,yfit=yfit2,parinfo=parinfo2,$
                     functargs=fa2,status=status2,perror=perror2,bestnorm=chisq2,dof=dof2,/quiet)
    yfit2 = fit_pix2wave(lines2.x,pars2,chipnum=lines2.chipnum,xb=xb2)
  
    ; Final values
    xb = xb2
    yfit = yfit2
    fpars = pars2
    perror = perror2 * sqrt(chisq2/dof2)
    rms = STDDEV(lines2.model_wave-yfit2)
    sig = MAD(lines2.model_wave-yfit2,/zero)
    rchisq = chisq2/dof2
    nlines = n_elements(lines2)
  
    ; Plotting
    ;ploterror,xb,lines2.model_wave-yfit3,err2*0,err2,ps=8,sym=1,xtit='X',ytit='Resid (Ang)',$
    ;          yr=[-0.4,0.4],ys=1,/nohat,tit='Final Wavelength Solution for Fiber '+strtrim(i+1,2)+' RMS='+stringize(rms,ndec=3)+' Ang'
    ;;plotc,xb,lines2.model_wave-yfit3,err2,ps=8
    ;oplot,[-4000,4000],[0,0],linestyle=2
    ;;oplot,xb2,poly(xb3,coef_resid),ps=1,co=250
  
    ; Plug values into the structure
    ;-------------------------------
    coefstr2[i].fiber = ifiber
    coefstr2[i].nlines = nlines
    coefstr2[i].coef = fpars
    coefstr2[i].coeferr = perror
    coefstr2[i].rchisq = rchisq
    coefstr2[i].rms = rms
    coefstr2[i].sig = sig
  
    ;stop
  
    ; Sine ONLY part
    fpars_sineonly = fpars
    fpars_sineonly[7:*] = 0.0
    yfit_sineonly = fit_pix2wave(lines2.x,fpars_sineonly,chipnum=lines2.chipnum)
    lines2.wave_fit = yfit
    if tag_exist(lines2,'XB') eq 0 then ADD_TAG,lines2,'XB',0.0,lines2
    lines2.xb = xb
  
    if keyword_set(verbose) then $
      print,'Initial RMS = ',strtrim(rms,2)
  
  
    ; Keep track of the fits
    PUSH,flinestr,lines2
  
    ; Print information
    if not keyword_set(silent) then begin
      form = '(I3,I5,F6.3,15G9.2)'
      print,ifiber+1,nlines,rms,fpars,format=form
    endif
  
  
    xoffset = fpars[0]
    chipgap1 = fpars[5]
    chipgap2 = fpars[6]  
  
    ; Plotting the solutions
    ;-----------------------
    if keyword_set(pl) or keyword_set(save) then begin
  
      co = 255
      ;if keyword_set(save) then begin
      ;  plots_dir = wave_dir+'plots/'
      ;  if FILE_TEST(plots_dir,/directory) eq 0 then FILE_MKDIR,plots_dir
      ;  file = plots_dir+'apwavecal_'+lampframeid+'_fiber'+strtrim(ifiber+1,2)
      ;  ps_open,file,/color,thick=4
      ;  device,/inches,xsize=10,ysize=8
      ;  co = 0
      ;endif
   
      if keyword_set(save) and ifiber eq nfibers/2 then begin
        ;set_plot,'PS' & !p=ps_p & !x=ps_x & !y=ps_y & !z=ps_z
        psfile1 = plots_dir+dirs.prefix+'Wave_'+outname+'_fiberfit'+strtrim(ifiber+1,2)
        PUSH,psfiles,psfile1
        ps_open,psfile1,thick=4,/color,/encap
        loadct,39,/silent
        device,/inches,xsize=8,ysize=11
        co = 0
      endif
  
      ; Total fit
      ;xx = findgen(6432)
      xx = cgscalevector(findgen(6432),-2198,4246)
      ;xx = cgscalevector(findgen(6432),-3222,3222)
      ff = fit_pix2wave(xx,fpars,chipnum=xx*0+2,xb=xx2)  ; all chip=2 since we are inputting XB
      yr1 = [min([lines2.model_wave,ff]),max([lines2.model_wave,ff])]
      yr1 = [yr1[0]-0.1*range(yr1),yr1[1]+0.1*range(yr1)]/1e4
  
      !p.multi = [0,1,3]
      ;xr = [0,6432]
      xr = minmax(xx2)
      charsize = 2.3
      
      if (keyword_set(save) and ifiber eq nfibers/2) or keyword_set(pl) then begin
      plot,[0],[0],/nodata,ps=8,xr=xr,yr=yr1,xs=1,ys=1,xticklen=0.04,charsize=charsize,$
           xtit='X (pix)',ytit='Wavelength (Microns)',tit='Lamp '+strjoin(lampframeid,'-')+' Fiber '+strtrim(ifiber+1,2)
      ;oplot,lines.xb,yfit1/1e4,co=250
      oplot,xx2,ff/1e4,co=250
      oplot,xb,lines2.model_wave/1e4,ps=8,sym=0.8
  
      ; plot the chipgaps
      oplot,[0,0]-1023.5-chipgap1,[-10,10],linestyle=1
      oplot,[0,0]-1023.5,[-10,10],linestyle=1
      oplot,[0,0]+1023.5,[-10,10],linestyle=1
      oplot,[0,0]+1023.5+chipgap2,[-10,10],linestyle=1
  
      legend,['Data','Total Fit'],textcolor=[co,250],charsize=1.3,/top,/left
      endif
  
  
      ; Polynomial Component - Residuals 
      ff_sineonly = fit_pix2wave(xx,fpars_sineonly,chipnum=xx*0+2) ; all chip=2 since we are inputting XB
      ff_resid = ff-ff_sineonly
      wdiff = lines2.model_wave - yfit_sineonly
  
      yrange = max(wdiff)-min(wdiff)
      yr2 = [min(wdiff)-0.1*yrange,max(wdiff)+0.1*yrange]
      if (keyword_set(save) and ifiber eq nfibers/2) or keyword_set(pl) then begin
      plot,[0],[0],/nodata,ps=1,xr=xr,yr=yr2,xs=1,ys=1,xticklen=0.04,xtit='X (pix)',$
           ytit='Wavelength (Ang)',tit='Polynomial Component',charsize=charsize
      oplot,[-4000,4000],[0,0],linestyle=2
      oplot,xb,yfit-yfit_sineonly,ps=1,co=250
      oplot,xx2,ff_resid,co=250
      oplot,xb,wdiff,ps=8,sym=0.8
  
      ; plot the chipgaps
      oplot,[0,0]-1023.5-chipgap1,[-10,10],linestyle=1
      oplot,[0,0]-1023.5,[-10,10],linestyle=1
      oplot,[0,0]+1023.5,[-10,10],linestyle=1
      oplot,[0,0]+1023.5+chipgap2,[-10,10],linestyle=1
      endif
  
      ; True Residuals
      yrange = max(wdiff)-min(wdiff)
      ;yr2 = [min(wdiff)-0.1*yrange,max(wdiff)+0.1*yrange]
      ;yr3 = [-0.05,0.05]
      yr3 = minmax(lines2.model_wave-yfit)
      if (keyword_set(save) and ifiber eq nfibers/2) or keyword_set(pl) then begin
      plot,[0],[0],/nodata,ps=1,xr=xr,yr=yr3,xs=1,ys=1,xticklen=0.04,xtit='X (pix)',$
           ytit='Wavelength (Ang)',tit='Residuals  RMS = '+strtrim(string(rms,format='(F10.5)'),2)+$
           ' A',charsize=charsize
      oplot,[-4000,4000],[0,0],linestyle=2
      oplot,xb,lines2.model_wave-yfit,ps=8,sym=0.8
  
      ; plot the chipgaps
      oplot,[0,0]-1023.5-chipgap1,[-10,10],linestyle=1
      oplot,[0,0]-1023.5,[-10,10],linestyle=1
      oplot,[0,0]+1023.5,[-10,10],linestyle=1
      oplot,[0,0]+1023.5+chipgap2,[-10,10],linestyle=1
  
      !p.multi = 0
      endif
  
      if keyword_set(save) and ifiber eq nfibers/2 then begin
        ;set_plot,'X' & !p=orig_p & !x=orig_x & !y=orig_y & !z=orig_z
        ps_close
      endif
  
      ;if keyword_set(save) then begin
      ;  ps_close
      ;  ps2gif,file+'.ps',rot=-90
      ;endif
  
      ;wait,0.3
  
    endif
  
    ;wait,0.5  ;0.2
    ;stop
  
  ENDFOR
  
  ; Fixing fibers with bad coefficients
  medrms = median(coefstr2.rms)
  sigrms = mad(coefstr2.rms)
  bdfibers = where(coefstr2.rms gt 6*sigrms+medrms,nbdfibers,comp=gdfibers)
  
  if nbdfibers gt 0 then begin
    print,'Fixing coefficients for ',strtrim(nbdfibers,2),' fibers'
  
    flinestr_orig = flinestr
    coefstr2_orig = coefstr2
  
    ; Use median smoothing to fix bad coefficients
    coef0 = coefstr2.coef
    coef0[bdfibers,*] = !values.f_nan
    coef1 = MEDFILT2D(coef0,15,dim=2,/edge)  ;31
    coefstr2[bdfibers].coef = coef1[*,bdfibers]
  
    ; Fix flinestr and coefstr2
    x = findgen(2048)
    for i=0,nbdfibers-1 do begin
  
      ; Get new fitted values
      bdlines = where(flinestr.fiber eq bdfibers[i],nbdlines)
      yfit = fit_pix2wave(flinestr[bdlines].x,coefstr2[bdfibers[i]].coef,chipnum=flinestr[bdlines].chipnum)
      flinestr[bdlines].wave_fit = yfit
  
      rms = STDDEV(flinestr[bdlines].model_wave-flinestr[bdlines].wave_fit)
      coefstr2[bdfibers[i]].rms = rms
  
    end
  
  endif  ; fix bad fibers
  
  fsig = mad(flinestr.model_wave-flinestr.wave_fit,/zero)
  print,'Final SIG=',stringize(fsig,ndec=3)
  
  ; Plot final residuals
  if keyword_set(pl) then begin
    plotc,flinestr.xb,flinestr.model_wave-flinestr.wave_fit,flinestr.ypos,ps=1,/ysty,xr=[-3222,3222],xs=1
    oplot,[0,0]-1023.5-75.,[-100,100],linestyle=2
    oplot,[0,0]+1023.5+75.,[-100,100],linestyle=2
    oplot,[-4000,4000],[0,0],linestyle=2
  endif
  if keyword_set(save) then begin
    ;set_plot,'PS' & !p=ps_p & !x=ps_x & !y=ps_y & !z=ps_z
    psfile1 = plots_dir+dirs.prefix+'Wave_'+outname+'_finalresid'
    PUSH,psfiles,psfile1
    ps_open,psfile1,thick=4,/color,/encap
    loadct,39,/silent
    plotc,flinestr.xb,flinestr.model_wave-flinestr.wave_fit,flinestr.ypos,ps=1,/ysty,xr=[-3222,3222],xs=1,$
          xtit='X',ytit='Resid (Ang)',tit='Final Residuals SIG='+stringize(fsig,ndec=3)+' Ang',charsize=1.3
    oplot,[0,0]-1023.5-75.,[-100,100],linestyle=2
    oplot,[0,0]+1023.5+75.,[-100,100],linestyle=2
    oplot,[-4000,4000],[0,0],linestyle=2
    ps_close
    ;set_plot,'X' & !p=orig_p & !x=orig_x & !y=orig_y & !z=orig_z
  endif
  
  
  
  ; Final structure
  coefstr = coefstr2
  
  
  ; Get the wavelength coeffcients and arrays
  ;------------------------------------------
  wavestr = REPLICATE({coef:dblarr(nfibers,npar-2),wave:dblarr(npix,nfibers)},3)
  x1 = findgen(npix)
  for i=0,nfibers-1 do begin
  
    pars = coefstr[i].coef
  
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
    xoffset = pars[0]
    chipgap1 = pars[5]
    chipgap2 = pars[6]  
    coef_chip1 = [-1023.5-2048-chipgap1+xoffset, pars[1:4], pars[7:*] ]
    coef_chip2 = [-1023.5+xoffset, pars[1:4], pars[7:*] ]
    coef_chip3 = [-1023.5+2048+chipgap2+xoffset, pars[1:4], pars[7:*] ]
  
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
  mjd5 = strtrim(long(allinfo[0].mjd5),2)
  ;outname = mjd5+'-'+strjoin(lampframeid,'-')
  if keyword_set(name) then outname = name $
    else outname = strjoin(lampframeid,'-')
  if not keyword_set(silent) then $
    print,'Writing output files to = '+wave_dir+dirs.prefix+'Wave-[abc]-'+outname+'.fits'
  
  lampfiles0 = file_dirname(lampid[0])+'/'+dirs.prefix+'1D-'+chiptag+'-'+lampframeid[0]+'.fits'
  
  ; Loop through the chips
  for i=0,2 do begin
  
    outfile = wave_dir+dirs.prefix+'Wave-'+chiptag[i]+'-'+outname+'.fits'
    head0 = headfits(lampfiles0[i])
    ;sxaddpar,head0,'LAMPTYPE',lamptype
  
    ; Put LAMP filenames into header
    for j=0,n_elements(lampid)-1 do $
      sxaddpar,head0,'LAMPFIL'+strtrim(j+1,2),lampid[j]
  
  
    ; Initialize the Output header
    ;------------------------------
    leadstr = 'APWAVECAL: '
    sxaddhist,leadstr+systime(0),head0
    info = GET_LOGIN_INFO()
    sxaddhist,leadstr+info.user_name+' on '+info.machine_name,head0
    sxaddhist,leadstr+'IDL '+!version.release+' '+!version.os+' '+!version.arch,head0
    sxaddhist,leadstr+'Output File:',head0
    sxaddhist,leadstr+' HDU0 - Header only',head0
    sxaddhist,leadstr+' HDU1 - Chip-specific wavelength coefficients [Nfiber,Ncoef]',head0
    sxaddhist,leadstr+' HDU2 - Wavelength arrays [Npix,Nfiber]',head0
    sxaddhist,leadstr+' HDU3 - Full coefficients [Nfiber,Ncoef]',head0
    sxaddhist,leadstr+'',head0
    sxaddhist,leadstr+' Fitting Method 3: Sine+poly fit',head0
    sxaddpar,head0,'FITMETH',fitmethod,' Wavelength fitting method'  

    ; Put LAMP filenames into header
    for j=0,n_elements(lampid)-1 do begin
      lampfil = file_dirname(lampid[j])+'/'+dirs.prefix+'1D-a-'+lampframeid[j]+'.fits'
      hd = headfits(lampfil)
      lamptype = ''
      if sxpar(hd,'LAMPUNE') eq 1 then lamptype='URANIUM'
      if sxpar(hd,'LAMPTHAR') eq 1 then lamptype='THARNE'
      lamptype = strupcase(lamptype)
  
      sxaddhist,leadstr+' LAMPFILE'+strtrim(j+1,2)+'='+lampid[j],head0
      sxaddhist,leadstr+'  LAMPTYPE='+lamptype,head0
    endfor
  
    medrms = median(coefstr.rms)
    fsig = mad(flinestr.model_wave-flinestr.wave_fit,/zero)
    sxaddhist,leadstr+' Median RMS='+stringize(medrms,ndec=4),head0
    sxaddhist,leadstr+' Final SIG='+stringize(fsig,ndec=4),head0
  
  
    ; HDU0 - header only
    FITS_WRITE,outfile,0,head0  ; write the coeffient array
    ; HDU1 - chip-specific coefficients
    MKHDR,head1,wavestr[i].coef,/image
    sxaddpar,head1,'CTYPE1','Fiber'
    sxaddpar,head1,'CTYPE2','Coefficient'
    ;sxaddpar,head1,'BUNIT',''
    MWRFITS,wavestr[i].coef,outfile,head1,/silent  ; write the coeffient array
    ; HDU2 - wavelength arrays
    MKHDR,head2,wavestr[i].wave,/image
    sxaddpar,head2,'CTYPE1','Pixel'
    sxaddpar,head2,'CTYPE2','Fiber'
    sxaddpar,head2,'BUNIT','Wavelength (Angstroms)'
    MWRFITS,wavestr[i].wave,outfile,head2,/silent   ; add the wavelength array to the first extension
    ; HDU3 - full coefficients
    MKHDR,head3,transpose(coefstr.coef),/image
    sxaddpar,head3,'CTYPE1','Fiber'
    sxaddpar,head3,'CTYPE2','Coefficient'
    ;sxaddpar,head3,'BUNIT',''
    MWRFITS,transpose(coefstr.coef),outfile,head3,/silent  ; the actual coefficients
  
  endfor ; chip loop
  
  ; Close the psfile
  ;if keyword_set(save) then begin
  ;  set_plot,'PS'
  ;  ps_close
  ;  spawn,['ps2pdf',psfile+'.ps',psfile+'.pdf'],/noshell
  ;  print,'Converting ',psfile,'.ps to ',psfile,'.pdf'
  ;endif
  
  ; Convert the figures
  if keyword_set(save) then begin
    print,'Converting figures'
    for i=0,n_elements(psfiles)-1 do begin
      FILE_DELETE,psfiles[i]+['.jpg','.pdf','.eps.gz'],/allow
      spawn,['convert',psfiles[i]+'.eps',psfiles[i]+'.jpg'],out,errout,/noshell
      spawn,['convert',psfiles[i]+'.eps',psfiles[i]+'.pdf'],out,errout,/noshell
      spawn,['gzip',psfiles[i]+'.eps'],out,errout,/noshell   ; compress the files
    endfor
  
    ; Combine them into one PDF
    comb_file = plots_dir+dirs.prefix+'Wave_'+outname+'.pdf'
    FILE_DELETE,comb_file,/allow
    print,'Writing combined PDF plots file to ',comb_file
    cmd = 'gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile='+comb_file+' '
    cmd += strjoin(psfiles+'.pdf',' ')
    spawn,cmd,out,errout
  
  endif
  
  
  ; Save the structures
  savefile = wave_dir+dirs.prefix+'Wave-'+outname+'.dat'
  print,'Saving fitting information to ',savefile
  SAVE,linestr,mlinestr,flinestr,dispstr,foutstr,coefstr,wavestr,fil=savefile
  
  medrms = median(coefstr.rms)
  print,'Median RMS = ',strtrim(medrms,2)

  ;stop

End ; poly+sine fit method

Else: begin
  print,'Fitmethod = ',fitmethod,' Not supported'
  return
End

ENDCASE

  
;stop
  
if keyword_set(stp) then stop

end
  
