;+
;
; APWAVECAL
;
; Do the wavelength calibration on a ThAr frame.  Fit all
; three chips together.  Fit multiple nights worth of lines
; and allow the chips to move slightly with time and dither shift.
;
; INPUTS:
;  lampid       The directory and ID8 string for the ThAr 1D frame concatenated
;  lsfid        The direcotry and ID8 string for the LSF calibration file.
;  =psfid       The directory and ID8 string for the PSF calibration file.
;  /multi       Use input .sav files if they exist for a multi-night solution
;  /pl          Show plots.
;  /save        Save diagnostic plots.
;  /verbose     Verbose output to the screen.
; 
; OUTPUTS:
;  The wavelength solutions and wavelength arrays will be written to
;  output FITS files.
;
; USAGE:
;  IDL>apwavecal,lampid8,lsfid,psfid=psfid
;
; By D.Nidever  Feb. 2010
;-

pro apmultiwavecal,lampid,lsfid,coefstr,psfid=psfid,verbose=verbose,pl=pl,multi=multi,$
              save=save,stp=stp,name=name,medfilt=medfilt,clobber=clobber


; Error Handling
;------------------
; Establish error handler. When errors occur, the index of the  
; error is returned in the variable Error_status:  
;CATCH, Error_status 

;This statement begins the error handler:  
;if (Error_status ne 0) then begin 
;   error = !ERROR_STATE.MSG  
;   if not keyword_set(silent) then print,error
;   CATCH, /CANCEL 
;stop
;   return
;endif

; Get APOGEE directories
dirs = getdir(apogee_dir,cal_dir,spectro_dir,apogee_vers,lib_dir)
psf_dir = spectro_dir+'cal/psf/'
wave_dir = spectro_dir+'cal/wave/'
lsf_dir = spectro_dir+'cal/lsf/'
linelist_dir = lib_dir+'arclines/'

nlampid = n_elements(lampid)

; Not enough inputs
if nlampid eq 0 then begin
  print,'Syntax - apwavecal,lampid,lsfid,psfid=psfid,multi=multi'
  return
endif

chiptag = ['a','b','c']
npix = 2048L

; Defaults
if n_elements(pl) eq 0 then pl=0                ; no plotting by default
if n_elements(save) eq 0 then save=1            ; save by default

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
info = APFILEINFO(file_dirname(lampid)+'/'+dirs.prefix+'1D-a-'+lampframeid+'.fits',/silent)
mjd5 = strtrim(long(info[0].mjd5),2)

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
endif else if ~keyword_set(multi) then  begin
  lampfiles0 = apogee_filename('1D',chip=chiptag,num=lampframeid[0])
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
endif


; Load the lamp linelists from Save files
;----------------------------------------
info_schema = {frame:'',lamptype:'',mjd:0L,dateobs:'',jd:0.0d0,ditherpix:0.0,group:0L,lo:0LL,hi:0LL}
if keyword_set(multi) then begin
  apgundef,linestr_all,mlinestr_all,linestr,mlinestr,allinfo,alldispstr
  ;; Loop over the LAMPIDs
  for i=0,nlampid-1 do begin
    savefile = wave_dir+dirs.prefix+'Wave-'+lampframeid[i]+'.dat'
    if file_test(savefile) eq 1 then begin
      print,'Loading ',savefile
      restore,savefile
      push,linestr_all,linestr
      push,mlinestr_all,mlinestr
      ui = uniq(mlinestr.frame,sort(mlinestr.frame))
      nuframe = n_elements(ui)
      uframe = mlinestr[ui].frame
      info1 = replicate(info_schema,nuframe)
      info1.frame = uframe
      info1.lamptype = mlinestr[ui].lamptype
      cmdj = getcmjd(long(info1.frame),mjd=mjd)
      info1.mjd = mjd
      for j=0,nuframe-1 do begin
        lampfiles = apogee_filename('2D',chip=chiptag,num=uframe[j])
        if file_test(lampfiles[0]) eq 0 then lampfiles+='.fz'
        if file_test(lampfiles[0]) eq 1 then begin
          head0 = headfits(lampfiles[0],exten=0)
          info1[j].ditherpix = sxpar(head0,'dithpix')
          info1[j].dateobs = sxpar(head0,'date-obs')
          info1[j].jd = date2jd(info1[j].dateobs)
        endif else print,'CANOT get dither shift for '+uframe[j]
      endfor
      push,alldispstr,dispstr
      push,allinfo,info1
    endif
  endfor
  linestr = linestr_all
  mlinestr = mlinestr_all
  apgundef,mlinestr_all

endif else begin

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
  
  ; Loop through the lamp frames and fit the lines
  ;-----------------------------------------------
  apgundef,dispstr
  For i=0,n_elements(lampframeid)-1 do begin
    apgundef,lampfiles,info
    apgundef,linestr1,linestr2,linestr3
    apgundef,dispstr1,dispstr2,dispstr3
  
    ; Construct chip names
    lampfiles = apogee_filename('1D',chip=chiptag,num=lampframeid[i])
  
    ; Get file info
    info = APFILEINFO(lampfiles,/silent)
  
    ; Check the files
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
      bdmask = where(frame0.(k).mask and badmask(),nbdmask)
      frame0.(k).flux[bdmask] = 0.0
      frame0.(k).err[bdmask] = baderr()
  
      medfilt=5
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
      endfor
  
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
  
    nfibers1 = n_elements(dispstr1)
    nfibers2 = n_elements(dispstr2)
    nfibers3 = n_elements(dispstr3)
  
    if nfibers2 ne nfibers1 or nfibers3 ne nfibers1 then begin
      print,'Number of fibers in three chips are NOT the same'
      return
    endif
  
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
  
    ; Get Y-positions from the trace in the middle of the GREEN chip
    ;  ADD TO LINESTR
    ypos1 = fltarr(nfibers)
    for j=0,nfibers-1 do ypos1[j] = poly(0,tstr2[j].coef)
    ADD_TAG,frame_linestr,'FIBER_YPOS1',0.0,frame_linestr
    frame_linestr.fiber_ypos1 = ypos1[frame_linestr.fiber]
    ypos2 = fltarr(nfibers)
    for j=0,nfibers-1 do ypos2[j] = poly(2047,tstr2[j].coef)
    ADD_TAG,frame_linestr,'FIBER_YPOS2',0.0,frame_linestr
    frame_linestr.fiber_ypos2 = ypos2[frame_linestr.fiber]
    ymid = fltarr(nfibers)
    for j=0,nfibers-1 do ymid[j] = poly(1023.5,tstr2[j].coef)
    ADD_TAG,frame_linestr,'FIBER_YMID',0.0,frame_linestr
    frame_linestr.fiber_ymid = ymid[frame_linestr.fiber]
    ADD_TAG,frame_linestr,'YPOS',0.0,frame_linestr
    for j=0L,n_elements(frame_linestr)-1 do $
        frame_linestr[j].ypos = poly(frame_linestr[j].x,(tracestr.(frame_linestr[j].chipnum-1))[frame_linestr[j].fiber].coef)
    add_tag,frame_linestr,'XB',0.0,frame_linestr
  
    ;; Only want lines with matches
    indmatch = where(frame_linestr.model_match eq 1 and frame_linestr.model_usewave eq 1,nindmatch)
    frame_mlinestr = frame_linestr[indmatch]  ; lines with matches
  
    ; Combine with previous linelist
    PUSH,linestr,frame_linestr
    PUSH,mlinestr,frame_mlinestr
  
    ;; Get info for this frame
    head0 = headfits(lampfiles[0],exten=0)
    info1 = replicate(info_schema,1)
    info1.frame = lampframeid[i]
    info1.lamptype = lamptype
    cmdj = getcmjd(long(info1.frame),mjd=mjd)
    info1.mjd = mjd
    info1.ditherpix = sxpar(head0,'dithpix')
    info1.dateobs = spxar(head0,'date-obs')
    info1.jd = date2jd(info1.dateobs)
    push,allinfo,info1
  
  
    ; Add extra tags to the dispersion structures
    ADD_TAG,dispstr1,'CHIPNUM',1L,dispstr1
    ADD_TAG,dispstr2,'CHIPNUM',2L,dispstr2
    ADD_TAG,dispstr3,'CHIPNUM',3L,dispstr3
    frame_dispstr = [dispstr1,dispstr2,dispstr3]
    ADD_TAG,frame_dispstr,'LAMPTYPE',lamptype,frame_dispstr
    ADD_TAG,frame_dispstr,'FRAME',lampframeid[i],frame_dispstr
  
    ; Combine with previous dispersion structures
    PUSH,alldispstr,frame_dispstr
  
  Endfor
endelse

nallinfo = n_elements(allinfo)
dispstr = alldispstr
apgundef,alldispstr

;; Group the frames
;;  subsequent frames at the same dither shift can be in the same "group"
si = sort(allinfo.jd)
allinfo = allinfo[si]
si = sort(long(mlinestr.frame))
mlinestr = mlinestr[si]
add_tag,mlinestr,'group',0L,mlinestr
index = create_index(mlinestr.frame)
allinfo[0].group = 1
group = 1
for i=0,nallinfo-1 do begin
  ;; Not the same grouop
  if (i gt 0) and ((long(allinfo[i].frame) ne long(allinfo[i-1].frame)+1) or $
     abs(allinfo[i].ditherpix-allinfo[i-1].ditherpix) gt 0.01) then group++
  allinfo[i].group = group
  MATCH,index.value,allinfo[i].frame,ind1,ind2,/sort,count=nmatch
  allinfo[i].lo = index.lo[ind1[0]]
  allinfo[i].hi = index.hi[ind1[0]]
  mlinestr[index.lo[ind1[0]]:index.hi[ind1[0]]].group = group
endfor
ngroups = max(allinfo.group)
print,strtrim(ngroups,2),' groups'

fibindex = create_index(mlinestr.fiber)


;; Check if the output already exists
mjd5 = strtrim(long(allinfo[0].mjd),2)
outname = allinfo[0].frame
; Multi-night/group solution, use day number of 1st frame
if ngroups gt 1 then outname=strmid(allinfo[0].frame,0,4)+'0000'
if file_test(wave_dir+dirs.prefix+'Wave-'+outname+'.dat') eq 1 and not keyword_set(clobber) then begin
  print,wave_dir+dirs.prefix+'Wave-'+outname+'.dat already exists and /clobber not set'
  return
endif


;################################################
;# Polynomial fits to each Fiber separately
;################################################

;; Loop over the groups and fit wavelength solutions to them
flinestr = mlinestr
for i=0,ngroups-1 do begin
  print,''
  print,'FITTING WAVELENGTH SOLUTIONS TO GROUP=',strtrim(i+1,2)
  MATCH,allinfo.group,i+1,ind1,ind2,/sort,count=nmatch
  lo = min(allinfo[ind1].lo)
  hi = max(allinfo[ind1].hi)
  mlinestr1 = mlinestr[lo:hi]
  APWAVECAL_GROUP,mlinestr1,flinestr1,parstr1,coefstr1,npoly=npoly,nfibers=nfibers
  ; chipgap1_coef, chipgap2_coef
  flinestr[lo:hi] = flinestr1
  if i eq 0 then begin
     parstr_schema = parstr1[0]
     struct_assign,{dumdum:''},parstr_schema
     parstr = replicate(parstr_schema,nfibers,ngroups)
     coefstr_schema = coefstr1[0]
     struct_assign,{dumdum:''},coefstr_schema
     coefstr = replicate(coefstr_schema,nfibers,ngroups)
  endif
  parstr[*,i] = parstr1
  coefstr[*,i] = coefstr1
endfor
;; Remove second dimension for single group
if ngroups eq 1 then begin
  parstr = reform(parstr)
  coefstr = reform(coefstr)
endif

;; Now fit all lines simultaneously for a fiber with small shifts
;; per chip/group
if ngroups gt 1 then begin
  fparstr = replicate({fiber:0,pars:fltarr(2+ngroups*3+npoly),perror:dblarr(2+ngroups*3+npoly),$
                       polycoef:dblarr(npoly),xshifts:fltarr(3,ngroups),ypos:0.0,sig:0.0,rms:0.0},nfibers)
  for i=0,nfibers-1 do begin
    MATCH,fibindex.value,i,ind1,ind2,/sort,count=nmatch
    ind = fibindex.index[fibindex.lo[ind1[0]]:fibindex.hi[ind1[0]]]    
    nind = n_elements(ind)    
    flinestr1 = flinestr[ind] 

    ;; Initial coefficients
    coefstr1 = reform(coefstr[i,*])
    chipgap1 = median([coefstr1.chipgap1])
    chipgap2 = median([coefstr1.chipgap2])
    ;; Group dither shifts
    ditherpix = fltarr(ngroups)
    ditherpix[allinfo.group-1] = allinfo.ditherpix
    ditherpix -= ditherpix[0]

    ;; Parameters
    ;; require the first group green chip to be "anchored"
    ;; small x-shift for each chip and group, wavelength coefficients

    ;; Not all groups might be represented for each fiber

    if nind gt 0 then begin
      xx = flinestr1.x
      yy = flinestr1.model_wave
      err = flinestr1.gfit_perror[1] > 0.05  ;0.1
      ;; Use initial chipgap solutions and dither shifts to
      ;;  to get initial small offsets
      grpchipxoff = fltarr(ngroups*3)
      ;; Add dither shift offsets relative to 1st group
      ;;  only for green chip (offset for all chips of group)
      grpchipxoff[indgen(ngroups)*3+1] -= ditherpix
      ;; Add offsets to red chip using chipgap1 relative to median
      grpchipxoff[indgen(ngroups)*3] += -(coefstr1.chipgap1-chipgap1)
      ;; Add offsets to blue chip using chipgap2 relative to median
      grpchipxoff[indgen(ngroups)*3+2] += coefstr1.chipgap2-chipgap2
      initpars = [chipgap1, chipgap2, grpchipxoff, coefstr1[0].coef]
      fa = {group:flinestr1.group,chipnum:flinestr1.chipnum}
      parinfo = replicate({limited:[0,0],limits:[0.0,0.0],fixed:0},n_elements(initpars))
      parinfo[0:1].fixed = 1  ; chipgaps are fixed
      parinfo[2:ngroups*3+1].limited = 1
      for j=0,ngroups*3-1 do parinfo[j+2].limits = [-0.5,0.5]+grpchipxoff[j]
      parinfo[indgen(ngroups)*3+3].limits[0] = -15   ; large limits for pixel shifts
      parinfo[indgen(ngroups)*3+3].limits[1] = 15
      parinfo[3].fixed = 1  ; force group 1 green chip xoff=0, ANCHOR
      pars1 = MPFITFUN('func_multi_poly',xx,yy,err,initpars,status=status,dof=dof,$
                       parinfo=parinfo,functargs=fa,bestnorm=chisq,perror=perror,yfit=yfit,/quiet)
      yfit1 = func_multi_poly(xx,pars1,group=flinestr1.group,chipnum=flinestr1.chipnum,xb=xb)
  
      ;plotc,xb,yy-yfit1,flinestr1.group,ps=1 

      ; Remove outliers and refit
      diff = yy-yfit1
      sig = mad(diff)
      gd = where(abs(diff) lt 2.5*sig,nbd)

      ;; Add median red/blue xoff corrections to chipgaps
      initpars2 = pars1
      medxoff_red = median(pars1[indgen(ngroups)*3+2],/even)
      initpars2[0] -= medxoff_red
      initpars2[indgen(ngroups)*3+2] -= medxoff_red
      medxoff_blue = median(pars1[indgen(ngroups)*3+4],/even)
      initpars2[1] += medxoff_blue
      initpars2[indgen(ngroups)*3+4] -= medxoff_blue
     
      fa2 = {group:flinestr1[gd].group,chipnum:flinestr1[gd].chipnum}
      pars2 = MPFITFUN('func_multi_poly',xx[gd],yy[gd],err[gd],initpars2,status=status2,dof=dof2,$
                       parinfo=parinfo,functargs=fa2,bestnorm=chisq2,perror=perror2,yfit=yfit2,/quiet)
      yfit2 = func_multi_poly(xx,pars2,group=flinestr1.group,chipnum=flinestr1.chipnum,xb=xb)
      rchisq2 = chisq2/dof2
      diff2 = flinestr1.model_wave-yfit2
      sig2 = mad(diff2)

      ;; Add median red/blue xoff corrections to chipgaps
      fpars = pars2
      medxoff_red = median(fpars[indgen(ngroups)*3+2],/even)
      fpars[0] -= medxoff_red
      fpars[indgen(ngroups)*3+2] -= medxoff_red
      medxoff_blue = median(fpars[indgen(ngroups)*3+4],/even)
      fpars[1] += medxoff_blue
      fpars[indgen(ngroups)*3+4] -= medxoff_blue

      fparstr[i].fiber = i
      fparstr[i].pars = fpars
      fparstr[i].perror = perror2
      fparstr[i].polycoef = fpars[2+ngroups*3:*]
      fparstr[i].xshifts = reform(fpars[2:ngroups*3+1],3,ngroups)
      fparstr[i].ypos = median(flinestr1.ypos)
      fparstr[i].sig = mad(yy-yfit2)
      fparstr[i].rms = stddev(yy[gd]-yfit2[gd])
      flinestr[ind].wave_fit = yfit2
    endif else begin
    ;; no lines for this fiber
      print,'No lines for fiber=',strtrim(i+1,2)
      fparstr[i].fiber = i
      fparstr[i].pars = 999999.
      fparstr[i].perror = 999999.
      ;fparstr[i].ypos = median(mlinestr[ind].ypos)
      fparstr[i].sig = 999999.
    endelse
  endfor  ; fiber loop

;; Only ONE group
endif else begin
  fparstr = replicate({fiber:0,pars:fltarr(2+npoly),perror:dblarr(2+npoly),$
                       polycoef:dblarr(npoly),xshifts:fltarr(3),ypos:0.0,sig:0.0,rms:0.0},nfibers)
  struct_assign,coefstr,fparstr
  ;; The COEFSTR has the best parameteres
  fparstr.pars[0] = coefstr.chipgap1
  fparstr.pars[1] = coefstr.chipgap2
  fparstr.pars[2:*] = coefstr.coef
  fparstr.perror[2:*] = coefstr.perror
  fparstr.polycoef = coefstr.coef
endelse

print,''  
print,'Final poly fits to each fiber. Sig = ',stringize(mad(flinestr.model_wave-flinestr.wave_fit),ndec=4),' A'


; Get the wavelength coeffcients and arrays
;------------------------------------------
wavestr = REPLICATE({coef:dblarr(nfibers,6+npoly),wave:dblarr(npix,nfibers)},3)
x1 = findgen(npix)
for i=0,nfibers-1 do begin
  polypars = fparstr[i].polycoef
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
  xoffset = 0.0
  chipgap1 = fparstr[i].pars[0]
  chipgap2 = fparstr[i].pars[1]
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
  print,'Writing output files to = '+wave_dir+dirs.prefix+'Wave-[abc]-'+outname+'.fits'

lampfiles0 = apogee_filename('1D',chip=chiptag,num=lampframeid[0])

; Loop through the chips
for i=0,2 do begin  
  outfile = wave_dir+dirs.prefix+'Wave-'+chiptag[i]+'-'+outname+'.fits'
  head0 = headfits(lampfiles0[i])
  ; Put LAMP filenames into header
  for j=0,n_elements(lampframeid)-1 do $
    sxaddpar,head0,'LAMPFIL'+strtrim(j+1,2),lampframeid[j]
  
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
  sxaddhist,leadstr+' HDU4 - Structure of fitting results',head0
  sxaddhist,leadstr+'',head0
  
  ; Put LAMP filenames into header
  for j=0,n_elements(lampid)-1 do begin
    ;lampfil = file_dirname(lampid[j])+'/'+dirs.prefix+'1D-a-'+lampframeid[j]+'.fits'
    lampfil = apogee_filename('1D',chip='a',num=lampframeid[j])
    hd = headfits(lampfil)
    lamptype = ''
    if sxpar(hd,'LAMPUNE') eq 1 then lamptype='URANIUM'
    if sxpar(hd,'LAMPTHAR') eq 1 then lamptype='THARNE'
    lamptype = strupcase(lamptype)
    sxaddhist,leadstr+' LAMPFILE'+strtrim(j+1,2)+'='+lampframeid[j],head0
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
  ; HDU4 - fitting structure
  MWRFITS,fparstr,outfile,/silent
endfor ; chip loop
    
; Save the structures
savefile = wave_dir+dirs.prefix+'Wave-'+outname+'.dat'
print,'Saving fitting information to ',savefile
SAVE,allinfo,linestr,mlinestr,flinestr,parstr,coefstr,fparstr,dispstr,wavestr,fil=savefile
  
medrms = median(coefstr.rms)
print,'Median RMS = ',strtrim(medrms,2)
  
if keyword_set(stp) then stop

end
  
