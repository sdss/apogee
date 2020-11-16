pro apskysub,frame,plugmap,outframe,suboption=suboption,nearest=nearest,$
             error=error,silent=silent,verbose=verbose,pl=pl,stp=stp

;+
;
; APSKYSUB
;
; This subtracts the airglow lines for the spectra.
;
; INPUTS:
;  frame     A structure with the header/data information for an
;                undersampled frame that has been wavelength calibrated
;  plugmap   The Plug Map structure for this plate
;  =suboption   The sky subtraction option: (1) nearest fibers, or
;                (2) model line fitting.
;  /nearest  Use the nearest sky fibers to subtract the sky background
;  /silent   Don't print anything to the screen.
;  /verbose  Print lots of information to the screen
;  /pl       Make some plots
;  /stp      Stop at the end of the program.
;
; OUTPUTS:
;  outframe  The same frame but with airglow lines subtracted and
;              sky spectrum (and error) added to the arrays.
;  =error    The error message if one occurred.
;
; USAGE:
;  IDL>apskysub,frame,lsfid,plutmap,outframe
;
; By D. Nidever  April 2010
;-

apgundef,outframe
;setdisp,/silent

; Not enough inputs
if n_elements(frame) eq 0 or n_elements(plugmap) eq 0 then begin
  print,'Syntax - apskysub,frame,plugmap,outframe,suboption=suboption,nearest=nearest,error=error,silent=silent,verbose=verbose,pl=pl,stp=stp'
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
needtags1 = ['CHIPA','CHIPB','CHIPC','SHIFT']
for i=0,n_elements(needtags1)-1 do begin
  if (where(tags eq needtags1[i]))[0] eq -1 then begin
    print,'TAG ',needtags1[i],' NOT FOUND in input structure'
    return
  end
end
needtags2 = ['HEADER','FLUX','ERR','MASK','WAVELENGTH','LSFCOEF','WCOEF']
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


chiptag = ['a','b','c']

; Initialize outframe, add SKY, SKYERR
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

    ; Add SKY/SKYERR after WAVELENGTH
    if tags[j] eq 'WAVELENGTH' then begin
      chstr = CREATE_STRUCT(chstr,'SKY',fltarr(npix,nfibers))
      chstr = CREATE_STRUCT(chstr,'SKYERR',fltarr(npix,nfibers))
    end
  end

  ; Add to the final OUTFRAME
  if i eq 0 then begin
    outframe = CREATE_STRUCT('chip'+chiptag[i],chstr)
  endif else begin
    outframe = CREATE_STRUCT(outframe,'chip'+chiptag[i],chstr)
  endelse

end
outframe = CREATE_STRUCT(outframe,'shift',frame.shift)


; Load the AIRGLOW linelist
;---------------------------
airstr = IMPORTASCII(linelist_dir+'airglow.txt',/header,/silent)
;airstr = IMPORTASCII(linelist_dir+'airglow_test.txt',/header,/silent)
;print,'USING AIRGLOW_TEST.TXT FOR DEBUGGING!!!'
nairstr = n_elements(airstr)


; Now identify the sky fibers with the plug map structure
;  type: 0=target, 1=sky, 2=telluric
;skyind = where(plugmap.fiberdata.type eq 1,nsky)
nplug=n_elements(plugmap.fiberdata.target1)
sky=intarr(nplug)
sci=intarr(nplug)
for i=0,nplug-1 do begin
  sky[i]=issky(plugmap.fiberdata[i].target1,plugmap.fiberdata[i].target2) 
  sci[i]=isscience(plugmap.fiberdata[i].target1,plugmap.fiberdata[i].target2) 
endfor
skyplugind = where(plugmap.fiberdata.spectrographid eq 2 and $
                   plugmap.fiberdata.holetype eq 'OBJECT' and $
                   sky and not sci,nsky)
;                   plugmap.fiberdata.objtype eq 'SKY',nsky)

if nsky eq 0 then begin
  error = 'halt: NO SKY fibers.  CANNOT do sky subtraction.'
  if not keyword_set(silent) then print,error
  return
endif
skyfiberid = plugmap.fiberdata[skyplugind].fiberid
;; need to subtract 1 to get IDL base-0 indices
;skyindex = skyfiberid-1
; fiberid=1 is at the top of the detector or index=299
; index = 300-fiberid
skyindex = 300-skyfiberid

; Make SKYFRAME
For i=0,2 do begin
  chstr0 = frame.(i)
  tags = tag_names(chstr0)

  ; Make the new chip structure
  for j=0,n_elements(tags)-1 do begin
    arr = chstr0.(j)
    ; Data arrays. These are [NPix,Nfiber]
    dum = where(stregex(['FLUX','ERR','MASK','WAVELENGTH','SKY','SKYERR'],tags[j],/boolean) eq 1,ndata)
    if ndata gt 0 then arr=arr[*,skyindex]
    ; Calib arrays. These are [Nfiber,Npar]
    dum = where(stregex(['LSFCOEF','WCOEF'],tags[j],/boolean) eq 1,ncalib)
    if ncalib gt 0 then arr=arr[skyindex,*]

    if j eq 0 then begin
      chstr = CREATE_STRUCT(tags[j],arr)
    endif else begin
      chstr = CREATE_STRUCT(chstr,tags[j],arr)
    endelse
  end

  ; Add to the final OUTFRAME
  if i eq 0 then begin
    skyframe = CREATE_STRUCT('chip'+chiptag[i],chstr)
  endif else begin
    skyframe = CREATE_STRUCT(skyframe,'chip'+chiptag[i],chstr)
  endelse

end

; Remove the continuum from the spectra
;print,'Removing continuum from the sky spectra'
;skyframe0 = skyframe
;for i=0,2 do begin
;  flux = skyframe.(i).flux
;  cont100 = MEDFILT2D(flux,101,dim=1,/edge_copy)
;  sig = MAD(flux-cont100,/zero,dim=1)
;  sig2d = replicate(1,npix)#sig
;  bd = where( (flux-cont100)/(sig2d>0.5) gt 4,nbd)
;  temp = flux
;  if nbd gt 0 then temp[bd]=!values.f_nan
;
;  cont100_2 = MEDFILT2D(temp,101,dim=1,/edge_copy)
;
;  skyframe.(i).flux = flux-cont100_2
;
;end
;
x = dindgen(npix)

; Defaults
if n_elements(suboption) eq 0 and keyword_set(suboption) then suboption=1
if n_elements(suboption) eq 0 then suboption=1  ; nearest sky fibers by default


;===============================
; DIFFERENT SUBTRACTION METHODS
;===============================
CASE suboption of

;################################################
;# 1.) Nearest sky fibers
;################################################

1: Begin  


  ;stop
  ;iplug=intarr(300)
  ;for i=0,299 do iplug[i]=where(plugmap.fiberdata.spectrographid eq 2 and plugmap.fiberdata.fiberid eq 300-i)
  ;isky=where(plugmap.fiberdata[iplug].objtype eq 'SKY')


  if not keyword_set(silent) then print,'Using nearest sky fibers for sky subtraction'

  print,'Subtracting Sky Lines from Spectra'

  ; check for persistence in fibers in blue chip
  persist=intarr(nfibers)
  For i=0,nfibers-1 do begin
    mask = frame.(2).mask[*,i]
    ; do we have a lot of persistence pixels?
    junk=where((mask and maskval('PERSIST_LOW')) or $
                (mask and maskval('PERSIST_MED')) or (mask and maskval('PERSIST_HIGH')),nbad)
    if float(nbad)/n_elements(mask) gt 0.5 then persist[i]=1 
  endfor


  ; Loop through the fibers
  For i=0,nfibers-1 do begin

    if (i+1) mod 50 eq 0 then print,strtrim(i+1,2),'/',strtrim(nfibers,2)

    ; The plugmap index for this fiber
    ; fiberid=1 is at the top of the detector or index=299
    ; index = 300-fiberid
    iplugind = where(plugmap.fiberdata.spectrographid eq 2 and $
                     plugmap.fiberdata.fiberid eq 300-i,niplugind)
    telstarplugind = where(plugmap.fiberdata.spectrographid eq 2 and $
                     plugmap.fiberdata.holetype eq 'OBJECT' and $
                     plugmap.fiberdata.objtype eq 'HOT_STD',nstar)


    ; No information for this fiber
    if niplugind eq 0 then begin
      print,'No information for Fiber=',strtrim(i+1,2),' in the plugmap file'
      goto,BOMB0
    endif
    objtype  = plugmap.fiberdata[iplugind].objtype
    ptarg  = plugmap.fiberdata[iplugind].target1
    starg  = plugmap.fiberdata[iplugind].target2

    ; Only correct non-sky fibers
;    If objtype ne 'SKY' then begin

      ; Loop through the chips
      For j=0,2 do begin

        ; Get object fiber data
        fiber = frame.(j).flux[*,i]
        fibererr = frame.(j).err[*,i]
        mask = frame.(j).mask[*,i]
        wave = frame.(j).wavelength[*,i]

        ; Find the closest sky fiber
        ira = plugmap.fiberdata[iplugind].ra
        idec = plugmap.fiberdata[iplugind].dec

        ;pskyplugind=skyplugind
        ; only take fibers that are "nearby" on chip, for LSF matching
        sgood=where(abs(plugmap.fiberdata[skyplugind].fiberid - (300-i)) lt 75,nsgood)
        if nsgood le 0 then begin
          print,'not halted: no sky fibers found within 75 fibers, fiber: ', 300-i
          goto,BOMB0
        endif

        ;oldpskyplugind=skyplugind[sgood]
        ; only include sky fibers for which adjacent fiber has H>9.5, but require 2, so adjust limit as needed
        pskyplugind=[]
        skymax=9.5
        while n_elements(pskyplugind) le 1 and skymax gt 5 do begin
         for isky=0,nsgood-1 do begin
          ii=plugmap.fiberdata[skyplugind[sgood[isky]]].fiberid 
          iplus = where(plugmap.fiberdata.spectrographid eq 2 and $
                        plugmap.fiberdata.fiberid  eq ii+1, nadj)
          if nadj eq 0 or plugmap.fiberdata[iplus].mag[1] gt skymax then  begin
            iplus = where(plugmap.fiberdata.spectrographid eq 2 and $
                        plugmap.fiberdata.fiberid  eq ii-1, nadj)
            if nadj eq 0 or plugmap.fiberdata[iplus].mag[1] gt skymax then pskyplugind=[pskyplugind,skyplugind[sgood[isky]]]
          endif
         endfor
         skymax-=0.25
        endwhile
        if skymax le 5 then begin
          print,'not halted: problem with enough sky fibers, fiber: ', 300-i
          goto,BOMB0
        endif

        ;if blue chip, consider persisnce
        if j eq 2 then begin
          ; only use "matching" persistence sky fibers
          ;sgood=where(persist[skyindex] eq persist[i],ngood)
          ; only use non-persistence sky fibers
          sgood=where(persist[skyindex] eq 0,ngood)
          if ngood gt 0 then pskyplugind=skyplugind[sgood] 
        endif 
        npsky=n_elements(pskyplugind)
        
        dist = sphdist(plugmap.fiberdata[pskyplugind].ra,plugmap.fiberdata[pskyplugind].dec,ira,idec,/deg)

        si = sort(dist)
        bestind = si[0:4<(npsky-1)]

        ; if we are subtracting from sky fiber, don't use it to determine the sky!
        if issky(ptarg,starg) and not isscience(ptarg,starg) then bestind = si[1:5<(npsky-1)]
        ;if objtype eq 'SKY' then bestind = si[1:5<(npsky-1)]
        ;bestind = first_el(minloc(dist))
        best_skyplugind = pskyplugind[bestind]
        best_skyplugindex = 300-plugmap.fiberdata[best_skyplugind].fiberid
        nusesky = n_elements(best_skyplugind)
        skydist = dist[bestind]

        ; print indices of selected sky fibers
        ;if j eq 0 then print,'fiber index: ', i, 'sky fibers: ', best_skyplugindex,skydist

        ; Get object fiber data
        fiber = frame.(j).flux[*,i]
        fibererr = frame.(j).err[*,i]
        wave = frame.(j).wavelength[*,i]
        ; Get fiber spectrum with no continuum
        fiberlines = fiber - MEDFILT1D(fiber,150,/edge)
        bd = where((mask and badmask()) gt 0,nbd)
        if nbd gt 0 then fiberlines[bd]=0.

        ; Initialize sky fiber arrays
        skyfiber = fltarr(npix,nusesky)
        skyfibererr = fltarr(npix,nusesky)
        skywave = dblarr(npix,nusesky)
        xsharr = fltarr(nusesky)

        totskyspec = fiber*0.0
        totskyerr = fiber*0.0
        ;totwt = 0.0
        totwt = fiber*0.0

        ;i Loop through neighboring sky fibers
        For k=0,nusesky-1 do begin
          ; Get sky fiber data
          skyfiber1 = frame.(j).flux[*,best_skyplugindex[k]]
          skyfibererr1 = frame.(j).err[*,best_skyplugindex[k]]
          skywave1 = frame.(j).wavelength[*,best_skyplugindex[k]]
          skymask = frame.(j).mask[*,best_skyplugindex[k]] and badmask()
          bd=where(skymask gt 0,nbd)
          skymask1 = fiber*0.0
          if nbd gt 0 then skymask1[bd]=1.

          ;Sort by wavelength
          wsi = sort(skywave1)
          skyfiber2 = skyfiber1[wsi]
          skywave2 = skywave1[wsi]
          skyfibererr2 = skyfibererr1[wsi]
          skymask2=skymask1[wsi]
          
          ; Get pixels with sky fiber wavelength range and sort
          wind = where(wave ge min(skywave2) and wave le max(skywave2),nwind)
          wsi2 = sort(wave[wind])
          wind2 = wind[wsi2]
          
          ; Need to interpolate the sky fiber to the object wavelength scale
          skyfiber3 = fiber*0.0
          skyfibererr3 = fiber*0.0 + median(skyfibererr2)
          ; initialize with all pixels bad in mask
          skymask3 = skymask2 * 0. + 1.
          skyfiber3[wind2] = SPLINE(skywave2,skyfiber2,wave[wind2],/double)
          ; Fixed following line 2/13/12, Holtz.
          ;skyfibererr3[wind2] = SPLINE(skyfibererr2,skyfiber2,wave[wind2],/double)
          skyfibererr3[wind2] = SPLINE(skywave2,skyfibererr2,wave[wind2],/double)
          skymask3[wind2] = SPLINE(skywave2,skymask2,wave[wind2],/double)

          ; Check for shifts: maybe not, bad pixels can throw things off!
          ;XCORLB,fiberlines,skyfiber3,5,xsh
          xsh=0
          xsharr[k] = xsh

          ; Shift the sky spectrum
          if abs(xsharr[k]) gt 0.01 then begin
            skyfiber4 = fiber*0.0
            skyfibererr4 = fiber*0.0 + median(skyfibererr1)
            skymask4 = skymask3 * 0. + 1.
            x2 = x-xsharr[k]
            ind = where(x2 ge 0.0 and x2 le npix-1,nind)
            skyfiber4[ind] = SPLINE(x,skyfiber3,x2[ind],/double)
            skyfibererr4[ind] = SPLINE(x,skyfibererr3,x2[ind],/double)
            skymask4[ind] = SPLINE(x,skymask3,x2[ind],/double)
          endif else begin
            skyfiber4 = skyfiber3
            skyfibererr4 = skyfibererr3
            skymask4 = skymask3
          endelse

          ; instead of scaling each individual sky frame, scale the combined one?
          ; or maybe avoid scaling altogether?
          doscale=0
          if doscale then begin
            ; Get the sky continuum
            skyfiber4_cont1 = MEDFILT1D(skyfiber4,251,/edge)
            skysig = MAD(skyfiber4-skyfiber4_cont1)
            bdpix = where(skyfiber4-skyfiber4_cont1 gt 3*skysig,nbdpix)
            temp = skyfiber4
            if nbdpix gt 0 then temp[bdpix]=!values.f_nan
            skyfiber4_cont = MEDFILT1D(temp,201,/edge)
            skyfiber4_lines = skyfiber4 - skyfiber4_cont
  
            ; Scale the sky fiber using one line

            ; Get the lines to use for scaling
            ;case j of
            ;0: refwave=16502.342
            ;1: refwave=16502.342
            ;2: refwave=16502.342
            ;endcase

            ;refpix = where(abs(wave-refwave) lt 1.5 and skyfiber4_lines gt 5*skysig and $
            ;               fiberlines gt 5*fibererr,nrefpix)
            refpix = where(skyfiber4_lines gt 10*skysig,nrefpix)
            if nrefpix eq 0 then refpix = where(skyfiber4_lines gt 5*skysig,nrefpix)
            if nrefpix eq 0 then refpix = where(skyfiber4_lines gt 2*skysig,nrefpix)
            if nrefpix eq 0 then goto,BOMB1
  
            scal = MEDIAN(fiberlines[refpix]/skyfiber4_lines[refpix])
            ;scal = 1.0
            ;print,scal
            ; corrected following line Holtz 110922
            ;fskyfiber = scal*skyfiber4 + skyfiber4_cont
            fskyfiber = scal*skyfiber4_lines + skyfiber4_cont
            fskyfibererr = scal*skyfibererr4         ; this is probably not quite right!!!

          endif else begin
            fskyfiber = skyfiber4
            fskyfibererr = skyfibererr4         ; this is probably not quite right!!!
          endelse


          ; The scaling doesn't seem to make a huge difference.

          ; save sky spectra, but these are not used?
          skyfiber[*,k] = fskyfiber
          skyfibererr[*,k] = fskyfibererr

          ; raw interpolated sky fibers
          ;skyfiber[*,k] = skyfiber4
          ;skyfibererr[*,k] = skyfibererr4

          ; Get weighted average
          ;totskyspec += fskyfiber * (1.0/skydist[k])
          ;totskyerr += fskyfibererr * (1.0/skydist[k])
          ;totskyspec += fskyfiber * (1.0/skydist[k]^2)
          ;totskyerr += fskyfibererr * (1.0/skydist[k]^2)
          ;totwt += 1.0/skydist[k]^2

          ; only include "good" pixels, i.e. where there is corresponding sky
          ;gd=where((skymask and badmask()) eq 0,ngd)
          gd=where(skyfiber4 ne 0. and abs(skymask4) le 0.01,ngd)
          if ngd gt 0 then begin
            ;totskyspec[gd] += fskyfiber[gd] * (1.0/skydist[k]^2)
            ;totskyerr[gd] += fskyfibererr[gd] * (1.0/skydist[k]^2)
            ;totwt[gd] += 1.0/skydist[k]^2
            totskyspec[gd] += fskyfiber[gd] * (1.0/fskyfibererr[gd]^2)
            totwt[gd] += 1.0/fskyfibererr[gd]^2
          endif
          BOMB1:
        Endfor ; sky fiber loop

        ; Subtract sky spectrum from object fiber
        ; Use a weighted average based on the distance
        ;skyspec = totskyspec/totwt
        ;skyerr = totskyerr/totwt
        ;fibersub = fiber - skyspec
        ; use error weighting
        skyspec = totskyspec/totwt
        skyerr = sqrt(1./totwt)
        ; use sttdev of sky fibers for sky error ??? definitely lowers S/N!!
        ; if nusesky gt 1 then skyerr=stddev(skyfiber,dim=2) 

        skycont = MEDFILT1D(skyspec,251,/edge)
        skylines = skyspec-skycont
        skysig = MAD(skylines)
        refpix = where(skylines gt 10*skysig,nrefpix)
        if nrefpix eq 0 then refpix = where(skylines gt 5*skysig,nrefpix)
        if nrefpix eq 0 then refpix = where(skylines gt 2*skysig,nrefpix)
        if nrefpix eq 0 then goto,BOMB2
        scalpix=fiberlines[refpix]/skylines[refpix]
        gd=where(scalpix gt 0.5 and scalpix lt 1.5,ngd)
        if ngd gt 0 then scal = MEDIAN(scalpix[gd]) else scal=MEDIAN(scalpix)

        skyspec = skycont + scal*skylines
        fibersub = fiber - skyspec
        fibererr = sqrt(fibererr^2 + skyerr^2)

;jj=where(plugmap.fiberdata[telstarplugind].fiberid eq 300-i,nj)
;if j eq 2 and nj gt 0 then begin
;print,'scale: ', scal
;  @s
;  stop
;endif


;set_plot,'X'
;plot,skyspec,yrange=[0,500]
;oplot,skyspec,color=2
;oplot,skyerr,color=3
;oplot,skyerr,color=4
;if j eq 2 then begin
;plot,fiber,xr=[0,500]
;oplot,fiber-skyspec,color=2
;if j eq 2 and i eq 300-4 then stop
;endif

        ; Plotting
        ;pl = 1 ;0
        if keyword_set(pl) then begin
          yr = [-1000,median(fiber)*4]
          yr = [yr[0],yr[1]+range(yr)*0.2]
          plot,fiber,yr=yr,xs=1,ys=1,tit='Fiber '+strtrim(i+1,2)+' Chip '+chiptag[j]
          ;oplot,fiber,co=250
          oplot,fibersub,co=200
          ;oplot,skyspec,co=150
          ;legend,['Before Subtraction','After Subtraction','Sky Lines'],textcolor=[250,200,150],charsize=1.2,$
          ;       pos=[30,yr[1]-range(yr)*0.02]
          legend,['Before Subtraction','After Subtraction'],textcolor=[255,250],charsize=1.2,$
                 pos=[30,yr[1]-range(yr)*0.02]
          ;stop
        endif

        ; Save the results
        outframe.(j).flux[*,i] = fibersub
        outframe.(j).err[*,i] = fibererr
        outframe.(j).sky[*,i] = skyspec
        outframe.(j).skyerr[*,i] = skyerr

        bd = where(totwt le 0,nbd)
        if nbd gt 0 then begin
          outframe.(j).flux[bd,i] = 0.
          outframe.(j).sky[bd,i] = 0.
          outframe.(j).err[bd,i] = baderr()
          outframe.(j).skyerr[bd,i] = baderr()
          outframe.(j).mask[bd,i] = outframe.(j).mask[bd,i] or maskval('NOSKY')
        endif

        ; find and set bit for locations where (smoothed) sky line is significant contributor
        ;  to spectrum
        if not issky(ptarg,starg) then begin
          highsky=skylinemask(outframe.(j).flux[*,i],outframe.(j).sky[*,i],nhigh,width=5,thresh=2)
          if nhigh gt 0 then outframe.(j).mask[highsky,i] = $
            outframe.(j).mask[highsky,i] or maskval('SIG_SKYLINE')
        endif
        ;stop

        BOMB2:
      Endfor ; chip loop

      ;stop

;    endif ; not SKY type

    BOMB0:

    ;stop

  Endfor  ; fiber loop
End  ; nearest sky fibers




;##############################################################
;# 2.) Model line fitting
;##############################################################


2: Begin

  chiptag = ['a','b','c']

  ;print,'RESTORING apskysub_03700023.dat'
  ;restore,'apskysub_03700023.dat'
  ;goto,skiptohere

  ;-------------------------------------
  ; Step 1. Fit sky fiber line fluxes 
  ;-------------------------------------
  print,'Fitting LSF to the lines'
  APSKYSUB_SKYLINEFIT,frame,skyindex,airstr,linestr,/verbose

  ;gdlines = where(linestr.model_match eq 1 and linestr.lsffit_status gt 0,ngdlines)
  ;linestr0 = linestr
  ;linestr = linestr[gdlines]
  
  ADD_TAG,linestr,'NUMFIBERS',0L,linestr
  ADD_TAG,linestr,'MEDFLUX',0.0,linestr
  ADD_TAG,linestr,'FLUXFRAC',0.0,linestr
  ADD_TAG,linestr,'FLUXFRACERR',0.0,linestr
  ADD_TAG,linestr,'FIBERSPECINORM',1.0,linestr
  
  skiptohere:
  
  ;stop
  
  ; Initialize the "Model" line structure
  ;---------------------------------------
  ; This keeps track of all the important fitting information for
  ; each unique sky line
  modlinestr = REPLICATE({ID:0L,model_wave:0.0d0,type:'',species:0,doublet:0,dbl_wsep:0.0,$
                          numfibers:0L,use:0,wave:0.0d0,wavesig:0.0d0,$
                          medflux:0.0d0,sigflux:0.0,coef:dblarr(11)},n_elements(airstr))
  modlinestr.id = airstr.id
  modlinestr.model_wave = airstr.wave
  modlinestr.type = airstr.name ;airstr.type
  modlinestr.doublet = airstr.doublet
  modlinestr.dbl_wsep = airstr.dbl_wsep
  species = ['OH','O2']
  specind = where(modlinestr.type eq 'OH',nspecind)
  if nspecind gt 0 then modlinestr[specind].species = 1
  specind = where(modlinestr.type eq 'O2',nspecind)
  if nspecind gt 0 then modlinestr[specind].species = 2
  
  numfibers_thresh = nsky/2  ; minimum number of fibers the lines was detected in
                             ;  to make it okay to use
  species = ['OH','O2']
  nspecies = n_elements(species)
  
  
  ;-----------------------------------------------------------------
  ; Step 2. Iteratively Measure the median flux value for each line
  ;           and the normalization for each fiber/species
  ;-----------------------------------------------------------------
  
  ; This structure is used to calculate the MEDIAN FLUXFRAC for each
  ; FIBER and SPECIES
  fiberstr = REPLICATE({skyfiber:0,zeta:0.0d0,eta:0.0d0,fiberspecinorm:dblarr(nspecies)+1.0,$
                        sigfiberspecinorm:dblarr(nspecies)},nsky)
  fiberstr.skyfiber = lindgen(nsky)
  fiberstr.zeta = plugmap.fiberdata[skyplugind].zeta
  fiberstr.eta = plugmap.fiberdata[skyplugind].eta
  
  
  ; Iterate
  count = 0
  endflag = 0
  lastmedflux = modlinestr.medflux
  lastnorm = fiberstr.fiberspecinorm
  WHILE (endflag eq 0) do begin
  
    ; Getting median flux values for each line
    ;-----------------------------------------
  
    ; Get median flux for each unique line
    ui = uniq(linestr.model_id,sort(linestr.model_id))
    modelid = linestr[ui].model_id
    nmodelid = n_elements(ui)
  
    ; Loop through the unique measured lines
    for i=0,nmodelid-1 do begin
      gdlines = where(linestr.model_id eq modelid[i],ngdlines)
      gdmod = where(modlinestr.id eq modelid[i],ngdmod)
      modlinestr[gdmod].numfibers = ngdlines
      if ngdlines gt 1 then begin
        ; Divide by the Fiber/Species normalization for this line
        ;  to get a better measurement of the median flux
        ;medflux = MEDIAN(linestr[gdlines].lsffit_flux / linestr[gdlines].fiberspecinorm,/even)
   
        ; First estimate
        ratioflux = linestr[gdlines].lsffit_flux / linestr[gdlines].fiberspecinorm
        ratioflux_error = (linestr[gdlines].lsffit_perror[0] > 1) / linestr[gdlines].fiberspecinorm
        medflux0 = MEDIAN(ratioflux,/even)
        sigflux0 = MAD(ratioflux) > 1
        stdflux = STDDEV(ratioflux)
  
        ; Get good values and weight by errors
        gpts = where(abs(ratioflux-medflux0) lt 3*sigflux0,ngpts)
        ROBUST_MEAN,ratioflux[gpts],medflux,sigflux,sig=ratioflux_error[gpts],error=robust_error
        ; there was an error
        if n_elements(robust_error) then begin
          medflux = medflux0
          sigflux = sigflux0
        endif
  
        ; Plug the values into the structures
        modlinestr[gdmod].medflux = medflux
        modlinestr[gdmod].sigflux = sigflux
        ;modlinestr[gdmod].wave = MEDIAN(linestr[gdlines].wave,/even)
        ;modlinestr[gdmod].wavesig = MAD(linestr[gdlines].wave)
        modlinestr[gdmod].wave = MEDIAN(linestr[gdlines].lsffit_wave,/even)
        modlinestr[gdmod].wavesig = MAD(linestr[gdlines].lsffit_wave)
        linestr[gdlines].numfibers = ngdlines
        linestr[gdlines].medflux = medflux
        linestr[gdlines].fluxfrac = linestr[gdlines].lsffit_flux/medflux
        linestr[gdlines].fluxfracerr = linestr[gdlines].lsffit_perror[0]/medflux
       ; print,i
       ; print,'WAVE:',median(linestr[gdlines].wave-linestr[gdlines].model_wave),stddev(linestr[gdlines].wave)  ; debugging
       ; print,'WAVE:',median(linestr[gdlines].lsffit_wave-linestr[gdlines].model_wave),stddev(linestr[gdlines].lsffit_wave)  ; debugging
       ; print,'FLUX: ',medflux,sigflux,stdflux
        ;if medflux lt 5000 then stop
        ;stop
      endif
    endfor

  
    ; Lines to subtract
    uselines = where(modlinestr.numfibers gt numfibers_thresh,nuselines)
    modlinestr[uselines].use = 1
  
    ; Calculate an average FLUXFRAC for each species/fiber
    ;------------------------------------------------------
  
    ; Species loop
    For i=0,nspecies-1 do begin
      ; Loop through the sky fibers
      For j=0,nsky-1 do begin
  
        ; Exclude "superpersistence" lines
        if skyindex[j] gt 200 then begin
          gd = where(linestr.skyfiber eq j and linestr.model_type eq species[i] and $
                     linestr.numfibers gt numfibers_thresh and linestr.chip lt 2,ngd)
        endif else begin
          gd = where(linestr.skyfiber eq j and linestr.model_type eq species[i] and $
                     linestr.numfibers gt numfibers_thresh,ngd)
        endelse

        ; Calculate the median fluxfrac for all the lines of the same species
        ;  This gives a fiber/species normalization
  
        ; We have this species in this fiber
        if ngd gt 0 then begin
          if ngd gt 1 then begin
  
            fluxfrac = linestr[gd].fluxfrac
            fluxfracerr = linestr[gd].fluxfracerr
            fiberspecinorm0 = MEDIAN(fluxfrac,/even)
            sigfiberspecinorm0 = MAD(fluxfrac)
            stdfiberspecinorm = STDDEV(fluxfrac)
  
            ; Get good values and weight by errors
            gpts = where(abs(fluxfrac-fiberspecinorm0) lt 3*sigfiberspecinorm0,ngpts)
            ROBUST_MEAN,fluxfrac[gpts],fiberspecinorm,sigfiberspecinorm,sig=fluxfracerr[gpts],error=robust_error
            ; there was an error
            if n_elements(robust_error) then begin
              fiberspecinorm = fiberspecinorm0
              sigfiberspecinorm = sigfiberspecinorm0
            endif
  
            ;plot,fluxfrac,ps=8,yr=[0,3]
            ;stop

            ;print,i,j,fiberspecinorm,sigfiberspecinorm,stdfiberspecinorm
          endif else begin
            fiberspecinorm = linestr[gd[0]].fluxfrac
            sigfiberspecinorm = 0.0
          endelse
  
          ; Put it in the structures
          fiberstr[j].fiberspecinorm[i] = fiberspecinorm
          fiberstr[j].sigfiberspecinorm[i] = sigfiberspecinorm
          linestr[gd].fiberspecinorm = fiberspecinorm
        end
  
        ;stop
  
      Endfor ; fiber loop
    Endfor ; species looop
  
    ; Should we stop?
    tocheck = where(modlinestr.use eq 1,ntocheck)
    medchange = abs(modlinestr[tocheck].medflux-lastmedflux[tocheck])/abs(modlinestr[tocheck].medflux)*100
    normchange = abs(fiberstr.fiberspecinorm-lastnorm)/abs(fiberstr.fiberspecinorm)*100
    maxchange = max([medchange,(normchange)(*)])
    if maxchange lt 0.01 or count ge 10 then endflag=1
    lastmedflux = modlinestr.medflux
    lastnorm = fiberstr.fiberspecinorm
  
    ;print,maxchange
  
    count++
  
    ;stop
  
  ENDWHILE

  ;gd = where(linestr.fluxfrac lt 100 and finite(linestr.fluxfrac) eq 1)
  ;linestr2 = linestr[gd]
  ;plotc,linestr2.model_wave,linestr2.fiber,linestr2.fluxfrac/linestr2.fiberspecinorm,ps=1,min=0.8,max=1.2
  ;
  ;stop
  
  headstr = 'APSKYSUB: '


  ;---------------------------------------------------------
  ; Step 3. Fit the species normalization across the plate
  ;---------------------------------------------------------
  npars = 11
  
  species = ['OH','O2']
  nspecies = n_elements(species)
  
  ; This structure is to keep track of the spatial fits to MEDFLUXFRAC
  ; for each species
  speciesfitstr = REPLICATE({npts:0L,pars:dblarr(npars),perror:dblarr(npars),sig:0.0,$
                             chisq:0.0,rchisq:0.0,dof:0.0,status:0},nspecies)
  
  For i=0,nspecies-1 do begin
  
    ; Put the coefficients into MODLINESTR
    gdmod = where(modlinestr.type eq species[i],ngdmod)

    ; No lines for this species
    if ngdmod eq 0 then begin
      goto,BOMB3
    endif

    ; Now fit species normaliation as a function of ZETA/ETA
    xx = fiberstr.zeta
    yy = fiberstr.eta
    zz = fiberstr.fiberspecinorm[i]
    err = zz*0+1
    npts = n_elements(xx)
  
    apgundef,status,dof,chisq,rchisq,yfit
    initpars = dblarr(11)
    initpars[0] = 1.0
    parinfo = REPLICATE({limited:[0,0],limits:[0.0,0.0],fixed:0},11)
    ;  we can change the order by fixing certain parameters
    ; 0-constant, 1-x, 2-x^2, 3-x^3
    ; 4-x*y, 5-x^2*y, 6-x*y^2, 7-x^2*y^2
    ; 8-y, 9-y^2, 10-y^3
    pars = MPFIT2DFUN('func_poly2d',xx,yy,zz,err,initpars,status=status,dof=dof,$
                    bestnorm=chisq,parinfo=parinfo,perror=perror,yfit=yfit,/quiet)      
    if status lt 1 then begin
      print,'Error in the fitting.  Using median value'
      pars = [median([zz]), 0.0, 0.0]
      perror = -1
      chisq = -1
      dof = 0
      rchisq = -1
    endif else begin
      rchisq = chisq/dof
    endelse

    ; Remove outliers and refit
    diff = zz-yfit
    sig = mad(diff)
    bd = where(abs(diff) gt 2.5*sig,nbd)
    if nbd gt 0 and status gt 0 then begin
      xx_orig = xx
      yy_orig = yy
      zz_orig = zz
      pars1 = pars
      remove,bd,xx,yy,zz
      initpars = pars1

      pars = MPFIT2DFUN('func_poly2d',xx,yy,zz,err,initpars,status=status,dof=dof,$
                      bestnorm=chisq,parinfo=parinfo,perror=perror,yfit=yfit,/quiet)      
      if status lt 1 then begin
        print,'Error in the fitting.  Using median value'
        pars = [median([zz]), 0.0, 0.0]
        perror = -1
        chisq = -1
        dof = 0
        rchisq = -1
      endif else begin
        rchisq = chisq/dof
      endelse

    endif ; refit

    ; Measure the scatter and get parameter errors
    if status gt 0 then begin
      sig = MAD(zz-yfit,/zero)
      err = zz*0+sig
 
      pars2 = MPFIT2DFUN('func_poly2d',xx,yy,zz,err,pars,status=status2,dof=dof2,$
                      bestnorm=chisq2,parinfo=parinfo,perror=perror2,yfit=yfit2,/quiet)
      pcerror = perror2*sqrt(chisq2/dof2)

    endif else begin
      sig = MAD(zz-median([zz]))
      pcerror = perror*0
      pcerror[0] = sig
    endelse

    ; KLUDGE, no fitting for now
    ;print,'NO spatial fitting for now'
    ;pars = [median(zz),fltarr(10)]

    ; Save the values
    speciesfitstr[i].npts = npts
    speciesfitstr[i].pars = pars
    speciesfitstr[i].perror = pcerror ;perror
    speciesfitstr[i].sig = sig
    speciesfitstr[i].chisq = chisq
    speciesfitstr[i].dof = dof
    speciesfitstr[i].rchisq = rchisq
    speciesfitstr[i].status = status

    ; Put the coefficients in the header
    APADDPAR,outframe,headstr+'AIRGLOW SPECIES 2D SPATIAL POLYNOMIAL FIT',/history
    APADDPAR,outframe,headstr+'SPECIES '+species[i],/history
    APADDPAR,outframe,headstr+'NPARS = '+strtrim(n_elements(speciesfitstr[i].pars),2),/history
    for k=0,n_elements(speciesfitstr[i].pars)-1 do begin
      parname = 'SKPR'+strtrim(i+1,2)+'_'+strtrim(k+1,2)
      APADDPAR,outframe,parname,speciesfitstr[i].pars[k]
      errname = 'SKER'+strtrim(i+1,2)+'_'+strtrim(k+1,2)
      APADDPAR,outframe,errname,speciesfitstr[i].perror[k]
    endfor
    APADDPAR,outframe,headstr+'SIG = '+stringize(speciesfitstr[i].sig,ndec=5),/history


    ; Put the coefficients into MODLINESTR
    gdmod = where(modlinestr.type eq species[i],ngdmod)
    modlinestr[gdmod].coef = pars
  
    print,'Variation across the plate = ',strtrim(100*mad(zz),2),' percent'
  
    ; Plotting
    ;pl = 1
    if keyword_set(pl) then begin
      erase
      ; Restore color table
      if file_test(spectro_dir+'lib/colors/coltable1.fits') eq 1 then begin
        fits_read,spectro_dir+'lib/colors/coltable1.fits',rgb,/no_abort
        tvlct,rgb[0,*],rgb[1,*],rgb[2,*]
      endif
      ;psfile = 'skysub2'
      ;ps_open,psfile,/color,thick=3,/encap
      ;device,/inches,xsize=11,ysize=11
      if !d.name eq 'X' then co1=255 else co1=0
  
      yr = minmax(zz)
      yr = [yr[0]-range(yr)*0.1,yr[1]+range(yr)*0.1]
  
      ; Med Fluxfrac vs. Zeta
      xr = minmax(xx)
      xr = [xr[0]-range(xr)*0.1,xr[1]+range(xr)*0.1]
      pos = [0.08,0.58,0.50,0.95]
      plot,xx,zz,ps=1,xtit='Zeta (deg)',ytit='Median Flux Fraction',tit='Species '+species[i],$
           xr=xr,yr=yr,xs=1,ys=1,position=pos
      oplot,xx,yfit,ps=4,co=250
      legend,['Data','Model'],textcolor=[co1,250],/top,/left
  
      ; Med Fluxfrac vs. eta
      xr = minmax(yy)
      xr = [xr[0]-range(xr)*0.1,xr[1]+range(xr)*0.1]
      pos = [0.58,0.58,0.99,0.95]
      plot,yy,zz,ps=1,xtit='Eta (deg)',ytit='Median Flux Fraction',tit='Species '+species[i],$
           xr=xr,yr=yr,xs=1,ys=1,position=pos,/noerase
      oplot,yy,yfit,ps=4,co=250
      legend,['Data','Model'],textcolor=[co1,250],/top,/left
  
      ; Colored points on the sky
      pos = [0.08,0.08,0.50,0.42]
      colpos = [0.08,0.49,0.50,0.51]
      arr1 = cgscalevector(findgen(400),-1.5,1.5)
      xarr = arr1#replicate(1.0,400)
      yarr = replicate(1.0,400)#arr1
      zarr = FUNC_POLY2D(xarr,yarr,pars)
      bd = where(sqrt(xarr^2+yarr^2) gt 1.5,nbd)
      zarr[bd] = 1000
      ;zmin = min([(zarr)(*),zz])
      ;zmax = max([(zarr)(*),zz])
      zmin = min(zz)
      zmax = max(zz)
      dln_display,zarr,arr1,arr1,xtit='Zeta (deg)',ytit='Eta (deg)',tit='Species '+species[i],min=zmin,max=zmax,$
              position=pos,/noerase,xminor=5,xticklen=0.03,yminor=5,yticklen=0.02,maskv=1000
      oplot,xx,yy,ps=8,sym=1.5,co=0
      plotc,xx,yy,zz,ps=8,xtit='Zeta (deg)',ytit='Eta (deg)',tit='Species '+species[i],position=pos,$
            colpos=colpos,/noerase,/over,min=zmin,max=zmax,bottom=1,ncolors=253
  
      ; Overplot the circle
      phi = cgscalevector(findgen(100),0.0,2*!dpi)
      oplot,1.5*sin(phi),1.5*cos(phi),co=255,thick=1.5
  
      ;ps_close
      ;ps2jpg,psfile+'.eps',/eps
  
      ;stop
      wait,1
    endif
  
    ;stop
  
    BOMB3:

  Endfor

  
  ;stop
  
  
  ;---------------------------------------------
  ; Remove the sky lines from all of the fibers
  ;---------------------------------------------
  print,'Subtracting Sky Lines from Spectra'
  ; Loop through the fibers
  For i=0,nfibers-1 do begin
  
    if (i+1) mod 25 eq 0 then print,strtrim(i+1,2),'/',strtrim(nfibers,2)
  
    ; The plugmap index for this fiber
    ; fiberid=1 is at the top of the detector or index=299
    ; index = 300-fiberid
    ;iplugind = where(plugmap.fiberdata.fiberid eq i+1,niplugind)
    iplugind = where(plugmap.fiberdata.spectrographid eq 2 and $
                     plugmap.fiberdata.fiberid eq 300-i,niplugind)
  
    ; No information for this fiber
    if niplugind eq 0 then begin
      print,'No information for Fiber=',strtrim(i+1,2),' in the plugmap file'
      goto,BOMB
    endif
    objtype  = plugmap.fiberdata[iplugind].objtype
  
    ; Loop through the chips
    For j=0,2 do begin
  
      ;ieta = plugmap.fiberdata[i].eta
      ;izeta = plugmap.fiberdata[i].zeta
      izeta = plugmap.fiberdata[iplugind].zeta
      ieta = plugmap.fiberdata[iplugind].eta
  
      ;fiber = frame.(j).data[*,i,0]
      ;wave = frame.(j).data[*,i,1]
      ;var = frame.(j).data[*,i,2]
      fiber = frame.(j).flux[*,i]
      wave = frame.(j).wavelength[*,i]
      err = frame.(j).err[*,i]
      wcoef = reform(frame.(j).wcoef[i,*])
      lsfcoef = reform(frame.(j).lsfcoef[i,*])
      ;lsfcoef = reform(lsfdata.(j).data[i,*])
  
      wsi = sort(wave)
      wr = [min(wave),max(wave)]
      dw = abs(MEDIAN(wave[1:*]-wave[0:npix-2]))
      pix = lindgen(npix)
  
      ; Initializing output arrays
      fibersub = fiber
      skyspec = fiber*0.0
      skyerr = fiber*0.0
  
      ; Get sky lines to subtract for this spectrum
      uselines = where(modlinestr.use eq 1 and modlinestr.wave ge wr[0]-2 and $
                       modlinestr.wave le wr[1]+2,nuselines)
  
      ; Loop through the lines
      For k=0,nuselines-1 do begin
  
        imodline = modlinestr[uselines[k]]
        medflux = imodline.medflux * FUNC_POLY2D(izeta,ieta,imodline.coef)
        ;lwave = imodline.wave             ; This is the median of the MEASURED wavelengths
        lwave = imodline.model_wave        ; model wavelength
  
        usepix = where(abs(wave-lwave) le 50.0*dw,nusepix)  ; 5.0
        linepix = pix[usepix]*xscale
  
        ; Get the Pixel center
        cenpix = spline(wave[wsi],pix[wsi],lwave,/double)*xscale
 
        ; Make the line
        if imodline.doublet eq 0 then begin
          line = medflux * LSF_GH(linepix,cenpix,lsfcoef)
        endif else begin  ; doublet
          lwave1 = lwave - 0.5*imodline.dbl_wsep
          lwave2 = lwave + 0.5*imodline.dbl_wsep
          cenpix1 = spline(wave[wsi],pix[wsi],lwave1,/double)*xscale
          cenpix2 = spline(wave[wsi],pix[wsi],lwave2,/double)*xscale          
          line1 = medflux * LSF_GH(linepix,cenpix1,lsfcoef)
          line2 = medflux * LSF_GH(linepix,cenpix2,lsfcoef)
          line = line1 + line2
        endelse  

        ;line2 = medflux * LSF_GH(linepix,cenpix+0.1,lsfcoef)
        ;ggg = where(linestr.fiber eq i and linestr.chip eq j+1 and abs(linestr.lsffit_pars[1]-cenpix) lt 1,nggg)
        ;if nggg gt 0 then begin
        ;  line2 = medflux * LSF_GH(linepix,linestr[ggg[0]].lsffit_pars[1],lsfcoef)
        ;  diff = linestr[ggg[0]].lsffit_pars[1]-cenpix
        ;  print,cenpix,linestr[ggg[0]].lsffit_pars[1],diff      
        ;endif
  
        ; Add to skyspec
        skyspec[usepix] += line
        skyerr[usepix] += sqrt(line)   ; sky error is just sqrt(flux)
  
        ; Subtract from fiber spectrum
        fibersub[usepix] -= line
  
        ; Plotting
        ;if i ge 99 and i le 108 then pl=1 else pl=0
        ;if i gt 30 then pl=1 else pl=0
        ;pl = 1
        ;pl = 0
        ;if objtype eq 'SKY' then pl=1
        ;pl = 1
        pl = 0
        if keyword_set(pl) then begin
          ;yr = [0,max(fiber[usepix])]
          yr = [min(fibersub[usepix]),max(fiber[usepix])]
          yr = [yr[0],yr[1]+range(yr)*0.1]
          plot,linepix/xscale,fiber[usepix],/nodata,yr=yr,xs=1,ys=1,$
               tit='Fiber '+strtrim(i+1,2)+' Chip '+chiptag[j]+' Line '+strtrim(k+1,2)
          oplot,linepix/xscale,fiber[usepix],co=250
          oplot,linepix/xscale,fibersub[usepix],co=200
          oplot,linepix/xscale,line,co=150
          ;oplot,linepix/2,line2,co=80
          legend,['Before Subtraction','After Subtraction','Sky Line'],textcolor=[250,200,150],/top,/left,charsize=1.2
          ;wait,0.3
          ;stop
        endif
  
        ;stop
  
      Endfor ; line loop
  
      ; Plotting
      ;if i ge 99 and i le 108 then pl=1 else pl=0
      ;if i gt 30 then pl=1 else pl=0
      ;pl = 1
      pl = 0 ;1 ;0
      if keyword_set(pl) then begin
        yr = [-8000,max(fiber)]
        yr = [yr[0],yr[1]+range(yr)*0.2]
        plot,fiber,/nodata,yr=yr,xs=1,ys=1,tit='Fiber '+strtrim(i+1,2)+' Chip '+chiptag[j]
        oplot,fiber
        oplot,skyspec,co=250,linestyle=2
        oplot,fibersub,co=150
        legend,['Before Subtraction','After Subtraction','Sky Lines'],textcolor=[255,150,250],charsize=1.2,$
               pos=[30,yr[1]-range(yr)*0.02]
        ;stop
      endif
  
  
      ; Save the results
      ; The data planes: [spec, wave, error, flag, sky, errsky]
      ;outframe.(j).data[*,i,0] = fibersub
      ;outframe.(j).data[*,i,4] = skyspec
      ;outframe.(j).data[*,i,5] = skyvar
      outframe.(j).flux[*,i] = fibersub
      outframe.(j).sky[*,i] = skyspec
      outframe.(j).skyerr[*,i] = skyerr
  
      ; Do I need to add anything to FLAGS???
  
      ;stop
  
    Endfor ; chip loop
  
    BOMB:
    
    ;if objtype eq 'SKY' then stop
  
    ;stop
  
  Endfor  ; fiber loop
  
  
  ; If there is still some "continuum" sky left then we can take
  ; a median of the "sky line subtracted" sky fiber spectra and
  ; then remove that from each fiber (maybe use the closest one).
  ;
  ; We could use a similar procedure where we normalize the continua
  ; spectra and then fit the normalization with ZETA/ETA across the
  ; plate.
  ;
  ; Interpolate over the regions with emission lines in the
  ; line-subtracted sky fibers
  
  ; Will the sky continuum screw up the determination of the lines???????
  ; Maybe I should remove a median filtered "background" from the sky
  ; spectra before fitting them with Gaussians, etc.
  
  
End ; model line fitting

Else: begin
  print,'suboption = ',suboption,' Not supported'
  return
End

ENDCASE

;stop

if keyword_set(stp) then stop

end
  
