;+
;
; APSKYSUB_SKYLINEFIT
;
; This fits airglow line fluxes to sky fibers using the LSF
;
; INPUTS:
;  frame     The frame of 1D extracted, wavelength-calibrated spectra.
;  skyindex  An array with sky fiber indices.
;  airstr    The airglow linelist.
;  /pl       Plot the fits.
;  /verbose  Verbose output.
;  /stp      Stop at the end of the program.
;
; OUTPUTS:
;  linestr   The linelist of airglow lines and the best-fit parameters.
;
; USAGE:
;  IDL>apskysub_skylinefit,frame,skyindex,airstr,linestr
;
; By D. Nidever  Feb. 2012
;-
pro apskysub_skylinefit,frame,skyindex,airstr,linestr,pl=pl,stp=stp,verbose=verbose

apgundef,linestr

; Not enough inputs
if n_elements(frame) eq 0 or n_elements(skyindex) eq 0 or n_elements(airstr) eq 0 then begin
  print,'Syntax - apskysub_skylinefit,frame,skyindex,airstr,linestr,pl=pl,verbose=verbose,stp=stp'
  return
endif

chiptag = ['a','b','c']

sz = size(frame.(0).flux)
npix = sz[1]
nsky = n_elements(skyindex)

; Loop through the sky fibers
For i=0,nsky-1 do begin

  ifiber = skyindex[i]

  ; Loop through the chips
  For j=0,2 do begin

    ichip = j

    lsfcoef = reform(frame.(ichip).lsfcoef[ifiber,*])
    wcoef = reform(frame.(ichip).wcoef[ifiber,*])
    fiberflux = reform(frame.(j).flux[*,ifiber])
    fibererr = reform(frame.(j).err[*,ifiber])
    fibermask = reform(frame.(j).mask[*,ifiber])
    fiberwave = reform(frame.(j).wavelength[*,ifiber])
    modelflux = fiberflux*0.0
    fibersnr = fiberflux/(fibererr>1)  

     ; Increase the noise in the error spectrum
    ;   NOT SURE THIS IS NEEDED ANYMORE!!!
    diff = shift(fiberflux,1)-fiberflux
    sig = MAD(diff,/zero)
    mederr = median(fibererr)
    adderr = sqrt((sig^2 - mederr^2) > 0)  ; want total median err to be SIG
    if adderr gt 5 then $
      fibererr = sqrt( fibererr^2 + adderr^2 )
    ;print,'ADDERR=',adderr  

     ; Do a better job of removing the continuum
    med0 = MEDFILT1D(fiberflux,101,/edge)
    temp = fiberflux
    sig1 = MAD(temp)
    mask = long(abs(temp-med0) gt 2*sig1)
    mask = convol(mask,indgen(5))
    mask = mask/(mask>1)
    bd1 = where(mask eq 1,nbd1)
    if nbd1 gt 0 then temp[bd1]=!values.f_nan
    med1 = MEDFILT1D(temp,101,/edge)

    temp2 = fiberflux
    mask2 = long(abs(temp2-med1) gt 2*sig1)
    mask2 = convol(mask2,indgen(5))
    mask2 = mask2/(mask2>1)
    bd2 = where(mask2 eq 1,nbd2)
    if nbd2 gt 0 then temp2[bd2]=!values.f_nan
    med2 = MEDFILT1D(temp2,101,/edge)

    fiberflux -= med2  ; remove the continuum

    ; Getting the LSF parameters
    binsize = lsfcoef[0]
    Xglobalcenter = lsfcoef[1]
    Horder = lsfcoef[2]
    Porder = lsfcoef[3:3+Horder]
    nGHcoefs = total(Porder+1)
    GHcoefs = lsfcoef[Horder+4:Horder+4+nGHcoefs-1]

    ; Wing parameters
    wproftype = lsfcoef[Horder+4+nGHcoefs]
    nWpar = lsfcoef[Horder+4+nGHcoefs+1]
    WPorder = lsfcoef[Horder+4+nGHcoefs+2:Horder+4+nGHcoefs+2+nWpar-1]
    nWcoefs = total(WPorder+1)
    Wcoefs = lsfcoef[Horder+4+nGHcoefs+2+nWpar:Horder+4+nGHcoefs+2+nWpar+nWcoefs-1]


    ;; Select lines for this fiber
    ;fibchind = where(linestr.skyfiber eq ifiber and linestr.chip eq ichip+1,nind)
    ;chiplinestr = linestr[fibchind]
    ;
    ;; Match the detected lines with known Airglow lines
    ;dcr = MEDIAN([chiplinestr.gpar[2]])*3
    ;SRCOR2,airstr.wave,airstr.wave*0,chiplinestr.gwave,chiplinestr.gwave*0,dcr,mind1,mind2,opt=1,/silent
    ;dum = where(mind1 gt -1,nmatch)
    ;
    ;; We have some matches
    ;if nmatch gt 0 then begin
    ;  chiplinestr[mind2].model_match = 1  
    ;  chiplinestr[mind2].model_id = airstr[mind1].id
    ;  chiplinestr[mind2].model_wave = airstr[mind1].wave
    ;  chiplinestr[mind2].model_type = airstr[mind1].name  ;airstr[mind1].type
    ;endif
    ;
    ;gdlines = where(chiplinestr.model_match eq 1,ngdlines)

    ; Get all AIRGLOW lines for this chip
    ;gdlines = where(airstr.wave ge min(fiberwave) and airstr.wave le max(fiberwave),ngdlines)
    ;gdlines = where(airstr.wave ge min(fiberwave) and airstr.wave le max(fiberwave) and airstr.emission gt 50,ngdlines)
    gdlines = where(airstr.wave ge min(fiberwave) and airstr.wave le max(fiberwave),ngdlines)
    chiplinestr = airstr[gdlines]



    ; Start the CHIPLINESTR structure
    chiplinestr = REPLICATE({skyfiber:0L,fiber:0L,chip:0L,x:0.0,model_id:0L,model_wave:0.0d0,model_emission:0.0,$
                             model_doublet:0L,model_dbl_wsep:0.0d0,$
                             model_dbl_xsep:0.0d0,model_type:'',lsffit_pars:dblarr(3),$
                             lsffit_perror:dblarr(3),lsffit_wave:0.0d0,lsffit_flux:0.0,lsffit_chisq:0.0,$
                             lsffit_status:0L,lsffit_ghcoefs:dblarr(nghcoefs+nwcoefs)},ngdlines)
    chiplinestr.skyfiber = i
    chiplinestr.fiber = ifiber
    chiplinestr.chip = ichip
    chiplinestr.model_id = airstr[gdlines].id
    chiplinestr.model_wave = airstr[gdlines].wave
    chiplinestr.model_type = airstr[gdlines].name
    if tag_exist(airstr,'EMISSION') then chiplinestr.model_emission=airstr[gdlines].emission
    if tag_exist(airstr,'DOUBLET') then chiplinestr.model_doublet=airstr[gdlines].doublet
    if tag_exist(airstr,'DBL_WSEP') then chiplinestr.model_dbl_wsep=airstr[gdlines].dbl_wsep

    si = sort(fiberwave)
    x = findgen(npix)

    ; Get X-values for the airglow lines
    si = sort(fiberwave)
    x = findgen(npix)
    for k=0,ngdlines-1 do chiplinestr[k].x = spline(fiberwave[si],x[si],chiplinestr[k].model_wave)

    ; Convert DOUBLET WSEP to pixels
    dwall = abs(fiberwave[1:*]-fiberwave[0:npix-2])
    dblind = where(chiplinestr.model_doublet eq 1,ndblind)
    for k=0,ndblind-1 do begin
      ind1 = dblind[k]
      dw1 = dwall[round(chiplinestr[ind1].x)]
      chiplinestr[ind1].model_dbl_xsep = chiplinestr[ind1].model_dbl_wsep / dw1
    end

    ; Need to fit all lines together like in APLSF
    ;----------------------------------------------
    nsigfit = 8


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
      if ngd gt 0 then begin
        xlo = min(gd)
        xhi = max(gd)
        mask[xlo:xhi] = 1
        xloarr[k] = xlo
        xhiarr[k] = xhi
      endif else xloarr[k]=-1
    endfor
    useind = where(mask eq 1,nuseind)

    ; Input values
    xin = x[useind]
    specin = fiberflux[useind]
    errspecin = fibererr[useind]

    ; "Downweight" pixel with low S/N
    snr_thresh = 5 ; 1
    bdpix = where(specin/(errspecin>1) lt snr_thresh,nbdpix)
    if nbdpix gt 0 then errspecin[bdpix] = ( errspecin[bdpix]*5 < 1e30 )  ; 10

    ; lo/hi indices for the "input" arrays
    npixarr = xhiarr-xloarr+1
    loarr = lonarr(ngdlines)
    hiarr = lonarr(ngdlines)
    ;hiarr[0] = npixarr[0]-1
    for k=0,ngdlines-1 do begin
      loarr[k] = where(useind eq xloarr[k])
      hiarr[k] = where(useind eq xhiarr[k])
    endfor

     ; Initial parameters
    initpar = dblarr(3*ngdlines+1+nGHcoefs+nWcoefs)
    ; We need to multiply the heights by sig*sqrt(2*pi) because the LSFs
    ;  are normalized but the Gaussians that APPEAKFIT.PRO have an
    ;  area under the curve of ht*sig*sqrt(2*pi).
    initpar[0:3*ngdlines-3:3] = (fiberflux[chiplinestr.x]>1) * 1.0 * sqrt(2*!dpi)   ; flux
    initpar[1:3*ngdlines-2:3] = chiplinestr.x
    initpar[2:3*ngdlines-1:3] = 0.0      ;            constant Y-offset
    initpar[3*ngdlines] = Xglobalcenter
    initpar[3*ngdlines+1:3*ngdlines+nGHcoefs] = GHcoefs
    if nWcoefs gt 0 then $
      initpar[3*ngdlines+nGHcoefs+1:3*ngdlines+nGHcoefs+nWcoefs] = Wcoefs

    ; Initial parameter constraints
    npars = n_elements(initpar)
    parinfo = REPLICATE({limited:[0,0],limits:[0.0,0.0],fixed:0},npars)
    parinfo[3*ngdlines:*].fixed = 1  ; Keep all LSF params fixed
    parinfo[0:3*ngdlines-3:3].limited = [1,0] ; flux must be positive
    parinfo[0:3*ngdlines-3:3].limits[0] = 0.0
    parinfo[1:3*ngdlines-2:3].fixed = 1       ; keep x-values fixed
    parinfo[2:3*ngdlines-1:3].fixed = 1       ; keep y-offsets fixed
    parinfo[3*ngdlines:*].fixed = 1     ; keep all LSF params fixed    

    fa = {x:xin,y:specin,err:errspecin,binsize:binsize,nlines:ngdlines,loarr:loarr,$
          hiarr:hiarr,porder:porder,wproftype:wproftype,wporder:wporder,$
          doublet:chiplinestr.model_doublet,dbl_sep:chiplinestr.model_dbl_xsep}

    ; Initial first guess
    yfit0 = SKYFIT_LSF_GH(xin,initpar,_extra=fa)


    ; First fit, only heights allowed to vary
    ;-----------------------------------------
    initpar1 = initpar
    parinfo1 = parinfo
    t0 = systime(1)
    par1 = MPFIT('skyfit_lsf_gh_dev',initpar1,$
                  parinfo=parinfo1,dof=dof,status=status1,bestnorm=chisq1,$
                  perror=perror1,functargs=fa,niter=niter1,autoderivative=0,maxiter=15,/quiet)
    yfit1 = skyfit_lsf_gh(xin,par1,_extra=fa)
    t1 = systime(1)


    ; Fit each line separately
    ;--------------------------
    ;  This is essentially what is done below (but fitting
    ;  all fiber simultaneously), but sometimes that code
    ;  craps out and doesn't get the right answer
    ;  This simple code is fairly safe and gets very close.
    allpar = par1
    lastallpar = allpar
    endflag = 0
    count = 0
    t0 = systime(1)     
    percheight = fltarr(ngdlines)+1e5
    diffcenter = fltarr(ngdlines)+1e5
    WHILE (endflag eq 0) do begin

      ; Loop through the lines
      for k=0,ngdlines-1 do begin
         if percheight[k] lt 0.01 and abs(diffcenter[k]) lt 0.01 then goto,skipfitline
         lo = loarr[k]
        hi = hiarr[k]
        xmid = chiplinestr[k].x

        ; Remove all other lines
        tpar = allpar
        tpar[k*3:k*3+2] = 0
        yfit_others = skyfit_lsf_gh(xin,tpar,_extra=fa)
        resid = specin-yfit_others

        fa2 = {x:xin[lo:hi],y:resid[lo:hi],err:errspecin[lo:hi],binsize:binsize,nlines:1,$
               loarr:[0],hiarr:[hi-lo],porder:porder,wproftype:wproftype,wporder:wporder,$
               doublet:[chiplinestr[k].model_doublet],dbl_sep:[chiplinestr[k].model_dbl_xsep]}

        ; Fake a doublet, skyfit_lsf_gh_dev can't handle single non-doublets
        if chiplinestr[k].model_doublet eq 0 then begin
          fa2.doublet = 1
          fa2.dbl_sep = 0.0
        endif

        ; Initial parameter guesses and constraints
        initpar2 = [allpar[k*3:k*3+2], allpar[3*ngdlines:*]]
        parinfo2 = [parinfo[k*3:k*3+2], parinfo[3*ngdlines:*] ]
        parinfo2[0].fixed = 0
        parinfo2[0].limited[0] = 1
        parinfo2[0].limits[0] = 0
        parinfo2[1].fixed = 0
        parinfo2[1].limited = 1
        parinfo2[1].limits = [-3,3]+initpar2[1]
        parinfo2[2].fixed = 1
        parinfo2[3:*].fixed = 1
        ;parinfo2[3].fixed = 0
        par2 = MPFIT('skyfit_lsf_gh_dev',initpar2,$
                         parinfo=parinfo2,dof=dof2,status=status2,bestnorm=chisq2,$
                         perror=perror2,functargs=fa2,niter=niter2,autoderivative=0,$
                         nfev=nfev2,ftol=1e-10,gtol=1e-10,xtol=1e-10,maxiter=50,/quiet)
        yfit2 = skyfit_lsf_gh(xin[lo:hi],par2,_extra=fa2)

        ; fix flux if it's not a doublet
        if chiplinestr[k].model_doublet eq 0 then par2[0]*=2

        ; Stuff back into the large parameter array
        allpar[k*3:k*3+2] = par2[0:2]

        skipfitline:

      endfor

      allyfit = skyfit_lsf_gh(xin,allpar,_extra=fa)

      ; Are we ending
      diffheight = allpar[0:3*ngdlines-3:3]-lastallpar[0:3*ngdlines-3:3]
      percheight = abs(diffheight)/(lastallpar[0:3*ngdlines-3:3]>1)
      diffcenter = allpar[1:3*ngdlines-2:3]-lastallpar[1:3*ngdlines-2:3]
      if (max(percheight) lt 0.01 and max(abs(diffcenter)) lt 0.01) or count gt 20 then endflag=1
if keyword_set(verbose) then print,max(percheight),max(abs(diffcenter)),count

      ;plot,specin,/xsty,tit='Iter='+strtrim(count,2)
      ;oplot,allyfit,co=250
      ;oplot,specin-allyfit-4000

       lastallpar = allpar
       count++

     ENDWHILE
    t1 = systime(1)
    ;plot,specin,/xsty
    ;oplot,allyfit,co=250
    ;oplot,specin-allyfit-4000
    chisq = total(((specin-allyfit)/errspecin)^2)


    ; Let the heights AND centers vary a little bit
    ;----------------------------------------------
    initpar3 = allpar
    parinfo3 = parinfo1
    parinfo3[0:3*ngdlines-3:3].fixed = 0
    parinfo3[0:3*ngdlines-3:3].limited[0] = 1
    parinfo3[0:3*ngdlines-3:3].limits[0] = 0
    parinfo3[1:3*ngdlines-2:3].fixed = 0
    parinfo3[1:3*ngdlines-2:3].limited = 1
    parinfo3[1:3*ngdlines-2:3].limits[0] = -3 + initpar3[1:3*ngdlines-2:3]
    parinfo3[1:3*ngdlines-2:3].limits[1] = 3 + initpar3[1:3*ngdlines-2:3]
    parinfo3[2:3*ngdlines-1:3].fixed = 1       ; keep y-offsets fixed
    t0 = systime(1)
    par3 = MPFIT('skyfit_lsf_gh_dev',initpar3,$
                     parinfo=parinfo3,dof=dof3,status=status3,bestnorm=chisq3,$
                     perror=perror3,functargs=fa,niter=niter3,autoderivative=0,$
                     nfev=nfev3,ftol=1e-10,gtol=1e-10,xtol=1e-10,maxiter=50,/quiet)
    t1 = systime(1)
    yfit3 = skyfit_lsf_gh(xin,par3,_extra=fa)
 

    ; Final parameters
    yfit = yfit3
    fpar = par3
    fperror = perror3 * sqrt(chisq3/dof3) 
    fchisq = chisq3/dof
    status = status3

    ; Plotting
    ;----------
    if keyword_set(pl) then begin
      xr = [0,n_elements(xin)-1]
      plot,specin,xs=1,xr=xr,ytit='Counts',tit='Sky Fiber='+strtrim(i,2)+' Chip='+strtrim(ichip,2)+$
           ' Chisq='+stringize(fchisq,ndec=3)+' Nlines='+strtrim(ngdlines,2),charsize=1.3
      oplot,yfit,co=250,linestyle=2
      oplot,specin-yfit-4000   ; residuals
      al_legend,['Data','Fit'],textcolor=[255,250],/top,/left,charsize=1.2
    endif

    ; Print the parameters
    fGHcoefs = fpar[3*ngdlines+1:*]  ; Gauss-Hermite parameters
    if keyword_set(verbose) then print,i,ichip,ngdlines,fchisq,format='(I5,I5,I5,F9.2)'

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
;stop
    ; Save the best-fitting model spectrum
    modelflux[useind] = yfit

    BOMB2:

    ; Stuff the information back into the large structure
    PUSH,linestr,chiplinestr

    ; Plot the chip spectrum and the model
    ;pl = 0
    ;if keyword_set(pl) then begin
    ;  plot,fiberflux
    ;  oplot,modelflux,co=250,linestyle=2
    ;  oplot,fiberflux-modelflux-5000
    ;  ;stop
    ;endif


  Endfor ; chip loop

Endfor  ; sky fiber loop

if keyword_set(stp) then stop

end
