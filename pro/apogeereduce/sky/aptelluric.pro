pro aptelluric,frame,plugmap,outframe,tellstar,nearest=nearest,silent=silent,verbose=verbose,$
               starfitopt=starfitopt,specfitopt=specfitopt,error=error,pl=pl,pltelstarfit=pltelstarfit,stp=stp,$
               save=save,plots_dir=plots_dir,preconv=preconv,single=single,visitstr=visitstr,$
               usetelstarfit=usetelstarfit,maxtellstars=maxtellstars,tellzones=tellzones,test=test

;+
;
; APTELLURIC
;
; This corrects for the telluric absorption of the atmosphere.
;
; INPUTS:
;  frame       A structure with the header/data information for an
;                  undersampled frame that has been wavelength
;                  calibrated and airglow line subtracted (not essential)
;  plugmap     The Plug Map structure for this plate
;  /nearest    Use the nearest hot stars to do a telluric correction.
;  /silent     Don't print anything to the screen.
;  =starfitopt  What types of stars to use for the telluric model
;                fitting. 1-Hot stars, 2-All stars.  By default starfitopt=1.
;  =specfitopt What type of spectral fitting method to use when doing
;                the telluci model fitting. 1-Chi-square minimization,
;                2-RMS minimization.  3-RMS for initial species model, then Chi-square
;                Default specfitopt=2.
;  /preconv    Use pre-convolved model spectra.
;  /verbose    Print lots of information to the screen
;  /pl         Make some plots
;  /save       Save plots.
;  =plots_dir  The directory to write the plots to.
;  /stp        Stop at the end of the program.
;
; OUTPUTS:
;  outframe    The same frame but with airglow lines subtracted and
;                sky spectrum (and error) added to the arrays.
;  =error      The error message if one occurred.
;
; USAGE:
;  IDL>aptelluric,frame,plugmap,outframe
;
; By D. Nidever  April 2010
;-

common telluric,convolved_telluric

;setdisp,/silent
apgundef,outframe

; Not enough inputs
if n_elements(frame) eq 0 or n_elements(plugmap) eq 0 then begin
  print,'Syntax - aptelluric,frame,plugmap,outframe,starfitopt=starfitopt,specfitopt=specfitopt,'
  print,'                    preconv=preconv,error=error,silent=silent,verbose=verbose,pl=pl,stp=stp'
  return
endif

; Get APOGEE directories
dirs=getdir(apogee_dir,cal_dir,spectro_dir,apogee_vers,lib_dir)
telluric_dir = lib_dir+'telluric/'
if FILE_TEST(telluric_dir,/directory) eq 0 then begin
  print,'TELLURIC Directory ',telluric_dir,' NOT FOUND'
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
needtags2 = ['HEADER','FLUX','ERR','MASK','WAVELENGTH','SKY','SKYERR','LSFCOEF','WCOEF']
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

apgundef,psfiles

; Getting exposure name
basename = file_basename(frame.(0).filename,'.fits')
expname = first_el(strsplit(basename,'-',/extract),/last)  ; ap1D-a-01600097

; Initialize outframe, add TELLURIC, TELLURICERR
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

    ; Add TELLURIC/TELLURICERR after SKYERR
    if tags[j] eq 'SKYERR' then begin
      chstr = CREATE_STRUCT(chstr,'TELLURIC',fltarr(npix,nfibers))
      chstr = CREATE_STRUCT(chstr,'TELLURICERR',fltarr(npix,nfibers))
    endif
  endfor

  ; Add to the final OUTFRAME
  if i eq 0 then begin
    outframe = CREATE_STRUCT('chip'+chiptag[i],chstr)
  endif else begin
    outframe = CREATE_STRUCT(outframe,'chip'+chiptag[i],chstr)
  endelse
endfor
outframe = CREATE_STRUCT(outframe,'shift',frame.shift)

species = ['CH4','CO2','H2O']
nspecies = n_elements(species)
maxpars=6
; tellstar structure will be saved with summary information for ALL stars
tellstar={im: 0L, scale: fltarr(nfibers,nspecies), sig: fltarr(3), $
          nfit: intarr(3), bestmod: intarr(3), fitpars: fltarr(maxpars,3), $
          fitscale: fltarr(nfibers,nspecies), rchisq: fltarr(nfibers), status: intarr(nfibers), $
          zeta: fltarr(nfibers), eta: fltarr(nfibers), mag: fltarr(nfibers,5)}  ;, fitscale2: fltarr(nfibers,nspecies)}

; Pre-convolve
if keyword_set(preconv) then begin
  ; aptelluric_convolve will return the array of LSF-convolved telluric spectra appropriate
  ;   for the specific wavelength solutions of this frame
  ; There are 3 species, and there may be models computed with different "scale" factor, i.e.
  ;   columns and precipitable water values. If this is the case, we fit each star not
  ;   only for a scaling factor of the model spectrum, but also for which of the models
  ;   is the best fit. For self-consistency, we adopt the model for each species that provides
  ;   the best fit for the largest number of stars, and then go back and refit all stars
  ;   with this model
  ; The telluric array is thus 4D: [npixels,nfibers,nspecies,nmodels]
  iplugind = where(plugmap.fiberdata.spectrographid eq 2 and $
                   plugmap.fiberdata.objtype ne 'SKY')
  convolved_telluric = APTELLURIC_CONVOLVE(frame,fiber=300-plugmap.fiberdata[iplugind].fiberid)
  sz=size(convolved_telluric)
  nmodel=1
  if sz[0] eq 4 then if sz[4] gt 1 then nmodel=sz[4] else nmodel=1
  if nmodel gt 1 then niter=2 else niter=1
endif else begin
  ; Read in the Telluric model spectra
  FITS_READ,telluric_dir+'CH4.fits',telim1,telhead1,message=message1,/no_abort
  FITS_READ,telluric_dir+'CO2.fits',telim2,telhead2,message=message2,/no_abort
  FITS_READ,telluric_dir+'H2O.fits',telim3,telhead3,message=message3,/no_abort
  if message1+message2+message3 ne '' then begin
    print,'ERROR loading the Telluric spectra'
    return
  endif
  telmodelstr = {telim1:telim1,telhead1:telhead1,telim2:telim2,telhead2:telhead2,telim3:telim3,telhead3:telhead3}
endelse

; Default fitting values
if n_elements(starfitopt) eq 0 then starfitopt=1
if n_elements(specfitopt) eq 0 then specfitopt=2

; Get the stars to fit the telluric model spectrum
CASE starfitopt of
0: begin
  print,'skipping telluric correction'
  return
end
; Hot stars
1: begin
  starplugind = where(plugmap.fiberdata.spectrographid eq 2 and $
                     plugmap.fiberdata.holetype eq 'OBJECT' and $
                     plugmap.fiberdata.objtype eq 'HOT_STD',nstar)
end
; All stars
2: begin
  starplugind = where(plugmap.fiberdata.spectrographid eq 2 and $
                     plugmap.fiberdata.holetype eq 'OBJECT' and $
                     (plugmap.fiberdata.objtype eq 'HOT_STD' or plugmap.fiberdata.objtype eq 'STAR'),nstar)
end
else: begin
  print,'STARFITOPT=',starfitopt,' is not a supported option'
  return
end
ENDCASE

if nstar eq 0 then stop,'halt: no telluric stars!'

; special keyword to limits number of telluric stars
if keyword_set(maxtellstars) then begin
  if keyword_set(tellzones) then begin
    ; fill out stars in spatial zones first
    nrad=tellzones[0]
    naz=tellzones[1]
    daz=360./naz
    drad=1.5^2/nrad
    rad=plugmap.fiberdata.eta^2+plugmap.fiberdata.zeta^2
    az=atan(plugmap.fiberdata.eta,plugmap.fiberdata.zeta)*180./!pi+180
    nuse=0
    for irad=0,nrad-1 do begin
      rmin=irad*drad
      rmax=(irad+1)*drad
      for iaz=0,naz-1 do begin
        azmin=iaz*daz
        azmax=(iaz+1)*daz
        gd=where(az[starplugind] ge azmin and az[starplugind] lt azmax and $
                 rad[starplugind] ge rmin and rad[starplugind] lt rmax, ngd)
        if ngd gt 0 then begin
          junk=min(plugmap.fiberdata[starplugind[gd]].mag[1],imin)
          if nuse eq 0 then use=gd[imin] else use=[use,gd[imin]]
          nuse=n_elements(use)
        endif 
      endfor
    endfor
    ; fill out remainder of stars taking the brightest first
    is=sort(plugmap.fiberdata[starplugind].mag[1])
    for i=0,n_elements(is)-1 do begin
      junk=where(use eq is[i],nused)
      if nused eq 0 then use=[use,is[i]]
    endfor
    stmp=starplugind[use[0:maxtellstars-1]]
    starplugind=stmp
    nstar=n_elements(starplugind)
  endif else begin
    tmp=sort(randomu(seed,nstar))
    stmp=starplugind[tmp[0:maxtellstars-1]]
    starplugind=stmp
    nstar=n_elements(starplugind)
  endelse
  print,'using ',nstar,' tellurics: ',plugmap.fiberdata[starplugind].mag[1]
endif

starfiberid = plugmap.fiberdata[starplugind].fiberid
starind = 300-starfiberid
starobjtype = plugmap.fiberdata[starplugind].objtype

print,strtrim(nstar,2),' Stars to fit'

;-------------------------------------------
;  Fit the Telluric model to the stars
;-------------------------------------------
; telstr structure will contain information for the stars used for the fit solution
telstr = REPLICATE({specfitopt:0,fiber:-1,objtype:'',zeta:0.0d0,eta:0.0d0,snr:0.0,par:fltarr(3),$
                    cont_coef:fltarr(7),bestmod:[-1,-1,-1],$
                    status:0,rchisq:99.99,rchisq2:99.99,xshift:0.0},nstar)
telstr.zeta = plugmap.fiberdata[starplugind].zeta
telstr.eta = plugmap.fiberdata[starplugind].eta
telstr.objtype = starobjtype

headstr = 'APTELLURIC: '
APADDPAR,outframe,headstr+'Fitting tellucic spectrum to '+strtrim(nstar,2)+' Stars ('+strjoin(species,', ')+')',/history

; iterate if we have multiple models for a species, so we can choose the one that most stars report
;   is the best
for iter=0,niter-1 do begin

; Loop through the stars
 For i=0,nstar-1 do begin

  print,'Star ',strtrim(i+1,2),'/',strtrim(nstar,2),' Fiber Index ',strtrim(starind[i],2)

  if keyword_set(visitstr) then begin
    j=where(visitstr.fiberid eq starfiberid[i])
    undefine,grid
    wave=[frame.chipa.wavelength[*,starind[i]],frame.chipb.wavelength[*,starind[i]],frame.chipc.wavelength[*,starind[i]]]
    apgetgrid,'apg_rvsynthgrid.fits',grid=grid,wave=reverse(wave)/(1.+visitstr[j].vrel/3.e5)
    bestgrid=where(grid.teff eq visitstr[j].rv_teff and grid.logg eq visitstr[j].rv_logg and grid.metals eq visitstr[j].rv_feh)
    template=reverse(reform(grid.ndata[bestgrid,*]))
;stop
  endif

  ; Do the fitting
  if keyword_set(preconv) then begin
    ; first iteration will try all models, second iteration contrained to a single "best" model
    if iter eq 0 then $
    APTELLURIC_SPECFIT,frame,starind[i],telmodelstr,outstr,specfitopt=2,$
       convolved_telluric=convolved_telluric,template=template  else $
    APTELLURIC_SPECFIT,frame,starind[i],telmodelstr,outstr,specfitopt=specfitopt,$
       convolved_telluric=convolved_telluric,template=template,bestmod=bestmod,test=test
  endif else $
    APTELLURIC_SPECFIT,frame,starind[i],telmodelstr,outstr,specfitopt=specfitopt,template=template

  itelstr = telstr[i]
  STRUCT_ASSIGN,outstr,itelstr,/nozero
  telstr[i] = itelstr

  ; Put these values in the header
  APADDPAR,outframe,headstr+'Fiber='+strtrim(300-starind[i],2)+' Norm='+strjoin(strtrim(string(outstr.par,format='(F8.4)'),2),', '),/history
  tellstar.fitscale[starind[i],*]=outstr.par
  tellstar.rchisq[starind[i]]=outstr.rchisq
  tellstar.status[starind[i]]=itelstr.status

  ;if iter eq 1 then begin
  ;  APTELLURIC_SPECFIT,frame,starind[i],telmodelstr,outstr2,specfitopt=2,$
  ;     convolved_telluric=convolved_telluric,template=template,bestmod=bestmod
  ;  tellstar.fitscale2[starind[i],*]=outstr2.par
  ;endif

  ; Plotting
  if iter eq niter-1 and keyword_set(pltelstarfit) and $
        (keyword_set(pl) or keyword_set(save)) and itelstr.status gt 0 then begin

    ifiber = 300-starind[i]

    co = 255
    if keyword_set(save) then begin
      ;psfile1 = plots_dir+'aptelluric_'+expname+'_hotstarfit_fiber'+strtrim(ifiber,2)
      psfile1 = plots_dir+dirs.prefix+'telluric_'+expname+'_telstarfit_fiber'+strtrim(ifiber,2)
      PUSH,psfiles,psfile1
      ps_open,psfile1,thick=4,/color,/encap
      device,/inches,xsize=42,ysize=7
      loadct,39,/silent
      co = 0
    endif

    medspec = median(outstr.spec)
    yr = [(medspec-0.5*abs(medspec))<0, medspec*1.5]
    xr = minmax(outstr.x)

    !p.multi=[0,1,2]
    plot,outstr.x,outstr.spec,xtit='Pixels',ytit='Counts',xs=1,ys=1,xr=xr,yr=yr,tit='Fiber '+strtrim(ifiber,2)
    ;oplot,x[gd],yfit1,co=250,linestyle=2
    ;oplot,x,yfit1*smspec,co=250,linestyle=2
    oplot,outstr.x,outstr.telluric*outstr.cont,co=250,linestyle=2
    oplot,outstr.x,outstr.spec/outstr.telluric,co=150
    legend,['Original','Telluric','Corrected'],textcolor=[co,250,150],/bottom,/left
    xyouts,mean(xr),yr[1]-0.05*range(yr),'Normalization='+strjoin(stringize(outstr.par[0:2],ndec=4),' '),align=0.5,charsize=1.5,charthick=4

    plot,outstr.x,outstr.spec,xtit='Pixels',ytit='Counts',xs=1,ys=1,xr=[3000,5000],yr=yr,tit='Fiber '+strtrim(ifiber,2)
    ;oplot,x[gd],yfit1,co=250,linestyle=2
    ;oplot,x,yfit1*smspec,co=250,linestyle=2
    oplot,outstr.x,outstr.telluric*outstr.cont,co=250,linestyle=2
    oplot,outstr.x,outstr.spec/outstr.telluric,co=150
    legend,['Original','Telluric','Corrected'],textcolor=[co,250,150],/bottom,/left
    xyouts,mean(xr),yr[1]-0.05*range(yr),'Normalization='+strjoin(stringize(outstr.par[0:2],ndec=4),' '),align=0.5,charsize=1.5,charthick=4
    !p.multi=[0,0,0]

    if keyword_set(save) then ps_close

  endif

  ;stop

 endfor  ; star loop

 ; for first iteration, determine which model spectrum to adopt for each species
 if iter eq 0 then begin
   bestmod=[-1,-1,-1]
   for ispecies=0,2 do begin
     gd=where(telstr.bestmod[ispecies] ge 0,ngd)
     if ngd gt 1 then $
       bestmod[ispecies]=nint(median([telstr[gd].bestmod[ispecies]])) $
     else if ngd eq 1 then $
       bestmod[ispecies]=nint(telstr[gd].bestmod[ispecies])
   endfor
   if nmodel gt 1 and min(bestmod) lt 0 then begin
     error = 'Not enough good telluric spectrum fits'
     if not keyword_set(silent) then print,error
     return
   endif
   tellstar.bestmod=bestmod
 endif
 print,'bestmod: ', bestmod
endfor   ; iteration

; Get the "median" normalization values
gdtelfits = where(telstr.status gt 0,ngdtelfits)
if ngdtelfits gt 0 then begin
  if ngdtelfits ge 2 then mednormpar = MEDIAN(telstr[gdtelfits].par,dim=2) else $
     mednormpar = telstr[gdtelfits[0]].par
endif else begin
  error = 'No good telluric spectrum fits'
  if not keyword_set(silent) then print,error
  return
endelse
print,'Median normalization parameters: ',strtrim(mednormpar,2)

; Determine which fits are good: parameter fits must not be too deviant
med_rchisq = median([telstr[gdtelfits].rchisq])
sig_rchisq = mad([telstr[gdtelfits].rchisq]) > 1.0
if ngdtelfits ge 2 then medpar=median([telstr[gdtelfits].par],dim=2) else $
  medpar = telstr[gdtelfits[0]].par
if ngdtelfits ge 2 then sigpar = mad([telstr[gdtelfits].par],dim=2) else $
  sigpar = [2.0,2.0,2.0]
snr_thresh = 20 
sigpar_thresh = 4.0
; need to split out different species for sigpar criteria to get telstr index right!
gdfits = where(telstr.rchisq le med_rchisq+2.5*sig_rchisq and $
               telstr.status gt 0 and telstr.snr gt snr_thresh and $
               abs(telstr.par[0]-medpar[0]#replicate(1,nstar)) lt sigpar_thresh*(sigpar[0]#replicate(1,nstar)) and $
               abs(telstr.par[1]-medpar[1]#replicate(1,nstar)) lt sigpar_thresh*(sigpar[1]#replicate(1,nstar)) and $
               abs(telstr.par[2]-medpar[2]#replicate(1,nstar)) lt sigpar_thresh*(sigpar[2]#replicate(1,nstar)),$
               ngdfits,comp=bdfits,ncomp=nbdfits)
; Lower constraints if we don't have enough stars
if ngdfits lt 4 then gdfits = where(telstr.rchisq le med_rchisq+2.5*sig_rchisq and $
                                    telstr.status gt 0 and telstr.snr gt snr_thresh,$
                                    ngdfits,comp=bdfits,ncomp=nbdfits)
; Lower constraints again if we still don't have enough stars
if ngdfits lt 4 then gdfits = where(telstr.status gt 0 and telstr.snr gt 10,ngdfits,comp=bdfits,ncomp=nbdfits)

; save bad fit designation
if nbdfits gt 0 then tellstar.status[starind[bdfits]] = -1

if ngdfits lt 4 and not keyword_set(single) then begin
  error = 'Not enough good telluric spectrum fits: '+string(ngdfits)
  if not keyword_set(silent) then print,error
  return
endif
print,strtrim(ngdfits,2),' good telluric spectral fits to use for 2D spatial fitting'
if keyword_set(test) then stop

; Fit the species normalization across the plate
;------------------------------------------------
; This structure is to keep track of the spatial fits for each species
speciesfitstr = REPLICATE({npts:0L,pars:dblarr(maxpars),perror:dblarr(maxpars),$
                           sig:0.0,chisq:0.0,rchisq:0.0,dof:0.0,status:0},nspecies)

For i=0,nspecies-1 do begin

  ; Now fit species normaliation as a function of ZETA/ETA
  xx = [telstr[gdfits].zeta]
  yy = [telstr[gdfits].eta]
  zz = [telstr[gdfits].par[i]]
  err = zz*0+1
  npts = n_elements(xx)

  ; use first order for CH4 and CO2, second order for H2O
  if i eq 2 then npars = 6 else npars = 3;
  ; if not enough stars, use first order for everything
  if ngdfits lt 10 and ngdfits ge 4 then npars=3   ; use linear fit if not enough points
  ; if just one star, then no surface fit!
  if keyword_set(single) then npars=1

  apgundef,status,dof,chisq,rchisq,yfit
  initpars = dblarr(npars)
  initpars[0] = 1.0
  parinfo = REPLICATE({limited:[0,0],limits:[0.0,0.0],fixed:0},npars)
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
  if nbd gt ngdfits/2 then begin
    error = 'More than half the fits are outliers!'
    if not keyword_set(silent) then print,error
    return
  endif
  if ngdfits-nbd lt 4 and npars gt 1 then begin
    error = 'Not enough good telluric spectrum fits after rejection: '+string(ngdfits-nbd)
    if not keyword_set(silent) then print,error
    return
  endif

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
status=0
  if status gt 0 and not keyword_set(single) then begin
    sig = MAD(zz-yfit,/zero)
    err = zz*0+sig

    pars2 = MPFIT2DFUN('func_poly2d',xx,yy,zz,err,pars,status=status2,dof=dof2,$
                    bestnorm=chisq2,parinfo=parinfo2,perror=perror2,yfit=yfit2,/quiet)
    pcerror = perror2*sqrt(chisq2/dof2)

  endif else begin
    sig = MAD(zz-median([zz]))
    pcerror = perror*0
    pcerror[0] = sig
  endelse

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
  APADDPAR,outframe,headstr+'TELLURIC 2D SPATIAL POLYNOMIAL FIT',/history
  APADDPAR,outframe,headstr+'SPECIES '+species[i],/history
  APADDPAR,outframe,headstr+'NPARS = '+strtrim(n_elements(speciesfitstr[i].pars),2),/history
  for k=0,n_elements(speciesfitstr[i].pars)-1 do begin
    parname = 'TLPR'+strtrim(i+1,2)+'_'+strtrim(k+1,2)
    APADDPAR,outframe,parname,speciesfitstr[i].pars[k]
    errname = 'TLER'+strtrim(i+1,2)+'_'+strtrim(k+1,2)
    APADDPAR,outframe,errname,speciesfitstr[i].perror[k]
  endfor
  APADDPAR,outframe,headstr+'SIG = '+stringize(speciesfitstr[i].sig,ndec=5),/history


  ;; Put the coefficients into MODLINESTR
  ;gdmod = where(modlinestr.type eq species[i],ngdmod)
  ;modlinestr[gdmod].coef = pars

  ; Plotting
  psym8
  pl = 1
  if (keyword_set(pl) or keyword_set(save)) and not keyword_set(single) then begin
    ;erase

    if keyword_set(save) then begin
      psfile1 = plots_dir+dirs.prefix+'telluric_'+expname+'_skyfit_'+species[i]
      PUSH,psfiles,psfile1
      ps_open,psfile1,thick=4,/color,/encap
      device,/inches,xsize=10,ysize=10
      loadct,39,/silent
      psym8
    endif

    ; Restore color table
    if file_test(spectro_dir+'lib/colors/coltable1.fits') eq 1 then begin
      fits_read,spectro_dir+'lib/colors/coltable1.fits',rgb,/no_abort
      tvlct,rgb[0,*],rgb[1,*],rgb[2,*]
    endif

    if !d.name eq 'X' then co1=255 else co1=0

    yr = minmax(zz)
    yr = [yr[0]-range(yr)*0.1,yr[1]+range(yr)*0.1]
    if range(yr) eq 0 then yr=[median([zz])-1,median([zz])+1]

    ; Normalization vs. Zeta
    xr = minmax(xx)
    xr = [xr[0]-range(xr)*0.1,xr[1]+range(xr)*0.1]
    if range(xr) eq 0 then xr=[median([xx])-1,median([xx])+1]
    pos = [0.08,0.58,0.50,0.95]
    pos = [0.08,0.07,0.50,0.44]
    plot,xx,zz,ps=1,xtit='Zeta (deg)',ytit='Normalization',tit='Species '+species[i],$
         xr=xr,yr=yr,xs=1,ys=1,position=pos
    oplot,xx,yfit,ps=4,co=250
    legend,['Data','Model'],textcolor=[co1,250],/top,/left

    ; Normalization vs. Eta
    xr = minmax(yy)
    xr = [xr[0]-range(xr)*0.1,xr[1]+range(xr)*0.1]
    if range(xr) eq 0 then xr=[median([xx])-1,median([xx])+1]
    pos = [0.58,0.58,0.99,0.95]
    plot,yy,zz,ps=1,xtit='Eta (deg)',ytit='Normalization',tit='Species '+species[i],$
         xr=xr,yr=yr,xs=1,ys=1,position=pos,/noerase
    oplot,yy,yfit,ps=4,co=250
    legend,['Data','Model'],textcolor=[co1,250],/top,/left

    ; Colored points on the sky
    pos = [0.08,0.08,0.50,0.42]
    pos = [0.08,0.58,0.50,0.95]
    colpos = [0.08,0.49,0.50,0.51]
    if dirs.telescope eq 'lco25m' then maxrad=1.1 else maxrad=1.5
    arr1 = scale_vector(findgen(400),-maxrad,maxrad)
    xarr = arr1#replicate(1.0,400)
    yarr = replicate(1.0,400)#arr1
    zarr = FUNC_POLY2D(xarr,yarr,pars)
    bd = where(sqrt(xarr^2+yarr^2) gt maxrad,nbd)
    zarr[bd] = 1000
    ;zmin = min([(zarr)(*),zz])
    ;zmax = max([(zarr)(*),zz])
    ;zmin = min(zz) < median([zz])-1
    ;zmax = max(zz) > median([zz])+1
    zmin = min(zz) > median([zz])-2.5*(mad([zz])>0.1)
    zmax = max(zz) < median([zz])+2.5*(mad([zz])>0.1)
    dln_display,zarr,arr1,arr1,xtit='Zeta (deg)',ytit='Eta (deg)',tit='Species '+species[i],min=zmin,max=zmax,$
            position=pos,/noerase,xminor=5,xticklen=0.03,yminor=5,yticklen=0.02,maskv=1000
    oplot,xx,yy,ps=8,sym=1.5,co=0
    plotc,xx,yy,zz,ps=8,xtit='Zeta (deg)',ytit='Eta (deg)',tit='Species '+species[i],position=pos,$
          colpos=colpos,/noerase,/over,min=zmin,max=zmax,bottom=1,ncolors=253

    ; Overplot the circle
    phi = scale_vector(findgen(100),0.0,2*!dpi)
    oplot,3.5*sin(phi),3.5*cos(phi),co=255,thick=1.5

    ;ps_close
    ;ps2jpg,psfile+'.eps',/eps

    if keyword_set(save) then begin
      ps_close
    endif else wait,1

    ;stop
    ;wait,1
  endif ; plotting

  ;stop

Endfor  ; species loop

; Convert the figures
if keyword_set(save) then begin

  ; Put the skyfit ones first
  gdskyfit = where(stregex(psfiles,'_skyfit_',/boolean) eq 1,ngdskyfit)
  if ngdskyfit gt 0 then begin
    left = psfiles
    psfiles = psfiles[gdskyfit]
    if ngdskyfit lt n_elements(left) then REMOVE,gdskyfit,left else apgundef,left
    if n_elements(left) gt 0 then PUSH,psfiles,left    
  endif

  print,'Converting figures'
  for i=0,n_elements(psfiles)-1 do begin
    ps2jpg,psfiles[i]+'.eps',/eps,chmod='664'o,/delete
;    spawn,['convert',psfiles[i]+'.eps',psfiles[i]+'.jpg'],out,errout,/noshell
;    spawn,['convert',psfiles[i]+'.eps',psfiles[i]+'.pdf'],out,errout,/noshell
;    spawn,['gzip',psfiles[i]+'.eps'],out,errout,/noshell   ; compress the files
  endfor

  ; Combine them into one PDF
;  print,'Writing combined PDF plots file to ',plots_dir+'aptelluric_'+expname+'.pdf'
;  cmd = 'gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile='+plots_dir+'aptelluric_'+expname+'.pdf '
;  cmd += strjoin(psfiles+'.pdf',' ')
;  spawn,cmd,out,errout
;  file_chmod,plots_dir+'aptelluric_'+expname+'.pdf','664'o
endif

;stop

;goto,SKIPTOEND

;---------------------
; Correct all fibers
;---------------------

tellstar.sig = speciesfitstr.sig
tellstar.nfit = speciesfitstr.npts
tellstar.fitpars = speciesfitstr.pars
for i=0,nfibers-1 do begin

  ;print,'Fiber ',strtrim(i+1,2)
  if (i+1) mod 25 eq 0 or i eq 0 then print,strtrim(i+1,2),'/',strtrim(nfibers,2)

  ; The plugmap index for this fiber
  ; fiberid=1 is at the top of the detector or index=299
  ; index = 300-fiberid
  ;iplugind = where(plugmap.fiberdata.fiberid eq i+1,niplugind)
  iplugind = where(plugmap.fiberdata.spectrographid eq 2 and $
                   plugmap.fiberdata.objtype ne 'SKY' and $
                   plugmap.fiberdata.fiberid eq 300-i,niplugind)

  ; No information for this fiber
  if niplugind eq 0 then begin
    ;print,'No information for Fiber=',strtrim(i+1,2),' in the plugmap file'
    goto,BOMB
  endif

  ; Getting the spectrum, concatenate them
  spec = [frame.(0).flux[*,i], frame.(1).flux[*,i], $
          frame.(2).flux[*,i] ]
  spec = double(spec)
  errspec = [frame.(0).err[*,i], frame.(1).err[*,i], $
             frame.(2).err[*,i] ] > 1
  errspec = double(errspec)

  ; The coordinates for this fiber
  izeta = plugmap.fiberdata[iplugind].zeta
  ieta = plugmap.fiberdata[iplugind].eta

  ; Get the normalizations for this position
  normpar = fltarr(3)
  for k=0,2 do normpar[k] = FUNC_POLY2D(izeta,ieta,speciesfitstr[k].pars)
  if keyword_set(usetelstarfit) then begin
    j=where(300-starfiberid eq i,nj)
    if nj gt 0 then begin
      print,izeta,ieta,normpar
      print,telstr[j].zeta,telstr[j].eta,telstr[j].par
      normpar=telstr[j].par
    endif
  endif
  tellstar.mag[i,*]=plugmap.fiberdata[iplugind].mag
  tellstar.zeta[i]=izeta
  tellstar.eta[i]=ieta
  tellstar.scale[i,*] =normpar

  ; Get convolved model spectrum
  ;-----------------------------
  ; Pre-convolved
  if keyword_set(preconv) then begin

    x = dblarr(3*npix)
    ; create the model array with correct dimensions
    mspec=reform(convolved_telluric[*,i,*])
    ; use the best model as determined from the telluric stars 
    for j=0,2 do mspec[*,j]=reform(convolved_telluric[*,i,j,bestmod[j]])

  ; On-the-fly convolution
  endif else begin

    apgundef,lsfpars,wcoef
    for j=0,2 do PUSH,lsfpars,frame.(j).lsfcoef[i,*]
    for j=0,2 do PUSH,wcoef,frame.(j).wcoef[i,*]
    APTELLURIC_SINGLECONVOLVE,lsfpars,wcoef,telmodelstr,x,mspec

  endelse


  ; Now correct the spectrum
  telnorm = fit_telluric(x,[normpar,1.0,0.0],mspec=mspec)
  spec2 = spec/telnorm
  ;errspec2 = errspec/telnorm   ; this is now done properly below


  ; Measure telluric errors
  ;------------------------
  ; Calculate the error in the normalization values
  normpar_error = fltarr(3)
  step = 1d-5
  for k=0,2 do begin
    dpar_poly = speciesfitstr[k].perror*0d

    ; Calculate the partial derivaties of the 2D polynomial
    ; normalization fit with each poly coefficient.
    normpar1 = FUNC_POLY2D(izeta,ieta,speciesfitstr[k].pars)
    for l=0,n_elements(dpar_poly)-1 do begin
      tpars = speciesfitstr[k].pars
      tpars[l] += step
      normpar2 = FUNC_POLY2D(izeta,ieta,tpars)
      dpar_poly[l] = (normpar2-normpar1)/step
    endfor ; loop through poly coefficients

    ; Propagate the errors:
    ; sig^2 = (del f/ del A sig_A)^2 + (del f/ del B sig_B)^2 + ...
    ; (del f/del A) is the partial derivative we just calculated
    ; sig_A are the PERRORs we measured previously
    normpar_error[k] = sqrt( TOTAL( (dpar_poly * speciesfitstr[k].perror )^2 ) )

  endfor ; loop through species

  ; Now get the telluric error spectrum
  telnorm_error = telnorm*0
  dpar_norm = dblarr(3,3*npix)
  for k=0,2 do begin

    ; Calculate the partial derivative of the telluric spectrum
    ; with each normalization value    
    tnormpar = normpar
    tnormpar[k] += step
    telnorm2 = fit_telluric(x,[tnormpar,1.0,0.0],mspec=mspec)
    dpar_norm[k,*] = (telnorm2-telnorm)/step
  endfor

  ; Propagate the errors
  telnorm_error = sqrt( TOTAL( (dpar_norm * (normpar_error#replicate(1,3*npix)) )^2,1 ) )

  ; Now add this in quadrature to the error spectrum
  ; f=A/B, or final = orig/telluric
  ; (sig_f/f)^2 = (sig_A/A)^2 + (sig_B/B)^2 + cross-term
  ; sig_final^2 = f^2 * ( (sig_orig/orig)^2 + (sig_telluric/telluric)^2 )
  ferrspec = (spec2>1) * sqrt( (errspec/(spec>1))^2 + (telnorm_error/telnorm)^2 )


  ; Put in OUTFRAME
  ;---------------------
  ;  The planes are: [spec, wave, error, flag, sky, errsky, telluric, error_telluric]
  for j=0,2 do outframe.(j).flux[*,i] = spec2[j*npix:j*npix+npix-1]
  ; wavelength remains unchanged
  for j=0,2 do outframe.(j).err[*,i] = ferrspec[j*npix:j*npix+npix-1]      ; error
  for j=0,2 do outframe.(j).telluric[*,i] = telnorm[j*npix:j*npix+npix-1]  ; telluric
  for j=0,2 do outframe.(j).telluricerr[*,i] = telnorm_error[j*npix:j*npix+npix-1]  ; telluric error

  for j=0,2 do begin
    bd=telluricmask(outframe.(j).telluric[*,i],nmask,width=5)
    if nmask gt 0 then outframe.(j).mask[bd,i] = outframe.(j).mask[bd,i] or maskval('SIG_TELLURIC')
  endfor

  ; Plotting
  debugpl = 0 ;1
  if keyword_set(debugpl) then begin
    yr = [0,median([spec])*1.5]
    plot,x,spec,xtit='Pixels',ytit='Counts',xs=1,yr=yr,ys=1,tit='Fiber='+strtrim(i+1,2)
    oplot,x,spec2,co=250
    legend,['Original','Corrected'],textcolor=[255,250],/bottom,/left
    print,normpar
  endif

  BOMB:

  ;stop

endfor

;SKIPTOEND:

;stop

end
