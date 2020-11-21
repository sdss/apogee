;+
;
; APRV,visitfile,visitstr
;
; This program measures the RV of a single apVisit file
;
; INPUTS:
;  visitfiles  The name of the visit file
;  visitstr    The visitstr structure with information on the
;                 visit spectra for this star.
;
; OUTPUTS:
;  The RVs are added to the headers of the apVisit files
;  visitstr    The visitstr structure with information on the
;                 visit spectra for this star.
;
; USAGE:
;  IDL>aprv,visitfile,visitstr
;
; Assembled by J.Holtzman  October 2011 based on ap1dobject_star
;-

pro aprv,visitfile,visitstr,save=save,dir_plots=dir_plots,stp=stp,pl=pl,noupdate=noupdate, trimgrid = trimgrid

cspeed = 2.99792458d5  ; speed of light in km/s

; Get APOGEE directories
dirs=getdir(apogee_dir,cal_dir,spectro_dir,apogee_vers)
cal_dir=dirs.caldir
lsf_dir = cal_dir+'lsf/'

; Do we have enough inputs
nvisitfiles = n_elements(visitfile)
if nvisitfiles eq 0 then begin
  print,'Syntax - aprv,visitfile,visitstr'
  return
endif

plate_dir = file_dirname(visitfile)

print,''
print,'--------------------------------------------------------------------------------'
print,'  Processing Visit file ',visitfile
print,'--------------------------------------------------------------------------------'
print,''

; Load the visit file
;--------------------
apgundef,str
APLOADVISIT,visitfile,str

visitstr.file = file_basename(visitfile)
visitstr.ra = str.ra
visitstr.dec = str.dec
;visitstr.ra_targ = str.ra_targ
;visitstr.dec_targ = str.dec_targ
visitstr.glon = str.glon
visitstr.glat = str.glat
visitstr.j = str.j
visitstr.j_err = str.j_err
visitstr.h = str.h
visitstr.h_err = str.h_err
visitstr.k = str.k
visitstr.k_err = str.k_err
visitstr.ak_targ = str.ak_targ
visitstr.ak_targ_method = str.ak_targ_method
visitstr.ak_wise = str.ak_wise
visitstr.sfd_ebv = str.sfd_ebv
visitstr.dateobs = str.dateobs
if tag_exist(str,'JDMID') then visitstr.jd=str.jdmid else visitstr.jd=str.jd
if tag_exist(str,'JDMID') then aprvjd=str.jdmid else aprvjd=str.jd
;visitstr.jd = str.jd
visitstr.snr = str.snr
visitstr.starflag = str.starflag
visitstr.starflags = starflag(str.starflag)

; Check that we have a decent spectrum
gdpix = where(str.spec gt 0.0 and (str.mask and badmask()) eq 0,ngdpix)
; modify to allow faint stars with low S/N and/or very bright neighbors
;if ngdpix lt 500 or (str.starflag and badstarflag()) gt 0 then begin
if ngdpix lt 500 then begin
  print,'not enough good points for RV: ',ngdpix, str.starflag
  goto,BOMB1
endif

;----------------------------------------------
; STEP 1:  Measure the Barycentric Correction
;----------------------------------------------
print,'STEP 1: Measuring the Barycentric Correction'
bc = BARYCENTRIC(DOUBLE(str.ra),DOUBLE(str.dec),jd=DOUBLE(aprvjd))  ; apo by default
add_tag,str,'bc',bc,str
print,' BC = ',stringize(bc,ndec=3),' km/s'
visitstr.BC = str.bc

;-------------------------------------------------------------
; STEP 2:  Cross-Correlate with rough synthetic spectral grid
;-------------------------------------------------------------
print,'STEP 2: Cross-correlate with rough synthetic spectral grid'
;apgetgrid,'apg_rvsynthgrid.fits',grid=grid
apgetgrid,'apg_synthgrid.fits',grid=grid
;NEW 9/1: Trim grid for color 
IF keyword_set(trimgrid) THEN grid = aptrimgrid(grid, visitstr)

aprvprep,grid.wave,str,specout,errout,ncorder=3, pl=pl
; set chip limits
pixlim = fltarr(2,3)
specouttot=total(specout,2,/nan) & errouttot=total(errout,2,/nan)
for i=0,2 do begin
  pixlim[*,i]=minmax(where((finite(specout[*,i])) and (specout[*,i] ne 0.0)))
  ; give up on bad spectra
  lo=pixlim[0,i] & hi=pixlim[1,i]
  obserr1=errouttot[lo:hi]
  bderr = where(obserr1 gt 1 or obserr1 le 0.0,nbderr,comp=gderr,ncomp=ngderr)                       
  if nbderr gt ngderr/2 then goto,BOMB1
endfor
ADD_TAG,grid,'PIXLIM',pixlim,grid
ADD_TAG,str,'INTERP_NSPEC',specout,str
ADD_TAG,str,'INTERP_NERR',errout,str

apxcorr_newgrid,grid,specouttot,errouttot,xcorr,bestgrid,pl=pl  ;/sum
;apxcorr_template,grid.wave,grid.data,specout,errout,xcorr,bestgrid,/sum
if abs(xcorr.vrel) gt 1000 then begin
  xcorr.vrel=0.
  xcorr.vrelerr=999999.
endif
ADD_TAG,str,'vrel',xcorr.vrel,str
ADD_TAG,str,'vrelerr',xcorr.vrelerr,str
ADD_TAG,str,'vhelio',xcorr.vrel+str.bc,str
ADD_TAG,str,'syn_wave',grid.wave,str
ADD_TAG,str,'syn_spec',reform(grid.ndata[bestgrid,*]),str
str = CREATE_STRUCT(str,'xcorr',xcorr)
visitstr.synthfile = file_basename(grid.file)
visitstr.rv_feh = grid.metals[bestgrid]
visitstr.rv_teff = grid.teff[bestgrid]
visitstr.rv_logg = grid.logg[bestgrid]
if tag_exist(grid,'ALPHA') then visitstr.rv_alpha=float(grid.alpha[bestgrid])
if tag_exist(grid,'CARBON') then visitstr.rv_carb=float(grid.carbon[bestgrid])
print,' Best synth parameters: [Fe/H]=',stringize(visitstr.rv_feh,ndec=1),' Teff=',stringize(visitstr.rv_teff,ndec=0), ' K logg=',stringize(visitstr.rv_logg,ndec=1)

;-------------------------------------------------------------------
; STEP 3:  Get radial velocity of pieces with ChiSq minimization
;            Comparing to best-fit synthetic spectrum
;-------------------------------------------------------------------
print,'STEP 3: ChiSq Minimization of Spectrum Chunks with Best-fit Synthetic Spectrum'
APNORMSPEC,str,error=normerror,/sepchip,/fixbadpix,/noerrcorr
APRADVEL_PIECES,str,error=radvel_error,pl=pl

; THIS HAS PROBLEMS WITH THE HOT STARS BECAUSE THERE AREN'T ENOUGH
; SPECTRAL FEATURES PER CHUNK.  Use the XCORR value instead.
;if str.vrelerr gt 0.5 then stop

; Determine which Velocity to use
use_chisq = 1  ; use chisq by default  
if n_elements(radvel_error) eq 0 then begin
  if (str.xcorr.vrelerr lt str.radvel.vrelerr) then use_chisq=0  ; use xcorr, lower error
endif else use_chisq=0  ; radvel/chisq error, use xcorr

if (use_chisq eq 1) then begin
  print,'Using Chisq Velocity values'
  str.vrel = str.radvel.vrel
  str.vrelerr = str.radvel.vrelerr
  visitstr.vtype = 1
; Using the XCORR value
endif else begin
  print,'Using XCORR Velocity values'
  str.vrel = str.xcorr.vrel
  str.vrelerr = str.xcorr.vrelerr
  visitstr.vtype = 2
endelse
str.vhelio = str.vrel+str.bc

; Final chisq
;--------------
synwave = str.syn_wave*(1.0d0 + str.vrel/cspeed)
synspec = str.syn_spec
; Loop through the chips
sz = size(str.spec)
npix = sz[1]
diffsq = fltarr(npix,3)
numpts = 0
for k=0,2 do begin
  ; Get the spectrum
  wobs = str.wave[*,k]       ; wavelength
  si = sort(wobs)             ; sort by wavelength
  wobs = wobs[si]
  ;nspec = str.nspec[si,k]      ; normalized spectrum
  spec = str.spec[si,k]      ; normalized spectrum
  cont = str.continuum[si,k]  ; continuum
  nspec = spec/cont
  mask = str.mask[si,k]
  err = str.err[si,k]        ; error
  err = err/cont                                 ; normalized error spectrum
  sky = str.sky[si,k]         ; sky
  ; Interpolate the synthetic spectrum
  ind = where(wobs ge min(synwave) and wobs le max(synwave) and $
              nspec gt 0.0 and (mask and badmask()) eq 0 and sky lt 10*cont,nind)
  medcont = median(spec)         ; make sure the continuum is positive
  if nind gt 10 and medcont gt 0 then begin
    synspec1 = SPLINE(synwave,synspec,wobs[ind])
    ; Measure the difference
    diffsq[ind,k] = (synspec1-nspec[ind])^2/err[ind]^2
    numpts += nind    ; number of points used
  endif
endfor
if numpts gt 0 then chisq = sqrt( total(diffsq) / numpts ) else chisq=-1.
ADD_TAG,str,'CHISQ',chisq,str
print,'Final chisq = ',stringize(str.chisq,ndec=2)

; Mesure Vsini to convolving the best-fit synthetic spectrum with
; various vsini velocities and find best chisq match.

;stop

print,' Vhelio = ',stringize(str.vhelio,ndec=3),' km/s'

; Calculate VLSR, VGSR
VCONV,str.vhelio,str.glon,str.glat,vlsr,vgsr
add_tag,str,'VLSR',vlsr,str
add_tag,str,'VGSR',vgsr,str

; Update apVisit header with RV and other parameters
;---------------------------------------------------
head = str.head
sxaddpar,head,'SNR',visitstr.snr
sxaddpar,head,'JD',str.jd
; Get HJD
hjd = helio_jd(DOUBLE(aprvjd)-2400000.0D,DOUBLE(str.ra),DOUBLE(str.dec))
sxaddpar,head,'HJD',hjd,' Reduced Heliocentric JD'
sxaddpar,head,'GLON',str.glon,' degrees'
sxaddpar,head,'GLAT',str.glat,' degrees'
sxaddpar,head,'BC',str.bc,' barycentric correction (km/s)'
sxaddpar,head,'VTYPE',visitstr.vtype,' RV type (1=chisq, 2=xcorr)'
sxaddpar,head,'VRAD',str.vrel,' doppler shift (km/s)'
sxaddpar,head,'VRADERR',str.vrelerr,' VRAD error (km/s)'
;sxaddpar,head,'Z',str.vrel/cspeed,' doppler shift'
;sxaddpar,head,'Z_ERR',str.vrelerr/cspeed,' Z error'
if n_elements(radvel_error) eq 0 then $
  sxaddpar,head,'VSCATTER',str.radvel.sigvrel,' km/s' else $
  sxaddpar,head,'VSCATTER',999999.,' km/s'
sxaddpar,head,'VHELIO',str.vhelio,' heliocentric velocity (km/s)'
sxaddpar,head,'VLSR',str.vlsr,' LSR velocity (km/s)'
sxaddpar,head,'VGSR',str.vgsr,' GSR velocity (km/s)'
sxaddpar,head,'RVTEFF',visitstr.rv_teff,' Teff (K) from mini-grid x correlation'
sxaddpar,head,'RVLOGG',visitstr.rv_logg,' logg (dex) from mini-grid x correlation'
sxaddpar,head,'RVFEH',visitstr.rv_feh,' [Fe/H] from mini-grid x correlation'
sxaddpar,head,'RVALPH',visitstr.rv_alpha,' [alpha/H] from mini-grid x correlation'
sxaddpar,head,'RVCARB',visitstr.rv_carb,' [C/H] from mini-grid x correlation'
sxaddpar,head,'XSHIFT',str.xcorr.xshift,' pixel shift'
sxaddpar,head,'CCPEAK',str.xcorr.ccpeak,' cross-corr peak'
sxaddpar,head,'CCPFWHM',str.xcorr.ccpfwhm,' cross-corr peak FWHM'
fxaddpar,head,'SYNTHFIL',grid.file,' RV Synthetic spec grid file',before='JD'
sxaddpar,head,'CHISQ',str.chisq,' Chi Squared'
leadstr = 'APRV '
sxaddhist,leadstr+'APOGEE Reduction Pipeline Version: '+getvers(),head
sxaddhist,leadstr+'Update on '+systime(0),head
info = GET_LOGIN_INFO()
sxaddhist,leadstr+info.user_name+' on '+info.machine_name,head
if not keyword_set(noupdate) then begin
  print,'Updating header in ',visitfile
  MODFITS,visitfile,0,head
endif

; PUT IN THE DATE OF THE UDPATE!!!!

; Update VISITSTR
visitstr.vrel = str.vrel
visitstr.vrelerr = str.vrelerr
visitstr.vhelio = str.vhelio
visitstr.vlsr = str.vlsr
visitstr.vgsr = str.vgsr
if tag_exist(visitstr,'CHISQ') eq 0 then add_tag,visitstr,'CHISQ',0.0,visitstr
visitstr.chisq = str.chisq
;stop

; Make a plot of the spectrum and best-fitting synthetic spectrum
if keyword_set(pl) or keyword_set(save) then begin

    co = 255
    charsize = 1.0
    if keyword_set(save) and not keyword_set(noupdate) then begin
      psfile1 = dir_plots+file_basename(visitfile,'.fits')
      ;visitstr.psfile = psfile1
      ps_open,psfile1,thick=3,/color,/encap
      loadct,39,/silent
      device,/inches,xsize=14,ysize=7
      charsize = 1.3
      co = 0
    endif

    ;!p.multi=[0,1,3]

    x0 = 0.06
    x1 = 0.99
    y0 = 0.09
    y1 = 1.00
    dy = (y1-y0)/3.0
    yoff = 0.05

    ; Chip loop
    for k=0,2 do begin

      owave = str.wave[*,k]
      ospec = str.nspec[*,k]
      ;spec = str.spec[*,i]
      ;cont = str.continuum[*,i]
      ;ospec = spec/cont
      swave = str.syn_wave*(1.0+str.vrel/cspeed)
      ind = where(swave ge min(owave) and swave le max(owave),nind)

      position = [x0,y0+(2-k)*dy,x1,y0+(3-k)*dy-yoff]
      if k gt 0 then noerase=1 else noerase=0
      yr = [0.0,1.5]
      xr = minmax(owave)
      if k eq 2 then xtit = 'Wavelength (!7l!xm)' else xtit=''
      plot,owave/1e4,ospec,xr=xr/1e4,yr=yr,xs=1,ys=1,xtit=xtit,charsize=charsize,$
           ytit='Normalized Flux',position=position,noerase=noerase,$
           xticklen=0.06,yticklen=0.01
      if nind gt 0 then oplot,swave[ind]/1e4,str.syn_spec[ind],co=250,linestyle=2

      if k eq 0 then $
        xyouts,mean(xr)/1e4,yr[1]-0.15*range(yr),'S/N='+stringize(median(str.spec/str.err),ndec=1)+$
               ' H='+stringize(str.h,ndec=2)+' J-Ks='+stringize(str.j-str.k,ndec=2),align=0.5,$
               color=co,charsize=1.2
      if k eq 1 then $
        xyouts,mean(xr)/1e4,yr[1]-0.15*range(yr),'Vrel='+stringize(str.vrel,ndec=2)+' '+$
               ' Vhelio='+stringize(str.vhelio,ndec=2)+'+/-'+stringize(str.vrelerr,ndec=2)+' km/s',$
               align=0.5,color=co,charsize=1.2
      if k eq 2 then $
        xyouts,mean(xr)/1e4,yr[1]-0.15*range(yr),'Best RV-grid (not ASPCAP!) Model Spectrum: [Fe/H]='+stringize(visitstr.rv_feh,ndec=2)+$
               ' Teff='+stringize(visitstr.rv_teff,ndec=0)+' logg='+stringize(visitstr.rv_logg,ndec=2)+$
               '  '+textoidl('\chi^2')+'='+stringize(str.chisq,ndec=2),align=0.5,color=co,charsize=1.2

      if k eq 0 then $
        al_legend,['Observed','Model'],textcolor=[co,250],/bottom,/left,charsize=1.2

    endfor ; chip loop
    xyouts,0.5,0.97,file_basename(visitfile,'.fits')+'  '+visitstr.apogee_id,align=0.5,/norm,charsize=1.4

    if keyword_set(save) and not keyword_set(noupdate) then begin
      ps_close
      FILE_DELETE,psfile1+['.gif','.pdf','.eps.gz'],/allow
      spawn,['convert',psfile1+'.eps',psfile1+'.jpg'],out,errout,/noshell
      ;spawn,['convert',psfile1+'.eps',psfile1+'.pdf'],out,errout,/noshell
      file_delete,psfile1+'.eps',/allow_nonexistent
      ;spawn,['gzip',psfile1+'.eps'],out,errout,/noshell   ; compress the files
    endif

endif

;stop
BOMB1:

if keyword_set(stp) then stop

end
