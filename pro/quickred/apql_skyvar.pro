;+
;
; APQL_SKYVAR
;
; INPUTS:
;  str      The quicklook structure.
;  allstr   The quicklook structure of all previous reads of this exposure.
;  /silent  Don't print anything to the screen.
;
; OUTPUTS:
;  The sky variation values are updated in the STR structure.
;  =error   The error message if one occurred.
;
; USAGE:
;  IDL>apl_skyvar,str,allstr
;
; By D.Nidever  2010
;-
pro apql_skyvar,str,allstr,silent=silent,error=error

; Initialize all values to BAD
str.skyvar_meddev_perc = -1
str.skyvar_stddev_perc = -1
; thes are pointers
;str.skyfiberstr = 
;str.skylinestr = 
;str.medskylinestr = 
str.skyvar_contflux = -1
str.skyvar_contflux_rate = -1
str.skyvar_avglineflux = -1
str.skyvar_avglineflux_rate = -1
str.sky_status = 0

; Error Handling
;------------------
; Establish error handler. When errors occur, the index of the
; error is returned in the variable Error_status:  
; CATCH, Error_status 

;This statement begins the error handler:  
; if (Error_status ne 0) then begin 
;    error = !ERROR_STATE.MSG  
;    print,' WE GET THE FOLLOWING ERROR in APQL_SKYVAR:  '
;    if not keyword_set(silent) then print,error
;    CATCH, /CANCEL 
;    return
; endif

; Not enough inputs
if n_elements(str) eq 0 then begin
  error = 'Not enough inputs'
  if not keyword_set(silent) then print,'Syntax - apql_skyvar,str,allstr,silent=silent,error=error'
  return
endif

; make sure the system variable exists before we proceed
DEFSYSV,'!apql',exists=apql_exists
if not apql_exists then begin
  error = 'APQL_SKYVAR: system variable !apql is not defined'
  if not keyword_set(silent) then print,error
  return
endif

; Check that plugmap exists
if not PTR_VALID(!apql.plugmap.datastr) then begin
  error = 'APQL_SKYVAR: no plugmap data available'
  if not keyword_set(silent) then print,error
  return
endif

; Check that the extracted spectra are there
if not PTR_VALID(str.frame) then begin
  error = 'APQL_SKYCHECK: NO extracted spectra'
  if not keyword_set(silent) then print,error
  return
endif


; Plot up the sky fibers and only show the regions around where
; you expect the sky lines to be
; Also show the estimated positions of the sky lines


;spectro_dir = APGETDIR('APOSPECTRO_DIR',/exists,error=direrr)
;if n_elements(direrr) gt 0 then return
;airglow_dir = spectro_dir+'lib/airglow/'

npix = 2048L

;skyind = where(plugmap.fiberdata.objtype eq 'SKY',nskyind)
skyind = where( (*!apql.plugmap.datastr).fiberdata.objtype eq 'SKY' and $
                (*!apql.plugmap.datastr).fiberdata.holetype eq 'OBJECT' AND $
                (*!apql.plugmap.datastr).fiberdata.spectrographid eq 2 and $
                (*!apql.plugmap.datastr).fiberdata.fiberid ne -1,nskyind)
if nskyind eq 0 then begin
  error = 'APQL_SKYVAR: No Sky fibers'
  if not keyword_set(silent) then print,error
  return
endif 

;skyfiberind = plugmap.fiberdata[skyind].fiberid-1
;skyfiberind = (*!apql.plugmap.datastr).fiberdata[skyind].fiberid-1
; fiberid=1 is at the top of the detector or index=299
; index = 300-fiberid
skyfiberind = 300-(*!apql.plugmap.datastr).fiberdata[skyind].fiberid

skyframe2 = {flux:(*str.frame)[1].flux[*,skyfiberind],err:(*str.frame)[1].err[*,skyfiberind],$
             mask:(*str.frame)[1].mask[*,skyfiberind]}



; Cross-correlate all of the fibers against each other to get
; zero-point shifts and a "global" spectrum
;------------------------------------------
;xshift = fltarr(nskyind)
;fiber_arr = skyframe2.flux*0
;fiber_arr[*,0] = skyframe2.flux[*,0]
;x = lindgen(npix)
;for i=1,nskyind-1 do begin
;  fiber1 = skyframe2.flux[*,i]
;  fiber1 -= MEDFILT1D(fiber1,150,/edge)
;  ;;XCORLB,refspec,fiber1,20,xsh
;  ;XCORLB,skyframe2.flux[*,0],fiber1,20,xsh
;  XCORLB,str.medsky[1].flux,fiber1,20,xsh
;  fiber2 = spline(x,fiber1,x-xsh)
;  if xsh gt 0 then fiber2[0:ceil(xsh)]=0    ; fix the ends
;  if xsh lt 0 then fiber2[npix-floor(abs(xsh))-1:npix-1] = 0
;  fiber_arr[*,i] = fiber2
;  xshift[i] = xsh
;  ;stop
;end

;stop

; Create continuum-subtracted MEDIAN spectrum
;--------------------------------------------
;medfiber0 = MEDIAN(fiber_arr,dim=2)
medfiber0 = str.medsky[1].flux
; remove a smooth background component
cont_medfiber150 = MEDFILT1D(medfiber0,150,/edge)
cont_medfiber50 = MEDFILT1D(medfiber0,50,/edge)
cont_medfiber = cont_medfiber150
cont_medfiber[0:80] = cont_medfiber50[0:80]
cont_medfiber[npix-80:npix-1] = cont_medfiber50[npix-80:npix-1]
medfiber = medfiber0 - cont_medfiber
medfiber[0:4] = 0.0  ; reference pixels
medfiber[npix-4:npix-1] = 0.0

; Maybe just median sky spectrum
apgundef,medlinestr2
;medflux2 = str.medsky[1].flux
;mederr2 = str.medsky[1].err
;medframe2 = str.medsky[1]
medframe2 = str.medsky[1]
medframe2.flux = medfiber
medframe2.err = sqrt(medfiber>1)
APPEAKFIT,medframe2,medlinestr2,/nogauss,/silent,error=error,count=nmedlinestr2
if n_elements(error) gt 0 then begin 
  if not keyword_set(silent) then print,'APQL_SKYVAR: APPEAKFIT: ',error
  return
endif
if nmedlinestr2 eq 0 then begin
  error = 'APQL_SKYVAR: no airglow lines found in median sky spectrum'
  if not keyword_set(silent) then print,error
  return
endif

;chipgap1 = 147
;chipgap2 = 140

; Maybe just check how the entire sky spectrum varies compared to the
; median spectrum

APPEAKFIT,skyframe2,linestr2,/nogauss,/silent,error=error,count=nlinestr2
if n_elements(error) gt 0 then begin
  if not keyword_set(silent) then print,'APQL_SKYVAR: APPEAKFIT: ',error
  return
endif
if nlinestr2 eq 0 then begin
  error = 'APQL_SKYVAR: no airglow lines found'
  if not keyword_set(silent) then print,error
  return
endif


ADD_TAG,linestr2,'MATCH',0,linestr2
ADD_TAG,linestr2,'MEDIND',0L,linestr2
ADD_TAG,linestr2,'MEDFLUX',0.0,linestr2
ADD_TAG,linestr2,'FLUXFRAC',0.0,linestr2

; Only keep good lines
;  this should remove most of the CRs
gdlines = where(linestr2.par0[2] gt 0.5,ngdlines)
if ngdlines eq 0 then begin
  error = 'APQL_SKYVAR: No good lines'
  if not keyword_set(silent) then print,error
  return
endif
linestr2 = linestr2[gdlines]

; Just use the middle chip

x = findgen(npix)


fiberstr = REPLICATE({fiber:-1,nmatch:0,medfrac:0.0,sigfrac:0.0},nskyind)
; print,'n_elemenst(medlinestr2)=',n_elements(medlinestr2)
; print,'n_elemenst(linestr2)=',n_elements(linestr2)

; Loop through the sky fibers
For i=0,nskyind-1 do begin

  gd = where(linestr2.fiber eq i,count)
  if count eq 0 then continue
  fiberlines = linestr2[gd]

  ; Cross-correlate the whole spectrum
  fiber1 = skyframe2.flux[*,i]
  fiber1 -= MEDFILT1D(fiber1,150,/edge)
  XCORLB,fiber1,medframe2.flux,20,xsh

  ; Match the lines 
  SRCOR2,medlinestr2.gaussx,medlinestr2.gaussx*0,fiberlines.gaussx-xsh,fiberlines.gaussx*0,5,ind1,ind2,opt=1,/silent
  dum = where(ind1 ne -1,nmatch)

  if nmatch eq 0 then continue
  ;if n_elements(ind2) le 1 or n_elements(ind1) le 1 then continue
  linestr2[gd[ind2]].match = 1
  linestr2[gd[ind2]].medind = ind1
  linestr2[gd[ind2]].medflux = medlinestr2[ind1].sumflux
  linestr2[gd[ind2]].fluxfrac = fiberlines[ind2].sumflux/medlinestr2[ind1].sumflux

  frac = fiberlines[ind2].sumflux/medlinestr2[ind1].sumflux
  ; print,'n_elements(ind1)=',n_elements(ind1)
  ; print,'n_elements(ind2)=',n_elements(ind2)
  ; print,ind1,ind2,frac
  ; print,'fiberlines[ind2].sumflux=',fiberlines[ind2].sumflux
  ; print,'medlinestr2[ind1].sumflux=',medlinestr2[ind1].sumflux
  ; print,'n_elements(frac)=',n_elements(frac)
  if n_elements(frac) gt 1 then medfrac = median([frac]) else medfrac=frac

  fiberstr[i].fiber = i
  fiberstr[i].nmatch = n_elements(ind1)
  fiberstr[i].medfrac = medfrac
  fiberstr[i].sigfrac = MAD([frac])
  ;fiberstr[i].sigfrac = stddev(frac)

  ;if long(str.readnum) eq 5 then stop

  ;stop

Endfor

; remove the mis-identified fibers
;p = where(fiberstr.fiber gt 0,count)
;if count eq 0 then begin
;endif else begin
;   fiberstr = fiberstr[p]
;endelse
gdfibers = where(fiberstr.fiber gt -1,ngdfibers)
if ngdfibers eq 0 then return
fiberstr = fiberstr[gdfibers]

; Now let's see how the lines vary
;if n_elements(fiberstr.medfrac) gt 1 then med = median([fiberstr.medfrac]) else med=fiberstr.medfrac
med = median([fiberstr.medfrac])
diffperc = abs(1-med)*100
;sig = stddev(fiberstr.medfrac)
sig = mad(fiberstr.medfrac)
sigperc = sig*100
print,'Median flux deviation = ',strtrim(diffperc,2),' %'
print,'Stddev flux deviation = ',strtrim(sigperc,2),' %'
str.skyvar_meddev_perc = diffperc
str.skyvar_stddev_perc = sigperc

; A lot of this variation is actually from the fiber-to-fiber
; throughput variation.  We need to remove this somehow.
; Maybe with cartridge-specific throughput values on disk.


; Save the current linelist
if PTR_VALID(str.skyfiberstr) then PTR_FREE,str.skyfiberstr
str.skyfiberstr = PTR_NEW(fiberstr,/no_copy)
if PTR_VALID(str.skylinestr) then PTR_FREE,str.skylinestr
str.skylinestr = PTR_NEW(linestr2,/no_copy)
str.medskylinestr = PTR_NEW(medlinestr2,/no_copy)
medlinestr2 = *str.medskylinestr

fiberstr = *str.skyfiberstr
linestr2 = *str.skylinestr

; How are the sky lines/flux changing with time
;------------------------------------------------

; Check previous reads
nallstr = n_elements(allstr)
if nallstr gt 0 then begin

  ; Load the line structures
  readstr = REPLICATE({sreadnum:'',readnum:0L,contflux:0.0,avglineflux:0.0,totlineflux:0.0},nallstr)
  for i=0,nallstr-1 do begin

    readnum = allstr[i].readnum
    ; Compare to airglow line structures from previous reads
    if PTR_VALID(allstr[i].medskylinestr) then begin
      imedlinestr = (*allstr[i].medskylinestr)
      imedflux = allstr[i].medsky[1].flux

      readstr[i].sreadnum = readnum
      readstr[i].readnum = long(readnum)
      readstr[i].contflux = median([imedflux])
      readstr[i].avglineflux = mean(imedlinestr.sumflux)
      ;readstr[i].avglinefluxrate = mean(imedlinestr.sumflux)/exptime

      ; Sum up the flux
      dum = imedflux
      dum -= medfilt1d(dum,150,/edge)
      readstr[i].totlineflux = total(dum)

      ; Match to current list
      ;SRCOR2,ilinestr.fiber,ilinestr.gaussx,linestr2.fiber,linestr2.gaussx,0.3,ind1,ind2,opt=1,/silent
      SRCOR2,imedlinestr.fiber,imedlinestr.gaussx,medlinestr2.fiber,medlinestr2.gaussx,0.3,ind1,ind2,opt=1,/silent
      dum = where(ind1 ne -1,nmatch)
      if nmatch eq 0 then begin
        SRCOR2,imedlinestr.fiber,imedlinestr.gaussx,medlinestr2.fiber,medlinestr2.gaussx,0.8,ind1,ind2,opt=1,/silent
        dum = where(ind1 ne -1,nmatch)
      endif

      ADD_TAG,imedlinestr,'READNUM',long(readnum),imedlinestr
      ADD_TAG,imedlinestr,'LINEMATCH',0,imedlinestr
      ADD_TAG,imedlinestr,'LASTIND',0L,imedlinestr

      if nmatch gt 0 then begin
         imedlinestr[ind1].linematch = 1
         imedlinestr[ind1].lastind = ind2
      endif

      PUSH,biglinestr,imedlinestr

    endif ; allstr[i].medskylinestr exists

  endfor  ; loop through previous reads

  ; gd = where(biglinestr.linematch eq 1,ngd)
  ; biglinestr2 = biglinestr[gd]
  ; ui = uniq(biglinestr2.lastind,sort(biglinestr2.lastind))

  ;if n_elements(readstr) ge 2 then begin
  ;  ;plot,[readstr.readnum],[readstr.avglineflux],ps=-1
  ;  plot,[readstr.readnum],[readstr.totlineflux]/32.,ps=-1
  ;  ;plot,[readstr.readnum],[readstr.avglineflux],ps=-1
  ;endif

endif  ; allstr exists


; Get the current avglineflux
;curind = where(readstr.readnum eq long(str.readnum),ncurind)
;medskyflux = fltarr(3*npix)
;for i=0,2 do medskyflux[i*npix:i*npix+npix-1]=str.medsky[i].flux
;str.skyvar_contflux = median(medskyflux)                            ; counts
str.skyvar_contflux = median([medframe2.flux])                        ; counts
str.skyvar_contflux_rate = str.skyvar_contflux/str.exptime          ; counts/sec
str.skyvar_avglineflux = mean( (*str.medskylinestr).sumflux)        ; counts
str.skyvar_avglineflux_rate = str.skyvar_avglineflux/str.exptime    ; counts/sec

; the avglineflux is varying a lot, too much
; there are ~32 airglow lines in the green


; Fit avglineflux as a function of time
;if nallstr gt 0 then begin
;  readnum = [readstr.readnum,long(str.readnum)]
;  avglineflux = [readstr.avglineflux,str.skyvar_avglineflux]
;  coef = poly_fit(readnum,avglineflux,1)
;  ;x = scale_vector(findgen(100),xr[0],xr[1])
;
;  fit = poly(readnum,coef)
;  resid = avglineflux-fit
;  fluxsig = stddev(resid)
;  sigperc = stddev(resid/fit)*100
;
;endif

; Median continuum flux
;med = median([readstr.contflux])
;sig = mad([readstr.contflux])>1
;yr = [med-3*sig,med+3*sig]

; Fit contflux as a function of time
;if nallstr gt 0 then begin
;  readnum = [readstr.readnum,long(str.readnum)]
;  contflux = [readstr.contflux,str.skyvar_contflux]
;  coef2 = poly_fit(readnum,contflux,1)
;
;  fit2 = poly(readnum,coef2)
;  resid2 = contflux-fit2
;  contsig = stddev(resid2)
;  sigperc2 = stddev(resid2/fit2)*100
;  if finite(sigperc2) eq 0 then sigperc2=0.0 ; divide by zero
;
;endif


; SKY STATUS
;--------------
; 0-good, 1-bad
str.sky_status = 0  ; okay for now

; Average line flux
;  don't want it to saturate
if str.skyvar_avglineflux gt 60000L then str.sky_status = 1

; Average line flux RATE
;  don't want the lines to saturate in 600sec
if str.skyvar_avglineflux_rate gt 100L then str.sky_status = 1

; Continuum Flux
if str.skyvar_contflux gt 10000 then str.sky_status = 1

; Continuum Flux RATE
if str.skyvar_contflux_rate gt 20 then str.sky_status = 1

; Median Flux Deviation
if str.skyvar_meddev_perc gt 5 then str.sky_status = 1

; Stddev Flux Deviation
if str.skyvar_stddev_perc gt 5 then str.sky_status = 1

;if long(str.readnum) eq 5 then stop

;stop

end
