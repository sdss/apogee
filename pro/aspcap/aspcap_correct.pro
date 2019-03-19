;+
; aspcap_correct: routine to correct/calibrate parameters/abundances, i.e.
;    populate param and elem from fparam, felem
;    Also, set flags!
;-
pro aspcap_correct,str,elem,caldir,dr10=dr10,dr13=dr13,noelem=noelem

giant_tepar=mrdfits(caldir+'/giant_tecal.fits',1)
dwarf_tepar=mrdfits(caldir+'/dwarf_tecal.fits',1)
loggpar=mrdfits(caldir+'/loggcal.fits',1)
if keyword_set(dr13) then begin
  a=mrdfits(caldir+'/.fits',1)
  b=mrdfits(caldir+'/giantcal.fits',2)
  c=mrdfits(caldir+'/giantcal.fits',3)
  d=mrdfits(caldir+'/giantcal.fits',4)
  giantcal = { elem: a.elem, elempar: b, errpar: c, errpar_h: d}
  a=mrdfits(caldir+'/dwarfcal.fits',1)
  b=mrdfits(caldir+'/dwarfcal.fits',2)
  c=mrdfits(caldir+'/dwarfcal.fits',3)
  d=mrdfits(caldir+'/dwarfcal.fits',4)
  dwarfcal = { elem: a.elem, elempar: b, errpar: c, errpar_h: d}
endif
giantcal=mrdfits(caldir+'/giantcal.fits',1)
dwarfcal=mrdfits(caldir+'/dwarfcal.fits',1)
;kludge
;errpar = {const: giantcal.errpar[0], tpar: giantcal.errpar[1],  fepar: giantcal.errpar[2], snpar: giantcal.errpar[3]}
;giantcal = { elem: giantcal.elem, elempar: giantcal, errpar: errpar, errpar_h: errpar }
;errpar = {const: dwarfcal.errpar[0], tpar: dwarfcal.errpar[1],  fepar: dwarfcal.errpar[2], snpar: dwarfcal.errpar[3]}
;dwarfcal = { elem: dwarfcal.elem, elempar: dwarfcal, errpar: errpar, errpar_h: errpar }


nelem=n_elements(elem)
if keyword_set(noelem) then nelem=0

mask='FF'XL
help,mask
j=where((str.aspcapflag and mask) gt 0, nj)
print,'limit warning: ', nj

if tag_exist(str,'param_chi2') then param_chi2=str.param_chi2 else param_chi2=str.aspcap_chi2
if tag_exist(str,'snrev') then strsnr=str.snrev else strsnr=str.snr
if tag_exist(str,'class') then class=str.class else class=str.aspcap_class
if tag_exist(str,'x_m') then x_m_exists=1 else x_m_exists=0

; get median chi2 for chi2 flagging criterion
dt=100 & medt=3550+indgen(70)*dt
medchi=fltarr(n_elements(medt)) & medn=intarr(n_elements(medt))
for i=0,n_elements(medt)-1 do begin
  j=where(str.fparam[0] ge medt[i]-dt/2. and str.fparam[0] lt medt[i]+dt/2.,nj)
  medn[i]=nj
  if nj gt 1 then medchi[i]=median(param_chi2[j])
endfor

; do initial flags star by star
for i=0L,n_elements(str)-1 do begin
  if i mod 100 eq 0 then print,'flag: ',i, ' of ',n_elements(str)
  teff=str[i].fparam[0]
  logg=str[i].fparam[1]
  feh=str[i].fparam[3]
  chi2=param_chi2[i]
  snr=strsnr[i]
  starflag=long(str[i].starflag)
  flag=long(str[i].aspcapflag)
  pflag=long(str[i].paramflag)
  jk0=(str[i].j-str[i].k)-1.5*str[i].ak_targ

  ; chi^2
  ;j=min(abs(teff-medt),it)
  ;chicrit=10+fltarr(n_elements(medt))
  ;chicrit=max([[medchi*2],[chicrit]],dim=2)
  ;if chi2 gt chicrit[it] then flag=flag or aspcapflagval('CHI2_WARN')
;
;  chicrit=30+fltarr(n_elements(medt))
;  chicrit=max([[medchi*3],[chicrit]],dim=2)
;  if chi2 gt chicrit[it] then flag=flag or aspcapflagval('CHI2_BAD')
  if chi2/(snr^2/100.^2) gt 30 then flag=flag or aspcapflagval('CHI2_WARN')
  if chi2/(snr^2/100.^2) gt 50 then flag=flag or aspcapflagval('CHI2_BAD')

  ; broad lines
  ;if ((starflag and starflagval('SUSPECT_BROAD_LINES')) gt 0) then flag=flag or aspcapflagval('ROTATION_WARN')
  ; only flag ASPCAP flag if result is from a giant grid without rotation
  if strpos(class[i],'GKg') ge 0 or strpos(class[i],'Mg') ge 0 then begin
    if str[i].rv_ccfwhm/str[i].rv_autofwhm gt 1.5 then flag=flag or aspcapflagval('ROTATION_WARN')
    if str[i].rv_ccfwhm/str[i].rv_autofwhm gt 2.0 then flag=flag or aspcapflagval('ROTATION_BAD')
  endif

  ; S/N 
  if snr lt 70 then flag=flag or aspcapflagval('SN_WARN')
  if snr lt 30 then flag=flag or aspcapflagval('SN_BAD')

  ; color=Teff
  if str[i].h lt 90 and str[i].ak_targ gt -1 and str[i].ak_targ ne 0. then begin
    if abs(teff-aspcap_colorte(jk0,feh)) gt  500 then flag=flag or aspcapflagval('COLORTE_WARN')
    if abs(teff-aspcap_colorte(jk0,feh)) gt 1000 then flag=flag or aspcapflagval('COLORTE_BAD')
  endif

  ; star level flag
  if ((flag and warnaspcapflag()) gt 0) then flag=flag or aspcapflagval('STAR_WARN')
  if ((flag and badaspcapflag()) gt 0) then flag=flag or aspcapflagval('STAR_BAD')
  ; if any parameter is near grid edge, set STAR_BAD
  for ipar=0,n_elements(pflag)-1 do if (pflag[ipar] and paramflagval('GRIDEDGE_BAD')) gt 0 then flag=flag or aspcapflagval('STAR_BAD')

  ; bad targets (for ASPCAP)
  if strpos(str[i].targflags,'EMBEDDED') ge 0 or (strpos(str[i].targflags,'EXTENDED') ge 0 and strpos(str[i].field,'198+08') lt 0) or strpos(str[i].targflags,'M31_CLUSTER') ge 0 then $
     flag=flag or aspcapflagval('STAR_BAD')

  ; load new flags
  str[i].aspcapflag=flag
  str[i].paramflag=pflag

endfor

; we will only calibrate stars that are not flagged as bad by criteria above
gd=where((str.aspcapflag and aspcapflagval('STAR_BAD')) eq 0,ngd)
if ngd le 0 then stop,'NO good stars!! something wrong?'

; define giants and dwarfs, only for stars that are not STAR_BAD
giants=where(str[gd].fparam[1] lt 2./1300.*(str[gd].fparam[0]-3500.)+2. and str[gd].fparam[1] lt 4 and str[gd].fparam[0] lt 7000,ngiants,comp=dwarfs,ncomp=ndwarfs)
giants=gd[giants]
dwarfs=gd[dwarfs]

; cap S/N at 200 for errors
snrerr=strsnr
bd=where(snrerr gt 200.,nbd)
if nbd gt 0 then snrerr[bd]=200.

; param array first
; calibrate Teff: note DR13 used teff, while DR14 uses mh
; giants
flag=str.paramflag[0]
teffclip=aspcap_clip(str.fparam[0],giant_tepar.temin,giant_tepar.temax,flag)
mhclip=aspcap_clip(str.fparam[3],giant_tepar.mhmin,giant_tepar.mhmax,flag)
str[giants].paramflag[0]=flag[giants]
tcal=where(str[giants].fparam[0] ge giant_tepar.caltemin and str[giants].fparam[0] le giant_tepar.caltemax,ngd,comp=tbd,ncomp=nbd)
if ngd gt 0 then begin
  str[giants[tcal]].param[0]=str[giants[tcal]].fparam[0]-(giant_tepar.par[0]+giant_tepar.par[1]*mhclip[giants[tcal]]+giant_tepar.par[2]*mhclip[giants[tcal]]^2)
  str[giants[tcal]].param_cov[0,0]=giant_tepar.rms^2
  str[giants[tcal]].param_cov[0,0]=(exp(giant_tepar.errpar[0]+giant_tepar.errpar[1]*(str[giants[tcal]].fparam[0]-4500)+giant_tepar.errpar[3]*str[giants[tcal]].fparam[3]+giant_tepar.errpar[2]*(snrerr[giants[tcal]]-100.)))^2
endif
if nbd gt 0 then begin
  for ibd=0,nbd-1 do str[giants[tbd[ibd]]].paramflag[0]=str[giants[tbd[ibd]]].paramflag[0] or paramflagval('CALRANGE_BAD')
endif
; dwarfs
flag=str.paramflag[0]
teffclip=aspcap_clip(str.fparam[0],dwarf_tepar.temin,dwarf_tepar.temax,flag)
mhclip=aspcap_clip(str.fparam[3],dwarf_tepar.mhmin,dwarf_tepar.mhmax,flag)
str[dwarfs].paramflag[0]=flag[dwarfs]
tcal=where(str[dwarfs].fparam[0] ge dwarf_tepar.caltemin and str[dwarfs].fparam[0] le dwarf_tepar.caltemax,ngd,comp=tbd,ncomp=nbd)
if ngd gt 0 then begin
  str[dwarfs[tcal]].param[0]=str[dwarfs[tcal]].fparam[0]-(dwarf_tepar.par[0]+dwarf_tepar.par[1]*mhclip[dwarfs[tcal]]+dwarf_tepar.par[2]*mhclip[dwarfs[tcal]]^2)
  str[dwarfs[tcal]].param_cov[0,0]=dwarf_tepar.rms^2
  str[dwarfs[tcal]].param_cov[0,0]=(exp(dwarf_tepar.errpar[0]+dwarf_tepar.errpar[1]*(str[dwarfs[tcal]].fparam[0]-4500)+dwarf_tepar.errpar[3]*str[dwarfs[tcal]].fparam[3]+dwarf_tepar.errpar[2]*(snrerr[dwarfs[tcal]]-100.)))^2
endif
if nbd gt 0 then begin
  for ibd=0,nbd-1 do str[dwarfs[tbd[ibd]]].paramflag[0]=str[dwarfs[tbd[ibd]]].paramflag[0] or paramflagval('CALRANGE_BAD')
endif

; calibrate logg, giants only
teff=str.fparam[0]
; clip logg and mh correction term
flag=str.paramflag[1]
logg=aspcap_clip(str.fparam[1],loggpar.loggmin,loggpar.loggmax,flag)
mh=aspcap_clip(str.fparam[3],loggpar.mhmin,loggpar.mhmax,flag)
str.paramflag[1]=flag

; DR14
dt=str.fparam[0]-(loggpar.rgbsep[0]+(str.fparam[1]-2.5)*loggpar.rgbsep[1]+str.fparam[3]*loggpar.rgbsep[2])
cn=str.fparam[4]-str.fparam[5]
rgb=where(str[gd].fparam[1] gt loggpar.calloggmin and str[gd].fparam[1] lt loggpar.calloggmax and $
          str[gd].fparam[0] gt loggpar.calteffmin and str[gd].fparam[0] lt loggpar.calteffmax and $
          (str[gd].fparam[1] lt loggpar.rclim[0]  or str[gd].fparam[1] gt loggpar.rclim[1] or  $
           cn[gd] lt loggpar.cnsep[0]+str[gd].fparam[3]*loggpar.cnsep[1]+dt[gd]*loggpar.cnsep[2]),nrgb )
if  nrgb gt 0 then begin
  rgb=gd[rgb]
  corr=(loggpar.rgbfit[0]+loggpar.rgbfit[1]*logg[rgb]+loggpar.rgbfit[2]*mh[rgb])
  ; ramp to zero at calloggmax (not needed if we don't calibrate the dwarfs!)
  ;j=where(str[rgb].fparam[1] gt 3.2,nj)
  ;if nj gt 0 then corr[j]*=(loggpar.calloggmax-str[rgb[j]].fparam[1])/(loggpar.calloggmax-3.2)
  str[rgb].param[1]=str[rgb].fparam[1]-corr
  str[rgb].param_cov[1,1]=loggpar.rgbrms^2
  str[rgb].param_cov[1,1]=(exp(loggpar.rgberrpar[0]+loggpar.rgberrpar[1]*(str[rgb].fparam[0]-4500)+loggpar.rgberrpar[3]*str[rgb].fparam[3]+loggpar.rgberrpar[2]*(snrerr[rgb]-100.)))^2
endif
rc=where(str[gd].fparam[1] gt loggpar.calloggmin and str[gd].fparam[1] lt loggpar.calloggmax and $
         str[gd].fparam[0] gt loggpar.calteffmin and str[gd].fparam[0] lt loggpar.calteffmax and $
         str[gd].fparam[1] gt loggpar.rclim[0] and str[gd].fparam[1] lt loggpar.rclim[1] and  $
         cn[gd] gt (loggpar.cnsep[0]+str[gd].fparam[3]*loggpar.cnsep[1]+dt[gd]*loggpar.cnsep[2]),nrc)
if  nrc gt 0 then begin
  rc=gd[rc]
  corr=(loggpar.rcfit[0]+loggpar.rcfit[1]*logg[rc]+loggpar.rcfit[2]*mh[rc])
  str[rc].param[1]=str[rc].fparam[1]-corr
  str[rc].param_cov[1,1]=loggpar.rcrms^2
  str[rc].param_cov[1,1]=(exp(loggpar.rcerrpar[0]+loggpar.rcerrpar[1]*(str[rc].fparam[0]-4500)+loggpar.rcerrpar[3]*str[rc].fparam[3]+loggpar.rcerrpar[2]*(snrerr[rc]-100.)))^2
endif
tbd=where(logg[gd] lt loggpar.calloggmin or logg[gd] gt loggpar.calloggmax,nbd)
if nbd gt 0 then begin
  for ibd=0,nbd-1 do str[tbd[ibd]].paramflag[0]=str[tbd[ibd]].paramflag[0] or paramflagval('CALRANGE_BAD')
endif

; logg for DR13
if keyword_set(dr13) then begin
 ; use dt to distinguish RGB/RC
 dt=str.fparam[0]-(4468+(str.fparam[1]-2.5)/0.0018 - 382.5*str.fparam[3])
 ; clip logg and mh

 cn=str.fparam[4]-str.fparam[5]
 ; interpolate in intermediate regime first, inflate uncertainty. Some of these may get recalculated using C/N
 int=where(dt[gd] gt 0 and dt[gd] lt 100 and logg[gd] gt 2.385 and str[gd].fparam[1] lt loggpar.calloggmax and str[gd].fparam[1] gt loggpar.calloggmin,nint)
 if nint gt 0 then begin
  int=gd[int]
  corr=((loggpar.rcpar[0]+loggpar.rcpar[1]*logg[int]+loggpar.rcpar[2]*mh[int]+loggpar.rcpar[3]*logg[int]^2)*(dt[int]-0.)/100.+$
      (loggpar.rgbpar[0]+loggpar.rgbpar[1]*logg[int]+loggpar.rgbpar[2]*mh[int])*(100.-dt[int])/100.)
  ; ramp to 0 correction?
  j=where(str[int].fparam[1] gt 3.2,nj)
  if nj gt 0 then corr[j]*=(loggpar.calloggmax-str[int[j]].fparam[1])/(loggpar.calloggmax-3.2)
  str[int].param[1]=str[int].fparam[1]-corr
  str[int].param_cov[1,1]=loggpar.rms^2*sqrt(2.)
 endif
 ; RGB
 rgb=where(dt[gd] lt 0 or logg[gd] lt 2.385 or (dt[gd] lt 100 and cn[gd] lt -0.113-0.0043*dt[gd])  $
     and str[gd].fparam[1] lt loggpar.calloggmax and str[gd].fparam[1] gt loggpar.calloggmin,nrgb)
 if  nrgb gt 0 then begin
  rgb=gd[rgb]
  corr=(loggpar.rgbpar[0]+loggpar.rgbpar[1]*logg[rgb]+loggpar.rgbpar[2]*mh[rgb])
  j=where(str[rgb].fparam[1] gt 3.2,nj)
  if nj gt 0 then corr[j]*=(loggpar.calloggmax-str[rgb[j]].fparam[1])/(loggpar.calloggmax-3.2)
  str[rgb].param[1]=str[rgb].fparam[1]-corr
  str[rgb].param_cov[1,1]=loggpar.rms^2
 endif
 ; RC
 rc=where(logg[gd] gt 2.385 and (dt[gd] gt 100 or (dt[gd] ge 0 and cn[gd] gt -0.088-0.0018*dt[gd]))  $
      and str[gd].fparam[1] lt loggpar.calloggmax and str[gd].fparam[1] gt loggpar.calloggmin,nrc)
 if nrc gt 0 then begin
  rc=gd[rc]
  corr=(loggpar.rcpar[0]+loggpar.rcpar[1]*logg[rc]+loggpar.rcpar[2]*mh[rc]+loggpar.rcpar[3]*logg[rc]^2)
  j=where(str[rc].fparam[1] gt 3.2,nj)
  if nj gt 0 then corr[j]*=(loggpar.calloggmax-str[rc[j]].fparam[1])/(loggpar.calloggmax-3.2)
  str[rc].param[1]=str[rc].fparam[1]-corr
  str[rc].param_cov[1,1]=loggpar.rms^2
 endif
endif

; calibrate other parameters: [M/H], [C/M], [N/M], [alpha/M]
; need to create tmp structure so aspcap_elemfit can modify flags if outside calibration range

; giants
tmp=str[giants]
; [M/H]
i=where(strtrim(giantcal.elem,2) eq 'M')
i=i[0]
gd=where( tmp.fparam[0] gt giantcal[i].caltemin and tmp.fparam[0] lt giantcal[i].caltemax,ngd,comp=bd,ncomp=nbd) 
flag=tmp[gd].paramflag[3] 
tmp[gd].param[3]=tmp[gd].fparam[3]-aspcap_elemfit(tmp[gd],giantcal,'M',flag)
tmp[gd].paramflag[3]=flag
; "empirical" parameter errors: [M/H] 
tmp[gd].param_cov[3,3]=(exp(giantcal[i].errpar[0]+giantcal[i].errpar[1]*(tmp[gd].fparam[0]-4500)+giantcal[i].errpar[3]*tmp[gd].fparam[3]+giantcal[i].errpar[2]*(snrerr[giants[gd]]-100.)))^2
; [C/M], [N/M]
tmp.param[4]=str[giants].fparam[4]
tmp.param[5]=str[giants].fparam[5]
; [alpha/M]
i=where(strtrim(giantcal.elem,2) eq 'alpha')
i=i[0]
gd=where( tmp.fparam[0] gt giantcal[i].caltemin and tmp.fparam[0] lt giantcal[i].caltemax,ngd,comp=bd,ncomp=nbd) 
flag=tmp[gd].paramflag[6] 
tmp[gd].param[6]=tmp[gd].fparam[6]-aspcap_elemfit(tmp[gd],giantcal,'alpha',flag)
tmp[gd].paramflag[6]=flag
; "empirical" parameter errors: [alpha/M] 
tmp[gd].param_cov[6,6]=(exp(giantcal[i].errpar[0]+giantcal[i].errpar[1]*(tmp[gd].fparam[0]-4500)+giantcal[i].errpar[3]*tmp[gd].fparam[3]+giantcal[i].errpar[2]*(snrerr[giants[gd]]-100.)))^2
str[giants]=tmp

; dwarfs
tmp=str[dwarfs]
i=where(strtrim(dwarfcal.elem,2) eq 'M')
i=i[0]
gd=where( tmp.fparam[0] gt dwarfcal[i].caltemin and tmp.fparam[0] lt dwarfcal[i].caltemax,ngd,comp=bd,ncomp=nbd) 
flag=tmp[gd].paramflag[3] 
tmp[gd].param[3]=tmp[gd].fparam[3]-aspcap_elemfit(tmp[gd],dwarfcal,'M',flag)
tmp[gd].paramflag[3]=flag
tmp[gd].param_cov[3,3]=(exp(dwarfcal[i].errpar[0]+dwarfcal[i].errpar[1]*(tmp[gd].fparam[0]-4500)+dwarfcal[i].errpar[3]*tmp[gd].fparam[3]+dwarfcal[i].errpar[2]*(snrerr[dwarfs[gd]]-100.)))^2
; [C/M], [N/M]
tmp.param[4]=str[dwarfs].fparam[4]
tmp.param[5]=str[dwarfs].fparam[5]
; [alpha/M]
i=where(strtrim(dwarfcal.elem,2) eq 'alpha')
i=i[0]
gd=where( tmp.fparam[0] gt dwarfcal[i].caltemin and tmp.fparam[0] lt dwarfcal[i].caltemax,ngd,comp=bd,ncomp=nbd) 
flag=tmp[gd].paramflag[6] 
tmp[gd].param[6]=tmp[gd].fparam[6]-aspcap_elemfit(tmp[gd],dwarfcal,'alpha',flag)
tmp[gd].paramflag[6]=flag
tmp[gd].param_cov[6,6]=(exp(dwarfcal[i].errpar[0]+dwarfcal[i].errpar[1]*(tmp[gd].fparam[0]-4500)+dwarfcal[i].errpar[3]*tmp[gd].fparam[3]+dwarfcal[i].errpar[2]*(snrerr[dwarfs[gd]]-100.)))^2
str[dwarfs]=tmp

; tranfer vmicro, vsin/vmacro (they are already transferred if locked, but do it anyway)
str.param[2]=str.fparam[2]
str.param[7]=str.fparam[7]

; now calibrate elements
print,'calibrating elements...'
elems=aspcap_elems(tagnames,elemtoh)
if nelem gt 0 then begin
 for i=0,n_elements(elem)-1 do begin
  el=strtrim(elem[i],2)
  ii=where(strtrim(elems,2) eq el)
  ii=ii[0]
  etoh=elemtoh[ii]
  ii=where(strtrim(giantcal.elem,2) eq el)
  ii=ii[0]

  ; reset CALRANGE flags in case we are redoing
  str[gd].elemflag[i]=str[gd].elemflag[i] and (not paramflagval('CALRANGE_BAD'))
  str[gd].elemflag[i]=str[gd].elemflag[i] and (not paramflagval('CALRANGE_WARN'))

  if ngiants gt 0 then begin
    ; flag mismatches between parameter and window fits
    if el eq 'C' then begin
      j=where(abs(str[giants].felem[i]-str[giants].fparam[4]) gt 0.25,nj)
      if nj gt 0 then str[giants[j]].elemflag[i]=str[giants[j]].elemflag[i] or paramflagval('PARAM_MISMATCH_WARN')
      j=where(abs(str[giants].felem[i]-str[giants].fparam[4]) gt 0.5,nj)
      if nj gt 0 then str[giants[j]].elemflag[i]=str[giants[j]].elemflag[i] or paramflagval('PARAM_MISMATCH_BAD')
    endif
    if el eq 'N' then begin
      j=where(abs(str[giants].felem[i]-str[giants].fparam[5]) gt 0.35,nj)
      if nj gt 0 then str[giants[j]].elemflag[i]=str[giants[j]].elemflag[i] or paramflagval('PARAM_MISMATCH_WARN')
    endif
    if el eq 'Fe' then begin
      j=where(abs(str[giants].felem[i]-str[giants].fparam[3]) gt 0.25,nj)
      if nj gt 0 then str[giants[j]].elemflag[i]=str[giants[j]].elemflag[i] or paramflagval('PARAM_MISMATCH_WARN')
      j=where(abs(str[giants].felem[i]-str[giants].fparam[3]) gt 0.5,nj)
      if nj gt 0 then str[giants[j]].elemflag[i]=str[giants[j]].elemflag[i] or paramflagval('PARAM_MISMATCH_BAD')
    endif

    ; calibrate the abundances for giants that are not flagged as BAD for this element
    if giantcal[i].elemfit ge 0 then begin
      gd=where( (str[giants].elemflag[i] and badparamflag()) eq 0 and str[giants].fparam[0] gt giantcal[i].caltemin and str[giants].fparam[0] lt giantcal[i].caltemax,ngd,comp=bd,ncomp=nbd) 
    endif else begin
      ngd=0
      bd=indgen(n_elements(giants))
      nbd=n_elements(bd)
    endelse
    ; need to create tmp structure so aspcap_elemfit can modify flags
    if ngd gt 0 then begin
      tmp=str[giants[gd]]
      tmpclass=class[giants[gd]]
      flag=str[giants[gd]].elemflag[i]
      ; calibration fits are done in [X/M] space, so include fitting out [M/H] Teff trends for elements fit using M dimension (etoh=1),
      ;   so these elements will still have the resisual [M/H] Teff trends in them in the elem array
      ;tmp.elem[i]=str[giants[gd]].felem[i]-aspcap_elemfit(tmp,giantcal,el,flag)
      tmpelem=str[giants[gd]].felem[i]-aspcap_elemfit(tmp,giantcal,el,flag)
      tmp.elemflag[i]=flag
      ; populate x_h and x_m arrays
      ; if this is C or N and we are in dwarf grid, then parameter is already [X/H]
      if el eq 'C' or el eq 'CI' or el eq 'N' then begin
         jj=where(strpos(tmpclass,'GKd') ge 0  or strpos(tmpclass,'Fd') ge 0 or strpos(tmpclass,'Md') ge 0,nj,comp=jg,ncomp=njg) 
         ;if nj gt 0 then tmp[jj].x_m[i]=tmp[jj].elem[i]-tmp[jj].fparam[3]
         ;if njg gt 0 then tmp[jg].x_m[i]=tmp[jg].elem[i]
         if nj gt 0 then tmp[jj].x_m[i]=tmpelem[jj]-tmp[jj].fparam[3]
         if njg gt 0 then tmp[jg].x_m[i]=tmpelem[jg]
      endif else begin
         ;if etoh then tmp.x_m[i]=tmp.elem[i]-tmp.fparam[3] else tmp.x_m[i]=tmp.elem[i]
         if etoh then tmp.x_m[i]=tmpelem-tmp.fparam[3] else tmp.x_m[i]=tmpelem
      endelse
      ; for [X/H], calculate using all calibrated quantities
      tmp.x_h[i]=tmp.x_m[i]+tmp.param[3]
      str[giants[gd]]=tmp

      if x_m_exists then str[giants[gd]].x_m_err[i]= $
         max([[exp(giantcal[ii].errpar[0]+giantcal[ii].errpar[1]*(str[giants[gd]].fparam[0]-4500)+giantcal[ii].errpar[3]*str[giants[gd]].fparam[3]+giantcal[ii].errpar[2]*(snrerr[giants[gd]]-100.))],[str[giants[gd]].felem_err[i]]],dim=2)
      if x_m_exists then str[giants[gd]].x_h_err[i]= $
         max([[exp(giantcal[ii].errpar[0]+giantcal[ii].errpar[1]*(str[giants[gd]].fparam[0]-4500)+giantcal[ii].errpar[3]*str[giants[gd]].fparam[3]+giantcal[ii].errpar[2]*(snrerr[giants[gd]]-100.))],[str[giants[gd]].felem_err[i]]],dim=2)
    endif
    if nbd gt 0 then begin
      str[giants[bd]].elemflag[i]=str[giants[bd]].elemflag[i] or paramflagval('CALRANGE_BAD')
      ; calibrate the errors
    endif
  endif

  ; calibrate the abundances for dwarfs that are not flagged as BAD for this element, but adjust to match giant correction at 5500
  if ndwarfs gt 0 then begin
    if dwarfcal[i].elemfit ge 0 then begin
      gd=where( (str[dwarfs].elemflag[i] and badparamflag()) eq 0 and str[dwarfs].fparam[0] gt dwarfcal[i].caltemin and str[dwarfs].fparam[0] lt dwarfcal[i].caltemax,ngd,comp=bd,ncomp=nbd)
    endif else begin
      ngd=0
      bd=indgen(n_elements(dwarfs))
      nbd=n_elements(bd)
    endelse
    delta=0.
;    if giantcal.elempar.temax[i] lt 5500 then gtmp=giantcal.elempar.temax[i] else gtmp=5500.
;    if dwarfcal.elempar.temax[i] lt 5500 then dtmp=dwarfcal.elempar.temax[i] else dtmp=5500.
;    delta=dwarfcal.elempar.par[0,i]*(dtmp-4500.)-giantcal.elempar.par[0,i]*(gtmp-4500.)
    if ngd gt 0 then begin
      tmp=str[dwarfs[gd]]
      tmpclass=class[dwarfs[gd]]
      flag=str[dwarfs[gd]].elemflag[i]
      ;tmp.elem[i]=str[dwarfs[gd]].felem[i]-aspcap_elemfit(tmp,dwarfcal,el,flag)+delta
      tmpelem=str[dwarfs[gd]].felem[i]-aspcap_elemfit(tmp,dwarfcal,el,flag)+delta
      tmp.elemflag[i]=flag
      ; populate x_h and x_m arrays
      ; if this is C or N and we are in dwarf grid, then parameter is already [X/H]
      if el eq 'C' or el eq 'CI' or el eq 'N' then begin
         jj=where(strpos(tmpclass,'GKd') ge 0  or strpos(tmpclass,'Fd') ge 0 or strpos(tmpclass,'Md') ge 0,nj,comp=jg,ncomp=njg) 
         ;if nj gt 0 then tmp[jj].x_m[i]=tmp[jj].elem[i]-tmp[jj].fparam[3]
         ;if njg gt 0 then tmp[jg].x_m[i]=tmp[jg].elem[i]
         if nj gt 0 then tmp[jj].x_m[i]=tmpelem[jj]-tmp[jj].fparam[3]
         if njg gt 0 then tmp[jg].x_m[i]=tmpelem[jg]
      endif else begin
         ;if etoh then tmp.x_m[i]=tmp.elem[i]-tmp.fparam[3] else tmp.x_m[i]=tmp.elem[i]
         if etoh then tmp.x_m[i]=tmpelem-tmp.fparam[3] else tmp.x_m[i]=tmpelem
      endelse
      tmp.x_h[i]=tmp.x_m[i]+tmp.param[3]
      str[dwarfs[gd]]=tmp
      ; flag cool dwarfs
      cool = where(str[dwarfs[gd]].fparam[0] lt 3500, ncool)
      if ncool gt 0 then str[dwarfs[gd[cool]]].elemflag[i]=str[dwarfs[gd[cool]]].elemflag[i] or paramflagval('CALRANGE_WARN')
      ; calibrate the errors
      if x_m_exists then str[dwarfs[gd]].x_m_err[i]= $
         max([[exp(dwarfcal[ii].errpar[0]+dwarfcal[ii].errpar[1]*(str[dwarfs[gd]].fparam[0]-4500)+dwarfcal[ii].errpar[3]*str[dwarfs[gd]].fparam[3]+dwarfcal[ii].errpar[2]*(snrerr[dwarfs[gd]]-100.))],[str[dwarfs[gd]].felem_err[i]]],dim=2)
      if x_m_exists then str[dwarfs[gd]].x_h_err[i]= $
         max([[exp(dwarfcal[ii].errpar[0]+dwarfcal[ii].errpar[1]*(str[dwarfs[gd]].fparam[0]-4500)+dwarfcal[ii].errpar[3]*str[dwarfs[gd]].fparam[3]+dwarfcal[ii].errpar[2]*(snrerr[dwarfs[gd]]-100.))],[str[dwarfs[gd]].felem_err[i]]],dim=2)
    endif
    if nbd gt 0 then begin
      str[dwarfs[bd]].elemflag[i]=str[dwarfs[bd]].elemflag[i] or paramflagval('CALRANGE_BAD')
    endif
  endif
  

 endfor

endif


; DR10 correction was done star by star, with built-in relations, old flag names
if keyword_set(dr10) then begin
  ; get median chi2
  dt=100 & medt=3550+indgen(70)*dt
  medchi=fltarr(n_elements(medt)) & medn=intarr(n_elements(medt))
  for i=0,n_elements(medt)-1 do begin
    j=where(str.fparam[0] ge medt[i]-dt/2. and str.fparam[0] lt medt[i]+dt/2.,nj)
    medn[i]=nj
    if nj gt 1 then medchi[i]=median(param_chi2[j])
  endfor

; do flags star by star
 for i=0L,n_elements(str)-1 do begin

  teff=str[i].fparam[0]
  logg=str[i].fparam[1]
  feh=str[i].fparam[3]
  chi2=param_chi2[i]
  snr=snrstr[i].snr
  starflag=long(str[i].starflag)
  flag=long(str[i].aspcapflag)
  pflag=long(str[i].paramflag)
  jk0=(str[i].j-str[i].k)-1.5*str[i].ak_targ

  ; turn off parameter warning flags from GRIDEDGE_WARN
  flag=flag and not aspcapflagval('TEFF_WARN')
  flag=flag and not aspcapflagval('LOGG_WARN')
  flag=flag and not aspcapflagval('M_H_WARN')
  flag=flag and not aspcapflagval('ALPHA_M_WARN')
  flag=flag and not aspcapflagval('C_M_WARN')
  flag=flag and not aspcapflagval('N_M_WARN')

  ; zero out covariances
    str[i].param_cov=0.
    errscale=15
  ; Teff
   if logg lt 3.8 then begin
     if teff gt 4600  and teff lt 5500 then str[i].param[0]=str[i].fparam[0]-(0.3968*teff-1938.3) 
     if teff lt 4600 then str[i].param[0]=str[i].fparam[0]+113
     str[i].param_cov[0,0]=max([str[i].fparam_cov[0,0]*errscale^2,(83.8-39.8*feh)^2])
   endif
   if teff gt 5000 then begin
     flag=flag or aspcapflagval('TEFF_WARN')
     pflag[0]=pflag[0] or paramflagval('OTHER_WARN')
   endif 

  ; log g
  ; correct log g
  ;if logg lt 3.5 and teff gt 3800 and teff lt 5300 and feh gt -2.5 and feh lt 0.2 then begin
   if logg lt 3.8 and teff gt 3500 and teff lt 5300 and feh gt -2.5 and feh le 0.5 then begin
     ;str[i].param[1]=str[i].fparam[1]-(-0.1269*feh^2-5.59e-4*teff-8.128e-5*teff*feh+2.839)
     ;str[i].param[1]=str[i].fparam[1]-(0.06286*feh^2+2.941e-5*teff+0.1362)
     ;str[i].param[1]=str[i].fparam[1]-(-0.09736*feh-4.634e-6*feh*teff+0.2397)   ; "test4"
     str[i].param[1]=str[i].fparam[1]-(-0.1222*feh+0.2396)   ; "test_final4"
     str[i].param_cov[1,1]=max([str[i].fparam_cov[1,1]*errscale^2,0.2^2])
 
   endif else begin
     flag=flag or aspcapflagval('LOGG_WARN')
     pflag[1]=pflag[1] or paramflagval('CALRANGE_WARN')
     ; extrapolate for giants at "edge" values
     ;if logg lt 3.8 then begin 
     ;  tmpteff=max([3500,min([teff,5300])])
     ;  tmpfeh=max([-2.5,min([0.5,feh])])
      ; str[i].param[1]=str[i].fparam[1]-(-0.1269*tmpfeh^2-5.59e-4*tmpteff-8.128e-5*tmpteff*tmpfeh+2.839)
      ;  str[i].param[1]=str[i].fparam[1]-(0.06286*feh^2+2.941e-5*tmpteff+0.1362)
      ; str[i].param[1]=str[i].fparam[1]-(-0.09736*tmpfeh-4.634e-6*tmpfeh*tmpteff+0.2397)
     ;  str[i].param[1]=str[i].fparam[1]-(-0.122*tmpfeh+0.2367)
     ;endif
   endelse

  ; [M/H]
   ; correct [M/H]
   ;if teff le 5500 and logg lt 3.5 then begin
   if teff le 5500 and logg lt 3.8 then begin
     str[i].param[3]=str[i].fparam[3]-(-0.06199*feh^2+1.125e-4*teff-4.734e-5*teff*feh-0.544)
     str[i].param_cov[3,3]=max([str[i].fparam_cov[3,3]*errscale^2,(0.0548-0.0361*feh)^2])
   endif else begin
     flag=flag or aspcapflagval('METALS_WARN')
     pflag[3]=pflag[3] or paramflagval('CALRANGE_WARN')
     ; extrapolate for giants on the hot end
     ;if logg lt 3.8 and teff lt 6500 then begin
     ;  tmpteff=min([teff,5500])
     ;  str[i].param[3]=str[i].fparam[3]-(-0.06199*feh^2+1.125e-4*tmpteff-4.734e-5*tmpteff*feh-0.544)
     ;endif
   endelse
   if teff gt 6500 then begin
     flag=flag or aspcapflagval('METALS_BAD')
     pflag[3]=pflag[3] or paramflagval('OTHER_BAD')
     flag=flag or aspcapflagval('METALS_WARN')
     pflag[3]=pflag[3] or paramflagval('OTHER_WARN')
   endif

  ; [C/M]
   flag=flag or aspcapflagval('CFE_WARN')
   pflag[4]=pflag[4] or paramflagval('OTHER_WARN')
   str[i].param_cov[4,4]=max([str[i].fparam_cov[4,4]*errscale^2,(0.3)^2])

  ; [N/M]
   flag=flag or aspcapflagval('NFE_WARN')
   pflag[5]=pflag[5] or paramflagval('OTHER_WARN')
   str[i].param_cov[5,5]=max([str[i].fparam_cov[5,5]*errscale^2,(0.3)^2])

  ; [alpha/M]
   if logg lt 3.8 and teff lt 5300 then begin
     str[i].param[6]=str[i].fparam[6]
     str[i].param_cov[6,6]=max([str[i].fparam_cov[6,6]*errscale^2,(0.1)^2])
   endif
   if feh lt -0.5 or feh gt 0.2 or logg gt 3.5 or teff lt 4200 then begin
     flag=flag or aspcapflagval('ALPHAFE_WARN')
     pflag[6]=pflag[6] or paramflagval('OTHER_WARN')
   endif

  ; chi^2
  j=min(abs(teff-medt),it)
  chicrit=10+fltarr(n_elements(medt))
  chicrit=max([[medchi*2],[chicrit]],dim=2)
  if chi2 gt chicrit[it] then flag=flag or aspcapflagval('CHI2_WARN')

  chicrit=30+fltarr(n_elements(medt))
  chicrit=max([[medchi*2],[chicrit]],dim=2)
  if chi2 gt chicrit[it] then flag=flag or aspcapflagval('CHI2_BAD')

  ; broad lines
  ;if ((starflag and starflagval('SUSPECT_BROAD_LINES')) gt 0) then flag=flag or aspcapflagval('ROTATION_WARN')
  if str[i].rv_ccfwhm/str[i].rv_autofwhm gt 1.5 then flag=flag or aspcapflagval('ROTATION_WARN')
  if str[i].rv_ccfwhm/str[i].rv_autofwhm gt 2.0 then flag=flag or aspcapflagval('ROTATION_BAD')

  ; S/N 
  if snr lt 70 then flag=flag or aspcapflagval('SN_WARN')
  if snr lt 50 then flag=flag or aspcapflagval('SN_BAD')

  ; color=Teff
  if str[i].h lt 90 then begin
    if abs(teff-aspcap_colorte(jk0,feh)) gt 500 then flag=flag or aspcapflagval('COLORTE_WARN')
    if abs(teff-aspcap_colorte(jk0,feh)) gt 1000 then flag=flag or aspcapflagval('COLORTE_BAD')
  endif

  ; star level flag
  if ((flag and warnaspcapflag()) gt 0) then flag=flag or aspcapflagval('STAR_WARN')
  if ((flag and badaspcapflag()) gt 0) then flag=flag or aspcapflagval('STAR_BAD')
  ; if any parameter is near grid edge, set STAR_BAD
  for ipar=0,6 do if (pflag[ipar] and paramflagval('GRIDEDGE_BAD')) gt 0 then flag=flag or aspcapflagval('STAR_BAD')

  ; load new flags
  str[i].aspcapflag=flag
  str[i].paramflag=pflag

 endfor
endif  ; DR10

mask=0L
mask='FF'XL
help,mask
j=where((str.aspcapflag and mask) gt 0, nj)
print,'limit warning: ', nj
mask='FF00'XL 
help,mask
j=where((str.aspcapflag and mask) gt 0, nj)
print,'other warnings: ', nj
mask='FFFF0000'XL
help,mask
j=where((str.aspcapflag and mask) gt 0, nj)
print,'all bad: ', nj
j=where(str.aspcapflag  eq 0, nj)
print,'all good: ', nj

end
