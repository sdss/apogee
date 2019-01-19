function aspcap_loadferre,name,libpar,libfile,npar=npar,extra=extra,elemfit=elemfit

; routine to read output ferre files and load into structures to write
; into an output FITS table
rdlibhead,libfile,libhead0,libhead
nchip=n_elements(libhead)
for ichip=0,nchip-1 do begin
  w=libhead[ichip].wave[0]+indgen(libhead[ichip].npix)*libhead[ichip].wave[1]
  if ichip eq 0 then wave=w else wave=[wave,w]
endfor

; clean up empty error file (if it is empty!)
info=file_info(name+'.err.out')
if info.exists and (info.size eq 0) then file_delete,name+'.err.out'

; read the FERRE output files
aspcap_load,name+'.frd',data
flux=float(data)
aspcap_load,name+'.err',data
err=float(data)
aspcap_load,name+'.mdl',data
mdl=float(data)
if file_test(name+'.con') then aspcap_load,name+'.con',data else data*=0.
aspcap_load,name+'.spm',data
sz=size(data,/dim)
; if we only have one object, make 2D arrays
if n_elements(sz) eq 1 then begin
  flux=reform(flux,n_elements(flux),1)
  err=reform(err,n_elements(err),1)
  mdl=reform(mdl,n_elements(mdl),1)
  data=reform(data,n_elements(data),1)
  sz=size(data,/dim)
endif
if n_elements(data) eq n_elements(flux) then pseudo=float(data) else pseudo=0.*flux
if file_test(name+'.ipf') then readcol,name+'.ipf',inputstars,format='(a)'
;if file_test(name+'.vrd') then readcol,name+'.vrd',stars,format='(a)'
stars=strarr(sz[1])
for i=0,sz[1]-1 do stars[i]=data[0,i]
spm=float(data[1:sz[0]-1,*])

;get input parameters for FERRE
readcol,name+'.nml',delim='=',format='(a,a)',inputpar,inputval
inputpar=strupcase(strtrim(inputpar,2))

j=where(inputpar eq 'NDIM',nj)
if nj gt 0 then ndim=fix(inputval[j]) else ndim=0
j=where(inputpar eq 'NOV',nj)
if nj gt 0 then nov=fix(inputval[j]) else nov=0
j=where(inputpar eq 'INDV',nj)
if nj gt 0 then sindv=inputval[j] else sindv=''
; convert indv to integer array
s=strsplit(sindv,/ext)
indv=intarr(n_elements(s))
for i=0,n_elements(indv)-1 do indv[i]=s[i]
j=where(inputpar eq 'NRUNS',nj)
if nj gt 0 then nruns=fix(inputval[j]) else nruns=0
j=where(inputpar eq 'NTHREADS',nj)
if nj gt 0 then nthreads=fix(inputval[j]) else nthreads=0
j=where(inputpar eq 'INDINI',nj)
if nj gt 0 then sindini=inputval[j] else sindini=''
; convert indini to integer array
s=strsplit(sindini,/ext)
indini=intarr(n_elements(s))
for i=0,n_elements(indini)-1 do indini[i]=s[i]
j=where(inputpar eq 'INIT',nj)
if nj gt 0 then init=fix(inputval[j]) else init=0
j=where(inputpar eq 'PCAPROJECT',nj)
if nj gt 0 then pcaproject=fix(inputval[j]) else pcaproject=0
j=where(inputpar eq 'PCACHI',nj)
if nj gt 0 then pcachi=fix(inputval[j]) else pcachi=0
j=where(inputpar eq 'INTER',nj)
if nj gt 0 then inter=fix(inputval[j]) else inter=0
j=where(inputpar eq 'F_FORMAT',nj)
if nj gt 0 then f_format=fix(inputval[j]) else f_format=0
j=where(inputpar eq 'F_ACCESS',nj)
if nj gt 0 then f_access=fix(inputval[j]) else f_access=0
j=where(inputpar eq 'ALGOR',nj)
if nj gt 0 then algor=fix(inputval[j]) else algor=1
j=where(inputpar eq 'MONO',nj)
if nj gt 0 then mono=fix(inputval[j]) else mono=0
j=where(inputpar eq 'LSF',nj)
if nj gt 0 then lsf=fix(inputval[j]) else lsf=0
j=where(inputpar eq 'ERRBAR',nj)
if nj gt 0 then errbar=fix(inputval[j]) else errbar=0
j=where(inputpar eq 'FILTERFILE',nj)
if nj gt 0 then filterfile=inputval[j] else filterfile=''

nparam=libhead[0].n_of_dim

; define index order used by structure for named parameters
params=aspcap_params(tagnames,flagnames,npar=npar,extra=extra)
ntotparam=n_elements(params)
index=intarr(ntotparam)
for ipar=0,ntotparam-1 do index[ipar]=where(strtrim(libhead[0].label,2) eq params[ipar])
outparam=fltarr(ntotparam)
outcovar=fltarr(ntotparam,ntotparam)

; new format libpar files have lock information for missing dimensions
lock=libpar.plock.lock
;ilock=where(lock gt 0,nlock)
nlock=n_elements(libpar.plock)
for i=0,nlock-1 do begin
 ; if no PLOCKs, then we will have a single blank entry
 if strlen(libpar.plock[i].lock) gt 0 then begin
  j=where(strtrim(params,2) eq strtrim(libpar.plock[i].lock,2),nj)
  if nj eq 0 then if strtrim(libpar.plock[i].lock,2) eq 'alpha' then j=where(strtrim(params,2) eq 'O Mg Si S Ca Ti',nj)
  if nj eq 0 then stop,'Unknown locked parameter ...'
  if i eq 0 then ilock=j else ilock=[ilock,j]
 endif
endfor
;;if nlock gt 0 then paramflag[ilock]=0
;c=transpose(libpar.plock.c)

if nparam+nlock ne ntotparam then begin
  print,'halt: Number of parameters + number of locked parameters in .param file do not total ', ntotparam
  stop
endif

if tag_exist(libpar,'holefile') then begin
   holedata=mrdfits(getenv('APOGEE_SPECLIB')+'/atmos/'+libpar.holefile,0,holehdr) 
   hdim=intarr(sxpar(holehdr,'NAXIS'))
   for idim=0,sxpar(holehdr,'NAXIS')-1 do begin
     card='CTYPE'+string(format='(i1)',idim+1)
     hdim[idim]= where(libhead0.label eq strtrim(sxpar(holehdr,card)))
   endfor
endif else holedata=0
;sz=size(holedata)
;print,sz[0],nparam,ndim,nlock
;if tag_exist(libpar,'holefile') and sz[0] ne nparam then stop

for i=0,n_elements(stars)-1 do begin
  ; note that order in data file may not be the same as order in output file! Get correspondance
  istar=where(strtrim(inputstars,2) eq strtrim(stars[i],2))
  if n_elements(istar) gt 1 then begin
   print,'more than one ID match found in load ferre, not halting'
   istar=istar[0]
  endif

  ; ferre results
  param=spm[0:nparam-1,i]
  paramerr=spm[nparam:2*nparam-1,i]
  chi2=10^spm[2*nparam+2,i]
  covar=spm[2*nparam+3:2*nparam+3+nparam^2-1,i]
  covar=reform(covar,nparam,nparam)

  ; set limit flag as needed
  paramflag=lonarr(ntotparam)
  aspcapflag=0L
  if nlock gt 0 then paramflag[ilock] = paramflag[ilock] or paramflagval('PARAM_FIXED')
  bit=1
  outparam*=0.
  lparam=outparam
  outcovar*=0.

  for ipar=0,ntotparam-1 do begin
    if index[ipar] ge 0 then begin
      val=param[index[ipar]]
      if val lt libhead0.llimits[index[ipar]]+libhead0.steps[index[ipar]]/8. or $
         val gt libhead0.llimits[index[ipar]]+$
                libhead0.steps[index[ipar]]*(libhead0.n_p[index[ipar]]-1)-$
                libhead0.steps[index[ipar]]/8. then begin
        ; [N/M] edge only will give warn, not bad, to prevent star from being called bad
        ; [C/M], LOG10VDOP, LGVSINI lower edge will give warn, not bad, to prevent star from being called bad
        if ~keyword_set(elemfit) and params[ipar] eq 'N' then begin
          aspcapflag=aspcapflag OR aspcapflagval(flagnames[ipar]+'_WARN')
          paramflag[ipar]=paramflag[ipar] or paramflagval('GRIDEDGE_WARN')
        endif else if ~keyword_set(elemfit) and val lt libhead0.llimits[index[ipar]]+libhead0.steps[index[ipar]]/8. and $
                   (params[ipar] eq 'C' or params[ipar] eq 'LOG10VDOP' or params[ipar] eq 'LGVSINI') then begin
          aspcapflag=aspcapflag OR aspcapflagval(flagnames[ipar]+'_WARN')
          paramflag[ipar]=paramflag[ipar] or paramflagval('GRIDEDGE_WARN')
        endif else begin
          aspcapflag=aspcapflag OR aspcapflagval(flagnames[ipar]+'_BAD')
          paramflag[ipar]=paramflag[ipar] or paramflagval('GRIDEDGE_BAD')
        endelse
      endif
      if val lt libhead0.llimits[index[ipar]]+libhead0.steps[index[ipar]]/2. or $
         val gt libhead0.llimits[index[ipar]]+$
                libhead0.steps[index[ipar]]*(libhead0.n_p[index[ipar]]-1)-$
                libhead0.steps[index[ipar]]/2. then begin
        aspcapflag=aspcapflag OR aspcapflagval(flagnames[ipar]+'_WARN')
        paramflag[ipar]=paramflag[ipar] or paramflagval('GRIDEDGE_WARN')
      endif
      ; load full <NTOTPARAM>D parameter arrays
      if param[index[ipar]] lt -999 then begin
        outparam[ipar]=-9999. 
        paramflag[ipar]=paramflag[ipar] or paramflagval('FERRE_BAD')
        aspcapflag=aspcapflag OR aspcapflagval(flagnames[ipar]+'_BAD')
      endif else outparam[ipar]=param[index[ipar]]
      if covar[index[ipar],index[ipar]] lt 0 then paramflag[ipar]=paramflag[ipar] or paramflagval('FERRE_WARN')
      for jpar=0,ntotparam-1 do begin
        if index[jpar] ge 0 then outcovar[ipar,jpar]=covar[index[ipar],index[jpar]]
      endfor
    endif
    bit*=2
  endfor
  if tag_exist(libpar,'femin') then if outparam[3] lt libpar.femin then chi2+=1000.

  ; check for grid holes potentially affecting results and set flags if so
  if n_elements(holedata) gt 1 then begin
    gridloc=fix((param-libhead0.llimits)/libhead0.steps)
    ; if we are off grid, set to grid edge
    for ipar=0,n_elements(gridloc)-1 do begin
      ; for Carbon, the grid starts at -1.5, but the holefile starts at -1
      if libhead0.label[ipar] eq 'C' then gridloc[ipar]-=2
      if gridloc[ipar] lt 0 then gridloc[ipar] = 0
      if gridloc[ipar] gt libhead0.n_p[ipar]-1 then gridloc[ipar] = 0
    endfor
    if sxpar(holehdr,'NAXIS') eq 5 then begin
      sz=size(holedata,/dim)
      h=intarr(5)
      ;if we're not fitting this parameter, set to center of grid (not guaranteed to be right arbitrarily!)
      for ii=0,4 do if hdim[ii] ge 0 then h[ii]=gridloc[hdim[ii]] else h[ii]=sz[ii]/2+1
      hmax=max(holedata[max([0,h[0]-0]):min([sz[0]-1,h[0]+1]),$
                     max([0,h[1]-0]):min([sz[1]-1,h[1]+1]),$
                     max([0,h[2]-0]):min([sz[2]-1,h[2]+1]),$
                     max([0,h[3]-0]):min([sz[3]-1,h[3]+1]),$
                     max([0,h[4]-0]):min([sz[4]-1,h[4]+1])])
    ;  print, param
    ;  print, gridloc[hdim[0]],gridloc[hdim[1]],gridloc[hdim[2]], gridloc[hdim[3]],gridloc[hdim[4]]
    ;  print, h
    ;  stop
      if hmax gt 0 then begin
        ;print,param
        ;print,fix(gridloc)
        aspcapflag=aspcapflag OR aspcapflagval('ATMOS_HOLE_BAD')
        ;stop
      endif
      hmax=max(holedata[max([0,h[0]-1]):min([sz[0]-1,h[0]+2]),$
                     max([0,h[1]-1]):min([sz[1]-1,h[1]+2]),$
                     max([0,h[2]-1]):min([sz[2]-1,h[2]+2]),$
                     max([0,h[3]-1]):min([sz[3]-1,h[3]+2]),$
                     max([0,h[4]-1]):min([sz[4]-1,h[4]+2])])
      if hmax gt 0 then aspcapflag=aspcapflag OR aspcapflagval('ATMOS_HOLE_WARN')
    endif else stop,'naxis is not 5'
    ;status=get_hole_status(holedata,gridloc,[3],dist=dist,hole=bad)
    ;if status eq 2 then aspcapflag=aspcapflag OR aspcapflagval('SPEC_HOLE_BAD') else $
    ;if status eq 1 then aspcapflag=aspcapflag OR aspcapflagval('SPEC_HOLE_WARN')
    ;status=get_hole_status(holedata,gridloc,[2],dist=dist,hole=bad)
    ;if status eq 2 then aspcapflag=aspcapflag OR aspcapflagval('ATMOS_HOLE_BAD') else $
    ;if status eq 1 then aspcapflag=aspcapflag OR aspcapflagval('ATMOS_HOLE_WARN')
  endif

  ; set values of locked parameters if we can
  if nlock gt 0 then begin
   for jlock=0,nlock-1 do begin
    ; calculate value of locked parameter given input formula from PLOCK array
    jteff = where(params eq 'TEFF',nj)  
    if nj gt 0 then teff=outparam[jteff] else stop,'no TEFF!'
    jlogg = where(params eq 'LOGG',nj) 
    if nj gt 0 then logg=outparam[jlogg] else stop,'no LOGG!'
    jmh = where(params eq 'METALS',nj) 
    if nj gt 0 then mh=outparam[jmh] else stop,'no METALS!'
    lparam[ilock[jlock]] = libpar.plock[jlock].const + libpar.plock[jlock].te_coef*teff + $
       libpar.plock[jlock].logg_coef[0]*logg + libpar.plock[jlock].logg_coef[1]*logg^2 + libpar.plock[jlock].logg_coef[2]*logg^3 + $
       libpar.plock[jlock].mh_coef * mh
    outparam[ilock[jlock]]=lparam[ilock[jlock]]
    ;lparam[ilock[jlock]]=c[jlock,0]
    ;sz=size(c,/dim)
    ;;for j=0,ntotparam-1 do lparam[ilock]+=outparam[j]*c[ilock,j+1]
    ;for j=0,sz[1]-2 do lparam[ilock[jlock]]+=outparam[j]*c[jlock,j+1]
    ;; special case for log used in FERRE
    ;;loglock=where(lock eq 1 and lparam gt 0,nloglock)
    ;if strpos(libpar.plock[jlock].lock,'LOG') ge 0 or strpos(libpar.plock[jlock].lock,'LG') ge 0 then lparam[ilock[jlock]]=alog10(lparam[ilock[jlock]])
    ;; if lparam[2] gt 0 then lparam[2]=alog10(lparam[2])  ; vmicro is in log
    ;outparam[ilock[jlock]]=lparam[ilock[jlock]]
   endfor
  endif

  ; param array with systematic corrections applied
  sparam=outparam*0.-9999.00
  sparam_covar=outcovar*0.-9999.00
  if nlock gt 0 then sparam[ilock] = outparam[ilock]
  if nlock gt 0 then sparam_covar[ilock,*] = outcovar[ilock,*]
  if nlock gt 0 then sparam_covar[*,ilock] = outcovar[*,ilock]
  
  ; output structure: must include items that are unique to a class, even
  ; if they are the same for all members of a class, so when different classes
  ; are combined, info will be preserved. For items that are the same for
  ; all objects in all classes, put in HDU3
 
  allindv=intarr(ntotparam)-1 
  allindini=intarr(ntotparam)-1 
  gd=where(index ge 0,ngd)
  ;allindv[gd]=indv
  ;allindini[gd]=indini

  ; following to allow for fixed parameters that appear in library
  libindv=intarr(ngd)
  libindini=intarr(ngd)
  for j=0,n_elements(indv)-1 do begin
    libindv[indv[j]-1]=indv[j] 
    if n_elements(indini) gt 1 then libindini[indv[j]-1]=indini[j] 
    ;libindini[indini[j]-1]=indini[j] 
  endfor
  allindv[gd]=libindv
  if n_elements(indini) gt 1 then allindini[gd]=libindini

  gdmdl=where(mdl[*,i] gt 0.,nmdl)
  junk=where(err[gdmdl,i] gt 1.,nbd)
  a={apogee_id: stars[i], class: libpar.class, param: sparam, param_cov: sparam_covar, $
     fparam: outparam, fparam_cov: outcovar, paramflag: paramflag, $
     param_chi2: chi2[0], aspcapflag: aspcapflag, $
     ndim: ndim, nov: nov, indv: allindv, nruns: nruns, init: init, $
     indini: allindini, inter: inter,$
     nthreads: nthreads, pcaproject: pcaproject, pcachi: pcachi, $,
     f_format: f_format, f_access: f_access, algor: algor, mono: mono, $
     lsf: lsf, errbar: errbar, filterfile: filterfile, badfrac: float(nbd)/nmdl}

; fix for possible bug in library length, to bring to 7212 elements
  if tag_exist(libhead[0],'file_data19') then atomlinelist=file_basename(libhead[0].file_data19) else atomlinelist='unknown'
  if tag_exist(libhead[0],'file_data20') then moleclinelist=file_basename(libhead[0].file_data19) else moleclinelist='unknown'

;  if n_elements(flux[*,i]) eq 7214 or n_elements(flux[*,i]) eq 7213 then begin
;    sz=size(libhead.npix,/dim)
;    libhead.npix=reform([2920,2399,1893],sz[0],sz[1])
;    nwave=[reform(wave[0:5319]),reform(wave[5321:7212])]
;    npseudo=[reform(pseudo[0:5319,istar]),reform(pseudo[5321:7212,istar])]
;    nflux=[reform(flux[0:5319,istar]),reform(flux[5321:7212,istar])]
;    nerr=[reform(err[0:5319,istar]),reform(err[5321:7212,istar])]
;    nmdl=[reform(mdl[0:5319,i]),reform(mdl[5321:7212,i])]
;    b={spec: nflux, err: nerr, spec_bestfit: nmdl, pseudo_continuum: npseudo, grid: file_basename(libfile), $
;      atomlinelist: atomlinelist, moleclinelist: moleclinelist}
;  endif else b={spec: flux[*,istar], err: err[*,istar], spec_bestfit: mdl[*,i], pseudo_continuum: pseudo[*,istar], $
;    grid: file_basename(libfile), atomlinelist: atomlinelist,  moleclinelist: moleclinelist}

   ;if n_elements(flux[*,i]) ne 7214 then stop,'n_elements(flux) ne 7214!'

   b={spec: flux[*,istar], err: err[*,istar], spec_bestfit: mdl[*,i], $
      pseudo_continuum: pseudo[*,istar], grid: file_basename(libfile), $
      atomlinelist: atomlinelist,  moleclinelist: moleclinelist}

  if i eq 0 then all_a=a else all_a=[all_a,a]
  if i eq 0 then all_b=b else all_b=[all_b,b]
endfor

; final structure with parameters fixed for all objects
;if n_elements(wave) eq 7214 or n_elements(wave) eq 7213 then wave=nwave

c={npix: n_elements(wave), param_symbol: params, wave:wave, $
   wavemin: reform(libhead.wave[0,*]), $
   wavemax: reform(libhead.wave[0,*])+(libhead.npix-1)*reform(libhead.wave[1,*])}

; output structure
str={param: all_a, spec: all_b, lib: c}

; return structure
return,str

end
