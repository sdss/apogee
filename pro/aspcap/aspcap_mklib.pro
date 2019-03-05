;+
;  aspcap_mklib
;    makes an ASPCAP configuration directory and file, given a set of libraries, etc.
;    assumes library order is vmicro,[C/M],[N/M],[alpha/M],[M/H],logg,Teff unless carlos option is given
;-
pro aspcap_mklib,config,suffix=suffix,prefix=prefix,libdir=libdir,classes=classes,elem=elem,vpar=vpar,carlos=carlos,maskdir=maskdir,$
    femin=femin,rotpar=rotpar,opar=opar,vlock=vlock,vmlock=vmlock,outdir=outdir,nmlock=nmlock,cmlock=cmlock,apred=apred,inst=inst,$
    renorm=renorm,obscont=obscont,maxwind=maxwind

if n_elements(apred) gt 0 then apsetver,vers=apred,instrument=inst
;if ~keyword_set(elem) then elem=['C','Al','Ca','Fe','K','Mg','Mn','Na','Ni','N','O','Si','S','Ti','V']
;if ~keyword_set(maskdir) then maskdir='filters_05062014'
;if ~keyword_set(elem) then elem=['C','CI','Al','Ca','CI','Co','Cr','Cu','Fe','Ge','K','Mg','Mn','Na','Nd','Ni','N','O','P','Rb','Si','S','Ti','V','Y']
if ~keyword_set(elem) then elem=aspcap_elems()
;if ~keyword_set(maskdir) then stop,'Must specify a maskdir!'
if n_elements(npar) eq 0 then npar = 0.
if n_elements(opar) eq 0 then opar = 0.
if n_elements(cpar) eq 0 then cpar = 0.
if n_elements(apar) eq 0 then apar = 0.
if n_elements(rotpar) eq 0 then rotpar = [1., 0., 0., 0.]
if n_elements(vpar) eq 0 then vpar = [1.,0.,0.,0.]

dirs=getdir()
indir=getenv('APOGEE_DIR')+'/config/aspcap/'+config+'/'
if ~keyword_set(outdir) then outdir=getenv('APOGEE_DIR')+'/config/aspcap/'+config+'/'
file_mkdir,outdir
;if ~file_test(indir+'class-'+dirs.instrument+'.list') then stop,'no class-'+dirs.instrument+'.list file in '+indir
;readcol,indir+'class-'+dirs.instrument+'.list',classes,temin,temax,loggmin,loggmax,femin,femax,fibermin,fibermax,libs,holefile,format='(a,f,f,f,f,f,f,i,i,a,a)',comment='#'
;mask='global.mask'
aploadplan,indir+'class-'+dirs.instrument+'.par',p,str='CLASS'
mask=p.mask
maskdir=getenv('APOGEE_DIR')+'/data/windows/'+p.maskdir+'/'
if tag_exist(p,'obscont') then obscont=p.obscont else obscont=0
classes=p.class.class
nclasses=n_elements(classes)
libs=p.class.lib
temin=p.class.temin
temax=p.class.temax
loggmin=p.class.loggmin
loggmax=p.class.loggmax
femin=p.class.femin
femax=p.class.femax
fibermin=p.class.fibermin
fibermax=p.class.fibermax
holefile=p.class.holefile
vmacrofit=p.class.vmacrofit
rotpar=p.class.vmacro
inter=p.class.inter
print,'key: ', keyword_set(renorm)
print,'p: ', p.class.renorm
if keyword_set(renorm) then renorm=intarr(nclasses)+renorm else renorm=p.class.renorm 

nclasses=n_elements(classes)
cfit=intarr(nclasses)-1
nfit=intarr(nclasses)-1
afit=intarr(nclasses)-1
efit=intarr(nclasses)-1
havecoarse=0
for i=0,nclasses-1 do begin
  class=strtrim(classes[i],2)
  rdlibhead,getenv('APOGEE_SPECLIB')+'/synth/'+libs[i],libhead0
  for ipar=0,n_elements(libhead0.n_p)-1 do begin
    if strtrim(libhead0.label[ipar],2) eq 'C' then cfit[i]=ipar+1
    if strtrim(libhead0.label[ipar],2) eq 'N' then nfit[i]=ipar+1
    if strtrim(libhead0.label[ipar],2) eq 'O Mg Si S Ca Ti' then afit[i]=ipar+1
    if strtrim(libhead0.label[ipar],2) eq 'METALS' then efit[i]=ipar+1
  endfor
  if cfit[i] lt 0 then cfit[i] = efit[i]
  if nfit[i] lt 0 then nfit[i] = efit[i]
  if afit[i] lt 0 then afit[i] = efit[i]
  openw,out,outdir+class+'.par',/get_lun
  printf,out,'class '+class
  printf,out,'lib '+libs[i]
  ;if n_elements(holefile) gt 0 then printf,out,'holefile '+holefile[i]
  if strtrim(holefile[i],2) ne 'None' then printf,out,'holefile '+holefile[i]
  printf,out,'nov '+string(libhead0.n_of_dim)

  if tag_exist(p.class,'cmlock') then if p.class[i].cmlock gt -9 then cmlock=p.class[i].cmlock else undefine,cmlock
  if tag_exist(p.class,'nmlock') then if p.class[i].nmlock gt -9 then nmlock=p.class[i].nmlock else undefine,nmlock
  if tag_exist(p.class,'vlock') then if p.class[i].vlock gt -9 then vlock=p.class[i].vlock else undefine,vlock
  if ~keyword_set(maxwind) then if tag_exist(p.class,'maxwind') then maxwind=p.class[i].maxwind else maxwind=0

  ; if this is a coarse grid, remove C, N, and VMICRO dimensions, and adjust INDV and INDINI accordingly
  initpar=fltarr(libhead0.n_of_dim)
  if strpos(class,'coarse') ge 0 then begin
    havecoarse=1
    indv=[]
    indini=[]
    for ipar=0,n_elements(libhead0.n_p)-1 do begin
      if strtrim(libhead0.label[ipar],2) eq 'LOG10VDOP' and n_elements(vlock) gt 0 then begin
        initpar[ipar]=alog10(vlock)
      endif else if strtrim(libhead0.label[ipar],2) eq 'LGVSINI' and n_elements(vmlock) gt 0 and n_elements(libhead0.n_p) eq 7 then begin
        initpar[ipar]=alog10(vmlock)
      endif else if strtrim(libhead0.label[ipar],2) eq 'C' and n_elements(cmlock) gt 0 then begin
        initpar[ipar]=cmlock
      endif else if strtrim(libhead0.label[ipar],2) eq 'N' and n_elements(nmlock) gt 0 then begin
        initpar[ipar]=nmlock
      endif else begin
        indv=[indv,ipar]
        tmp=1
        ;if strtrim(libhead0.label[ipar],2) eq 'TEFF' then tmp=3
        ;if strtrim(libhead0.label[ipar],2) eq 'LOGG' then tmp=2
        ;if strtrim(libhead0.label[ipar],2) eq 'METALS' then tmp=2
        indini=[indini,tmp]
      endelse
    endfor
    printf,out,'indv '+string(indv+1,format='(9i2)')
    printf,out,'indini '+string(indini,format='(9i3)')
    printf,out,'initpar '+string(initpar,format='(9f12.4)')
    printf,out,'coarse 1 '
    printf,out,'noplot 1 '
  endif else begin
    params=aspcap_params()
    ;printf,out,'indv '+string(indgen(libhead0.n_of_dim)+1,format='(7i2)')
    pmin=[]
    pmax=[]
    indv=[]
    indini=[]
    for ipar=0,n_elements(libhead0.n_p)-1 do begin
      ; if we have vlock, then don't include vmicro
      amin=-99999 & amax=99999
      if strtrim(libhead0.label[ipar],2) eq 'LOG10VDOP' and n_elements(vlock) gt 0 then begin
        initpar[ipar]=alog10(vlock)
      endif else if strtrim(libhead0.label[ipar],2) eq 'LGVSINI' and n_elements(vmlock) gt 0 and n_elements(libhead0.n_p) eq 7 then begin
        initpar[ipar]=alog10(vmlock)
      endif else if strtrim(libhead0.label[ipar],2) eq 'C' and n_elements(cmlock) gt 0 then begin
        initpar[ipar]=cmlock
      endif else if strtrim(libhead0.label[ipar],2) eq 'N' and n_elements(nmlock) gt 0 then begin
        initpar[ipar]=nmlock
      endif else begin
        indv=[indv,ipar]
        ;tmp=-1
        ;if strtrim(libhead0.label[ipar],2) eq 'TEFF' then tmp=3
        ;if strtrim(libhead0.label[ipar],2) eq 'LOGG' then tmp=2
        ;if strtrim(libhead0.label[ipar],2) eq 'METALS' then tmp=2
        ;if strtrim(libhead0.label[ipar],2) eq 'C' then tmp=2
        ;if strtrim(libhead0.label[ipar],2) eq 'N' then tmp=2
        ;if strtrim(libhead0.label[ipar],2) eq 'O Mg Si S Ca Ti' then tmp=2
        j=where(params eq strtrim(libhead0.label[ipar],2))
        tmp=p.class[i].indini[j]
        indini=[indini,tmp]
        if strtrim(libhead0.label[ipar],2) eq 'TEFF' then begin
          amin=temin[i] & amax=temax[i]
        endif
        if strtrim(libhead0.label[ipar],2) eq 'LOGG' then begin
          amin=loggmin[i] & amax=loggmax[i]
        endif
        if strtrim(libhead0.label[ipar],2) eq 'METALS' then begin
          amin=femin[i] & amax=femax[i]
        endif
      endelse
      pmin=[pmin,amin]
      pmax=[pmax,amax]
    endfor
    printf,out,'indv '+string(indv+1,format='(9i2)')
    printf,out,'init  1 '
    printf,out,'indini '+string(indini,format='(9i3)')
    printf,out,'initpar '+string(initpar,format='(9f12.4)')
    printf,out,'pmin '+string(pmin,format='(9f12.2)')
    printf,out,'pmax '+string(pmax,format='(9f12.2)')
    printf,out,'fibermin '+string(fibermin[i],format='(i12)')
    printf,out,'fibermax '+string(fibermax[i],format='(i12)')
    if havecoarse eq 0 then printf,out,'coarse  -1 '
  endelse
  printf,out,'inter '+string(inter[i],format='(i12)')
  printf,out,'renorm '+string(renorm[i],format='(i12)')
  printf,out,'obscont '+string(obscont,format='(f8.2)')
  if file_test(maskdir+'/'+mask) then begin
    printf,out,'mask  '+mask
    file_delete,outdir+mask,/allow
    file_copy,maskdir+mask,outdir+mask
  endif else stop, 'No global mask, OK?',maskdir+mask
  ; information about locked parameters. If these are changed, must be in conjunction with changes in aspcap_loadferre
  printf,out,'typedef struct {'
  printf,out,'  char lock[16];'
  ;printf,out,'  float c[8];'
  printf,out,'  float const;'
  printf,out,'  float te_coef;'
  printf,out,'  float logg_coef[3];'
  printf,out,'  float mh_coef;'
  printf,out,'} PLOCK;'
  j=where(strtrim(libhead0.label,2) eq 'LOG10VDOP',nj)
  if nj eq 0 then begin
    if n_elements(vpar) eq 4 then $
    printf,out,'PLOCK LOG10VDOP '+string(format='(f8.4,f8.4," {",3f8.4," } ",f8.4)',vpar[0],0.,vpar[1:3],0.)+'      #  vmicro' $
    else if n_elements(vpar) eq 5 then $
    printf,out,'PLOCK LOG10VDOP '+string(format='(f8.4,f8.4," {",3f8.4," } ",f8.4)',vpar[0],0.,vpar[1:3],vpar[4])+'      #  vmicro' 
  endif
  j=where(strtrim(libhead0.label,2) eq 'LGVSINI',nj)
  if nj eq 0 then  $
    printf,out,'PLOCK LGVSINI '+string(format='(f8.4,f8.4," {",3f8.4," } ",f8.4)',rotpar[0],rotpar[1],[rotpar[2],0.,0.],rotpar[3])+'      #  rotation/vmacro' 
  j=where(strtrim(libhead0.label,2) eq 'O Mg Si S Ca Ti',nj)
  if nj eq 0 then  $
    printf,out,'PLOCK alpha '+string(format='(f8.4,f8.4," {",3f8.4," } ",f8.4)',apar[0],0.,[0.,0.,0.],0.)+'      #  alpha'
  j=where(strtrim(libhead0.label,2) eq 'O',nj)
  if nj eq 0 then  $
    printf,out,'PLOCK O '+string(format='(f8.4,f8.4," {",3f8.4," } ",f8.4)',opar[0],0.,[0.,0.,0.],0.)+'      #  oxygen'
  j=where(strtrim(libhead0.label,2) eq 'C',nj)
  if nj eq 0 then  $
    printf,out,'PLOCK C '+string(format='(f8.4,f8.4," {",3f8.4," } ",f8.4)',cpar[0],0.,[0.,0.,0.],0.)+'      #  carbon' 
  ;else if n_elements(cmlock) gt 0 then $
  ;  printf,out,'PLOCK C '+string(format='(f8.4,f8.4," {",3f8.4," } ",f8.4)',cmlock,0.,[0.,0.,0.],0.)+'      #  carbon'
  j=where(strtrim(libhead0.label,2) eq 'N',nj)
  if nj eq 0 then  $
    printf,out,'PLOCK N '+string(format='(f8.4,f8.4," {",3f8.4," } ",f8.4)',npar[0],0.,[0.,0.,0.],0.)+'      #  nitrogen' 
  ;else if n_elements(nmlock) gt 0 then $
  ;  printf,out,'PLOCK N '+string(format='(f8.4,f8.4," {",3f8.4," } ",f8.4)',nmlock,0.,[0.,0.,0.],0.)+'      #  nitrogen'
  free_lun,out 
endfor

;Now do the individual elements

openw,list,outdir+'elem.list',/get_lun
printf,list,'#'
printf,list,'#'
for i=0,n_elements(elem)-1 do begin
 el=strtrim(elem[i])
 if file_test(maskdir+'/'+el+'.filt') then begin
  printf,list,el
  file_delete,outdir+el+'.mask',/allow
  file_copy,maskdir+'/'+el+'.filt',outdir+el+'.mask'
  fit=efit
  if el eq 'C' or el eq 'CI' then fit=cfit 
  if el eq 'N' then fit=nfit 
  if el eq 'Ca' or el eq 'Mg' or el eq 'O' or el eq 'Si' or el eq 'S' or el eq 'Ti' or el eq 'TiII' then fit=afit

  ; old style puts everything in header with different columns for different grids
  ;openw,out,outdir+el+'.elem.par',/get_lun
  ;printf,out,'elem '+el
  ;line=''
  ;for j=0,n_elements(classes)-1 do line=line+classes[j]+' '
  ;printf,out,'class '+line
  ;line=''
  ;for j=0,n_elements(classes)-1 do line=line+libs[j]+' '
  ;printf,out,'lib '+line
  ;line=''
  ;for j=0,n_elements(classes)-1 do line=line+string(1)+' '
  ;printf,out,'nov '+line
  ;line=''
  ;for j=0,n_elements(classes)-1 do line=line+string(fit[j])+' '
  ;printf,out,'indv '+line
  ;line=''
  ;for j=0,n_elements(classes)-1 do line=line+string(1)+' '
  ;printf,out,'indini '+line
  ;line=''
  ;for j=0,n_elements(classes)-1 do line=line+string(inter[j])+' '
  ;printf,out,'inter '+line
  ;line=''
  ;for j=0,n_elements(classes)-1 do line=line+string(renorm[j])+' '
  ;printf,out,'renorm '+line
  ;line=''
  ;for j=0,n_elements(classes)-1 do line=line+el+' '
  ;printf,out,'mask '+line
  ;free_lun,out

  ; new style uses an INFO structure with one line per library
  openw,new,outdir+el+'.elem.par',/get_lun
  printf,new,'elem '+el
  printf,new,' typedef struct {'
  printf,new,'char class[16];'
  printf,new,'char libs[132];'
  printf,new,'int nov;'
  printf,new,'int indv;'
  printf,new,'int ttie[3];'
  printf,new,'int indini;'
  printf,new,'int inter;'
  printf,new,'int renorm;'
  printf,new,'char mask[8];'
  printf,new,'} INFO;'
  print,el
  for j=0,n_elements(classes) -1 do begin
    ttie=intarr(3)-1
    if fit[j] eq efit[j] then begin
      if cfit[j] ne efit[j] then ttie[0]=cfit[j]
      if nfit[j] ne efit[j] then ttie[1]=nfit[j]
      if afit[j] ne efit[j] then ttie[2]=afit[j]
    endif
    printf,new,'INFO ',classes[j],libs[j],1,fit[j],' { ',ttie,' } ',1,inter[j],renorm[j],el,format='(a,a,1x,a,i3,i3,a,3i3,a,i,i,i,1x,a)'
    print,j,fit[j],efit[j],cfit[j],nfit[j],afit[j],ttie
  endfor

  free_lun,new
 endif
endfor
free_lun,list
if keyword_set(maxwind) then  begin
  readcol,outdir+'/elem.list',format='(a)',els,stringskip='#',/silent
  for i=0,n_elements(els)-1 do n=filtsplit(els[i],maxwind=maxwind,outdir=outdir)
endif

end
