function speclib_mkturbospec,teff,logg,mh,am,cm,nm,wrange=wrange,dw=dw,atmod=atmod,elem=elem,welem=welem,save=save,linelist=linelist,vmicro=vmicro,solarisotopes=solarisotopes,specdir=specdir,nskip=nskip,endskip=endskip,kurucz=kurucz,h2o=h2o,norun=norun,split=split,linedir=linedir,marcsdir=marcsdir

if not keyword_set(wrange) then wrange=[15100.,17000]
if not keyword_set(split) then split=200
if not keyword_set(dw) then dw=0.1
if keyword_set(elem) then suffix=elem else suffix=''
if not keyword_set(vmicro) then vmicro=2.0
if not keyword_set(endskip) then endskip=0
if not keyword_set(nskip) then nskip=0
if keyword_set(save) then stdout=' ' else stdout = ' >& /dev/null'
if ~keyword_set(marcsdir) then marcsdir='marcs/edvarsson/'

; linelist
if ~keyword_set(linelist) then linelist='turbospec.201312161124.new.vac'
if ~keyword_set(linedir) then linelistdir=getenv('APOGEE_SPECLIB')+'/linelists/' else linelistdir=linedir

; output file name
if keyword_set(kurucz) then atmos='k' else atmos='m'
if keyword_set(solarisotopes) then prefix=atmos+'s' else prefix=atmos+'g'
root=prefix+'t'+strtrim(string(nint(teff)),2)+'g'+cval(logg)+$ 
  'm'+cval(mh)+'a'+cval(am)+'c'+cval(cm)+'n'+cval(nm)+'v'+cval(vmicro)+suffix

; spherical or plane-parallel?
if logg le 3 then geo='s' else geo='p'

; with specdir= keyword, look for already existing fluxes
if keyword_set(specdir) then begin
  if geo eq 's' then subdir='SPH' else subdir='PP'
  modfile=geo+strtrim(string(nint(teff)),2)+'_g'+cval(logg,/turbo,/one)+'*'+'_z'+cval(mh,/turbo)+'_a'+cval(am,/turbo)+'_c'+cval(cm,/turbo)+'_n'+cval(nm,/turbo)+'_o'+cval(am,/turbo)
print,modfile
  files=file_search(getenv('APOGEE_SPECLIB')+'/'+specdir+'/'+subdir+'/'+modfile+'*_xit'+string(format='(f3.1)',vmicro)+'*.spec')
print,files
 specfile=files[0]
 ; file=getenv('APOGEE_SPECLIB')+'/'+specdir+'/'+subdir+'/'+files[0]+'_xit'+string(format='(f3.1)',vmicro)+'.spec'
  ;spawn,'pwd'
  print,specfile
  info=file_info(specfile)
  if info.size gt 0 then readcol,specfile,w,data,flux,count=count else count=0
  if count gt 0 then return,flux
  ;if count le 0 then return,0. else return,flux
endif

; atmosphere: note that all [N/M] is solar for atmospheres (since it doesn't have a big effect)

if keyword_set(kurucz) then begin
  geo = 'p'
  ; atmosphere: need to convert to MOOG format
  modfile='m'+cval(mh)+'c'+cval(cm)+'o'+cval(am)+'t'+strtrim(string(nint(teff)),2)+$
    'g'+string(format='(i2.2)',nint(logg*10.))+'v20'
  moddir=getenv('APOGEE_SPECLIB')+'/kurucz/'+'m'+cval(mh)+'c'+cval(cm)+'o'+cval(am)+'/'
  print,' searching: '+moddir+'a'+modfile+'*'
  files=file_search(moddir+'a'+modfile+'*',count=count)
  if count eq 0 then begin
    stop,'no atmosphere found: ',modfile
  endif else if count gt 1 then begin
    stop,'more than one atmosphere found',modfile
  endif
  script=getenv('SPECLIB_DIR')+'/scripts/maketurbospecmodel.awk'
  print,'converting atmosphere...'
  ; trim top 7 layers off Kurucz atmopshere for log(tau)>-6
  if nskip eq 0 then trim=0
  if nskip eq 1 then trim=7
  if nskip eq 2 then trim=15
  if nskip gt 2 then return,0.
  print,'trimming: ', trim
  spawn,['awk','-f ',script,'trim='+string(trim,format='(i2.2)'),files[0]],out,/noshell
  openw,lun,file_basename(files[0]),/get_lun
  for i=0,n_elements(out)-1 do printf,lun,out[i]
  free_lun,lun
endif else begin
  modfile=geo+strtrim(string(nint(teff)),2)+'_g'+cval(logg,/turbo,/one)+'*'+'_z'+cval(mh,/turbo)+'_a'+cval(am,/turbo)+'_c'+cval(cm,/turbo)+'_n'+cval(0.,/turbo)     ;+'_o'+cval(am,/turbo)
  ;moddir=getenv('APOGEE_SPECLIB')+'/marcs/MARCS_M_models/Mgrid/'
  moddir=getenv('APOGEE_SPECLIB')+'/'+marcsdir+'/M/'
  print,' searching: '+moddir+modfile+'*'
  files=file_search(moddir+modfile+'*',count=count)
  if count eq 0 then begin
    moddir=getenv('APOGEE_SPECLIB')+'/'+marcsdir+'/GK/'
    print,' searching: '+moddir+modfile+'*'
    files=file_search(moddir+modfile+'*',count=count)
  endif
  if count eq 0 then begin
    moddir=getenv('APOGEE_SPECLIB')+'/'+marcsdir+'/K/'
    print,' searching: '+moddir+modfile+'*'
    files=file_search(moddir+modfile+'*',count=count)
  endif
  if count eq 0 then begin
    stop,'no atmosphere found: ',modfile
  endif else if count gt 1 then begin
    stop,'more than one atmosphere found',modfile
  endif
  print,' found: ', files

  ; if we have nskip= keyword, trim layers from the model atmosphere
  if keyword_set(nskip) || keyword_set(endskip) then begin
   ; count header lines
   com="awk '$1=="+'"k"'+" {print NR}' "+files[0]+" | head -1 "
   spawn,com,nhead
   com='grep "Number of depth" ' +files[0]+' | awk '+"'"+'{printf("%d",$1-'+string(nskip)+')}'+"'"
   spawn,com,nlayers
   com="awk 'NR<="+nhead+"|| ($1>"+string(nskip)+" && $1<="+string(fix(nlayers)+nskip-endskip)+")' "+files[0]+' | sed '+"'"+'/Number of depth/c\'+nlayers+'  Number of depth points'+"'"
   print,'nskip: ', nskip
   spawn,com,modout
   file_delete,file_basename(files[0]),/allow
   openw,lun,file_basename(files[0]),/get_lun
   for i=0,n_elements(modout)-1 do printf,lun,modout[i]
   free_lun,lun
  endif else begin
    file_delete,file_basename(files[0]),/allow
    file_link,files[0],file_basename(files[0]),/allow
  endelse
endelse
if not keyword_set(atmod) then atmod=file_basename(files[0])

; run calculations in a subdirectory to make it easy to clean
dir=prefix+'m'+cval(mh)+'a'+cval(am)+'c'+cval(cm)+'n'+cval(nm)+'v'+cval(vmicro)+suffix
file_mkdir,dir
file_move,atmod,dir,/over
file_link,getenv('SPECLIB_DIR')+'/src/turbospec/DATA',dir+'/DATA',/allow

root=dir+'/'+root

; set up abundance array for minigrid if requested
if keyword_set(elem) then begin
  elemnum=speclib_elem(elem,abun=abun)
  eabun=abun-0.75+indgen(10)*0.25 
  nelem=1
  linelistdir=linelistdir+elem+'/'
endif else begin
  eabun=0.
  nelem=0
endelse

; welem only computes in windows, but it is easier/faster to compute the whole range with a linelist that only has lines in windows!
if ~keyword_set(welem) then welem=wrange
sz=size(welem)
if sz[0] eq 2 then nrange=sz[2] else nrange=1

; loop over element abundances if we want minigrid
spec=[]
for ielem=0,n_elements(eabun)-1 do begin

file=root+string(format='(i2.2)',ielem)
bsynfile=file+'bsyn'+string(format='(i2.2)',ielem)+'.inp'

; only compute opacities for a single nominal abundance
if ielem eq 0 then begin
  openw,lun,root+'_babsma.csh',/get_lun
  printf,lun,"#!/bin/csh -f"
  printf,lun,getenv('SPECLIB_DIR')+"/bin/babsma_lu " +stdout+" << EOF"
  printf,lun,"'LAMBDA_MIN:'  '"+string(min(welem)-dw,format='(f12.3)')+"'"
  printf,lun,"'LAMBDA_MAX:'  '"+string(max(welem)+dw,format='(f12.3)')+"'"
  printf,lun,"'LAMBDA_STEP:'  '"+string(dw)+"'"
  printf,lun,"'MODELINPUT:'  '"+atmod+"'"
  if keyword_set(kurucz) then printf,lun,"'MARCS-FILE:'  '.false.'"
  printf,lun,"'MODELOPAC:'  '"+file_basename(root)+'opac'+"'"
  printf,lun,"'METALLICITY:'  '"+string(mh)+"'"
  printf,lun,"'ALPHA/Fe:'  '"+string(am)+"'"
  printf,lun,"'HELIUM:'  '"+string(0.00)+"'"
  printf,lun,"'R-PROCESS:'  '"+string(0.00)+"'"
  printf,lun,"'S-PROCESS:'  '"+string(0.00)+"'"
  printf,lun,"'INDIVIDUAL ABUNDANCES:'  '"+string(3)+"'"
  printf,lun,6,' ',string(8.39+mh+cm)
  printf,lun,7,' ',string(7.78+mh+nm)
  printf,lun,8,' ',string(8.66+mh+am)
  ;if keyword_set(elem) then begin
  ;  printf,lun,elem,' ',string(eabun[ielem]+mh)
  ;endif
  if ~keyword_set(solarisotopes) then begin
    printf,lun,"'ISOTOPES:'  '"+string(2)+"'"
    ; Arcturus?
    ;printf,lun,6.012,0.888
    ;printf,lun,6.013,0.112
    ; adopt ratio of 12C/13C=15
    printf,lun,6.012,0.9375
    printf,lun,6.013,0.0625
  endif 
  printf,lun,"'XIFIX:'  '"+'T'+"'"
  printf,lun,string(vmicro)
  printf,lun,"EOF"
  if ~keyword_set(norun) then begin
    free_lun,lun
    file_chmod,root+'_babsma.csh','770'o
    cd,dir
    spawn,['time','./'+file_basename(root)+'_babsma.csh'],/noshell
    cd,'..'
  endif
endif

; split the wavelength ranges into chunks: if they are separated by <split, just do a single turbospec call
iw=0
ntot=0
for irange=0,nrange-1 do begin
  ; make sure welem is aligned with wrange
  welem[0,irange]=nint((welem[0,irange]-wrange[0])/dw)*dw+wrange[0]
  welem[1,irange]=nint((welem[1,irange]-wrange[0])/dw)*dw+wrange[0]
  i1=nint((welem[0,irange]-wrange[0])/dw)
  i2=nint((welem[1,irange]-wrange[0])/dw)
  ntot+=(i2-i1+1)
  if irange gt 0 then if welem[0,irange]-welem[1,irange-1] gt split and split gt 0 then iw=[iw,irange]
endfor
iw=[iw,nrange]
print,'ntot: ', ntot

allflux=fltarr(nint((wrange[1]-wrange[0])/dw)+1)
;for irange=0,nrange-1 do begin

if split lt 0 then nloop=0 else nloop=n_elements(iw)-2

for irange=0,nloop do begin
w1=welem[0,iw[irange]]
w2=welem[1,iw[irange+1]-1]

; create bsyn control file
openw,lbsyn,bsynfile,/get_lun
printf,lbsyn,"'LAMBDA_STEP:'  '"+string(dw)+"'"
if split lt 0 then begin
  for ii=0,nrange-1 do begin
    printf,lbsyn,"'LAMBDA_MIN:'  '"+string(welem[0,ii],format='(f12.3)')+"'"
    printf,lbsyn,"'LAMBDA_MAX:'  '"+string(welem[1,ii],format='(f12.3)')+"'"
  endfor
endif else begin
  printf,lbsyn,"'LAMBDA_MIN:'  '"+string(w1,format='(f12.3)')+"'"
  printf,lbsyn,"'LAMBDA_MAX:'  '"+string(w2,format='(f12.3)')+"'"
endelse
printf,lbsyn,"'INTENSITY/FLUX:'  '"+'Flux'+"'"
printf,lbsyn,"'COS(THETA):'  '"+string(1.00)+"'"
printf,lbsyn,"'ABFIND:'  '"+'.false'+"'"
printf,lbsyn,"'MODELINPUT:'  '"+atmod+"'"
if keyword_set(kurucz) then printf,lbsyn,"'MARCS-FILE:'  '.false.'"
printf,lbsyn,"'MODELOPAC:'  '"+file_basename(root)+'opac'+"'"
printf,lbsyn,"'RESULTFILE:'  '"+file_basename(file)+"'"
printf,lbsyn,"'METALLICITY:'  '"+string(mh)+"'"
printf,lbsyn,"'ALPHA/Fe:'  '"+string(am)+"'"
printf,lbsyn,"'HELIUM:'  '"+string(0.00)+"'"
printf,lbsyn,"'R-PROCESS:'  '"+string(0.00)+"'"
printf,lbsyn,"'S-PROCESS:'  '"+string(0.00)+"'"
printf,lbsyn,"'INDIVIDUAL ABUNDANCES:'  '"+string(3+nelem)+"'"
printf,lbsyn,6,' ',string(8.39+mh+cm)
printf,lbsyn,7,' ',string(7.78+mh+nm)
printf,lbsyn,8,' ',string(8.66+mh+am)
if keyword_set(elem) then begin
  printf,lbsyn,speclib_elem(elem),' ',string(eabun[ielem]+mh)
endif
if ~keyword_set(solarisotopes) then begin
  printf,lbsyn,"'ISOTOPES:'  '"+string(2)+"'"
  ; Arcturus?
  ;printf,lbsyn,6.012,0.888
  ;printf,lbsyn,6.013,0.112
  ; adopt ratio of 12C/13C=15
  printf,lbsyn,6.012,0.9375
  printf,lbsyn,6.013,0.0625
endif 
; do we need to add the H2O linelist?
if ~keyword_set(h2o) then begin
  if teff lt 4000 then begin
    if mh+am lt -1.5 or teff gt 3250 then h2o=1 else h2o=2
  endif else h2o=0
endif
nlists=3
; if no HI lines, don't use that list: it takes a while to read
n_HI=file_lines(linelistdir+'turbospec.'+strtrim(linelist,2)+'.Hlinedata')
if n_HI lt 3 then nlists-=1
; if we are using H2O, add that list
if h2o gt 0 then nlists+=1
printf,lbsyn,"'NFILES:'  '"+string(nlists)+"'"
if n_HI ge 3 then printf,lbsyn,linelistdir+'turbospec.'+strtrim(linelist,2)+'.Hlinedata'
printf,lbsyn,linelistdir+'turbospec.'+strtrim(linelist,2)+'.atoms'
printf,lbsyn,linelistdir+'turbospec.'+strtrim(linelist,2)+'.molec'
if h2o eq 1 then printf,lbsyn,linelistdir+'turbospec.h2o-BC8.5V'+'.molec'
if h2o eq 2 then printf,lbsyn,linelistdir+'turbospec.h2o-BC9.5V'+'.molec'
if geo eq 's' then printf,lbsyn,"'SPHERICAL:'  '"+'T'+"'" else $
  printf,lbsyn,"'SPHERICAL:'  '"+'F'+"'" 
printf,lbsyn,30
printf,lbsyn,300.00
printf,lbsyn,15
printf,lbsyn,1.3
free_lun,lbsyn

; control file, with special handling in case bsyn goes into infinite loop ...
openw,lun,root+'_bsyn.csh',/get_lun
;if keyword_set(elem) then begin
;printf,lun,"#!/bin/csh -f"
;printf,lun,getenv('SPECLIB_DIR')+"/bin/bsyn_lu " +stdout+" < "+file_basename(bsynfile)
;
;endif else begin
printf,lun,"#!/bin/csh -f"
printf,lun,getenv('SPECLIB_DIR')+"/bin/bsyn_lu " +stdout+" < "+file_basename(bsynfile)+' &'
printf,lun,'set bsynjob = $!'
printf,lun,"set ok = 0"
printf,lun,"set runtime = `ps -q $bsynjob -o cputime | tail -1 | awk -F: '{print ($1*3600)+($2*60)+$3}'`"
tmax=120*fix(.05/min([0.05,dw]))
if h2o eq 1 then tmax=240*fix(.05/min([0.05,dw]))
if h2o eq 2 then tmax=600*fix(.05/min([0.05,dw]))
printf,lun,'while ( $runtime < '+string(tmax)+' && $ok == 0 )'
printf,lun,'  usleep 200000'
printf,lun,"  set runtime = `ps -q $bsynjob -o cputime | tail -1 | awk -F: '{print ($1*3600)+($2*60)+$3}'`"
printf,lun,'  if ( `ps -p $bsynjob -o comm=` == "" ) then'
printf,lun,'    echo process done, exiting!'
printf,lun,'    set ok = 1'
;printf,lun,'    exit'
printf,lun,'  endif'
printf,lun,'end'
printf,lun,'if ( $ok == 0 ) then'
printf,lun,'  echo expired, killing job'
printf,lun,'  kill $bsynjob'
printf,lun,'endif'
;endelse

; run bsyn and read output
if ~keyword_set(norun) then begin
  free_lun,lun
  file_chmod,root+'_bsyn.csh','770'o
  cd,dir
  spawn,['time','./'+file_basename(root)+'_bsyn.csh'],/noshell
  cd,'..'
  print,'reading file ', file, file_test(file)
  if ~file_test(file) then return,0.
  if file_lines(file) eq 0 then return,0.
  readcol,file,w,data,flux,count=count
  print,count
  ;if count lt nint((welem[1,irange]-welem[0,irange])/dw+1) then stop,'wrong length'
  ;if count lt nint((welem[1,irange]-welem[0,irange])/dw+1) then return,0.

  for ii=iw[irange],iw[irange+1]-1 do begin
    i1=nint((welem[0,ii]-wrange[0])/dw)
    i2=nint((welem[1,ii]-wrange[0])/dw)
    j1=nint((welem[0,ii]-w1)/dw)
    ; in case turbospec cuts off last wavelength
    j2=j1+i2-i1+1
    ;print,j1,j2,n_elements(flux)
    if j2 gt n_elements(flux)-1 then j2-=1
    ;print,j1,j2,n_elements(flux)
    ;print,ii,welem[0,ii],i1,j1,j2
    allflux[i1:i1+j2-j1]=flux[j1:j2]
  endfor
  ;print,'split: ', split
endif

endfor  ; loop over wavelength ranges

spec=[[spec],[allflux]]
;help,spec
endfor  ; loop over minigrid abundances if requested

if ~keyword_set(save) then file_delete,dir,/recursive

return,spec

end
