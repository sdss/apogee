pro aspcap_mkplan,files,apvisit=apvisit,single=single,aspcap_vers=aspcap_vers,apogee_vers=apogee_vers,queue=queue,apred_vers=apred_vers,apstar_vers=apstar_vers,aspcap_config=aspcap_config,results_vers=results_vers,nstars=nstars,ncpus=ncpus,noplot=noplot,noelem=noelem,commiss=commiss,nored=nored,visits=visits,caldir=caldir,npar=npar,renorm=renorm,maxwind=maxwind

; this routine makes the input parameter files for running ASPCAP

apsetver,vers=apred_vers
;if not keyword_set(aspcap_vers) then aspcap_vers = aspcap_version()
;if not keyword_set(aspcap_config) then aspcap_config='v6'
dirs=getdir(apogee_dir,calib_dir,spectro_dir)
;if keyword_set(apvisit) then  begin
;  outdir=apogee_dir+'/'+apred_vers+'/plates/'+'/'+aspcap_vers
;endif else begin
;  outdir=apogee_dir+'/'+apred_vers+'/'+apstar_vers+'/'+aspcap_vers
;endelse

outdir=getenv('APOGEE_ASPCAP')+'/'+apred_vers+'/'+aspcap_vers

if keyword_set(results_vers) then outdir=outdir+'/'+results_vers
print,'outdir: ', outdir
file_mkdir,outdir

if keyword_set(commiss) then aspcapstar='aspcapStarC' else aspcapstar='aspcapStar'
if keyword_set(single) then begin
  if keyword_set(visit) then openw,out,/get_lun,outdir+'aspcapPlate.par' else $
     openw,out,/get_lun,outdir+'/'+aspcapstar+'.par' 
endif

if not keyword_set(queue) then queue=0
qname='apogee'
qgroup='apogee'
;indini=[2,1,1,1,3,2]

iplate=0
if not keyword_set(apvisit) then begin
  for i=0,n_elements(files)-1 do begin
    dirname=file_basename(file_dirname(files[i]))
    if strpos(dirname,'apo1m') ge 0 then telescope='apo1m' $
    else if strpos(dirname,'apo25m') ge 0 then telescope='apo25m' $
    else if strpos(dirname,'lco25m') ge 0 then telescope='lco25m' $
    else stop,'No telescope name in directory name!'

    apsetver,vers=apred_vers,telescope=telescope
    dirs=getdir()
    field=file_basename(files[i])
    print,field
    file_mkdir,outdir+'/'+dirname+'/plan'
    if not keyword_set(single) then openw,out,/get_lun,outdir+'/'+dirname+'/plan/'+aspcapstar+'-'+field+'.par'
    if not keyword_set(single) or iplate eq 0 then begin
      ;printf,out,'idlwrap_version  ',idlwrap_version()
      ;printf,out,'speclib_version  ',speclib_version()
      ;printf,out,'ferre_version  ',ferre_version()
      printf,out,'apogee_version  ',getenv('APOGEE_VER')
      printf,out,'apvisit    0'
      printf,out,'apred_vers ''',apred_vers,'''
      printf,out,'telescope ''',telescope,'''
      printf,out,'instrument ''',dirs.instrument,'''
      printf,out,'apstar_vers ''',apstar_vers,'''
      printf,out,'aspcap_vers ''',aspcap_vers,'''
      if keyword_set(results_vers) then printf,out,'results_vers ''',results_vers,'''
      printf,out,'aspcap_config ''',aspcap_config,'''
      if keyword_set(ncpus) then printf,out,'ncpus       ',ncpus
      printf,out,'queue       ',queue
      printf,out,'qname       ''',qname,'''
      printf,out,'qgroup      ''',qgroup,'''
      ;printf,out,'indini      ',indini,format='(a,6i3)'
      if keyword_set(nstars) then printf,out,'nstars       ',nstars
      if keyword_set(visits) then printf,out,'visits       ',visits
      if keyword_set(noplot) then printf,out,'noplot       ',noplot
      if keyword_set(noelem) then printf,out,'noelem       ',noelem
      if keyword_set(caldir) then printf,out,'caldir       ''',caldir,'''
      if keyword_set(commiss) then printf,out,'commiss      ',commiss
      if keyword_set(nored) then printf,out,'nored      ',nored
      if keyword_set(npar) then printf,out,'npar      ',npar
      if keyword_set(renorm) then printf,out,'renorm      ',renorm
      if keyword_set(maxwind) then printf,out,'maxwind      ',maxwind
      printf,out,'typedef struct {'
      printf,out,'  char field[24];'
      printf,out,'  char outfield[24];'
      printf,out,'} ASPCAP;'
    endif
    printf,out,'ASPCAP '+dirname+'/'+field+' '+dirname+'/'+field
    if not keyword_set(single) then free_lun,out
    iplate+=1
  endfor
endif else begin
  for i=0,n_elements(files)-1 do begin
    plate=file_basename(files[i])
    print,plate
    aploadplan,files[i],planstr
    if not keyword_set(single) then openw,out,/get_lun,outdir+'/aspcapPlate-'+string(format='(i4.4)',planstr.plateid)+'-'+string(format='(i5.5)',planstr.mjd)+'.par'
    if not keyword_set(single) or iplate eq 0 then begin
      printf,out,'apvisit    1'
      printf,out,'apogee_vers ''',apogee_vers,'''
      printf,out,'aspcap_vers ''',aspcap_vers,'''
      printf,out,'telescope ''',telescope,'''
      printf,out,'instrument ''',dirs.instrument,'''
      if keyword_set(ncpus) then printf,out,'ncpus       ',ncpus
      printf,out,'queue       ',queue
      printf,out,'qname       ''',qname,'''
      printf,out,'qgroup      ''',qgroup,'''
      if keyword_set(noplot) then printf,out,'noplot       ',noplot
      if keyword_set(caldir) then printf,out,'caldir       ''',caldir,'''
      printf,out,'typedef struct {'
      printf,out,'  int plateid;'
      printf,out,'  int mjd;'
      printf,out,'} ASPCAP;'
    endif
    printf,out,'ASPCAP '+string(planstr.plateid)+' '+string(planstr.mjd)
    if not keyword_set(single) then free_lun,out
    iplate+=1
  endfor
endelse
if keyword_set(single) then free_lun,out
end

