; doaspcap is master routine for running ASPCAP pipline: pre-processing, FERRE, and post-processing

pro doaspcap,planfile,mjd=mjd,nruns=nruns,queue=queue,clobber=clobber,old=old,ncpus=ncpus,errbar=errbar,renorm=renorm,obscont=obscont,aspcap_vers=aspcap_vers,results_vers=results_vers,nstars=nstars,hmask=hmask,maskfile=maskfile,obspixmask=obspixmask,pixmask=pixmask,highbad=highbad,skyerr=skyerr,starlist=starlist,lowbad=lowbad,higherr=higherr,conthighbad=conthighbad,contlowbad=contlowbad,conthigherr=conthigherr,qaspcap=qaspcap,redux_root=redux_root,red_vers=red_vers,noplot=noplot,vacuum=vacuum,commiss=commiss,nored=nored,visits=visits,aspcap_config=aspcap_config,fits=fits,symlink=symlink,doelemplot=doelemplot,noelem=noelem,maxwind=maxwind,caldir=caldir,persist=persist,npar=npar,nelem=nelem,mask_telluric=mask_telluric,altmaskdir=altmaskdir,testmjd=testmjd,notie=notie,no_version_check=no_version_check

TIC
; read plan file and required fields: apvisit and plateid/mjd or field
; get appropriate file name templates and directory for input file
aploadplan,planfile,planstr,struct='ASPCAP'

if ~keyword_set(no_version_check) then begin
  if tag_exist(planstr,'apogee_version') then $
    if planstr.apogee_version ne getenv('APOGEE_VER') and planstr.apogee_version ne 'test' then  stop,'APOGEE version does not match!'
endif

; other parameters priority: 1. keyword, 2. planfile, 3. default
if not keyword_set(ncpus) then ncpus=apsetpar(planstr,'ncpus',4)
if not keyword_set(errbar) then errbar=apsetpar(planstr,'errbar',1)
if not keyword_set(renorm) then renorm=apsetpar(planstr,'renorm',0)
if keyword_set(renorm) then frdsuffix='.obs' else frdsuffix='.frd'
if not keyword_set(obscont) then obscont=apsetpar(planstr,'obscont',0)
if not keyword_set(nruns) then nruns=apsetpar(planstr,'nruns',1)
if n_elements(queue) eq 0 then queue=apsetpar(planstr,'queue',0)
if not keyword_set(apogee_vers) then apogee_vers=apsetpar(planstr,'apogee_vers','unknown')
if not keyword_set(aspcap_vers) then aspcap_vers=apsetpar(planstr,'aspcap_vers','aspcap_'+getenv('APOGEE_VER'))
if not keyword_set(results_vers) then results_vers=apsetpar(planstr,'results_vers','aspcap_'+getenv('APOGEE_VER'))
if not keyword_set(apred_vers) then apred_vers=apsetpar(planstr,'apred_vers','unknown')
if not keyword_set(apstar_vers) then apstar_vers=apsetpar(planstr,'apstar_vers','unknown')
if not keyword_set(aspcap_config) then aspcap_config=apsetpar(planstr,'aspcap_config','unknown')
if not keyword_set(telescope) then telescope=apsetpar(planstr,'telescope','apo25m')
if not keyword_set(hmask) then hmask=apsetpar(planstr,'hmask',0)
if not keyword_set(maskfile) then maskfile=apsetpar(planstr,'maskfile',0)
if not keyword_set(pixmask) then pixmask=apsetpar(planstr,'pixmask',0)
if not keyword_set(obspixmask) then obspixmask=apsetpar(planstr,'obspixmask',0)
if not keyword_set(conthighbad) then conthighbad=apsetpar(planstr,'conthighbad',1.1)
if not keyword_set(contlowbad) then contlowbad=apsetpar(planstr,'contlowbad',0.001)
if not keyword_set(conthigherr) then conthigherr=apsetpar(planstr,'conthigherr',1.e10)
if not keyword_set(highbad) then highbad=apsetpar(planstr,'highbad',1.1)
if not keyword_set(lowbad) then lowbad=apsetpar(planstr,'lowbad',0.001)
if not keyword_set(higherr) then higherr=apsetpar(planstr,'higherr',1.e10)
if not keyword_set(commiss) then commiss=apsetpar(planstr,'commiss',0)
if not keyword_set(nored) then nored=apsetpar(planstr,'nored',0)
if not keyword_set(minerr) then minerr=apsetpar(planstr,'minerr',0.005)
if not keyword_set(vacuum) then vacuum=apsetpar(planstr,'vacuum',0)
if n_elements(persist) eq 0 then persist=apsetpar(planstr,'persist',0)
if n_elements(skyerr) eq 0 then skyerr=apsetpar(planstr,'skyerr',3)
if n_elements(skyfact) eq 0 then skyfact=apsetpar(planstr,'skyfact',1000.)
if n_elements(visits) eq 0 then visits=apsetpar(planstr,'visits',0)
if n_elements(npar) eq 0 then npar=apsetpar(planstr,'npar',7)
if n_elements(nelem) eq 0 then nelem=apsetpar(planstr,'nelem',0)
if n_elements(notie) eq 0 then notie=apsetpar(planstr,'notie',0)
if tag_exist(planstr,'fits') then fits=planstr.fits
if not keyword_set(symlink) then symlink=apsetpar(planstr,'symlink',0)
if n_elements(indv) eq 0 then indv=apsetpar(planstr,'indv',0)
; following default only appropriate for 6 parameter library in "normal" parameter order!
;if n_elements(indini) eq 0 then if tag_exist(planstr,'INDINI') then indini=planstr.indini
;if n_elements(inter) eq 0 then if tag_exist(planstr,'INTER') then inter=planstr.inter
if not keyword_set(nstars) then nstars=apsetpar(planstr,'nstars',0)
if not keyword_set(starlist) then starlist=apsetpar(planstr,'starlist',0)
if n_elements(clobber) eq 0 then clobber=apsetpar(planstr,'clobber',0)
if keyword_set(qaspcap) then begin
  conthighbad=1.1
  contlowbad=0.001
  conthigherr=1.
  highbad=100.
  lowbad=0.001
  higherr=100.
  skyerr=0
  minerr=0
endif
if n_elements(noplot) eq 0 then noplot=apsetpar(planstr,'noplot',0)
if not keyword_set(doelemplot) then doelemplot=apsetpar(planstr,'doelemplot',0)
if keyword_set(altmaskdir) then altmaskdir=getenev('APOGEE_DIR')+'/data/windows/'+altmaskdir else undefine,altmaskdir
if n_elements(noelem) eq 0 then noelem=apsetpar(planstr,'noelem',0)
if not keyword_set(maxwind) then maxwind=apsetpar(planstr,'maxwind',0)
if not keyword_set(caldir) then caldir=apsetpar(planstr,'caldir',0)
; directories
apsetver,vers=apred_vers,telescope=telescope
dirs=getdir(apogee_dir,cal_dir,spectro_path,vers)
if not keyword_set(aspcap_root) then aspcap_root=apsetpar(planstr,'aspcap_root',dirs.aspcap+'/')
if not keyword_set(redux_root) then redux_root=apsetpar(planstr,'redux_root',dirs.redux)
if not keyword_set(red_vers) then red_vers=apsetpar(planstr,'red_vers',apogee_vers)
if not keyword_set(liblist_path) then liblist_path=apsetpar(planstr,'liblist_path',getenv('APOGEE_DIR')+'/config/aspcap/')
if not keyword_set(libr_path) then libr_path=apsetpar(planstr,'libr_path',dirs.speclib+'/synth/')
src_path=getenv('APOGEE_DIR') 
exec_path=src_path+'/bin/'

; set base names and file locations
if planstr.apvisit eq 1 then begin
  ; for apVisit files, still in development?
  apvisit=planstr.apvisit
  id=planstr.aspcap.plateid
  mjd=planstr.aspcap.mjd
  fits='apVisit-*-*-*.fits'
  datadir='plates/'+string(format='(i4.4)',id)+'/'+string(format='(i5.5)',mjd)+'/'
  outdir=datadir
  oname=string(format='(i4.4)',id)+'-'+string(format='(i5.5)',mjd)
  odir='plates/'
  ofile='aspcapPlate'
endif else begin
  ; for apStar files
  field=planstr.aspcap.field
  ;if not keyword_set(fits) then fits='apStar-*2M????????????????.fits'
  if not keyword_set(fits) then fits='a?Star-*.fits'
  if keyword_set(commiss) then fits='a?StarC-*2M????????????????.fits'
  if tag_exist(planstr.aspcap,'outfield') then outfield=planstr.aspcap.outfield  else outfield=field
  datadir=apstar_vers+'/'+strtrim(field,2)+'/'
  outdir=strtrim(outfield,2)+'/'
  file_mkdir,outdir
  oname=strtrim(file_basename(outfield),2)
  odir='stars/'
  ofile='aspcapField'
  if keyword_set(commiss) then begin
    oname=oname+'C'
    ofile='aspcapFieldC'
    outdir=strtrim(outfield,2)+'C/'
  endif
endelse

; set up ASPCAP parameter names
params=aspcap_params(npar=npar)

; get library name(s) 
;readcol,liblist_path+aspcap_config+'/class-'+dirs.instrument+'.list',format='(a)',class,stringskip='#'
aploadplan,liblist_path+aspcap_config+'/class-'+dirs.instrument+'.par',classplan,str='CLASS'
class=classplan.class.class

; create configuration files
configdir=aspcap_root+apred_vers+'/'+aspcap_vers+'/config/'+dirs.instrument+'/'
file_mkdir,configdir
while file_test(configdir+'lockfile') do apwait,configdir+'lockfile',10
if ~file_test(configdir+'/done') then begin
  openw,lock,/get_lun,configdir+'/lockfile'
  free_lun,lock
  aspcap_mklib,aspcap_config,outdir=configdir,maskdir=getenv('APOGEE_DIR')+'/data/windows/filters_26042016/',maxwind=maxwind
  ; create individual windows if requested
  if keyword_set(maxwind) then  begin
    readcol,configdir+'/elem.list',format='(a)',els,stringskip='#',/silent
    for i=0,n_elements(els)-1 do n=filtsplit(els[i],maxwind=maxwind,outdir=configdir)
  endif
  openw,done,/get_lun,configdir+'/done'
  free_lun,done
  file_delete,configdir+'/lockfile'
endif
  
if keyword_set(testmjd) then openw,done,'mjd.lis',/get_lun
; loop over the input plates
for idir=0,n_elements(datadir)-1 do begin
 ; input data for this plate
 indir=redux_root+'/'+apred_vers+'/'+datadir[idir]

 ; get stars to process
 mjddir=''
 if keyword_set(starlist) then begin
   readcol,starlist,names,te,logg,metals,format='(a,f,f,f)'
   files='apStar-'+apred_vers+'-'+names+'.fits'
   nstars=0
 endif else begin
   print,indir,field[idir]
   list =mrdfits(apogee_filename('Field',field=field[idir]),1,status=status)
   if keyword_set(commiss) then $
   list=mrdfits(indir+'apFieldC-'+strtrim(file_basename(field[idir]),2)+'.fits',1,status=status) 
   if status lt 0 or size(list,/type) eq 3 then $
     files=file_basename(file_search(indir+fits,test_symlink=symlink)) else begin
       files=strtrim(list.file,2)
       list=mrdfits(indir+'apFieldVisits-'+strtrim(file_basename(field[idir]),2)+'.fits',1)
       mjdlast=max(list[uniq(list.mjd,sort(list.mjd))].mjd)
       mjddir=strtrim(mjdlast,2)
     endelse
 endelse
 resultsdir=aspcap_root+apred_vers+'/'+aspcap_vers+'/'+outdir[idir]+'/'

 ; see if we've already done this up to latest MJD
 print,outdir[idir], mjddir,file_test(resultsdir+mjddir)
 if keyword_set(testmjd) then printf,done,outdir[idir], mjddir,file_test(resultsdir+mjddir)
 if mjddir ne '' and file_test(resultsdir+mjddir) then goto,nextdir
 if files[0] eq '' then goto,nextdir
 if keyword_set(testmjd) then goto,nextdir

 nclass=n_elements(class)
 first=1
 for iclass=0,nclass-1 do begin
  clock=TIC(class[iclass])
  print,iclass,class[iclass]
  aploadplan,configdir+'/'+class[iclass]+'.par',libpar,str='PLOCK'
  libfile=libr_path+libpar.lib

  ; get library parameters for this class
  rdlibhead,libfile,libhead0,libhead
  nlib=n_elements(libhead)
  npix=0
  for i=0,nlib-1 do npix+=libhead[i].npix
  info=file_info(libfile+'.unf')
  libsize=strcompress(string(round((info.size/1d9)+0.5),format='(I2)'),/remove_all)
  if tag_exist(libpar,'INITPAR') then par=libpar.initpar else $ 
    par=fltarr(libhead0.n_of_dim)
  parerr=fltarr(libhead0.n_of_dim)
  ; if this grid is marked with init=0, then get initial guess of parameters from previous coarse runs
  if tag_exist(libpar,'coarse') then coarse=libpar.coarse else coarse=0
  if tag_exist(libpar,'init') then init=libpar.init else init=1
  if first eq 1 and coarse eq 0 then begin
     beststr=aspcap_bestclass(allparam,allspec,alllib)
     first=0
  endif
  if tag_exist(libpar,'indini') then indini=libpar.indini else indini=0
  if tag_exist(libpar,'renorm') then renorm=libpar.renorm
  if tag_exist(libpar,'obscont') then obscont=libpar.obscont
  if keyword_set(renorm) then frdsuffix='.obs' else frdsuffix='.frd'

  ; output directory for this class
  specdir=aspcap_root+apred_vers+'/'+aspcap_vers+'/'+outdir[idir]+'/ferre/spectra/'
  file_mkdir,specdir
  workdir=aspcap_root+apred_vers+'/'+aspcap_vers+'/'+outdir[idir]+'/ferre/class_'+class[iclass]+'/'
  file_mkdir,workdir
  outname=class[iclass]+'-'+oname[idir]

  ; loop through the requested directory and construct input FERRE files
  print,workdir+outname+'.ipf  ',file_test(workdir+outname+'.ipf')
  nfit=0
  nobj=0
  if not file_test(workdir+outname+'.ipf') or keyword_set(clobber) then begin
   
    ; write the FERRE control file, and open the FERRE input files
    if tag_exist(libpar,'indi') then indi=libpar.indi else undefine,indi
    if tag_exist(libpar,'mask') then filterfile=configdir+'/'+libpar.mask else undefine,filterfile
    writeferre,workdir,outname,libhead0,nruns=nruns,ncpus=ncpus,indv=libpar.indv,indini=indini,$
            interord=libpar.inter,findi=indi,errbar=errbar,init=init,renorm=abs(renorm),obscont=obscont,$
            filterfile=filterfile
    openw,ipf,workdir+outname+'.ipf',/get_lun
    openw,labl,workdir+outname+'.labl',/get_lun
    openw,cfile,workdir+outname+'.con',/get_lun
    openw,frd,workdir+outname+frdsuffix,/get_lun 
    openw,err,workdir+outname+'.err',/get_lun
    openw,wav,workdir+outname+'.wav',/get_lun
    if keyword_set(nstars) then nobj=nstars else nobj=n_elements(files)
    for i=0,nobj-1 do begin
     if file_test(indir+files[i]) then begin
      print,indir+files[i]

      ; now process and write the spectral data 
      if keyword_set(persist) then persistfile=outdir+oname+'.persist'
      aspcap_rdobsnew,indir+files[i],obs,gaps,endgaps=endgaps,apvisit=apvisit,$
        skyerr=skyerr,skyfact=skyfact,visits=visits,persist=persistfile
      sz=size(obs.flux)
      if sz[0] eq 1 then nspec=1 else nspec=sz[2]
      if sxpar(obs.head,'VHELIO') lt 999990. and $
         (sxpar(obs.head,'ANDFLAG') and starflagval('BAD_PIXELS')) eq 0 and $
         (sxpar(obs.head,'ANDFLAG') and starflagval('LOW_SNR')) eq 0 and $
         (~keyword_set(persist) or (sxpar(obs.head,'ANDFLAG') and starflagval('PERSIST_HIGH')) gt 0) then begin
        if keyword_set(hmask) then hmask,obs,wid=hmask
        if keyword_set(maskfile) then aspcap_mask,obs,liblist_path+aspcap_config+'/'+maskfile
        vactoair,obs.wave,wout
        if tag_exist(libhead[0],'VACUUM') then if libhead[0].vacuum eq 1 then wout=obs.wave
        if keyword_set(vacuum) then wout=obs.wave
        wave0=wout*((1d)-sxpar(obs.head,'VRAD')*1.d5/2.99792458d10)
        nchips=n_elements(gaps)/2
        if keyword_set(obspixmask) then $
          aspcap_pixmask,obs,liblist_path+aspcap_config+'/'+obspixmask,/obs
        ; loop over summed spectrum and visits if requested
        for ispec=0,nspec-1 do begin
         ; make sure visit spectrum doesn't have too many bad pixels
         jbad=where(obs.flag[*,ispec] and badmask(),nbad)
         print,ispec,nbad,size(obs.flag[*,ispec],/dim)
         if ispec eq 0 or float(nbad)/n_elements(obs.flux[*,ispec]) lt 0.4 then begin
           if nbad gt 0 then obs.invar[jbad] = 0.01
           if keyword_set(mask_telluric) then begin
             jbad=where(obs.flag[*,ispec] and maskval('SIG_TELLURIC'),nbad)
             if nbad gt 0 then obs.invar[jbad] = 0.01
           endif
           ; loop over the chips
           for ichip=0,nchips-1 do begin
            libwave=libhead[ichip].wave[0]+indgen(libhead[ichip].npix)*libhead[ichip].wave[1]
            if libhead[ichip].logw eq 1 then libwave=10^libwave
            out=interpol(obs.flux[*,ispec],alog10(wave0),alog10(libwave))
            outinvar=1./interpol(1./obs.invar[*,ispec],alog10(wave0),alog10(libwave))
            ; outerr interpolates in error, not inverse variance
            outerr=interpol(sqrt(1./obs.invar[*,ispec]),alog10(wave0),alog10(libwave))
            cpars=libhead[ichip].continuum
            if keyword_set(renorm) then begin
              cpars[0]=abs(renorm)
              cpars[1]=1
              cpars[2]=0.
              cpars[3]=0.
              conthighbad=2.0
              highbad=2.0
              conthigherr=0.5
            endif
    
            ; set initial mask for use in continuum determination
            mask=intarr(n_elements(out))+1
            bad=where(out eq 0 or outinvar le 0 or finite(outinvar) eq 0,nbad)
            if nbad gt 0 then mask[bad]=0
            if float(nbad)/n_elements(out) gt 0.4 then begin
              print,'badfrac 1 rejected: ', ichip,indir+files[i]
              goto,badfrac
            endif
            speclib_continuum,cpars[0],cpars[1],cpars[2],cpars[3],out,cont,mask=mask
            ; add to mask based on first iteration continuum, and redo continuum
            if n_elements(conthighbad) gt 0 then begin
              bad = where(out/cont gt conthighbad,nbad)
              if nbad gt 0 then mask[bad]=0
            endif
            if n_elements(contlowbad) gt 0 then begin
              bad = where(out/cont le contlowbad,nbad)
              if nbad gt 0 then mask[bad]=0
            endif
            if n_elements(conthigherr) gt 0 then begin
              bad = where(outerr/out gt conthigherr,nbad)
              if nbad gt 0 then mask[bad]=0
            endif
            bad=where(mask eq 0,nbad)
            if float(nbad)/n_elements(out) gt 0.4 then begin
              print,'badfrac 2 rejected: ', ichip,indir+files[i]
              goto,badfrac
            endif
            speclib_continuum,cpars[0],cpars[1],cpars[2],cpars[3],out,cont,mask=mask ;,ivar=outinvar  ivar weights by throughput!
            if keyword_set(nored) and ichip eq 2 then outinvar=0.
            ; concatenate chips
            if ichip eq 0 then pseudo=cont else pseudo=[pseudo,cont]
            if ichip eq 0 then new=out/cont else new=[new,out/cont]
            if ichip eq 0 then newerr=outerr/cont else newerr=[newerr,outerr/cont]
            if ichip eq 0 then newinvar=outinvar*cont^2 else newinvar=[newinvar,outinvar*cont^2]
            if ichip eq 0 then libwaveall=libwave else libwaveall=[libwaveall,libwave]
           endfor
           ; if we have a previous run, get best coarse run for initial parameters
           ; if renorm<0, get adjustment to continuum, before flagging bad pixels
           starname=file_basename(files[i],'.fits')
           if ispec gt 0 then starname=starname+'_v'+string(format='(i3.3)',ispec)
           if coarse eq 0 then begin
             best=where(strpos(beststr.param.apogee_id,starname) ge 0) 
             for ipar=0,n_elements(par)-1 do begin
               jpar=where(strtrim(beststr.lib.param_symbol,2) eq strtrim(libhead0.label[ipar],2))
               par[ipar]=beststr.param[best[0]].fparam[jpar]
             endfor
             coarsename=aspcap_root+apred_vers+'/'+aspcap_vers+'/'+outdir[idir]+$
                        '/ferre/class_'+beststr.param[best[0]].class+'/'+beststr.param[best[0]].class+'-'+oname[idir]
             aspcap_load,coarsename+'.norm',data
             normspec=float(data)
             if renorm lt 0 then new/=normspec[*,best]
           endif
           ; mask if we have a pixmask
           if keyword_set(pixmask) then aspcap_pixmask,new,liblist_path+aspcap_config+'/'+pixmask
           ; set masks for analysis
           if n_elements(highbad) gt 0 then begin
             bad=where(new gt highbad,nbad) 
             if nbad gt 0 then new[bad]=0.
           endif
           if n_elements(lowbad) gt 0 then begin
             bad=where(new lt lowbad,nbad) 
             if nbad gt 0 then new[bad]=0.
           endif
           if n_elements(higherr) gt 0 then begin
             bad=where(newerr/new gt higherr,nbad) 
             if nbad gt 0 then new[bad]=0.
           endif
           ferr=sqrt(1./newinvar)
           bad=where(finite(ferr) eq 0 or ferr gt 100,nbad)
           if nbad gt 0 then new[bad]=0.
           if nbad gt 0 then ferr[bad]=10.
           bad=where(new le 0,nbad)
           if nbad gt 0 then ferr[bad]=10.
           if keyword_set(minerr) then begin
             j=where(ferr lt minerr and ferr gt 0,nj)
             if nj gt 0 then ferr[j]=minerr
           endif
  
           ; write the FERRE input files
           cformat="("+string(npix)+"(E14.6))"
           fformat="("+string(npix)+"(F12.6))"
           eformat="("+string(npix)+"(F12.6))"
           wformat="("+string(npix)+"(F12.4))"
           iformat="(a,"+string(2*libhead0.n_of_dim)+"(F10.3))"
           quote="'"
           lformat='('+string(libhead0.n_of_dim)+'('+'"'+quote+'",a,'+'"'+quote+'",1x))'

      ; write input star file
           if keyword_set(starlist) then begin
             ipar=where(libhead0.label eq 'TEFF',nipar)
             if nipar gt 0 then par[ipar] = te[i]
             ipar=where(libhead0.label eq 'LOGG',nipar)
             if nipar gt 0 then par[ipar] = logg[i]
             ipar=where(libhead0.label eq 'METALS',nipar)
             if nipar gt 0 then par[ipar] = metals[i]
           endif
           ; determine whether we should skip the star for this grid, depending on parameters and mean fiber
           skip=0
           if tag_exist(libpar,'pmin') then $
             for ipar=0,n_elements(par)-1 do $
               if par[ipar] lt libpar.pmin[ipar] or par[ipar] gt libpar.pmax[ipar] then begin
                 skip=1
                 reason='parameter '+string(ipar)+' out of range'
               endif
           if sxpar(obs.head,'MEANFIB') gt 0 then meanfib=sxpar(obs.head,'MEANFIB') else meanfib=150
           if tag_exist(libpar,'fibermin') then $
             if meanfib lt libpar.fibermin then begin
               skip=1
               reason='fiber out of range'
             endif
           if tag_exist(libpar,'fibermax') then $
             if meanfib gt libpar.fibermax then begin
               skip=1
               reason='fiber out of range'
             endif
           if ~skip then begin
             for ipar=0,n_elements(par)-1 do begin
               ; make sure starting value is within grid
               if par[ipar] lt libhead0.llimits[ipar] then $
                 par[ipar] = libhead0.llimits[ipar]+libhead0.steps[ipar]/2.
               if par[ipar] gt libhead0.llimits[ipar]+(libhead0.n_p[ipar]-1)*libhead0.steps[ipar] then  $
                 par[ipar] =libhead0.llimits[ipar]+(libhead0.n_p[ipar]-1)*libhead0.steps[ipar]-libhead0.steps[ipar]/2.
             endfor
             ; write the FERRE data files 
             printf,cfile,pseudo,format=cformat
             printf,frd,new,format=fformat
             printf,err,ferr,format=eformat
             printf,wav,libwaveall,format=wformat
             printf,ipf,starname,format=iformat,par,parerr
             printf,labl,libhead.label,format=lformat
             nfit+=1
           endif else print,'  skipping : ', reason
          badfrac: 
         endif else print,'  skipping (too many bad pixels): ' ; fraction of bad pixels too high
        endfor   ; loop over visit spectra
      endif else print,'skipping (bad RV, STARFLAG, S/N or persist option): ', sxpar(obs.head,'STARFLAGS'), sxpar(obs.head,'VHELIO')
     endif else print,'missing: ', indir+files[i]
    endfor   ; loop over objects
    free_lun,ipf
    free_lun,labl
    free_lun,frd
    free_lun,cfile
    free_lun,err
    free_lun,wav
  endif
  ; run FERRE !  
  ; with queue=0, wait for FERRE to be done, and run in foreground
  ; with queue=1, wait for FERRE to be done, but run through PBS
  ; with queue=2, put FERRE job in background and move on
  cd,workdir,current=cwd
  ; only run FERRE if .spm file doesn't already exist, or clobber set
  if nfit gt 0 and (not file_test(workdir+outname+'.spm') or clobber ne 0) then begin
    file_delete,workdir+outname+'.spm',/allow_nonexistent    
    aspcap_wrpbsscript,outname,exec_path,ncpus,jobsid,libsize,queue=queue,qname=qname,workdir=workdir
    spawn,['./'+outname+'.pbs'],/noshell
  endif
  classend:
  dt=TOC(clock)
  print,'class:',class[iclass],dt,nobj,nfit,dt/nfit
  TOC
  ; create output plots and output FITS file
  resultsdir=aspcap_root+apred_vers+'/'+aspcap_vers+'/'+outdir[idir]+'/param/class_'+class[iclass]+'/'
  file_mkdir,resultsdir
  nin=file_lines(workdir+outname+'.ipf')
  if file_test(workdir+outname+'.spm') then nout=file_lines(workdir+outname+'.spm') else nout=0
  if nin gt 0 then begin
   if nout eq nin then begin
     undefine,path
     print,'plateresults: ',oname[idir],' ',class[iclass]
     if tag_exist(libpar,'noplot') then cnoplot=libpar.noplot else cnoplot=0
     if not keyword_set(noplot) and cnoplot eq 0 then $
       aspcap_plateresults,oname[idir],class[iclass],results,$
         interord=interord,version=aspcap_vers,/plot,/ps,/single,$
         ext=ext,versdat=apred_vers,out_path=resultsdir,obs_path=redux_root,path=path,$
         apvisit=apvisit,libr_path=libr_path,queue=queue,ppath=workdir,oobspath=indir;,/zoom,/append
   endif else stop,'HALT: Number of SPM output not equal to IPF input'
   path=aspcap_root+apred_vers+'/'+aspcap_vers+'/'+outdir[idir]+'/'
   file_delete,resultsdir+'/'+ofile+'-'+strcompress(oname[idir],/remove_all)+'.fits',/allow_nonexist
   str=aspcap_loadferre(workdir+outname,libpar,libfile,npar=npar)
   if iclass eq 0 then begin
     allparam=str.param & allspec=str.spec & alllib=str.lib 
   endif else begin
     allparam=[allparam,str.param] & allspec=[allspec,str.spec] & alllib=[alllib,str.lib]
   endelse
   ;aspcap_data,str,indir+files,/getobject
   ; output FITS file and HTML summary
   aspcap_writefits,str,resultsdir+'/'+outname+'.fits'
   ;aspcap_writehtml,str,resultsdir+'/class_'+class[iclass],class[iclass]+'-'+oname[idir],npar=npar
   print,'correct'
   if coarse eq 1 then spawn,['continuum_correct',workdir+outname,libfile,'--write',workdir+outname+'.norm'],/noshell
   print,'done correct'
  endif
  cd,cwd
  TOC
  nextclass: 
 endfor ; nclass

 ; now get the best result for each star from among the different classes
 resultsdir=aspcap_root+apred_vers+'/'+aspcap_vers+'/'+outdir[idir]+'/'
 finalstr=aspcap_bestclass(allparam,allspec,alllib,/nocoarse,classes=class)

 ; correct parameters if we have caldir
 elemskip=2
 if keyword_set(caldir) then begin
   aspcap_data,finalstr,indir+files,/keepid
   paramstr=finalstr.param
   aspcap_correct,paramstr,0,aspcap_root+apred_vers+'/'+aspcap_vers+'/'+caldir+'/',/noelem
   finalstr.param=paramstr
   elemskip=1
 endif
 ; write aspcapField file (no aspcapStar files)
 aspcap_writefits,finalstr,resultsdir+'/'+ofile+'-'+strcompress(oname[idir],/remove_all)+'.fits',mjddir=mjddir

 ; redo CNO for whatver grids might be configured to do so
 if file_test(configdir+'/CN.elem.par') then  begin
   aploadplan,configdir+'/CN.elem.par',libpar,str='CNOINFO' 
   for iclass=0,n_elements(libpar.cnoinfo)-1 do begin
     print,'  class: ', libpar.cnoinfo[iclass].class
     libfile=libr_path+libpar.cnoinfo[iclass].libs
     rdlibhead,libfile,libhead0,libhead

     elemdir='elem_CN/'
     outname='CN-'+libpar.cnoinfo[iclass].class+'-'+oname[idir]
     workdir=aspcap_root+apred_vers+'/'+aspcap_vers+'/'+outdir[idir]+'/ferre/'+elemdir
     if (not file_test(workdir+outname+'.spm') or keyword_set(clobber)) then begin
       file_mkdir,workdir
       index=intarr(n_elements(libhead[0].label))
       for ipar=0,n_elements(libhead[0].label)-1 do index[ipar]=where(params eq strtrim(libhead[0].label[ipar],2))
       if tag_exist(libpar.cnoinfo[iclass],'mask') then begin
         filterfile=configdir+'/'+libpar.cnoinfo[iclass].mask+'.mask'
         readcol,filterfile,mask,/silent
         npix=n_elements(mask)
       endif else begin
         npix=99999
         undefine,filterfile
       endelse
       if tag_exist(libpar.cnoinfo,'indini') then indini=libpar.cnoinfo[iclass].indini else undefine,indini
       if tag_exist(libpar.cnoinfo,'renorm') then renorm=libpar.cnoinfo[iclass].renorm       
       if keyword_set(notie) then undefine,ttie else ttie=libpar.cnoinfo[iclass].ttie
       ; use init=0, with no indini, to use starting guess from parameter run
       writeferre,workdir,outname,libhead0,nruns=nruns,ncpus=ncpus,indv=libpar.cnoinfo[iclass].indv,$
         init=0,interord=libpar.cnoinfo[iclass].inter,$
         filterfile=configdir+'/'+libpar.cnoinfo[iclass].mask+'.mask',$
         findi=indi,errbar=errbar,renorm=abs(renorm),obscont=obscont,ttie=ttie
       cd,workdir,current=cwd
       file_delete,workdir+outname+frdsuffix,/allow_non
       file_delete,workdir+outname+'.err',/allow_non
       file_delete,workdir+outname+'.ipf',/allow_non
       file_link,'../spectra/'+libpar.cnoinfo[iclass].class+'-'+oname+frdsuffix,workdir+outname+frdsuffix
       file_link,'../spectra/'+libpar.cnoinfo[iclass].class+'-'+oname+'.err',workdir+outname+'.err'
       file_link,'../spectra/'+libpar.cnoinfo[iclass].class+'-'+oname+'.ipf',workdir+outname+'.ipf'
       nfit=0
       if file_test(workdir+outname+'.ipf') then nfit=file_lines(workdir+outname+'.ipf')
       if nfit gt 0 then begin
         spawn,['ferre.x',outname+'.nml'],result,/noshell,/stderr
         openw,foutput,outname+'.out',/get_lun
         for iline = 0,n_elements(result)-1 do printf,foutput,result[iline]
         free_lun,foutput
       endif
       cd,cwd
     endif
     ; to read the FERRE files, use the main class parameter file with the correct PLOCKs for this class
     nfit=0
     if file_test(workdir+outname+'.ipf') then nfit=file_lines(workdir+outname+'.ipf')
     if nfit gt 0 then begin
       aploadplan,configdir+'/'+libpar.cnoinfo[iclass].class+'.par',classpar,str='PLOCK' 
       str=aspcap_loadferre(workdir+outname,classpar,libfile,npar=npar,/elemfit) 
       for istar=0,n_elements(str.param)-1 do begin
         j=where(finalstr.param.apogee_id eq str.param[istar].apogee_id,nj)
         if nj gt 0 then begin
           for ii=0,n_elements(libpar.cnoinfo[iclass].indv)-1 do begin
             if libpar.cnoinfo[iclass].indv[ii] ge 0 then begin
               outindex=where(aspcap_params() eq libhead0.label[libpar.cnoinfo[iclass].indv[ii]-1])
               ;print,outindex,str.param[istar].fparam[outindex],finalstr.param[j].fparam[outindex]
               finalstr.param[j].fparam[outindex]=str.param[istar].fparam[outindex]
               finalstr.param[j].fparam_cov[outindex,outindex] = str.param[istar].fparam_cov[outindex,outindex]
               finalstr.param[j].paramflag=str.param[istar].paramflag[outindex]
             endif
           endfor
         endif
       endfor
     endif
   endfor
 endif   ; CNO


 ; individual elements using parameters from best class
 if ~keyword_set(noelem) and file_test(configdir+'/elem.list') then begin
   resultsdir=aspcap_root+apred_vers+'/'+aspcap_vers+'/'+outdir[idir]+'/'
   file_mkdir,resultsdir
   ; in preparation for elements, write out the best spectra into subdirectories for each class
   ; we will link to these in the element subdirectories
   if keyword_set(nstars) then nobj=nstars else nobj=n_elements(finalstr.param)
   for iclass=0,nclass-1 do begin
     print,'  class: ', class[iclass]
     aploadplan,configdir+'/'+class[iclass]+'.par',libpar,str='PLOCK'
     libfile=libr_path+libpar.lib
     ; get library parameters for this class
     rdlibhead,libfile,libhead0,libhead
     index=intarr(n_elements(libhead[0].label))
     for ipar=0,n_elements(libhead[0].label)-1 do $
       index[ipar]=where(params eq strtrim(libhead[0].label[ipar],2))
     iformat="(a,"+string(2*libhead0.n_of_dim)+"(F10.3))"
     fformat="("+string(npix)+"(F12.6))"
     eformat="("+string(npix)+"(F12.6))"
     openw,frd,specdir+class[iclass]+'-'+oname[idir]+frdsuffix,/get_lun
     openw,err,specdir+class[iclass]+'-'+oname[idir]+'.err',/get_lun
     openw,ipf,specdir+class[iclass]+'-'+oname[idir]+'.ipf',/get_lun
     nfit=0
     for i=0,nobj-1 do begin
       if finalstr.param[i].class eq class[iclass] then begin
         init=finalstr.param[i].fparam[index]
         printf,ipf,file_basename(finalstr.param[i].apogee_id,'.fits'),format=iformat,init,init*0.
         if n_elements(finalstr.spec[i].spec) eq 7212 then begin
           spec=[finalstr.spec[i].spec[0:5319],0.,finalstr.spec[i].spec[5320:7211],0.]
           nerr=[finalstr.spec[i].err[0:5319],0.,finalstr.spec[i].err[5320:7211],0.]
         endif else begin
           spec=finalstr.spec[i].spec
           nerr=finalstr.spec[i].err
         endelse
         printf,frd,spec,format=fformat
         ferr=nerr
         printf,err,ferr,format=eformat 
         nfit+=1
       endif
     endfor
     free_lun,ipf
     free_lun,frd
     free_lun,err
   endfor

   ; now extend the structures for elemental abundances
   elem_order=aspcap_elems(elemtagnames,elemtoh,elem_fitnames)
   nelem=n_elements(elem_order)
   readcol,configdir+'/elem.list',format='(a)',elem,stringskip='#',/silent
   if keyword_set(maxwind) then begin
     newparam=struct_rename_tags(finalstr.param,'','',addtag='FELEM',addval=fltarr(nelem,maxwind+1)-9999.)
     newparam=struct_rename_tags(newparam,'','',addtag='FELEM_ERR',addval=fltarr(nelem,maxwind+1)-9999.)
     newparam=struct_rename_tags(newparam,'','',addtag='FELEM_CAL',addval=fltarr(nelem,maxwind+1)-9999.)
     newparam=struct_rename_tags(newparam,'','',addtag='FELEM_CAL_ERR',addval=fltarr(nelem,maxwind+1)-9999.)
   endif else begin 
     newparam=struct_rename_tags(finalstr.param,'','',addtag='FELEM',addval=fltarr(nelem)-9999.)
     newparam=struct_rename_tags(newparam,'','',addtag='FELEM_ERR',addval=fltarr(nelem)-9999.)
     newparam=struct_rename_tags(newparam,'','',addtag='FELEM_CAL',addval=fltarr(nelem)-9999.)
     newparam=struct_rename_tags(newparam,'','',addtag='FELEM_CAL_ERR',addval=fltarr(nelem)-9999.)
   endelse
   newparam=struct_rename_tags(newparam,'','',addtag='ELEM',addval=fltarr(nelem)-9999.)
   newparam=struct_rename_tags(newparam,'','',addtag='ELEM_ERR',addval=fltarr(nelem)-9999.)
   newparam=struct_rename_tags(newparam,'','',addtag='X_H',addval=fltarr(nelem)-9999.)
   newparam=struct_rename_tags(newparam,'','',addtag='X_H_ERR',addval=fltarr(nelem)-9999.)
   newparam=struct_rename_tags(newparam,'','',addtag='X_M',addval=fltarr(nelem)-9999.)
   newparam=struct_rename_tags(newparam,'','',addtag='X_M_ERR',addval=fltarr(nelem)-9999.)
   newparam=struct_rename_tags(newparam,'','',addtag='ELEM_CHI2',addval=fltarr(nelem))
   newparam=struct_rename_tags(newparam,'','',addtag='ELEMFLAG',addval=intarr(nelem))
   newspec=struct_rename_tags(finalstr.spec,'','',addtag='ELEM_MASK',addval=long(finalstr.spec[0].spec*0.))
   ;tmpspec=struct_rename_tags(finalstr.spec,'','',addtag='ELEM_BESTFIT',addval=finalstr.spec[0].spec*0.)
   ;newspec=struct_rename_tags(tmpspec,'','',addtag='ELEM_MASK',addval=long(finalstr.spec[0].spec*0.))
   newlib=struct_rename_tags(finalstr.lib,'','',addtag='ELEM_SYMBOL',addval=elem_order)
   newlib=struct_rename_tags(newlib,'','',addtag='ELEM_VALUE',addval=elem_fitnames)
   newlib=struct_rename_tags(newlib,'','',addtag='ELEMTOH',addval=elemtoh)
   if keyword_set(maxwind) then $
     newlib=struct_rename_tags(newlib,'','',addtag='FELEM_WIND',addval=fltarr(nelem,maxwind,3)-9999.)
   ; create new nmlfiles (list of nml files)
   for ielem=0,n_elements(elem)-1 do begin
     aploadplan,configdir+'/'+elem[ielem]+'.elem.par',libpar,str='INFO' 
     for iclass=0,n_elements(libpar.info)-1 do begin
       openw,nmlfiles,aspcap_root+apred_vers+'/'+aspcap_vers+'/'+outdir[idir]+'/ferre/'+libpar.info[iclass].class+'.nmlfiles',/get_lun
       free_lun,nmlfiles
     endfor
   endfor
   ; loop over elements twice, first time to create nml files and nmlfiles file lists, then run FERRE
   ; also include run with calibrated parameters
   ; elemloop=0 setup uncalibrated run, 1 setup calibrated run, and run
   ; elemloop=2 read uncalibrated run, 3 read calibrated run
   ; second time, read the results
   for elemloop=0,3,elemskip do begin
    for ielem=0,n_elements(elem)-1 do begin
     jelem=where(strtrim(elem_order,2) eq strtrim(elem[ielem],2),norder)
     if norder eq 0 then stop,'cant match elem in aspcap_elems'
     clock=TIC(elem[ielem])
     TOC
     print,'Element:', elem[ielem], ' elemloop: ', elemloop
     if file_test(configdir+'/'+elem[ielem]+'.elem.par') then $
     aploadplan,configdir+'/'+elem[ielem]+'.elem.par',libpar,str='INFO' else $
     aploadplan,configdir+'/'+elem[ielem]+'.par',libpar,str='INFO' 

     ; loop over all classes for which this element can be derived (given in .par file)
     for iclass=0,n_elements(libpar.info)-1 do begin
       print,'  class: ', libpar.info[iclass].class
       libfile=libr_path+libpar.info[iclass].libs
       rdlibhead,libfile,libhead0,libhead

       index=intarr(n_elements(libhead[0].label))
       for ipar=0,n_elements(libhead[0].label)-1 do $
            index[ipar]=where(params eq strtrim(libhead[0].label[ipar],2))
       cindex=where(libhead[0].label eq 'C')
       cfit=where(elem_order eq 'C')
       if tag_exist(libpar.info[iclass],'mask') then begin
         filterfile=configdir+'/'+libpar.info[iclass].mask+'.mask'
         readcol,filterfile,mask,/silent
         npix=n_elements(mask)
       endif else begin
         npix=99999
         undefine,filterfile
       endelse

       ; with maxwind keyword, we will loop over the full set of windows plus each individual window
       if keyword_set(maxwind) and file_test(configdir+elem[ielem]+'.wind') and ~tag_exist(libpar.info,'minigrid') then begin
         readcol,configdir+elem[ielem]+'.wind',w1,w2,w3
         if tag_exist(libpar,'indi') then indi=libpar.indi else undefine,indi
         nwind=n_elements(w1)
       endif else nwind=0
       for iwind=0,nwind do begin
         if iwind eq 0 then suffix='' else suffix='_'+string(format='(i2.2)',iwind)
         elemdir='elem_'+elem[ielem]+suffix+'/'
         workdir=aspcap_root+apred_vers+'/'+aspcap_vers+'/'+outdir[idir]+'/ferre/'+elemdir
         iformat="(a,"+string(2*libhead0.n_of_dim)+"(F10.3))"
         fformat="("+string(npix)+"(F12.6))"
         eformat="("+string(npix)+"(F12.6))"
         outname=elem[ielem]+'-'+libpar.info[iclass].class+'-'+oname[idir]
         if elemloop eq 1 or elemloop eq 3 then outname=outname+'_cal'
         file_mkdir,workdir
         if tag_exist(libpar.info,'indini') then indini=libpar.info[iclass].indini else undefine,indini
         if tag_exist(libpar.info,'renorm') then renorm=libpar.info[iclass].renorm       
         if keyword_set(notie) then undefine,ttie else ttie=libpar.info[iclass].ttie
         ; use init=0, with no indini, to use starting guess from parameter run
         writeferre,workdir,outname,libhead0,nruns=nruns,ncpus=ncpus,indv=libpar.info[iclass].indv,$
            init=0,interord=libpar.info[iclass].inter,$
            filterfile=configdir+'/'+libpar.info[iclass].mask+suffix+'.mask',$
            findi=indi,errbar=errbar,renorm=abs(renorm),obscont=obscont,ttie=ttie,elemdir=elemdir,$
            libdir='../lib_'+libpar.info[iclass].class
         if tag_exist(libpar.info,'minigrid') then minigrid=libpar.info[iclass].minigrid else minigrid=''
         if keyword_set(nstars) then nobj=nstars else nobj=n_elements(finalstr.param)
         ; get parameters for all objects whose best fit was current class
         if minigrid ne '' then begin
           ; with minigrid, write flux and err files in workdir
           readcol,getenv('APOGEE_DIR')+'/data/windows/'+minigrid+'/'+elem[ielem]+'.wave',w1,w2,format='(d,d)'
           wvac=[]
           for iii=0,n_elements(w1)-1 do wvac=[[wvac],[alog10(w1[iii]),alog10(w2[iii])]]
           openw,ipf,workdir+outname+'.ipf',/get_lun
           openw,frd,workdir+outname+frdsuffix,/get_lun
           openw,err,workdir+outname+'.err',/get_lun
           for i=0,nobj-1 do begin
             if finalstr.param[i].class eq libpar.info[iclass].class then begin
               init=finalstr.param[i].fparam[index]
               missing=where(index lt 0,nmissing)
               if nmissing gt 0 then init[missing] = 0.
               printf,ipf,file_basename(finalstr.param[i].apogee_id,'.fits'),format=iformat,init,init*0.
               help,finalstr.lib.wave
               spec=speclib_welem(finalstr.lib.wave,finalstr.spec[i].spec,wvac)
               nerr=speclib_welem(finalstr.lib.wave,finalstr.spec[i].err,wvac)
               printf,frd,spec,format=fformat
               printf,err,nerr,format=eformat 
             endif
           endfor
           free_lun,frd
           free_lun,err
         endif else begin
           file_delete,workdir+outname+frdsuffix,/allow_non
           file_delete,workdir+outname+'.err',/allow_non
           file_delete,workdir+outname+'.ipf',/allow_non
           file_link,'../spectra/'+libpar.info[iclass].class+'-'+oname+frdsuffix,workdir+outname+frdsuffix
           file_link,'../spectra/'+libpar.info[iclass].class+'-'+oname+'.err',workdir+outname+'.err'
           file_link,'../spectra/'+libpar.info[iclass].class+'-'+oname+'.ipf',workdir+outname+'.ipf'
         endelse
         nfit=0
         if file_test(workdir+outname+'.ipf') then nfit=file_lines(workdir+outname+'.ipf')
         if nfit gt 0 then begin
           if elemloop lt 2 and (not file_test(workdir+outname+'.spm') or keyword_set(clobber)) then begin
              file_delete,workdir+outname+'.spm',/allow_nonexistent    
              cd,workdir+'/../',current=cwd
              openu,nmlfiles,aspcap_root+apred_vers+'/'+aspcap_vers+'/'+outdir[idir]+'/ferre/'+libpar.info[iclass].class+'.nmlfiles',/get_lun,/append
              printf,nmlfiles,elemdir+'/'+outname+'.nml'
              free_lun,nmlfiles
              file_delete,'lib_'+libpar.info[iclass].class,/allow
              file_link,file_dirname(libhead.file),'lib_'+libpar.info[iclass].class,/allow
              cd,cwd
           endif
           TOC
           ; read the FERRE output and load into array
           if elemloop ge  2 then begin
             if minigrid ne '' then $
             outindex=where(aspcap_params(extra=elem[ielem]) eq libhead0.label[libpar.info[iclass].indv[0]-1]) else $
             outindex=where(aspcap_params() eq libhead0.label[libpar.info[iclass].indv[0]-1])
             ; to read the FERRE files, use the main class parameter file with the correct PLOCKs for this class
             aploadplan,configdir+'/'+libpar.info[iclass].class+'.par',classpar,str='PLOCK' 
             if minigrid eq '' then  begin
               str=aspcap_loadferre(workdir+outname,classpar,libfile,npar=npar,/elemfit) 
               readcol,configdir+'/'+libpar.info[iclass].mask+'.mask',mask,/silent
               gd=where(mask gt 0.,ngd)
             endif else begin
               str=aspcap_loadferre(workdir+outname,classpar,libfile,npar=npar,extra=elem[ielem],/elemfit)
               ngd=0
             endelse
             for istar=0,n_elements(str.param)-1 do begin
              j=where(newparam.apogee_id eq str.param[istar].apogee_id,nj)
              if nj gt 0 then begin
                if keyword_set(maxwind) then begin 
                  if iwind gt 0 then newlib.felem_wind[jelem,iwind-1,*] = [w1[iwind-1],w2[iwind-1],w3[iwind-1]]
                  if elemloop eq 2 then begin
                    newparam[j].felem[jelem,iwind] = str.param[istar].fparam[outindex] 
                    if str.param[istar].fparam_cov[outindex,outindex] gt 0 then $
                      newparam[j].felem_err[jelem,iwind] = sqrt(str.param[istar].fparam_cov[outindex,outindex])
                  endif else begin 
                    newparam[j].felem_cal[jelem,iwind] = str.param[istar].fparam[outindex] 
                    if str.param[istar].fparam_cov[outindex,outindex] gt 0 then $
                      newparam[j].felem_cal_err[jelem,iwind] = sqrt(str.param[istar].fparam_cov[outindex,outindex])
                  endelse
                endif else  begin
                  if elemloop eq 2 then begin
                    newparam[j].felem[jelem] = str.param[istar].fparam[outindex]
                    if str.param[istar].fparam_cov[outindex,outindex] gt 0 then $
                      newparam[j].felem_err[jelem] = sqrt(str.param[istar].fparam_cov[outindex,outindex])
                  endif else begin 
                    newparam[j].felem_cal[jelem] = str.param[istar].fparam[outindex]
                    if str.param[istar].fparam_cov[outindex,outindex] gt 0 then $
                      newparam[j].felem_cal_err[jelem] = sqrt(str.param[istar].fparam_cov[outindex,outindex])
                  endelse
                endelse
                if iwind eq 0 then begin
                  newparam[j].elem_chi2[jelem] = str.param[istar].param_chi2
                  newparam[j].elemflag[jelem] = str.param[istar].paramflag[outindex]
                endif
                if ngd gt 0 then begin
                  ;newspec[j].elem_bestfit[gd]+=str.spec[istar].spec_bestfit[gd]
                  newspec[j].elem_mask[gd]=newspec[j].elem_mask[gd] or 2L^jelem[0]
                endif
              endif
             endfor
           endif
         endif else begin
           ; clean up if we didn't do any fits
           file_delete,workdir+'input.nml',/allow_non
           file_delete,workdir+'lib',/allow_non
           file_delete,workdir+outname+'.nml',/allow_non
           file_delete,workdir+outname+'.ipf',/allow_non
           file_delete,workdir+outname+frdsuffix,/allow_non
           file_delete,workdir+outname+'.err',/allow_non
         endelse
       endfor ; nwind
     endfor  ; nclass
     ;if minigrid ne '' then stop
     dt=TOC(clock)
     print,'elem:',elem[ielem],dt
    endfor ; nelem
    ; run ferre for all nmlfiles
    if (elemloop eq 1 and elemskip eq 1) or (elemloop eq 0 and elemskip eq 2) then begin
      workdir=aspcap_root+apred_vers+'/'+aspcap_vers+'/'+outdir[idir]+'/ferre/'
      cd,workdir,current=cwd
      nmlfile=file_search('*.nmlfiles')
      for ifile=0,n_elements(nmlfile)-1 do begin
        if file_lines(nmlfile[ifile]) gt 0 then begin
          print,'running: ',nmlfile[ifile]
          spawn,['ferre.x','-l',nmlfile[ifile]],result,/noshell,/stderr
          openw,foutput,nmlfile[ifile]+'.out',/get_lun
          for iline = 0,n_elements(result)-1 do printf,foutput,result[iline]
          free_lun,foutput
        endif
      endfor
      cd,cwd
    endif
   endfor  ; elemloop
   finalstr={param: newparam, spec: newspec, lib: newlib}
   if not keyword_set(noplot) and doelemplot gt 0 then begin
     file_mkdir,resultsdir+'/elem'
     openw,script,resultsdir+'/elem/convert.csh',/get_lun
     printf,script,'cd '+resultsdir+'/elem/'
     for istar=0,n_elements(finalstr.param)-1 do elemplot,finalstr,istar,hard=resultsdir+'/elem/',script=script,maskdir=configdir,altmaskdir=altmaskdir
     free_lun,script
     file_chmod,resultsdir+'/elem/convert.csh','770'o
     spawn,[resultsdir+'/elem/convert.csh'],/noshell
     file_delete,resultsdir+'/elem/convert.csh'
     for iel=0,n_elements(finalstr.lib.elem_symbol)-1 do begin
        el=strtrim(finalstr.lib.elem_symbol[iel],2)
        openw,ehtml,resultsdir+'/elem/'+el+'.html',/get_lun
        printf,ehtml,'<HTML><BODY><TABLE BORDER=2>' 
        for istar=0,n_elements(finalstr.param)-1 do $
          printf,ehtml,'<TR><TD>'+finalstr.param[istar].apogee_id+$
                    '<TD><IMG SRC='+finalstr.param[istar].apogee_id+'.'+el+'.jpg height=400>'
        free_lun,ehtml
     endfor
   endif
 endif  ; elem
 aspcap_writefits,finalstr,resultsdir+'/'+ofile+'-'+strcompress(oname[idir],/remove_all)+'.fits',mjddir=mjddir
 TOC

 ; add information from apStar file
 aspcap_data,finalstr,indir+files,/getobject

 ; add corrections and flags
 paramstr=finalstr.param
 print,'caldir: ', caldir
 if keyword_set(caldir) then begin
   ;calstr=mrdfits(aspcap_root+aspcap_vers+'/'+calfile,1)
   ;restore,filename=aspcap_root+aspcap_vers+'/'+calfile
   ;aspcap_docal,paramstr,finalstr.lib.elem_symbol,calstr
   aspcap_correct,paramstr,finalstr.lib.elem_symbol,aspcap_root+apred_vers+'/'+aspcap_vers+'/'+caldir+'/'
 endif
 finalstr.param=paramstr

 resultsdir=aspcap_root+apred_vers+'/'+aspcap_vers+'/'+outdir[idir]+''
 file_mkdir,resultsdir
 ; output HTML pages
 aspcap_writehtml,finalstr,resultsdir,oname[idir],/bestclass,npar=npar
 aspcap_writehtml,finalstr,resultsdir,oname[idir],/bestclass,/spectra,npar=npar
 if ~keyword_set(noelem) then aspcap_writehtml,finalstr,resultsdir,oname[idir],/bestclass,/elem,npar=npar

 ; write it out
 aspcap_writefits,finalstr,resultsdir+'/'+ofile+'-'+strcompress(oname[idir],/remove_all)+'.fits',/aspcapstar,apred_vers=apred_vers,mjddir=mjddir

 nextdir:
endfor  ; datadir
if keyword_set(testmjd) then free_lun,done
TOC
end

