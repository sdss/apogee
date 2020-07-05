pro aspcap_calibrate, planfile, caldir=caldir, test=test

 aploadplan,planfile,planstr,struct='ASPCAP'

 apsetver,vers=planstr.apred_vers,telescope=planstr.telescope
 dirs=getdir(apogee_dir,cal_dir,spectro_path,vers)
 if not keyword_set(aspcap_root) then aspcap_root=apsetpar(planstr,'aspcap_root',dirs.aspcap+'/')
 if not keyword_set(caldir) then caldir=apsetpar(planstr,'caldir','cal/')

 field=planstr.aspcap.field
 if tag_exist(planstr.aspcap,'outfield') then outfield=planstr.aspcap.outfield  else outfield=field
 outdir=strtrim(outfield,2)+'/'
 ofile='aspcapField'

 ; get the classes
 if not keyword_set(liblist_path) then liblist_path=apsetpar(planstr,'liblist_path',getenv('APOGEE_DIR')+'/config/aspcap/')
 aploadplan,liblist_path+planstr.aspcap_config+'/class-'+dirs.instrument+'.par',classplan,str='CLASS'
 class=classplan.class.class

 for idir=0,n_elements(outdir)-1 do begin

   resultsdir=aspcap_root+planstr.apred_vers+'/'+planstr.aspcap_vers+'/'+outdir[idir]+'/'
   oname=strtrim(file_basename(outfield),2)

   param=mrdfits(resultsdir+'/'+ofile+'-'+strcompress(oname[idir],/remove_all)+'.fits',1)
   spec=mrdfits(resultsdir+'/'+ofile+'-'+strcompress(oname[idir],/remove_all)+'.fits',2)
   lib=mrdfits(resultsdir+'/'+ofile+'-'+strcompress(oname[idir],/remove_all)+'.fits',3)
   finalstr={param: param, spec: spec, lib: lib}

   ; repair paramflag[0,4,5] which got messed up by CN calculation for dwarfs
   j=where(strpos(param.class,'Fd') eq 0 or strpos(param.class,'GKd') eq 0 or strpos(param.class,'Md') eq 0, nj)
   ; do C and N flags for all dwarf classes
   if nj gt 0 then begin
     param[j].paramflag[4] = 0
     bd = where(param[j].fparam[4] lt (-0.5 + 0.25/8.) or param[j].fparam[4] gt (0.5 - 0.25/8.),nbd)
     if nbd gt 0 then param[j[bd]].paramflag[4] = param[j[bd]].paramflag[4] or paramflagval('GRIDEDGE_BAD')
     bd = where(param[j].fparam[4] lt (-0.5 + 0.25/2.) or param[j].fparam[4] gt (0.5 - 0.25/2.),nbd)
     if nbd gt 0 then param[j[bd]].paramflag[4] = param[j[bd]].paramflag[4] or paramflagval('GRIDEDGE_WARN')
     bd = where(param[j].fparam[4] lt -999,nbd)
     if nbd gt 0 then param[j[bd]].paramflag[4]=param[j[bd]].paramflag[4] or paramflagval('FERRE_BAD')
     param[j].paramflag[5] = 0
     bd = where(param[j].fparam[5] lt (-0.5 + 0.5/8.) or param[j].fparam[5] gt (1.5 - 0.5/8.),nbd)
     if nbd gt 0 then param[j[bd]].paramflag[5] = param[j[bd]].paramflag[5] or paramflagval('GRIDEDGE_BAD')
     bd = where(param[j].fparam[5] lt (-0.5 + 0.5/2.) or param[j].fparam[5] gt (1.5 - 0.5/2.),nbd)
     if nbd gt 0 then param[j[bd]].paramflag[5] = param[j[bd]].paramflag[5] or paramflagval('GRIDEDGE_WARN')
     bd = where(param[j].fparam[5] lt -999,nbd)
     if nbd gt 0 then param[j[bd]].paramflag[5]=param[j[bd]].paramflag[5] or paramflagval('FERRE_BAD')
   endif
   ; now Teff flags for each class
   j=where(strpos(param.class,'Fd') eq 0, nj)
   if nj gt 0 then begin
     param[j].paramflag[0] = 0
     bd = where(param[j].fparam[0] lt (5500 + 250/8.) or param[j].fparam[0] gt (8000 - 250/8.),nbd)
     if nbd gt 0 then param[j[bd]].paramflag[0] = param[j[bd]].paramflag[0] or paramflagval('GRIDEDGE_BAD')
     bd = where(param[j].fparam[0] lt (5500 + 250/2.) or param[j].fparam[0] gt (8000 - 250/2.),nbd)
     if nbd gt 0 then param[j[bd]].paramflag[0] = param[j[bd]].paramflag[0] or paramflagval('GRIDEDGE_WARN')
   endif
   j=where(strpos(param.class,'GKd') eq 0, nj)
   if nj gt 0 then begin
     param[j].paramflag[0] = 0
     bd = where(param[j].fparam[0] lt (3500 + 250/8.) or param[j].fparam[0] gt (6000 - 250/8.),nbd)
     if nbd gt 0 then param[j[bd]].paramflag[0] = param[j[bd]].paramflag[0] or paramflagval('GRIDEDGE_BAD')
     bd = where(param[j].fparam[0] lt (3500 + 250/2.) or param[j].fparam[0] gt (6000 - 250/2.),nbd)
     if nbd gt 0 then param[j[bd]].paramflag[0] = param[j[bd]].paramflag[0] or paramflagval('GRIDEDGE_WARN')
   endif
   j=where(strpos(param.class,'Md') eq 0, nj)
   if nj gt 0 then begin
     param[j].paramflag[0] = 0
     bd = where(param[j].fparam[0] lt (3000 + 100/8.) or param[j].fparam[0] gt (4000 - 100/8.),nbd)
     if nbd gt 0 then param[j[bd]].paramflag[0] = param[j[bd]].paramflag[0] or paramflagval('GRIDEDGE_BAD')
     bd = where(param[j].fparam[0] lt (3000 + 100/2.) or param[j].fparam[0] gt (4000 - 100/2.),nbd)
     if nbd gt 0 then param[j[bd]].paramflag[0] = param[j[bd]].paramflag[0] or paramflagval('GRIDEDGE_WARN')
   endif

   ; repair PARAM_FIXED bits that were corrupted by previous aspcap_correct
   param.paramflag[8] = param.paramflag[8] or paramflagval('PARAM_FIXED')
   j=where(param.class eq 'BA',nj)
   if nj gt 0 then begin
     param[j].paramflag[4] = param[j].paramflag[4] or paramflagval('PARAM_FIXED')
     param[j].paramflag[5] = param[j].paramflag[5] or paramflagval('PARAM_FIXED')
     param[j].paramflag[6] = param[j].paramflag[6] or paramflagval('PARAM_FIXED')
     param[j].paramflag[7] = param[j].paramflag[7] or paramflagval('PARAM_FIXED')
   endif
   j=where(strpos(param.class, 'GKg') eq 0 or strpos(param.class,'Mg') eq 0,nj)
   if nj gt 0 then param[j].paramflag[7] = param[j].paramflag[7] or paramflagval('PARAM_FIXED')

   print,'caldir: ', caldir
   aspcap_correct,param,finalstr.lib.elem_symbol,aspcap_root+planstr.apred_vers+'/'+planstr.aspcap_vers+'/'+caldir+'/'
   finalstr.param=param

   if keyword_set(test) then begin
     resultsdir=resultsdir+'test/'
     file_mkdir,resultsdir
   endif

   ; output HTML pages
   aspcap_writehtml,finalstr,resultsdir,oname[idir],/bestclass,npar=npar
   aspcap_writehtml,finalstr,resultsdir,oname[idir],/bestclass,/spectra,npar=npar
   aspcap_writehtml,finalstr,resultsdir,oname[idir],/bestclass,/elem,npar=npar
  
   ; write it out
   aspcap_writefits,finalstr,resultsdir+'/'+ofile+'-'+strcompress(oname[idir],/remove_all)+'.fits',/aspcapstar,apred_vers=planstr.apred_vers
 endfor

end
