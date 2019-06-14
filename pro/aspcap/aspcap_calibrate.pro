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

 for idir=0,n_elements(outdir)-1 do begin

   resultsdir=aspcap_root+planstr.apred_vers+'/'+planstr.aspcap_vers+'/'+outdir[idir]+'/'
   oname=strtrim(file_basename(outfield),2)
 
   param=mrdfits(resultsdir+'/'+ofile+'-'+strcompress(oname[idir],/remove_all)+'.fits',1)
   spec=mrdfits(resultsdir+'/'+ofile+'-'+strcompress(oname[idir],/remove_all)+'.fits',2)
   lib=mrdfits(resultsdir+'/'+ofile+'-'+strcompress(oname[idir],/remove_all)+'.fits',3)
   finalstr={param: param, spec: spec, lib: lib}

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
