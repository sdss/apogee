pro writeferre,outdir,file,libhead,nruns=nruns,interord=interord,direct=direct,pca=pca,ncpus=ncpus,indv=indv,findi=findi,indini=indini,init=init,filterfile=filterfile,errbar=errbar,renorm=renorm

; routine to write FERRE input file

; default FERRE parameters if not specified
if n_elements(ncpus) eq 0 then ncpus=2
if n_elements(nruns) eq 0 then nruns=1
if n_elements(interord) eq 0 then interord=3
if n_elements(direct) eq 0 then direct=1
if n_elements(pca) eq 0 then pca=1
if n_elements(errbar) eq 0 then errbar=-2

; open the output file
openw,nml,outdir+file+'.nml',/get_lun

printf,nml,' &LISTA'
ndim=libhead.n_of_dim
printf,nml,format='(A,1X,I1)',' NDIM =',ndim
if keyword_set(findi) then $
  printf,nml,format='(A,'+string(ndim)+'(I2))',' INDI =',findi 
if n_elements(indv) gt 0 then begin
  printf,nml,format='(A,1X,I1)',' NOV =',n_elements(indv)
  printf,nml,format='(A,'+string(ndim)+'(I2))',' INDV =',indv 
endif else begin
  printf,nml,format='(A,1X,I1)',' NOV =',ndim
  printf,nml,format='(A,'+string(ndim)+'(I2))',' INDV =',indgen(ndim)+1
endelse
print,'linking ...'+file_dirname(libhead.file)
file_delete,outdir+'/lib',/allow
file_link,file_dirname(libhead.file),outdir+'/lib',/allow
printf,nml,' SYNTHFILE(1) = '+"'"+'lib/'+file_basename(libhead.file)+'.hdr'+"'"
if keyword_set(filterfile) then printf,nml,' FILTERFILE ='+"'"+filterfile+"'"
printf,nml,' PFILE = '+"'./"+file+".ipf'"
printf,nml,' ERFILE = '+"'./"+file+".err'"
printf,nml,' OPFILE = '+"'./"+file+".spm'"
printf,nml,' OFFILE = '+"'./"+file+".mdl'"
printf,nml,' ERRBAR = '+string(errbar)
if keyword_set(renorm) then begin
  printf,nml,' CONT = 1'
  printf,nml,' NCONT = '+string(renorm)
  printf,nml,' OBSCONT =  0'
  printf,nml,' FFILE = '+"'./"+file+".obs'"
  printf,nml,' SFFILE = '+"'./"+file+".frd'"
endif else begin
  printf,nml,' FFILE = '+"'./"+file+".frd'"
endelse
if n_elements(init) gt 0 then printf,nml,' init='+strcompress(string(init),/remove_all)

if keyword_set(indini) then begin
  printf,nml,format='(A,'+string(ndim)+'(I2))',' INDINI =',indini 
  nruns=1
  for i=0,n_elements(indini)-1 do nruns*=indini[i]
  if nruns lt 1 then begin
    print,'Error in INDINI array!'
    stop
  endif
endif
printf,nml,' nruns='+strcompress(string(nruns),/remove_all)
printf,nml,' algor=3'
printf,nml,' nthreads='+string(ncpus)
printf,nml,' pcaproject=0'
printf,nml,' pcachi=0'
if n_elements(abundance) eq 0 then printf,nml,' covprint=1'
if n_elements(pca) eq 0 then printf,nml,' optimize=1'
if n_elements(interord) gt 0 then printf,nml,' inter='+strcompress(string(interord,format='(I3)'),/remove_all)
if n_elements(direct) gt 0 then begin
   printf,nml,' F_FORMAT= 1'
   printf,nml,' F_ACCESS= 0' ; TO CHANGE
endif
printf,nml,' /'
free_lun,nml
file_copy,outdir+file+'.nml',outdir+'input.nml',/overwrite

end
