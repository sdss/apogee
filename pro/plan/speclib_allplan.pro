pro speclib_allplan,planfile

aploadplan,planfile,planstr
root=strsplit(planfile,'.',/extract)

synthcode=planstr.synthcode
config=planstr.config
if tag_exist(planstr,'atmos') then atmos=planstr.atmos else atmos=0
if tag_exist(planstr,'synthdir') then synthdir=planstr.synthdir else synthdir=0
if tag_exist(planstr,'marcsdir') then marcsdir=planstr.marcsdir
if tag_exist(planstr,'specdir') then specdir=planstr.specdir else specdir=0
if tag_exist(planstr,'smooth') then smooth=planstr.smooth else smooth=0
if tag_exist(planstr,'vmsuffix') then vmsuffix=planstr.vmsuffix else vmsuffix=''
if tag_exist(planstr,'continuum') then continuum=planstr.continuum else continuum=0
if tag_exist(planstr,'vacuum') then if planstr.vacuum eq 0 then air=1 else air=0
if tag_exist(planstr,'resolution') then resolution=planstr.resolution else resolution=22500
if tag_exist(planstr,'apred_vers') then apred_vers=planstr.apred_vers 
if tag_exist(planstr,'lsfid') then lsfid=planstr.lsfid else lsfid=0
if tag_exist(planstr,'waveid') then waveid=planstr.waveid else waveid=0
if tag_exist(planstr,'lsffiber') then lsffiber=planstr.lsffiber else lsffiber=0
if tag_exist(planstr,'vmicro') then if planstr.vmicro[0] gt 0 then vmicro=planstr.vmicro else vmicro=0 else vmicro=0
if tag_exist(planstr,'vmacro') then vmacro=planstr.vmacro else vmacro=0
if tag_exist(planstr,'vmicrofit') then vmicrofit=planstr.vmicrofit else vmicrofit=0
if tag_exist(planstr,'vmacrofit') then vmacrofit=planstr.vmacrofit else vmacrofit=0
if tag_exist(planstr,'width') then width=planstr.width else width=0
if tag_exist(planstr,'twod') then twod=planstr.twod else twod=0
if tag_exist(planstr,'kernel') then kernel=planstr.kernel else kernel=0
linelist=planstr.linelist
solarisotopes=planstr.solarisotopes
dw=planstr.dw
wrange=planstr.wrange
if tag_exist(planstr,'elem') then nelem=n_elements(planstr.elem) else nelem=0

; write individual plan files for each am,cm,nm,vt tuple for paralle processing
for icm=0,planstr.ncm-1 do begin
 cm=planstr.cm0+icm*planstr.dcm
 for inm=0,planstr.nnm-1 do begin
  nm=planstr.nm0+inm*planstr.dnm
  for iam=0,planstr.nam-1 do begin
   am=planstr.am0+iam*planstr.dam
   for ivt=0,planstr.nvt-1 do begin
    vt=10.^(planstr.vt0+ivt*planstr.dvt)
    if twod then nmh = planstr.nmh else nmh = 1
    for imh=0,nmh-1 do begin
      if twod then begin
        mh=[planstr.mh0+imh*planstr.dmh,0.,1]
        name='m'+cval(mh[0])+'a'+cval(am)+'c'+cval(cm)+'n'+cval(nm)+'v'+cval(vt) 
      endif else begin
        name='a'+cval(am)+'c'+cval(cm)+'n'+cval(nm)+'v'+cval(vt)
        mh=[planstr.mh0,planstr.dmh,planstr.nmh]
      endelse
    if keyword_set(vmicrofit) then begin
     vt=vmicro
     if vmsuffix eq '' then vmsuffix='vfit'
     name='a'+cval(am)+'c'+cval(cm)+'n'+cval(nm)+'_'+vmsuffix
    endif
    ;file=root[0]+'_'+name+'.par'
    file=file_dirname(planfile)+'/'+planstr.name+'_'+name+'.par'
    speclib_mkplan,file,$
       teff=[planstr.teff0,planstr.dteff,planstr.nteff],$
       logg=[planstr.logg0,planstr.dlogg,planstr.nlogg],$
       mh=mh,$
       am=[am,0.,1],nm=[nm,0.,1],cm=[cm,0,1],$
       rot=[planstr.rot0,planstr.drot,planstr.nrot],$
       continuum=continuum,$
       name=name,wrange=wrange,dw=dw,width=width,linelist=linelist,$
       vmicro=vt,fitvmicro=vmicrofit,synthcode=synthcode,solarisotopes=solarisotopes,atmos=atmos,$
       synthdir=synthdir,marcsdir=marcsdir,specdir=specdir,smooth=smooth,air=air,resolution=resolution,lsfid=lsfid,$
       waveid=waveid,lsffiber=lsffiber,vmacro=vmacro,fitvmacro=vmacrofit,config=config,twod=twod,kernel=kernel,apred_vers=apred_vers

    if nelem gt 0 then begin
     for i=0,n_elements(planstr.elem)-1 do begin
       file='apPlan-'+name+'_'+strtrim(planstr.elem[i],2)+'.par'
       file=file_dirname(planfile)+'/'+planstr.name+'_'+name+strtrim(planstr.elem[i],2)+'.par'
       speclib_mkplan,file,$
         teff=[planstr.teff0,planstr.dteff,planstr.nteff],$
         logg=[planstr.logg0,planstr.dlogg,planstr.nlogg],$
         mh=mh,$
         am=[am,0.,1],nm=[nm,0.,1],cm=[cm,0,1],$
         rot=[planstr.rot0,planstr.drot,planstr.nrot],$
         continuum=continuum,$
         name=name,$
         specdir=specdir,marcsdir=marcsdir,smooth=smooth,$
         wrange=wrange,dw=dw,width=width,linelist=linelist,$
         vmicro=vt,fitvmicro=vmicrofit,vmacro=vmacro,fitvmacro=vmacrofit,$
         synthcode=synthcode,solarisotopes=solarisotopes,atmos=atmos,$
         elem=strtrim(planstr.elem[i],2),maskdir=planstr.maskdir,$
         air=air,resolution=resolution,lsfid=lsfid,waveid=waveid,lsffiber=lsffiber,config=config,twod=twod,kernel=kernel,apred_vers=apred_vers
     endfor
    endif
    endfor
   endfor
  endfor
 endfor
endfor

; break up file for bundling/PCA parallel processing

tags=tag_names(planstr)
spawn,"sed '/elem/d' "+root[0]+'.par >temp.par'
for i=0,planstr.npart do begin
 file_copy,'temp.par',root[0]+'_'+string(format='(i2.2)',i)+'.par',/over
 openw,1,root[0]+'_'+string(format='(i2.2)',i)+'.par',/append
 printf,1,'dopart ',i
 close,1
 if nelem gt 0 then begin
  for ielem=0,n_elements(planstr.elem)-1 do begin
   file_copy,root[0]+'_'+string(format='(i2.2)',i)+'.par',root[0]+'_'+strtrim(planstr.elem[ielem],2)+'_'+string(format='(i2.2)',i)+'.par',/over
   openw,1,root[0]+'_'+planstr.elem[ielem]+'_'+string(format='(i2.2)',i)+'.par',/append
   printf,1,'elem ',planstr.elem[ielem]
   close,1
  endfor
 endif
endfor
file_delete,'temp.par'

end
