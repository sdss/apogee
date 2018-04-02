pro speclib_mkplan,planfile,teff=teff,logg=logg,mh=mh,am=am,cm=cm,nm=nm,rot=rot,name=name,wrange=wrange,dw=dw,width=width,elem=elem,linelist=linelist,continuum=continuum,air=air,solarisotopes=solarisotopes,vmicro=vmicro,synthcode=synthcode,atmos=atmos,synthdir=synthdir,specdir=specdir,smooth=smooth,resolution=resolution,lsfid=lsfid,waveid=waveid,lsffiber=lsffiber,vmacro=vmacro,config=config,twod=twod,kernel=kernel,apred_vers=apred_vers,fitvmicro=vmicrofit,fitvmacro=vmacrofit,maskdir=maskdir,marcsdir=marcsdir

if n_elements(planfile) eq 0 then begin
  print, 'Usage: speclib_mkplan,planfile'
  return
endif

; set defaults
if not keyword_set(name) then name='test'
if not keyword_set(teff) then teff=[4000.,250.,1]
if not keyword_set(logg) then logg=[2.,0.5,1]
if not keyword_set(mh) then mh=[0.,0.5,1]
if not keyword_set(am) then am=[0.,0.25,1]
if not keyword_set(cm) then cm=[0.,0.25,1]
if not keyword_set(nm) then nm=[0.,alog10(2.0),1]
if not keyword_set(rot) then rot=[0.,0.5,1]
if not keyword_set(wrange) then wrange=[15150.,17000.]
if not keyword_set(dw) then dw=0.1
if not keyword_set(linelist) then linelist='moog.201312111200.vac'
if not keyword_set(resolution) then resolution=22500.d0
if not keyword_set(continuum) then continuum=[4.,10.,0.1,3.]


; write plan file
openw,lun,planfile,/get_lun
printf,lun,'name ',name
printf,lun,'config  '+config
printf,lun,'synthcode  '+synthcode
if keyword_set(atmos) then printf,lun,'atmos     '+atmos
if keyword_set(solarisotopes) then begin
  if abs(solarisotopes) eq 1 then isodir='solarisotopes' else isodir='giantisotopes'
  if solarisotopes lt 0 then isodir='tests/'+isodir
endif
if keyword_set(synthdir) then printf,lun,'synthdir     '+synthdir
if keyword_set(marcsdir) then printf,lun,'marcsdir     '+marcsdir
if keyword_set(specdir) then printf,lun,'specdir     '+synthcode+'/'+atmos+'/'+isodir+'/'+specdir
if keyword_set(smooth) then printf,lun,'smooth     '+smooth
if keyword_set(solarisotopes) then printf,lun,'solarisotopes   ',solarisotopes
printf,lun,'wrange ',wrange
if dw lt 0 and vmicrofit eq 0 then begin
  if vmicro lt 4 then printf,lun,'dw 0.05' else printf,lun,'dw 0.1'
endif else printf,lun,'dw ',dw
if keyword_set(air) then printf,lun,'vacuum   0'
if keyword_set(width) then printf,lun,'width ',width
if keyword_set(lsfid) then printf,lun,'lsfid'+'   '+string(lsfid)
if keyword_set(apred_vers) then printf,lun,'apred_vers'+'   '+apred_vers
if keyword_set(waveid) then printf,lun,'waveid'+'   '+string(waveid)
if keyword_set(lsffiber) then printf,lun,'lsffiber ',string(format='(i5)',lsffiber)
if keyword_set(twod) then printf,lun,'twod     ',twod
printf,lun,'resolution ',resolution
printf,lun,'continuum ',string(format='(4f8.2)',continuum)
printf,lun,'linelist ',linelist
printf,lun,'teff0 ',teff[0]
printf,lun,'dteff ',teff[1]
printf,lun,'nteff ',nint(teff[2])
printf,lun,'logg0 ',logg[0]
printf,lun,'dlogg ',logg[1]
printf,lun,'nlogg ',nint(logg[2])
printf,lun,'mh0 ',mh[0]
printf,lun,'dmh ',mh[1]
printf,lun,'nmh ',nint(mh[2])
printf,lun,'am0 ',am[0]
printf,lun,'dam ',am[1]
printf,lun,'nam ',nint(am[2])
printf,lun,'cm0 ',cm[0]
printf,lun,'dcm ',cm[1]
printf,lun,'ncm ',nint(cm[2])
printf,lun,'nm0 ',nm[0]
printf,lun,'dnm ',nm[1]
printf,lun,'nnm ',nint(nm[2])
printf,lun,'rot0 ',rot[0]
printf,lun,'drot ',rot[1]
printf,lun,'nrot ',nint(rot[2])
if keyword_set(kernel) then printf,lun,'kernel   ',kernel
if keyword_set(vmicro) then printf,lun,'vmicro   ',vmicro,format='(a,8f10.5)'
if keyword_set(vmacro) then printf,lun,'vmacro   ',vmacro,format='(a,8f10.5)'
if keyword_set(vmicrofit) then printf,lun,'vmicrofit  ',vmicrofit
if keyword_set(vmacrofit) then printf,lun,'vmacrofit  ',vmacrofit
if keyword_set(elem) then printf,lun,'elem ',elem
if keyword_set(maskdir) then printf,lun,'maskdir ',maskdir
free_lun,lun

end



