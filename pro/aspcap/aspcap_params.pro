function aspcap_params,tagnames,flagnames,npar=npar,extra=extra

if n_elements(npar) eq 0 then npar=0

tagnames=strupcase(['TEFF','LOGG','LOGVMICRO','M_H','C_M','N_M','ALPHA_M','LGVSINI','PARAM_O'])
flagnames=strupcase(['TEFF','LOGG','VMICRO','M_H','C_M','N_M','ALPHA_M','VSINI','O'])
params=['TEFF','LOGG','LOG10VDOP','METALS','C','N','O Mg Si S Ca Ti','LGVSINI','O']
if keyword_set(extra) then begin
  tagnames=[tagnames,strupcase(extra)]
  flagnames=[flagnames,strupcase(extra)+'_M']
  params=[params,extra]
endif
return,params

if npar eq 8 then begin
 tagnames=strupcase(['TEFF','LOGG','LOGVMICRO','PARAM_M_H','C_M','N_M','PARAM_ALPHA_M','LGVSINI'])
 flagnames=strupcase(['TEFF','LOGG','VMICRO','M_H','C_M','N_M','ALPHA_M','VSINI'])
 params=['TEFF','LOGG','LOG10VDOP','METALS','C','N','O Mg Si S Ca Ti','LGVSINI']
endif else if npar eq 0 or npar eq 7 then begin
 tagnames=strupcase(['TEFF','LOGG','LOGVMICRO','PARAM_M_H','C_M','N_M','PARAM_ALPHA_M'])
 flagnames=strupcase(['TEFF','LOGG','VMICRO','M_H','C_M','N_M','ALPHA_M'])
 ;tagnames=strupcase(['TEFF','LOGG','LOGVMICRO','METALS','CFE','NFE','ALPHAFE'])
 params=['TEFF','LOGG','LOG10VDOP','METALS','C','N','O Mg Si S Ca Ti']
endif else stop,'unknown value of npar in aspcap_params'
if keyword_set(elem) then begin
  params=[params,elem]
  tagnames=[tagnames,strupcase(elem)]
  flagnames=[flagnames,strupcase(elem)+'_H']
endif
return,params
;return,['TEFF','LOGG','LOGVT','METALS','C','N','O Mg Si S Ca Ti']
end
 
