function speclib_mkmoog,teff,logg,mh,am,cm,nm,wrange=wrange,dw=dw,width=width,atmod=atmod,elem=elem,save=save,linelist=linelist,vmicro=vmicro,solarisotopes=solarisotopes,marcs=marcs,nsyn=nsyn,noflux=noflux

if not keyword_set(wrange) then wrange=[15150.,17000]
if not keyword_set(dw) then dw=0.1
if not keyword_set(width) then width=7.0
if keyword_set(elem) then suffix=speclib_elem(abs(elem)) else suffix=''

; linelist
if not keyword_set(linelist) then linelist='moog.201312111200.vac'
linelistdir=getenv('APOGEE_REDUX')+'/speclib/linelists/'

; output file name
file='t'+strtrim(string(nint(teff)),2)+'g'+string(format='(i2.2)',nint(logg*10.))+$
  'm'+cval(mh)+'a'+cval(am)+'c'+cval(cm)+'n'+cval(nm)+'v'+cval(vmicro)+suffix


; atmosphere: note that all [N/M] is solar for atmospheres (since it doesn't have a big effect)
if keyword_set(marcs) then begin
  if logg le 3 then geo='s' else geo='p'
  modfile=geo+strtrim(string(nint(teff)),2)+'_g'+cval(logg,/turbo,/one)+'*'+'_z'+cval(mh,/turbo)+'_a'+cval(am,/turbo)+'_c'+cval(cm,/turbo)+'_n'+cval(0.,/turbo)+'_o'+cval(am,/turbo)

  moddir=getenv('APOGEE_REDUX')+'/speclib/marcs/MARCS_K_models/'
  print,' searching: '+moddir+modfile+'*'
  files=file_search(moddir+modfile+'*',count=count)
  script=getenv('SPECLIB_DIR')+'/scripts/marcs.awk'
endif else begin
  ; atmosphere: need to convert to MOOG format
  modfile='m'+cval(mh)+'c'+cval(cm)+'o'+cval(am)+'t'+strtrim(string(nint(teff)),2)+$
    'g'+string(format='(i2.2)',nint(logg*10.))+'v20'
  moddir=getenv('APOGEE_REDUX')+'/speclib/kurucz_filled/'+'m'+cval(mh)+'c'+cval(cm)+'o'+cval(am)+'/'
  print,' searching: '+moddir+'a'+modfile+'*'
  files=file_search(moddir+'a'+modfile+'*',count=count)
  script=getenv('SPECLIB_DIR')+'/scripts/makemoogmodel.awk'
endelse
if count eq 0 then begin
  stop,'no atmosphere found: ',modfile
endif else if count gt 1 then begin
  stop,'more than one atmosphere found',modfile
endif

print,' found: ', files

if keyword_set(vmicro) then $
spawn,'awk -f '+script+' vmicro='+string(format='(f3.1)',vmicro)+' '+files[0],out else $
spawn,'awk -f '+script+' '+files[0],out 
openw,lun,file+'.org',/get_lun
for i=0,n_elements(out)-1 do printf,lun,out[i]
free_lun,lun
if not keyword_set(atmod) then atmod=file+'.org'

; weed the line list
openw,lun,file+'.par',/get_lun
printf,lun,'weedout'
printf,lun,'terminal x11'
printf,lun,'plot 0'
printf,lun,'standard_out /dev/null'
printf,lun,'keeplines_out  '+file+'.lines.tmp'
printf,lun,'tosslines_out  /dev/null'
printf,lun,'keeplines_ratio 0.00001'
printf,lun,'summary_out /dev/null'
printf,lun,'smoothed_out /dev/null'
printf,lun,'damping 0'
printf,lun,'model_in '+atmod
printf,lun,'lines_in '+linelistdir+linelist
printf,lun,'atmosphere 1'
printf,lun,'molecules 2'
printf,lun,'lines 1'
printf,lun,'flux/int 0'
printf,lun,'synlimits'
printf,lun,string(wrange[0])+string(wrange[1])+string(dw)+string(width)
free_lun,lun
com=['time','MOOGSILENT',file+'.par']
spawn,com,/noshell
;spawn,'time MOOGSILENT '+file+'.par'
openw,lun,file+'.awk'
printf,lun,'BEGIN {'
printf,lun,'outfile = '+'"'+file+'.lines"'
printf,lun,'}'
printf,lun,'($1>'+string(wrange[0])+'&& $1<'+string(wrange[1])+') {print $0 > outfile}'
free_lun,lun
;spawn,'awk -f '+getenv('SPECLIB_DIR')+'/scripts/makemoogmodel.awk '+files[0],out
;com=['awk','-f',file+'.awk',file+'.lines.tmp','>',file+'.lines']
com=['awk','-f',file+'.awk',file+'.lines.tmp']
spawn,com,/noshell
;spawn,'awk -f'+file+'.awk '+file+'.lines.tmp '+'> '+file+'.lines'
print,'remove 1: ',file+'.awk',file+'.lines.tmp'
if not keyword_set(save) then file_delete,file+'.awk',file+'.lines.tmp'

; MOOG synth input file
for iflux=0,1 do begin
openw,lun,file+'.par',/get_lun
if iflux eq 0 then printf,lun,'doflux' else printf,lun,'synth'
printf,lun,'terminal x11'
printf,lun,'plot 0'
printf,lun,'standard_out /dev/null'
printf,lun,'summary_out '+file
printf,lun,'smoothed_out /dev/null'
printf,lun,'damping 0'
printf,lun,'strong 1'
printf,lun,'stronglines_in '+linelistdir+'/'+'stronglines.vac'
printf,lun,'model_in '+atmod
printf,lun,'lines_in '+file+'.lines'
printf,lun,'atmosphere 1'
printf,lun,'molecules 1'
printf,lun,'lines 1'
printf,lun,'flux/int 0'
isotopes=strarr(17,2)
isotopes[0,*]=['108.00116','1.001']
isotopes[1,*]=['606.01212','0.91']
isotopes[2,*]=['606.01213',8]
isotopes[3,*]=['606.01313','81']
isotopes[4,*]=['607.01214','0.91']
isotopes[5,*]=['607.01314','8']
isotopes[6,*]=['607.01215','273']
isotopes[7,*]=['608.01216','0.91']
isotopes[8,*]=['608.01316','8']
isotopes[9,*]=['608.01217','1101']
isotopes[10,*]=['608.01218','551']
isotopes[11,*]=['114.00128','1.011']
isotopes[12,*]=['114.00129','20']
isotopes[13,*]=['114.00130','30']
isotopes[14,*]=['101.00101','1.001']
isotopes[15,*]=['101.00102','1000']
isotopes[16,*]=['126.00156','1.00']
if keyword_set(solarisotopes) then begin
isotopes[0,*]=['108.00116','1.001']
isotopes[1,*]=['606.01212','1.01']
isotopes[2,*]=['606.01213','90']
isotopes[3,*]=['606.01313','180']
isotopes[4,*]=['607.01214','1.01']
isotopes[5,*]=['607.01314','90']
isotopes[6,*]=['607.01215','273']
isotopes[7,*]=['608.01216','1.01']
isotopes[8,*]=['608.01316','90']
isotopes[9,*]=['608.01217','1101']
isotopes[10,*]=['608.01218','551']
isotopes[11,*]=['114.00128','1.011']
isotopes[12,*]=['114.00129','20']
isotopes[13,*]=['114.00130','30']
isotopes[14,*]=['101.00101','1.001']
isotopes[15,*]=['101.00102','1000']
isotopes[16,*]=['126.00156','1.00']
endif
if keyword_set(elem) then begin
 if ~keyword_set(nsyn) then nsyn=6
 printf,lun,'isotopes 17 '+string(abs(nsyn))
 for i=0,16 do printf,lun,isotopes[i,0],' ',replicate(isotopes[i,1]+' ',abs(nsyn))
 ; need to specify N abundance since not in model!
 if elem eq 7 then begin
   printf,lun,'abundances    1      '+string(abs(nsyn))
   ;tmp=nm+[-0.75,-0.5,-0.25,0.25,0.5,0.75]
   if nsyn eq 6 then $
   tmp=nm+[-0.75,-0.5,-0.25,0.25,0.5,0.75] else if nsyn eq 7 then $
 ;  tmp=nm+[-0.75,-0.5,-0.25,0.00,0.25,0.5,0.75] else $
   tmp=nm+[-0.75,-0.5,-0.25,0.00,0.25,0.5,0.75] else if nsyn eq 2 then $
   tmp=nm+[0.00,0.1]  else $
   tmp=nm+[-0.3,-0.2,-0.1,0.00,0.1,0.2,0.3] 
   printf,lun,string(elem)+ string(format='(7f8.2)',tmp)
 endif else if elem gt 0 then begin 
   printf,lun,'abundances    2      '+string(abs(nsyn))
   if nsyn eq 6 then $
   printf,lun,string(elem)+ '  -0.75  -0.5 -0.25  0.25  0.5 0.75' else if nsyn eq 7 then $
   printf,lun,string(elem)+ '  -0.75  -0.5 -0.25  0.00  0.25  0.5 0.75' else if nsyn eq 2 then $
   printf,lun,string(elem)+ '   0.00  0.1' else $
   printf,lun,string(elem)+ '  -0.3  -0.2 -0.1  0.00  0.1  0.2 0.3' 
   printf,lun,'7 '+string(replicate(nm,abs(nsyn)),format='(7f8.2)')
 endif else begin
   alpha=[8,12,14,16,20,22]
   elem1=where(alpha ne abs(elem))
   printf,lun,'abundances         '+string(n_elements(elem1)+1)+'  '+string(abs(nsyn))
   for ii=0,n_elements(elem1)-1 do begin
     if nsyn eq 6 then $
     printf,lun,string(alpha[elem1[ii]])+ '   0.75   0.5  0.25  -0.25  -0.5 -0.75'  else if nsyn eq 7 then $
     printf,lun,string(alpha[elem1[ii]])+ '   0.75   0.5  0.25  0.00  -0.25  -0.5 -0.75' else $
     printf,lun,string(alpha[elem1[ii]])+ '   0.3   0.2  0.1  0.00  -0.1  -0.2 -0.3' 
   endfor
   printf,lun,'7 '+string(replicate(nm,abs(nsyn)),format='(7f8.2)')
 endelse
endif else begin
 nsyn=1
 printf,lun,'isotopes 17 1'
 for i=0,16 do printf,lun,isotopes[i,0]+' '+isotopes[i,1]
 ; need to specify N abundance since not in model!
 printf,lun,'abundances    1      1'
 printf,lun,'7 '+string(nm)
endelse
printf,lun,'obspectrum     5'
printf,lun,'plotpars      0'
printf,lun,'  15900.0  15940.0   0.40    1.02'
printf,lun,'    68     0.000   0.000    1.0000'
printf,lun,'   n        0.300    0.0     0.0     0.0     0.0'
printf,lun,'synlimits'
printf,lun,string(wrange[0])+string(wrange[1])+string(dw)+string(width)
free_lun,lun

; run MOOG!
com=['time','MOOGSILENT',file+'.par']
spawn,com,/noshell
;spawn,'time MOOGSILENT '+file+'.par'

; read in data and return it
if iflux eq 0 then readcol,file,w,flux else begin
  readmoog,file,w,data,nsyn=abs(nsyn)
  if ~keyword_set(noflux) then for isyn=0,abs(nsyn)-1 do data[*,isyn]*=flux
endelse

endfor
print,'cleaning: ',atmod,file,file+'.par',file+'.lines'
if not keyword_set(save) then file_delete,atmod,file,file+'.par',file+'.lines'
return,data

end
