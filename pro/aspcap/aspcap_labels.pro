pro aspcap_labels,filename,labels,slabels,npars,ps=ps,path=path

if n_elements(ps) gt 0 then begin
   ssmicro='!9' + String("170B) + '!X!lt!N'
   alpha='!9' + String("141B) + '!X'
endif else begin
   alpha='!4' + String("141B) + '!X'
   ssmicro='!4' + String("156B) + '!X!lt!N

endelse


value=file_test(filename+'.labl')

print,'Labels filename',filename+'.labl'
if value ne 0 then begin 
   load,filename+'.labl',rlabels
;   readcol,filename+'.labl',rlabels
   rlabels=''
   openr,num,filename+'.labl',/get_lun 
   readf,num,rlabels
   close,num & free_lun,num
   labels='l'
   num=strpos(rlabels,"'") & clabel=rlabels
   while num ge 0 do begin 
      clabel=strmid(clabel,num+1)
      num=strpos(clabel,"'")
      if num ge 0 then labels=[labels,"'"+strmid(clabel,0,num)+"'"]
      clabel=strmid(clabel,num+1) &  num=strpos(clabel,"'") 
   
   endwhile
   labels=labels[1:*]
;   npars= where(labels eq 'CPHOT')-1L
;   labels=rlabels[1:npars]
;   labels=rlabels[*,0]

   npars=n_elements(labels)
;   slabels=aspcap_strarr2str(rlabels[1:*])

endif else begin 
   pos=strpos(filename,'/',/reverse_search)
   input=strmid(filename,0,pos+1)+'input-'+strmid(filename,pos+1,strlen(filename)-pos-1)+'.nml'
   print,'input',input
   library=aspcap_getlibname(input)
   aspcap_getlabels,library,rlabels,path=path

   index=where(rlabels ne "")
   labels=rlabels[index]
   npars=n_elements(rlabels[index])
endelse

   index=where(labels eq "METALS" or labels eq "'METALS'",ncts)
   if ncts gt 0 then labels[index]="[Fe/H]"
   index=where(labels eq "TEFF" or labels eq "'TEFF'",ncts)
   if ncts gt 0 then labels[index]="Teff"
   index=where(labels eq "C" or labels eq "'C'",ncts)
   if ncts gt 0 then labels[index]='[C/Fe]'
   index=where(labels eq "N" or labels eq "'N'",ncts)
   if ncts gt 0 then labels[index]='[N/Fe]'

   index=where(labels eq "O Mg Si S Ca Ti" or labels eq "'O Mg Si S Ca Ti'",ncts)
   if ncts gt 0 then labels[index]='['+alpha+'/Fe]'
   index=where(labels eq "LOG10VDOP" or labels eq "'LOG10VDOP'",ncts)
   if ncts gt 0 then labels[index]=ssmicro
   index=where(labels eq "LOGG" or labels eq "'LOGG'",ncts)
   if ncts gt 0 then labels[index]='logg'
   llabels=labels
   print,'Npars',npars
   for i=0,npars-2 do begin 
     llabels[i]=llabels[i]+'/' 
   endfor
   slabels=aspcap_strarr2str(llabels)


end

