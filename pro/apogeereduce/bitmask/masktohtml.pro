pro masktohtml,masktype,file=file,tex=tex

 if keyword_set(file) then openw,out,file,/get_lun else out=-1
 if masktype eq 0 then $
 getmaskvals,flag,badflag,maskcontrib,descrip else $
 getstarflags,flag,badflag,descrip

 if keyword_set(tex) then begin
   for i=0L,n_elements(flag)-1 do begin
     printf,out,flag[i]+'&'+string(format='(i)',i)+ '&'+ descrip[i]+ '\\
   endfor

 endif else begin
 printf,out,'<TABLE>'
 printf,out,'<TR>'
 printf,out,'<TD>Bit name</TD>'
 printf,out,'<TD>Binary digit</TD>'
 printf,out,'<TD>Description</TD>'
 printf,out,'</TR>'
 for i=0L,n_elements(flag)-1 do begin
   printf,out,'<TR>'
   printf,out,'<TD> '+flag[i]+'</TD>'
   printf,out,'<TD>'+ string(format='(i)',i)+ '</TD>'
   printf,out,'<TD>'+ descrip[i]+ '</TD>'
   printf,out,'</TR>'
 endfor
 printf,out,'</TABLE>'
 endelse

 if keyword_set(file) then free_lun,out

end
