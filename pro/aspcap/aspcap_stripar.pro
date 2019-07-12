function aspcap_stripar,spm,label

; This function convert the scalar stellar parameter 
; into strings to add to the plots 


nobjects=n_elements(spm[0,*])
npar=n_elements(spm[*,0])
format=replicate('(F7.2)',npar)
sspm=transpose(replicate('S',nobjects))
for i=0,npar-1 do begin
  sspm=sspm+'/'+string(spm[i,*],format=format[i])
endfor
sspm=strmid(sspm,2,strlen(sspm)-1)

return,sspm
end
