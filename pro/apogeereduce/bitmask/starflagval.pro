function starflagval,flag

 getstarflags,flagvals
 j=where(strtrim(flagvals,2) eq strtrim(flag,2),nj)
 if nj gt 0 then val=2L^j[0] else val=0

 if val eq 0 then begin
   print, 'Undefined mask: ',flag
   stop
 endif
 return,val
end

