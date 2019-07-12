function aspcapflag,mask,type

if n_elements(type) eq 0 then type=0

flag=''
getaspcapflags,flags,badflag
for i=0,n_elements(flags)-1 do if is_bit_set(mask,i) eq 1 and $
  (type eq 0 or (badflag[i] eq type)) then flag=flag+flags[i]+','

lastcomma=strpos(flag,',',/reverse_search)
strput,flag,' ',lastcomma
return,strtrim(flag,2)

end

