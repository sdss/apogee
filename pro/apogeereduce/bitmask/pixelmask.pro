function pixelmask,mask

flag=''
getmaskvals,flags
for i=0,n_elements(flags)-1 do if is_bit_set(mask,i) eq 1 then flag=flag+flags[i]+','

lastcomma=strpos(flag,',',/reverse_search)
strput,flag,' ',lastcomma
return,flag

end

