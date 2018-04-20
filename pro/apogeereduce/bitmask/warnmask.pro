function warnmask

 getmaskvals,flags,badflags
 warn=0L
 for i=0,n_elements(badflags)-1 do if badflags[i] eq 0 then warn=warn or 2^i
 return,warn
end

