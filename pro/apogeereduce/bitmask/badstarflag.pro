function badstarflag
 getstarflags,flags,badflags
 bad=0L
 for i=0,n_elements(badflags)-1 do if badflags[i] then bad=bad or 2^i

 return,bad
end

