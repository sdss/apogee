function badparamflag
 getparamflags,flags,badflags
 bad=0L
 for i=0,n_elements(badflags)-1 do if badflags[i] eq 1 then bad=bad or 2L^i

 return,bad
end
