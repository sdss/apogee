function aspcap_colorte,jk0,feh,dwarf=dwarf

; (J-K)_0 - Teff relation from Gonzalez Hernandez & Bonifacio (2009)
if keyword_set(dwarf) then begin
  b0=0.6524 & b1=0.5813 & b2=0.1225 & b3=-0.0646 & b4=0.0370 & b5=0.0016 ; dwarfs
endif else begin
  b0=0.6517 & b1=0.6312 & b2=0.0168 & b3=-0.0381 & b4=0.0256 & b5=0.0013 ; giants
endelse
theta=b0+b1*jk0+b2*jk0^2+b3*jk0*feh+b4*feh+b5*feh^2

return,5040./theta
end
 
