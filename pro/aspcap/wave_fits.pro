function wave_fits,hdr,log=log,ln=ln

wave=sxpar(hdr,'CRVAL1')+lindgen(sxpar(hdr,'NAXIS1'))*sxpar(hdr,'CDELT1')
if keyword_set(log) then wave=10^wave else if keyword_set(ln) then wave=exp(wave)
return,wave

end
