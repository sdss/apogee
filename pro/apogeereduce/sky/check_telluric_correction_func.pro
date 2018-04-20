function check_telluric_correction_func,fspec,best_std,mspec=mspec,telpix=telpix,x=x,stp=stp

min_scales = 0.1
max_scales = 10.0 ; 4.0

; Loop through the scaling factors and measure rms
nscales = 100 ; 50
scales = scale_vector(findgen(nscales),min_scales,max_scales)
stdarr = fltarr(nscales)
sigarr = fltarr(nscales)
for i=0,nscales-1 do begin

  test_telspec = fit_telluric(x,[scales[i],0.0,0.0, 1.0,0.0],mspec=mspec)
  fspec_corr = fspec / (test_telspec > 0.01)

  ; Measure RMS in the telluric-dominates pixels
  std = STDDEV(fspec_corr[telpix])
  sig = MAD(fspec_corr[telpix])
  stdarr[i] = std
  sigarr[i] = sig

endfor

;bestind_sig = first_el(minloc(sigarr))
;best_scale_sig = scales[bestind_sig]
;monte_str[count].best_scale_sig = best_scale_sig
;monte_str[count].best_sig = min(sigarr)

scales2 = scale_vector(findgen(1000),min_scales,max_scales)
stdarr2 = spline(scales,stdarr,scales2)

bestind_std = first_el(minloc(stdarr2))
best_scale_std = scales2[bestind_std]
best_std = min(stdarr2)

if keyword_set(stp) then stop

return,best_scale_std

end
