pro smcolor,ps=ps

if keyword_set(ps) then begin
  ;tvlct,255,255,255,0
  tvlct,0,0,0,0
  tvlct,0,0,0,1
endif else begin
  tvlct,0,0,0,0
  tvlct,255,255,255,1
endelse
tvlct,255,0,0,2
tvlct,0,255,0,3
tvlct,0,0,255,4
tvlct,255,0,255,5
tvlct,0,255,255,6
tvlct,255,255,0,7
end
