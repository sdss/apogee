function telluricmask,tell,nmask,plot=plot,width=width

if not keyword_set(width) then width=5

sm=smooth(tell,width)
if keyword_set(plot) then begin
  set_plot,'X'
  plot,tell,color=2,psym=10
  oplot,sm,color=4,psym=10
endif
bd=where(sm lt 0.9,nmask,complement=gd)
if keyword_set(plot) then begin
  sm[gd]=0
  sm[bd]=1
  oplot,sm*10000,color=5,psym=10
  stop
endif
return,bd
end

