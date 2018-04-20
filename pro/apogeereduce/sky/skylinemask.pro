function skylinemask,flux,sky,nhigh,plot=plot,width=width,thresh=thresh

; procedure to identify pixels around significant sky lines 

if not keyword_set(width) then width=5
if not keyword_set(thresh) then thresh=0.5

; subtract median sky, and smooth
sm=smooth(sky-median(sky),width)
if keyword_set(plot) then begin
  set_plot,'X'
  plot,flux,xrange=[0,500],psym=10
  oplot,sky,color=2,psym=10
  oplot,sm,color=4,psym=10
endif

; identify all pixels with sky flux more than 0.5 of median filtered star flux
medflux=zap(reform(flux,n_elements(flux),1),[50,1])
bd=where(sm gt thresh*medflux,nhigh,complement=gd)
if keyword_set(plot) then begin
  sm[gd]=0
  sm[bd]=1
  oplot,sm*10000,color=5,psym=10
  stop
endif
return,bd
end

