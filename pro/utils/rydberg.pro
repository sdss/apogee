function rydberg,n1,n2
r=1.0973731568d7
me=9.109382d-31
mprot=1.672621e-27
rm=r/(1+me/mprot)

l=rm*(1.d0/n1^2-1.d0/n2^2)
w=1./l*1.e10
;  print,n2,1./l*1.e10
;  oplot,[w,w],[0.,2.],color=4
return,w
end
