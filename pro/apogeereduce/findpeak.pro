function findpeak,y,npk,level=level
d0 = y - shift(y,1)
bad=where(finite(d0) eq 0)
if bad[0] ge 0 then begin
  d0a=y-shift(y,2)
  d0[bad]=d0a[bad]
endif
d1 = y - shift(y,-1)
bad=where(finite(d1) eq 0)
if bad[0] ge 0 then begin
  d1a=y-shift(y,-2)
  d1[bad]=d1a[bad]
endif
pk = where(d0 gt 0 and d1 gt 0, npk)
;pk = where((d0 gt 0 or not finite(d0)) and (d1 gt 0 or not finite(d1)),npk)
if keyword_set(level) and npk gt 0 then begin
    bigind = where(y[pk] gt level, npk)
    if npk gt 0 then pk=pk[bigind] else pk=-1
    ;return,pk[bigind]
endif
return,pk
end
