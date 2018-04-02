function cval,x,turbo=turbo,one=one

if keyword_set(turbo) then begin
  if x lt 0 then prefix='-' else prefix='+' 
  if keyword_set(one) then return,prefix+string(format='(f3.1)',abs(x)) else $
  return,prefix+string(format='(f4.2)',abs(x))
endif else begin
  if x lt 0 then prefix='m' else prefix='p'
  return,prefix+string(format='(i2.2)',nint(abs(x)*10.))
endelse

end
