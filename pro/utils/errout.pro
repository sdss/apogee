; errout sets the value to output for the error for bad pixels
function errout,x

bd=where(x eq baderr() or x le 0 or finite(x) eq 0,nbd)
if nbd gt 0 then x[bd]=baderr()
return,x
end

