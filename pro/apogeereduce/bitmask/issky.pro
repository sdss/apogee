function issky,p,s
if is_bit_set(s,flagnum('APOGEE_SKY')) then return,1
return,0
end

