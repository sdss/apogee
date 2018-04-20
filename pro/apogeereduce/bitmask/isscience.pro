function isscience,p,s
if is_bit_set(p,flagnum('APOGEE_SHORT')) then return,1
if is_bit_set(p,flagnum('APOGEE_MEDIUM')) then return,1
if is_bit_set(p,flagnum('APOGEE_LONG')) then return,1
if is_bit_set(p,flagnum('APOGEE_ANCILLARY')) then return,1
if is_bit_set(p,flagnum('APOGEE_1MTARGET')) then return,1
return,0
end
