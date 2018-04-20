function is_bit_set, bitmasks, qbit
;===========
; Determines which given hexadecimal bitmasks have the Xth bit(s) set
; Input:
;  bitmask      Array of 32-bit integers, sum of 2^bits
;  qbit Bit # user wants to know the status of (yes/no)
; Output:
;  Returns array of 1's and 0's, corresponding to whether each bitmask has the given bit set
;   (if only querying 1 bitmask, result will be scalar)
; Example:
;  IDL> arr = is_bit_set([-2147483640, -2147352576], 17)
;  IDL> print, arr
;        0  1
;  To find, e.g., calibration cluster stars in data structure "data":
;  IDL> clust = where( is_bit_set(data.apogee_target2, 10) eq 1,nclust )
;
; GZ Jan 2012
;===========

nmasks = n_elements(bitmasks) & output=intarr(nmasks)

for i=0,nmasks-1 do begin
  bits=reverse(binary(long(bitmasks[i])))
  if bits[qbit] eq 1 then output[i] = 1
endfor

if nmasks eq 1 then return,(output)(0) else return,output
end
