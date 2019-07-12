function aspcap_strarr2str,var

; This procedure converts an string array in a scalar string with each
; element after the other, followed by a blank space. 
;
; INPUT:
;       var       -strarr     string array
;
; OUTPUT:
;      ovar       -string     scalar string
;
;
; BY ANA ELIA GARCIA PEREZ - JUNE 2011

; 


n=n_elements(var)

ovar=var[0]
print,n
for i=1,n-1 do begin 
   ovar=ovar+' '+var[i]
   ;stop
endfor
return,ovar
end
