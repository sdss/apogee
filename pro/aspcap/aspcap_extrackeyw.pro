function aspcap_extrackeyw,head,keyword
  
; This function find the value of keyword in the header
;
;
; INPUT:     
;   head          - strarr  image header
;   keyword       - string  the keyword to search
;
; OUTPUT:    
;           
;   value         - string  the value of keyword found in the header
;
;
; By Ana Elia Garcia Perez - March 2010
;
;


hobj=stregex(head,keyword,/SUBEXPR)
index=where(hobj ne -1,ncts)
if ncts gt 0 then begin
   value=STRTRIM(strmid(head(where(hobj ne -1)),10,71),2)
  dum=stregex(value(0),"/",/SUBEXPR)
  if dum ne -1 then value(0)=strtrim(strmid(value(0),0,dum))
endif else begin
   value="9d99"
endelse
return,reform(value(0))

end
