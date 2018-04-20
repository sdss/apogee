pro aploadcframe,files,outstr,plugmap,exthead=exthead,error=error,stp=stp

;+
;
; APLOADCFRAME
;
; This program loads the 3 chip apCframe files
;
; INPUTS:
;  files    The three apCframe chip files.  The chips are
;             assumed to be in order, i.e. [chip1,chip2,chip3].
;  /exthead Extend the header arrays to 5000 elements.
;  /stp     Stop at the end of the program
;
; OUTPUTS:
;  outstr   A structure of the chip data and headers.
;  plugmap  The plugmap structure
;  =error   The error message if one occurred.
;
; USAGE:
;  IDL>aploadcframe,files,outstr,plugmap
;
; By D.Nidever  July 2010
;-

apgundef,outstr,plugmap

; Not enough inputs
if n_elements(files) ne 3 then begin
  print,'Syntax - aploadcframe,files,outstr,plugmap,exthead=exthead,error=error'
  error = 'Not enough inputs'
  return
endif

if size(files,/type) ne 7 then begin
  error = 'FILES must be a 3-element STRING array'
  print,error
  return
endif

; Get file info
info = APFILEINFO(files,/silent)

bd = where(info.exists eq 0,nbd)
if nbd gt 0 then begin
  error = strjoin(files[bd],' ')+' NOT FOUND'
  print,error
  return
endif

chiptag = ['a','b','c']

; Loop through the chips
For i=0,2 do begin

  apgundef,head,flux,err,mask,wavelength,sky,skyerr
  apgundef,telluric,telerr,wcoef,lsfcoef

  ; Load the header
  head = headfits(files[i],errmsg=message0)

  ; Load the flux
  FITS_READ,files[i],flux,head1,message=message1,/no_abort,exten=1       ; flux
  FITS_READ,files[i],err,head2,message=message2,/no_abort,exten=2        ; error
  FITS_READ,files[i],mask,head3,message=message3,/no_abort,exten=3       ; mask
  FITS_READ,files[i],wavelength,head4,message=message4,/no_abort,exten=4 ; wavelength
  FITS_READ,files[i],sky,head5,message=message5,/no_abort,exten=5        ; sky
  FITS_READ,files[i],skyerr,head6,message=message6,/no_abort,exten=6     ; skyerr
  FITS_READ,files[i],telluric,head7,message=message7,/no_abort,exten=7   ; telluric
  FITS_READ,files[i],telerr,head8,message=message8,/no_abort,exten=8     ; telerr
  FITS_READ,files[i],wcoef,head9,message=message9,/no_abort,exten=9      ; wcoef
  FITS_READ,files[i],lsfcoef,head10,message=message10,/no_abort,exten=10 ; lsfcoef
  plugdata = MRDFITS(files[i],11,head11,status=status,/silent)            ; plugmap
  message11 = strtrim(status,2)
  if message11 eq '0' then message11=''  ; status=0 is success
  plughdr = MRDFITS(files[i],12,head12,status=status,/silent)            ; plugmap
  tellstar = MRDFITS(files[i],13,head13,status=status,/silent)            ; plugmap

  ; Errors occurred
  if message0+message1+message2+message3+message4+message5+message6+message7+$
       message8+message9+message10+message11 ne '' then begin
    error = 'ERROR loading file '+files[i]
    print,error
    return
  endif


  ; Extend the header to 5000 elements
  if keyword_set(exthead) then begin
    head0 = head
    head = strarr(5000)
    nhead = n_elements(head0)
    head[0:nhead-1] = head0
  endif

  ; Construct the original filename

  ; Make the chip structure
  (SCOPE_VARFETCH('chip'+chiptag[i],/enter)) = {filename:files[i],header:head,flux:flux,err:err,mask:mask,$
               wavelength:wavelength,sky:sky,skyerr:skyerr,telluric:telluric,telluricerr:telerr,wcoef:wcoef,$
               lsfcoef:lsfcoef}

End

; Combine the chip data
;-----------------------
outstr = {chipa:chipa,chipb:chipb,chipc:chipc,tellstar:tellstar}

;; Reconstruct the plugmap structure
;;-----------------------------------
;; plugmap = {HDR, fiberdata, other plate tags}
;; Get the original plPlugMap header lines
;plind1 = where(stregex(head11,'HISTORY PLMAPHD: ',/boolean) eq 1,nplind1)
;if nplind1 gt 0 then pluglines1 = head11[plind1] else pluglines1=['']
;plughist = strmid(pluglines1,17)  ; remove prefix
;plugmap = {hdr:plughist,fiberdata:plugdata}
;; Get the plate tags
;plind2 = where(stregex(head11,'HISTORY PLMAPTG: ',/boolean) eq 1,nplind2)
;if nplind2 gt 0 then pluglines2 = head11[plind2]
;for i=0,nplind2-1 do begin
;  line = strmid(pluglines2[i],17)                  ; remove prefix
;  arr = strsplit(line,'=',/extract)                ; get the keyword name
;  name = strtrim(arr[0],2)
;  ; Create dummy header for SXPAR
;  hdline = 'DUM     = '+arr[1]
;  len = strlen(hdline)
;  hdline = hdline+string(replicate(32B,80-len))     ; make it 80 chars long
;  dumhd = [hdline,'END'+string(replicate(32B,77))]  ; create dummy header array
;  value = sxpar(dumhd,'DUM')                        ; use sxpar to get the value
;  ; Add this tag to the plugmap structure
;  plugmap = CREATE_STRUCT(plugmap,name,value)  ; add to the plugmap
;end
plugmap = {hdr:plughdr,fiberdata:plugdata}

if keyword_set(stp) then stop

apgundef,head,flux,err,mask,wavelength,sky,skyerr
apgundef,telluric,telerr,wcoef,lsfcoef

end
