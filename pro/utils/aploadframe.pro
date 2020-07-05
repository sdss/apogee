pro aploadframe,input,outstr,exthead=exthead,error=error,stp=stp

;+
;
; APLOADFRAME
;
; This program loads the 3 chip files for one frame.
;
; INPUTS:
;  input    This can be either:
;             1.) The three chip files.  The chips are assumed
;                  to be in order, i.e. [chip1,chip2,chip3], OR
;             2.) A concatenate of the directory name, the appropriate
;                  suffix (i.e. ap1D or ap2D) and the frame ID.  For
;                  example, /net/stream/apogee/data/20110122/ap1D-00000303.
;  /exthead Extend the header arrays to 5000 elements.
;  /stp     Stop at the end of the program
;
; OUTPUTS:
;  outstr   A structure of the chip data and headers.
;  =error   The error message if one occurred.
;
; USAGE:
;  IDL>aploadframe,input,outstr
;
; By D.Nidever  March 2010
;-

apgundef,outstr
dirs=getdir()

; Not enough inputs
if n_elements(input) ne 1 and n_elements(input) ne 3 then begin
  print,'Syntax - aploadframe,input,outstr,exthead=exthead,error=error'
  error = 'Not enough inputs'
  return
endif

if size(input,/type) ne 7 then begin
  error = 'FILES be a 1- or 3-element STRING array'
  print,error
  return
endif

chiptag = ['a','b','c']

; Only one input, make chip files
;  Must be in format: /dirname/ap2D-00000303 or ap1D-00000303
if n_elements(input) eq 1 then begin

  dir = file_dirname(input)
  base = file_basename(input)
  arr = strsplit(base,'-',/extract)
  suffix = strtrim(arr[0],2)+'-'
  ;suffix = strmid(base,0,5)
  if suffix ne dirs.prefix+'1D-' and suffix ne dirs.prefix+'2D-' then begin
    print,'Single input must be in the format /dirname/ap2D-XXXXXXXX or /dirname/ap1D-XXXXXXXX'
    return
  endif
  ;frameid = strmid(base,5)
  frameid = strtrim(arr[1],2)
  files = dir+'/'+suffix+chiptag+'-'+frameid+'.fits'

; 3-element input.
endif else begin
  files = input
endelse

; Get file info
info = APFILEINFO(files,/silent)

bd = where(info.exists eq 0,nbd)
if nbd gt 0 then begin
  error = strjoin(files[bd],' ')+' NOT FOUND'
  print,error
  return
endif


; Loop through the chips
For i=0,2 do begin

  ; Load the header
  head = headfits(files[i],errmsg=message0)

  ; Load the flux
  if info[i].fpack then begin
    flux=mrdfits(files[i],1)
    message1=''
    err=mrdfits(files[i],2)
    message2=''
    mask=mrdfits(files[i],3)
    message3=''
  endif else begin
    FITS_READ,files[i],flux,head1,message=message1,/no_abort,exten=1  ; flux
    FITS_READ,files[i],err,head2,message=message2,/no_abort,exten=2   ; error
    FITS_READ,files[i],mask,head3,message=message3,/no_abort,exten=3  ; mask
  endelse

  ; set bad error values
  bderr=where(err le 0,nbd)
  if nbd gt 0 then err[bderr]=baderr()

  ; Load wavelengths
  if sxpar(head,'WAVEFILE') then begin
    if info[i].fpack then wave=mrdfits(files[i],4) $
    else FITS_READ,files[i],wave,head4,message=message4,/no_abort,exten=4  ; wave
    if sxpar(head,'WAVEHDU') then begin
      whdu=sxpar(head,'WAVEHDU')
      if info[i].fpack then wcoef=mrdfits(files[i],whdu) $
      else FITS_READ,files[i],wcoef,head4,message=message4,/no_abort,exten=whdu  ; wave
    endif
  endif else message4=''

  ; Errors occurred
  if message0+message1+message2+message3+message4 ne '' then begin
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

  ; Make the chip structure
  if n_elements(wave) gt 0 then begin
    if n_elements(wcoef) gt 0 then $
    (SCOPE_VARFETCH('chip'+chiptag[i],/enter)) = {filename:files[i],header:head,flux:flux,err:err,mask:mask,wavelength:wave,wcoef:wcoef} else $
    (SCOPE_VARFETCH('chip'+chiptag[i],/enter)) = {filename:files[i],header:head,flux:flux,err:err,mask:mask,wavelength:wave}
  endif else begin
    (SCOPE_VARFETCH('chip'+chiptag[i],/enter)) = {filename:files[i],header:head,flux:flux,err:err,mask:mask}
  endelse

Endfor

; Combine the data
outstr = {chipa:chipa,chipb:chipb,chipc:chipc}

;; Normal structure
;if not keyword_set(exthead) then begin
;  outstr = {chipa:{header:head1,data:data1,filename:files[0]},$
;            chipb:{header:head2,data:data2,filename:files[1]},$
;            chipc:{header:head3,data:data3,filename:files[2]}}
;
;; Extend the header arrays to 5000 elements
;endif else begin
;  outstr = {chipa:{header:strarr(5000),data:data1,filename:files[0]},$
;            chipb:{header:strarr(5000),data:data2,filename:files[1]},$
;            chipc:{header:strarr(5000),data:data3,filename:files[2]}}
;  nh1=n_elements(head1) & outstr.chipa.header[0:nh1-1]=head1
;  nh2=n_elements(head2) & outstr.chipb.header[0:nh2-1]=head2
;  nh3=n_elements(head3) & outstr.chipc.header[0:nh3-1]=head3
;endelse

if keyword_set(stp) then stop

;apgundef,data1,head1,data2,head2,data3,head3
apgundef,flux,head1,var,head2,mask,head3

end
