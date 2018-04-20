pro aploadraw,file,cube,head,silent=silent,error=error,real=real,stp=stp

;+
;
; APLOADRAW
;
; Load a bundled APOGEE datacube
; 
; INPUTS:
;  input   The bundles APOGEE datacube filename
;
; OUTPUTS:
;  cube    The datacube [2560,2048,Nreads]
;  head    The header from extension=0
;  =error  The error message if one occured.
;  /real   Change the array type to FLOAT instead of UINT.
;  /stp    Stop at the end of the program
;
; USAGE:
;  IDL>aploadraw,'apR-a-00000402.fits',cube,head
;
; By D.Nidever  Jan. 2011
;-

apgundef,cube,head

; Not enough inputs
if n_elements(file) eq 0 then begin
  error = 'Syntax - aploadraw,file,cube,head,error=error,real=real,stp=stp'
  print,error
  return
endif

; File doesn't exist
if file_test(file) eq 0 then begin
  error = file+' NOT FOUND'
  print,error
  return
endif

print,'Loading ',file

; Get header
head = headfits(file,errmsg=errmsg)
if errmsg ne '' then begin
  error = errmsg
  print,errmsg
  apgundef,cube,head
  return
end

; Figure out how many reads/extensions there are
;  the primary unit should be empty
nreads = sxpar(head,'NREAD',count=num_nread)
if num_nread eq 0 or nreads eq 0 then begin
  nreads = 0
  message = ''
  while (message eq '') do begin
    nreads++
    dum = headfits(file,exten=nreads,errmsg=message)
  end
  nreads--  ; removing the last one
end
print,'Nreads = ',strtrim(nreads,2)

head1 = headfits(file,exten=1)
nx = sxpar(head1,'NAXIS1')
ny = sxpar(head1,'NAXIS2')

; Initialize the cube
if keyword_set(real) then cube=fltarr(nx,ny,nreads) else $
  cube = uintarr(nx,ny,nreads)

; Reading in the data
for i=1,nreads do begin
  FITS_READ,file,im,exthead,exten=i,/no_abort,message=message

  if message ne '' then begin
    error = message
    apgundef,cube,head
    return
  endif

  if keyword_set(real) then im=float(im)   ; make it FLOAT
  cube[*,*,i-1] = im
end

if keyword_set(stp) then stop

end
