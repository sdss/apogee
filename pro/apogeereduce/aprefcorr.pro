;+
;
; APREFCORR
;
; This corrects a raw APOGEE datacube for the reference pixels
; and reference output
;
; INPUTS:
;  cube       The raw APOGEE datacube with reference array.  This
;               will be updated with the reference subtracted cube.
;  head       The header for CUBE.
;  indiv=n    Subtract the individual reference arrays after nxn median filter. If 
;             If <0, subtract mean reference array. If ==0, no reference array subtraction
;  /noflip    Do not flip the reference array.
;  /silent    Don't print anything to the screen.
;
; OUTPUTS:
;  cube is updated with the reference subtracted cube to save memory.
;  mask       The flag mask.
;  =readmask  Mask indicating if reads are bad (0-good, 1-bad)
;
; USAGE:
;  IDL>aprefcorr,cube,head,mask
;
; By J. Holtzman   2011
; Incorporated into ap3dproc.pro  D.Nidever May 2011
;-

; refsub subtracts the reference array from each quadrant with proper flipping
pro aprefcorr_sub,image,ref
 revref=reverse(ref)
 image[0:511,*]-=ref
 image[512:1023,*]-=revref
 image[1024:1535,*]-=ref
 image[1536:2047,*]-=revref
 return
end

; refcorr does the "bias" subtraction, using the reference array and
;    the reference pixels. Subtract a mean reference array (or individual
;    with /indiv), then subtract vertical ramps from each quadrant using
;    reference pixels, then subtract smoothed horizontal ramps
function aprefcorr,cube,head,mask,indiv=indiv,vert=vert,horz=horz,noflip=noflip,silent=silent,readmask=readmask,lastgood=lastgood,cds=cds,plot=plot,fix=fix,q3fix=q3fix,keepref=keepref

ncube = n_elements(cube)
nhead = n_elements(head)

apgundef,reduced,mask

; Not enough inputs
if ncube eq 0 or nhead eq 0 then begin
  print,'Syntax - aprefcorr,cube,head,mask,indiv=indiv,vert=vert,horz=horz,noflip=noflip,silent=silent,readmask=readmask'
  return,0
endif

; Number of reads
sz = size(cube)
nread = sz[3]

; create long output
out = lonarr(2048,2048,nread)
if keyword_set(keepref) then refout = lonarr(512,2048,nread)

; Ignore reference array by default
if n_elements(indiv) eq 0 then indiv=3
; Default is to do CDS, vertical, and horizontal correction
if n_elements(cds) eq 0 then cds=1
if n_elements(vert) eq 0 then vert=1
if n_elements(horz) eq 0 then horz=1
print,'in aprefcorr, indiv: ',indiv

satval=55000L

snmin = 10
if keyword_set(indiv) then hmax=1e10 else hmax=65530

if n_elements(mask) le 1 then mask=intarr(2048,2048)
readmask = intarr(nread)
if not keyword_set(silent) then $
  print,'Calculating mean reference'
meanref = fltarr(512,2048) & nref=intarr(512,2048)
for i=0,nread-1 do begin
  ref = cube[2048:2559,*,i]

  m = MEAN(ref[128:511-128,128:2047-128],/double)
  s = stddev(ref[128:511-128,128:2047-128],/double)
  h = MAX(ref[128:511-128,128:2047-256])
  sat = where(ref ge satval,nsat)
  if nsat gt 0 then ref[sat] = !values.f_nan
  ; SLICE business is just for special fast handling, ignored if
  ;   not in header
  card = string(format='(a,i3.3)','SLICE',i)
  iread = sxpar(head,card,count=count)
  if count eq 0 then iread=i+1
  if not keyword_set(silent) then $
    print,format='(%"reading ref: %3d %3d\r",$)',i,iread
  ; skip first read and any bad reads
  if iread gt 1 and m/s gt snmin and h lt hmax then begin
    good = where(finite(ref))
    meanref[good] += (ref[good]-m)
    nref[good] += 1
    readmask[i] = 0
  endif else begin
    print,'Rejecting: ',i,m,s,h
    readmask[i] = 1
  endelse
endfor
meanref /= nref

if not keyword_set(silent) then $
  print,'Reference processing '
;reduced = fltarr(2048,2048,nread)

; Create vertical and horizontal ramp images
rows = findgen(2048)
cols = intarr(512)
cols += 1
vramp = (cols#rows)/2048
vrramp = 1-vramp
cols = findgen(2048)
rows = intarr(2048)
rows += 1
hramp = (cols#rows)/2048
hrramp = 1-hramp
clo = fltarr(2048) & chi=fltarr(2048)

if keyword_set(cds) then cdsref=cube[0:2047,*,1]

; Loop over the reads
lastgood=nread-1
for iread=0,nread-1 do begin

  ; Subtract mean reference array
  red = long(cube[0:2047,*,iread])

  sat = where(red gt satval,nsat)
  if nsat gt 0 then begin
    if iread eq 0 then nsat0=nsat
    red[sat] = 65535
    mask[sat] = (mask[sat] or maskval('SATPIX'))
    ; if we have a lot of saturated pixels, note this read (but don't do anything)
    if nsat gt nsat0+2000 then begin
      if lastgood eq nread-1 then lastgood=iread-1
    endif
  endif else nsat0=0
  ; pixels that are identically zero are bad, see these in first few reads
  bad = where(red eq 0,nbad)
  if nbad gt 0 then mask[bad] = (mask[bad] or maskval('BADPIX'))
  if not keyword_set(silent) then $
    print,format='(%"Ref processing: %3d  nsat: %5d\r",$)',iread+1,n_elements(sat)
  if readmask[iread] gt 0 then begin
    red = !values.f_nan
    goto,nextread
  endif

  ; with cds keyword, subtract off first read before getting reference pixel values
  if keyword_set(cds) then red-=cdsref

  ref = cube[2048:2559,*,iread]
  if indiv eq 1 then begin
    APREFCORR_SUB,red,ref
    ref-=ref
  endif else if indiv gt 1 then begin
    APREFCORR_SUB,red,median(ref,indiv)
    ref-=median(ref,indiv)
  endif else if indiv lt 0 then begin
    APREFCORR_SUB,red,meanref
    ref-=meanref
  endif

  if keyword_set(vert) then begin
   ; Subtract vertical ramp
   for j=0,3 do begin
    ;rlo=MEAN(red[j*512:(j+1)*512-1,1:3],/nan)
    ;rhi=MEAN(red[j*512:(j+1)*512-1,2044:2046],/nan)
    rlo = MEAN(red[j*512:(j+1)*512-1,2:3],/nan)
    rhi = MEAN(red[j*512:(j+1)*512-1,2045:2046],/nan)
    red[j*512:(j+1)*512-1,*] -= rlo*vrramp
    red[j*512:(j+1)*512-1,*] -= rhi*vramp
    if keyword_set(plot) then begin
      plot,rlo*vrramp[0,*]+rhi*vramp[0,*]
      print,j,rlo,rhi 
      atv,cube[0:2047,*,iread]-cube[0:2047,*,1]
      atv,red
      stop
    endif
   endfor
  endif

  ; Subtract smoothed horizontal ramp
  if keyword_set(horz) then begin
   for i=0,2047 do begin
    clo[i] = MEAN(red[1:3,i],/nan)
    chi[i] = MEAN(red[2044:2046,i],/nan)
   endfor
   ;clo = total(red[1:3,*],1,/nan) / ( total(finite(red[1:3,*]),1) > 1)
   ;chi = total(red[2044:2046,*],1,/nan) / ( total(finite(red[2044:2046,*]),1) > 1)

   sm=7
   ;slo = SMOOTH(clo,sm,/edge_truncate,/nan)
   ;shi = SMOOTH(chi,sm,/edge_truncate,/nan)
   slo = medfilt1d(clo,sm,/edge)
   shi = medfilt1d(chi,sm,/edge)

   if keyword_set(plot) then begin
    plot,slo
    oplot,shi
    ;stop
   endif
   if keyword_set(noflip) then begin
    red -= (rows#slo)*hrramp
    red -= (rows#shi)*hramp
   endif else begin
    ;bias = (rows#slo)*hrramp+(rows#shi)*hramp
    ; just use single bias value of minimum of left and right to avoid bad regions in one
    bias = rows#min([[slo],[shi]],dim=2)
    fbias = bias
    fbias[512:1023,*] = reverse(bias[512:1023,*])
    fbias[1536:2047,*] = reverse(bias[1536:2047,*])
    red -= fbias
   endelse
   if keyword_set(plot) then begin
    atv,red,min=-50,max=50
    stop
   endif
  endif
  if keyword_set(q3fix) then begin
    ;fix=red
    q3offset=fltarr(2048)
    for irow=0,2047 do begin
      q2m=median(red[923:1023,irow])
      q3a=median(red[1024:1124,irow])
      q3b=median(red[1435:1535,irow])
      q4m=median(red[1536:1636,irow])
      ;fix[1024:1535,irow]+=((q2m-q3a)+(q4m-q3b))/2.
      q3offset[irow]=((q2m-q3a)+(q4m-q3b))/2.
    endfor
    ;plot,q3offset
    ;oplot,medfilt1d(q3offset,7,/edge),color=2
    ;red=fix
    red[1024:1535,*]+=(medfilt1d(q3offset,7,/edge)##(fltarr(512)+1))
    ;atv,red,min=-200,max=200,/linear
    ;stop
  endif

  ; Make sure saturated pixels are set to 65535
  ;  removing the reference values could have
  ;  bumped them lower
  if nsat gt 0 then red[sat]=65535

  nextread:
  ;reduced[*,*,iread] = red
  ;cube[0:2047,*,iread] = red  ; overwrite with the ref-subtracted image
  out[*,*,iread] = red
  if keyword_set(keepref) then refout[*,*,iread] = ref

endfor ; read loop

; Trim off the reference array
;cube = cube[0:2047,*,*]

; mask the reference pixels
mask[0:3,*] = (mask[0:3,*] or maskval('BADPIX'))
mask[2044:2047,*] = (mask[2044:2047,*] or maskval('BADPIX'))
mask[*,0:3] = (mask[*,0:3] or maskval('BADPIX'))
mask[*,2044:2047] = (mask[*,2044:2047] or maskval('BADPIX'))

if not keyword_set(silent) then begin
  print,''
  print,'lastgood: ',lastgood
endif

if keyword_set(keepref) then return,[out,refout] else return,out

end
