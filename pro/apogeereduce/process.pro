; Main routine for 3D and 2D processing of APOGEE data
; contains
;  maskval
;  getcmjd
;  approcess  : top level processing routine
;               processes list of frames by chip, reading in
;               cal products for each chip, then looping over
;               the frames
;    process  : processes a single frame from 3D -> 2D
;    refcorr  : does the reference pixel correction
;      refsub
;    gain
;    rn
;    crrej      : CR rejection routine
;    fsamp      : Fowler sampling (or up-the-ramp) calculation
;    extract    : 2D -> 1D extraction routine
;      pmul

function getnum,mjd
 return,mjd-55562L
end

function getid,mjd,frame
 return,getnum(mjd)*10000L+frame
; id=string(getnum(mjd)*10000+frame)
; return,strtrim(id,2)
end

; refsub subtracts the reference array from each quadrant with proper flipping
pro refsub,image,ref
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
pro refcorr,cube,head,nread,reduced,mask,indiv=indiv,noflip=noflip

snmin=10
if keyword_set(indiv) then hmax=1e10 else hmax=65530

readmask=intarr(nread)
print,'calculating mean reference'
meanref=fltarr(512,2048) & nref=intarr(512,2048)
for i=0,nread-1 do begin
 ref=cube[2048:2559,*,i]

 m=MEAN(ref[128:511-128,128:2047-128],/double)
 s=stddev(ref[128:511-128,128:2047-128],/double)
 h=MAX(ref[128:511-128,128:2047-256])
 sat=where(ref ge 65535)
 if sat[0] ge 0 then ref[sat]=!values.f_nan
 ; SLICE business is just for special fast handling, ignored if
 ;   not in header
 card=string(format='(a,i3.3)','SLICE',i)
 iread=sxpar(head,card,count=count)
 if count eq 0 then iread=i+1

 ;print,format='(%"reading ref: %3d %3d\r",$)',i,iread
; skip first read and any bad reads
 if iread gt 1 and m/s gt snmin and h lt hmax then begin
   good=where(finite(ref))
   meanref[good]+=(ref[good]-m)
   nref[good]+=1
   readmask[i]=0
 endif else begin
   print,'rejecting: ',i,m,s,h
   readmask[i]=1
 endelse
endfor
meanref/=nref

print,'reference processing '
reduced=fltarr(2048,2048,nread)

; create vertical and horizontal ramp images
rows=findgen(2048)
cols=intarr(512)
cols+=1
vramp=(cols#rows)/2048
vrramp=1-vramp
cols=findgen(2048)
rows=intarr(2048)
rows+=1
hramp=(cols#rows)/2048
hrramp=1-hramp
clo=fltarr(2048) & chi=fltarr(2048)

; loop over the reads
for iread=0,nread-1 do begin
 ; subtract mean reference array
 red=cube[0:2047,*,iread]
 sat=where(red gt 65534)
 if sat[0] ge 0 then begin
   red[sat]=!values.f_nan
   mask[sat]=(mask[sat] or maskval('SATPIX'))
 endif
 ; print,format='(%"ref processing: %3d  nsat: %5d\r",$)',iread,n_elements(sat)
 if readmask[iread] gt 0 then begin
   red=!values.f_nan
   goto,nextread
 endif
 if keyword_set(indiv) then begin
   ref=cube[2048:2559,*,iread]
   refsub,red,ref
 endif else refsub,red,meanref

 ; subtract vertical ramp
 for j=0,3 do begin
  ;rlo=MEAN(red[j*512:(j+1)*512-1,1:3],/nan)
  ;rhi=MEAN(red[j*512:(j+1)*512-1,2044:2046],/nan)
  rlo=MEAN(red[j*512:(j+1)*512-1,2:3],/nan)
  rhi=MEAN(red[j*512:(j+1)*512-1,2045:2046],/nan)
  red[j*512:(j+1)*512-1,*]-=rlo*vrramp
  red[j*512:(j+1)*512-1,*]-=rhi*vramp
 endfor

 ; subtract smoothed horizontal ramp
 for i=0,2047 do begin
   clo[i]=MEAN(red[1:3,i],/nan)
   chi[i]=MEAN(red[2044:2046,i],/nan)
 endfor

 slo=SMOOTH(clo,50,/edge_truncate,/nan)
 shi=SMOOTH(chi,50,/edge_truncate,/nan)

 if keyword_set(noflip) then begin
   red-=(rows#slo)*hrramp
   red-=(rows#shi)*hramp
 endif else begin
   bias=(rows#slo)*hrramp+(rows#shi)*hramp
   fbias=bias
   fbias[512:1023,*]=reverse(bias[512:1023,*])
   fbias[1536:2047,*]=reverse(bias[1536:2047,*])
   red-=fbias
 endelse

 nextread:
 reduced[*,*,iread]=red
endfor

end

; crrej does CR rejection from up the ramp sampling
pro crrej,red,nread,mask,gain,rn,nofix=nofix,mean=mean,ncr=ncr

gainavg=mean(gain)
rnavg=mean(rn)
if nread lt 4 then begin
  print,'cant do CR rejection with ', nread, ' reads'
endif
sz=size(red)
rate=reform(red,sz[1]*sz[2],nread)
red=reform(red,sz[1]*sz[2],nread,/overwrite)
reject=10^2
; calulate all of the rates between pairs of images
for i=nread-1,1,-1 do rate[*,i]-=rate[*,i-1] 
ncr=0
for i=1,nread-1 do begin
  ; get median rate of 7 reads
  if keyword_set(mean) then begin
    i1=i-3 lt 0?0:i-3
    i2=i-1 lt 0?0:i-1
    n=i2-i1+1
    if n gt 1 then mrate=total(rate[*,i1:i2],2) else mrate=rate[*,i1:i2]
    i1=i+1 gt nread-1?nread-1:i+1
    i2=i+3 gt nread-1?nread-1:i+3
    if (i2-i1+1) gt 1 then mrate+=total(rate[*,i1:i2],2) else mrate+=rate[*,i1:i2]
    n+=(i2-i1+1)
    mrate/=n
  endif else begin
    i1=i-3 lt 0?0:i-3
    i2=i+3 gt nread-1?nread-1:i+3
    window=rate[*,i1:i2]
    mrate=median(window,dim=2)
  endelse
  diff=(rate[*,i]-mrate)
  ; reject pixels where rate is large compared to expected variance
  cr=where(finite(diff) and diff gt 0 and diff^2 gt reject*apvariance(mrate,2,gainavg,rnavg))
  print,format='(%"CR rejection: %3d  %5d\r",$)',i,n_elements(cr)
  ; if nofix, set rejected pixels to NaN for this and all subsequent reads, i.e.
  ;   no CR fixing is done
  if cr[0] ge 0 then begin
    ncr+=n_elements(cr)
    mask[cr]=(mask[cr] or maskval('CRPIX'))
    for j=i,nread-1 do begin
      if keyword_set(nofix) then begin
        red[cr,j]=!values.f_nan
      endif else begin
        red[cr,j]-=diff[cr]
      endelse
    endfor
  endif
endfor
red=reform(red,sz[1],sz[2],nread,/overwrite)
window=0
rate=0

end

; fsamp does Fowler or up-the-ramp sampling
function fsamp,red,nfs,var,gain,rn

sz=size(red)
nread=sz[3]
if n_elements(gain) eq 0 then gain=fltarr(4)+2.
if n_elements(rn) eq 0 then rn=fltarr(4)+12.

if nfs eq 0 then begin
  ; up the ramp
  sumts = fltarr(sz[1],sz[2])   ; SUM t*s
  sums = fltarr(sz[1],sz[2])    ; SUM s
  sum = intarr(sz[1],sz[2])    ; SUM s
  sumt = fltarr(sz[1],sz[2])   ; SUM t*s
  sumt2 = fltarr(sz[1],sz[2])   ; SUM t*s
  for i=1L,nread-1 do begin
    slice=red[*,*,i]
    if size(slice,/type) eq 3 then $
      good=where(slice ne 2L^31) else $
      good=where(finite(slice))
    if good[0] ge 0 then begin
      sumts[good] += i*reform(slice[good])
      sums[good] += reform(slice[good])
      sum[good] += 1
      sumt[good]+=i
      sumt2[good]+=i^2
    endif
  end
  ;sumt = total(findgen(nread))     ; SUM t
  ;sumt2 = total(findgen(nread)^2)  ; SUM t^2
  ; The slope in Counts per read, similar to med_dCounts_im
  ;slope = (nread*sumts - sumt*sums)/(nread*sumt2 - sumt^2)
  slope = (sum*sumts - sumt*sums)/(sum*sumt2 - sumt^2)
  ; To get the total counts just multiply by nread
  im = slope * (nread-2L)
  ; the first read doesn't really add any signal, just a zero-point
  var=0.*im
  for iquad=0,3 do begin
    i1=iquad*512 & i2=i1+511
    var[i1:i2,*]=12.*(nread-1.)/(nread*(nread+1.))*rn(iquad)^2+$
                  6.*(nread^2+1)/(5*nread*(nread+1))*im[i1:i2,*]*gain(iquad)
    var[i1:i2,*]/=(gain(iquad)^2)
  endfor

  return,im
endif else begin
  ; fowler sampling
  if nfs gt nread/2 then nfs=nread/2
  sread=red[*,*,1]
  for i=2,nfs do sread+=red[*,*,i]
  eread=red[*,*,nread-nfs]
  for i=nread-nfs+1,nread-1 do eread+=red[*,*,i]
  var=0.*sread
  for iquad=0,3 do begin
    i1=iquad*512 & i2=i1+511
    var[i1:i2,*]=2*rn(iquad)^2/nfs+(eread[i1:i2,*]-sread[i1:i2,*])/nfs
    var[i1:i2,*]/=(gain(iquad)^2)
  endfor
  return,(eread-sread)/nfs
endelse

end

; process does all the processing for a single image (single chip)
;  does dark subtraction and flat field if given dark and flat frames
;  uses step= for sparse up-the-ramp sampling if desired
;  nfs=0 by default (up the ramp)
;  turn off bias processing with /noref
;  turn off CR rejection with /nocr
;  turn off 3D->2D with /nofs (returns DCS image)
function process,cmjd,num,chip,head,red,dark,flat,err,gain,rn,mask,step=step,$
     nfs=nfs,noref=noref,nocr=nocr,nofs=nofs,stp=stp,nofix=nofix,ncr=ncr,indiv=indiv,horz=horz,vert=vert,maxread=maxread,nread=nread

if not keyword_set(step) then step=1
step=0
if not keyword_set(nfs) then nfs=0
if not keyword_set(nofix) then nofix=0
if not keyword_set(ncr) then ncr=0
if n_elements(mask) le 1 then mask=bytarr(2048,2048)
if n_elements(gain) eq 0 then gain=2+fltarr(4)
if n_elements(rn) eq 0 then rn=12+fltarr(4)
if n_elements(indiv) eq 0 then indiv=3

dirs=getdir()
datadir=dirs.datadir
dir=datadir+cmjd+'/'
time0=systime(/seconds)
file=dirs.prefix+string(format='("R-",a,"-",i8.8)',chip,num)
; if another process is uncompressing, wait
while file_test(dir+file+'.apz.lock') do wait,10

fitsdir=getenv('APOGEE_LOCALDIR')+'/'
if fitsdir eq '' then fitsdir=dir+'/'
if not file_test(dir+file+'.fits') then apunzip,dir+file+'.apz',/no_checksum,fitsdir=fitsdir
print,'reading '+file
aploadraw,fitsdir+file+'.fits',cube,head,/real ;,step=step
file_delete,fitsdir+file+'.fits',/quiet

if keyword_set(maxread) then cube=cube[*,*,0:maxread-1]

sz=size(cube)
if sz[0] eq 0 then return,0
nread=sz[3]
print,'nread ',nread
print,'time: ',systime(/seconds)-time0

if not keyword_set(noref) then begin
  print,'bias/reference correcting...'
  ;refcorr,cube,head,nread,red,mask,/indiv
  red=aprefcorr(cube,head,mask,indiv=indiv,horz=horz,vert=vert)
;  red1=aprefcorr(cube,head,mask,indiv=1)
;  red1h=aprefcorr(cube,head,mask,indiv=1,horz=0)
;  red3=aprefcorr(cube,head,mask,indiv=3)
;  for iquad=0,3 do begin
;  print,robust_sigma(cube[iquad*512+5:(iquad+1)*512-5,iquad*512+5:(iquad+1)*512-5,nread-1]-cube[iquad*512+5:(iquad+1)*512-5,iquad*512+5:(iquad+1)*512-5,1])
;  print,robust_sigma(red1[iquad*512+5:(iquad+1)*512-5,iquad*512+5:(iquad+1)*512-5,nread-1]-red1[iquad*512+5:(iquad+1)*512-5,iquad*512+5:(iquad+1)*512-5,1])
;  print,robust_sigma(red1h[iquad*512+5:(iquad+1)*512-5,iquad*512+5:(iquad+1)*512-5,nread-1]-red1h[iquad*512+5:(iquad+1)*512-5,iquad*512+5:(iquad+1)*512-5,1])
;  print,robust_sigma(red3[iquad*512+5:(iquad+1)*512-5,iquad*512+5:(iquad+1)*512-5,nread-1]-red3[iquad*512+5:(iquad+1)*512-5,iquad*512+5:(iquad+1)*512-5,1])
;  endfor
;  atv,cube[0:2047,0:2047,nread-1]-cube[0:2047,0:2047,1],min=-50,max=50
;  atv,red1[0:2047,0:2047,nread-1]-red1[0:2047,0:2047,1],min=-50,max=50
;  atv,red1h[0:2047,0:2047,nread-1]-red1h[0:2047,0:2047,1],min=-50,max=50
;  atv,red3[0:2047,0:2047,nread-1]-red3[0:2047,0:2047,1],min=-50,max=50
;stop
endif else red=cube
print,'time: ',systime(/seconds)-time0

if n_elements(dark) gt 1 then begin
  print,'dark correction...'
  darksz=size(dark)
  if nread le darksz[3] then begin
    for i=0,nread-1 do begin
      card=string(format='(a,i3.3)','SLICE',i)
      islice=sxpar(head,card,count=ncard)
      if ncard eq 0 then islice=i+1
      red[*,*,i]-=dark[*,*,islice-1]
    endfor
    ; error from median dark counts
    err2=dark[*,*,nread-1]
  endif else begin
    print,'Not enough reads in dark frame!! .con to continue without dark subtraction'
    stop
    err2=fltarr(2048,2048)
  endelse
endif else err2=fltarr(2048,2048)
print,'time: ',systime(/seconds)-time0

if not keyword_set(nocr) then begin
  print,'cr rejection'
  crrej,red,nread,mask,gain,rn,nofix=nofix,ncr=ncr
endif
print,'time: ',systime(/seconds)-time0

if not keyword_set(nofs) then begin
  print,'fowler sample'
  d=fsamp(red,nfs,var,gain,rn)
  err2+=var
endif else begin
  d=red
  err2+=red[*,*,nread-1]-red[*,*,1]
endelse
print,'time: ',systime(/seconds)-time0

if n_elements(flat) gt 1 then begin
  print,'2D flat correction'
  if keyword_set(nofs) then $
  for i=0,nread-1 do red[*,*,i]/=flat $
  else d/=flat
  err2/=flat^2
endif
print,'time: ',systime(/seconds)-time0
err=sqrt(err2)
sz=size(d)

if keyword_set(stp) then stop

return,d

end

