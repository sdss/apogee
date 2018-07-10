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
