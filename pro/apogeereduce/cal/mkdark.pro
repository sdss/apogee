;======================================================================
;
; mkdark : makes superdarks
;
;  USAGE:  mkdark,ims,cmjd=cmjd,step=step,psfid=psfid
;
;   INPUT:  ims: list of image numbers to include in superdark
;           cmjd=cmjd : (optional,obsolete) gives MJD directory name if not encoded in file number
;           step=step : (optional,obsolete) process every step image in UTR
;           psfid=psfid : (optional) EPSF id to use to try to subtract off thermal traces
;
;   OUTPUT:  a set of apDark-[abc]-ID8.fits files
;
pro mkdark,ims,cmjd=cmjd,step=step,psfid=psfid,clobber=clobber

i1=ims[0]
nframes=n_elements(ims)

dirs=getdir()
caldir=dirs.caldir

; is another process already creating file
darkfile=caldir+'/darkcorr/'+dirs.prefix+string(format='("Dark-",i8.8)',i1) +'.tab'
while file_test(darkfile+'.lock') do apwait,darkfile,10

; does file already exist?
if file_test(darkfile) and ~keyword_set(clobber) then begin
  print,' Dark file: ', darkfile, ' already made'
  return
endif

; open lock file
openw,lock,/get_lun,darkfile+'.lock'
free_lun,lock

; initialize summary structure
sum={num: i1, nframes: 0, nreads: 0, nsat: 0L, nhot: 0L, nhotneigh: 0L, nbad: 0L, medrate: 0., psfid: 0L, nneg: 0L}
darklog=REPLICATE(sum,3)

if not keyword_set(step) then step=0
; loop over the chips
chips=['a','b','c']
for ichip=0,n_elements(chips)-1 do begin
 chip=chips[ichip]

 time0=systime(/seconds)
 ii=0
 for jj=0,n_elements(ims)-1 do begin
  i=ims[jj]
  if not keyword_set(cmjd) then cm=getcmjd(i) else cm=cmjd
  print,chip,i

  ; process (bias-only) each individual frame
  d=process(cm,i,chip,head,r,step=step,/nofs,/nofix,/nocr)
  print,'done process'
  sz=size(d)
  if sz[1] ne 2048 then stop,sz
  mask=bytarr(sz[1],sz[2])

  ; construct cube of reads minus second read
  if jj eq 0 then head0=head
  sz=size(r)
  if jj eq 0 then begin
    if ichip eq 0 then red=fltarr(2048,2048,sz[3],nframes) else red*=0.
  endif
  red[*,*,*,ii]=r
  apgundef,r
  for iread=sz[3]-1,1,-1 do begin
    red[*,*,iread,ii]-=red[*,*,1,ii]
  endfor
  ii=ii+1
  help,/mem
 endfor
 
 ; median them all
 print,'median...'
 dark=median(red,dimension=4)

 ; option to remove any trace of spectral traces
 if keyword_set(psfid) then begin
   darklog[ichip].psfid=psfid
   print,'reading epsf '
   file=dirs.prefix+string(format='("EPSF-",a,"-",i8.8)',chip,psfid)
   tmp=mrdfits(caldir+'/psf/'+file+'.fits',0,head)
   ntrace=sxpar(head,'NTRACE')
   img=ptrarr(ntrace,/allocate_heap)
   for i=0,ntrace-1 do begin
     ptmp=mrdfits(caldir+file+'.fits',i+1,/silent)
     *img[i] = ptmp.img
     p ={lo: ptmp.lo, hi: ptmp.hi, img: img[i]}
     if i eq 0 then psf=replicate(p,ntrace)
     psf[i] = p
   endfor
   nread=sz[3]
   for iread=2,nread-1 do begin
     var=dark[*,*,iread]
     ; want to subtract off mean background dark level before fitting traces
     ; iterate once for this
     back=median(dark[*,*,iread],10)
     niter=2
     for iter=0,niter-1 do begin
       print,iread,iter
       d=dark[*,*,iread]-back
       spec=extract(d,ntrace,psf,var)
       sspec=zap(spec,[200,1])
       d*=0
       for k=0,ntrace-1 do begin
         p1=psf[k]
         lo=psf[k].lo & hi=psf[k].hi
         img=*p1.img
         r=intarr(hi-lo+1)+1
         sub=sspec[*,k]#r
         bad=where(sub lt 0,nbad)
         if nbad gt 0 then sub[bad]=0
         d[*,lo:hi]+=sub*img
       endfor
       if iter lt niter-1 then back=median(dark[*,*,iread]-d,10)
     endfor
     dark[*,*,iread]-=d
   endfor
 endif

 ; flag "hot" pixels in mask image
 nread=sz[3]
 rate=(dark[*,*,nread-1]-dark[*,*,1])/(nread-2)

 ; create mask array
 ; NaN is bad!
 bad=where(finite(rate) eq 0,nsat) 
 if nsat gt 0 then mask[bad]=mask[bad] or 1

 ; flux accumulating very fast is bad!
 ;maxrate=25.
 maxrate=10.
 hot=where(rate gt maxrate,nhot)
 if nhot gt 0 then mask[hot]=mask[hot] or 2
 ; flag adjacent pixels to hot pixels as bad at 1/4 the maximum rate
 n=[-1,1,-2048,2048]
 nhotneigh=0
 for in=0,3 do begin
   ; only consider neighbors on the chip!
   neigh=hot+n[in]
   on=where(neigh ge 0 and neigh lt 2048L*2048L)
   nlow=where(rate[neigh[on]] gt maxrate/4.,nbad)
   if nbad ge 0 then mask[neigh[on[nlow]]]=mask[neigh[on[nlow]]] or 4
   nhotneigh+=n_elements(hot)
   ; same for bad
   neigh=bad+n[in]
   on=where(neigh ge 0 and neigh lt 2048L*2048L)
   nlow=where(rate[neigh[on]] gt maxrate/4.,nbad)
   if nbad gt 0 then mask[neigh[on[nlow]]]=mask[neigh[on[nlow]]] or 4
   nhotneigh+=n_elements(hot)
 endfor

 print,'Creating chi2 array ....'
 chi2=fltarr(2048L*2048*nread)
 n=intarr(2048L*2048*nread)
 dark=reform(dark,2048L*2048*nread,/overwrite)
 for ii=0,nframes-1 do begin
   tmp=reform(red[*,*,*,ii],2048L*2048*nread)
   good=where(finite(tmp),ngood)
   if ngood gt 0 then chi2[good]+=(tmp[good]-dark[good])^2/apvariance(dark[good],1)
   n[good]+=1
 ;  chi2+=(red[*,*,*,ii]-dark)^2/apvariance(dark,1)
 endfor
 chi2/=n
 dark=reform(dark,2048,2048,nread,/overwrite)
 chi2=reform(chi2,2048,2048,nread,/overwrite)

 ; set nans to 0 before writing
 bad=where(finite(dark) eq 0,nbad)
 if nbad gt 0 then dark[bad]=0.
 medrate=median(rate)

 ; clip pixels where rate is less than 3sigma of expected noise
 ;  for 2 readouts?
 ;low=where(dark[*,*,nread-1] lt 3*1.25*sqrt(apvariance(0,2)/nframes))

 ;median filter along reads dimenstion
 for i=0,2047 do begin
   slice=reform(dark[i,*,*])
   dark[i,*,*]=medfilt2d(slice,7,dim=2)
 endfor

 ; set negative pixels to zero
 neg=where(dark lt -10,nneg)
 if nneg gt 0 then dark[neg]=0.

 ; write them out
 if step gt 1 then  $
   file=dirs.prefix+string(format='("Dark",i1,"-",a,"-",i8.8)',step,chip,i1) $
 else $
   file=dirs.prefix+string(format='("Dark-",a,"-",i8.8)',chip,i1) 

 leadstr = 'APMKDARK: '
 sxaddhist,leadstr+systime(0),head0
 info = GET_LOGIN_INFO()
 sxaddhist,leadstr+info.user_name+' on '+info.machine_name,head0
 sxaddhist,leadstr+'IDL '+!version.release+' '+!version.os+' '+!version.arch,head0
 sxaddhist,leadstr+' APOGEE Reduction Pipeline Version: '+getvers(),head0
 mwrfits,0,caldir+'darkcorr/'+file+'.fits',head0,/create
 mwrfits,dark,caldir+'darkcorr/'+file+'.fits'
 mwrfits,chi2,caldir+'darkcorr/'+file+'.fits'
 mwrfits,mask,caldir+'darkcorr/'+file+'.fits'

 ; make some plots/images
 if not file_test(caldir+'darkcorr/plots',/dir) then file_mkdir,caldir+'darkcorr/plots'
 darkplot,dark,mask,caldir+'darkcorr/plots/'+file,/hard
 
 ; summary data table
 darklog[ichip].num=i1
 darklog[ichip].nframes=nframes
 darklog[ichip].nreads=nread
 darklog[ichip].nsat=nsat
 darklog[ichip].nhot=nhot
 darklog[ichip].nhotneigh=nhotneigh
 darklog[ichip].nbad=nbad
 darklog[ichip].medrate=medrate
 darklog[ichip].nneg=nneg

 ; save the rate file
 file=dirs.prefix+string(format='("DarkRate-",a,"-",i8.8)',chip,i1) 
 mwrfits,rate,caldir+'darkcorr/'+file+'.fits',/create

 dark=0
 time=systime(/seconds)
 print,'done '+chip,time-time0

endfor

apgundef,red

; write the summary log information
file=dirs.prefix+string(format='("Dark-",i8.8)',i1) 
mwrfits,darklog,caldir+'darkcorr/'+file+'.tab',/create

; remove lock file
file_delete,darkfile+'.lock'

; compile summary web page
darkhtml,caldir

end
