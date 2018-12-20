;======================================================================
;
; mkflat : makes superflats from dithered individual frames
;
;  USAGE:  mkdark,ims,cmjd=cmjd,darkid=darkid,/clobber,/kludge,nrep=nrep
;
;   INPUT:  ims: list of image numbers to include in superflat
;           cmjd=cmjd : (optional,obsolete) gives MJD directory name if not encoded in file number
;           darkid=darkid : dark frame to be used if images are reduced
;           /clobber : rereduce images even if they exist
;           /kludge : set bottom and top non-illuminated pixels to unity
;           nrep=nrep : median filters each batch of nrep frames before combining
;
;   OUTPUT:  a set of apFlat-[abc]-ID8.fits files
;
pro mkflat,ims,cmjd=cmjd,darkid=darkid,clobber=clobber,kludge=kludge,nrep=nrep,dithered=dithered

 i1=ims[0]
 nframes=n_elements(ims)
 if not keyword_set(nrep) then nrep=1

 dirs=getdir(apodir,caldir,specdir,apovers,libdir,datadir=datadir)

 flatfile=caldir+'/flatcorr/'+dirs.prefix+string(format='("Flat-",i8.8)',i1)+'.tab'
 ; is another process already creating file?
 while file_test(flatfile+'.lock') do apwait,flatfile,10

 ; does file already exist?
 if file_test(flatfile) and not keyword_set(clobber) then begin
   print,' Flat file: ', flatfile, ' already made'
   return
 endif

 ; open lock file
 openw,lock,/get_lun,flatfile+'.lock'
 free_lun,lock

 sum={name: '',num: i1, nframes: 0}
 flatlog=REPLICATE(sum,3)

 perclow=0.85              ; fraction for rejecting pixels
 nperclow=0.95              ; fraction for rejecting neighbor pixels 
 perchi=1.25              ; fraction for rejecting pixels
 nperchi=1.05              ; fraction for rejecting neighbor pixels 
 x1norm=800 & x2norm=1200  ; region for getting normalization 
 y1norm=800 & y2norm=1000
 filter=[50,1]            ; filter size for smoothing for large scale structure
 chip=['a','b','c']
 
; ; process frames if necessary
; if keyword_set(cmjd) then $
;   d=approcess(ims,cmjd=cmjd,darkid=darkid,nfs=1,/clobber) $
; else d=approcess(ims,darkid=darkid,nfs=1,/nocr,/clobber) 

  outdir=caldir+'/flatcorr/'
  if file_test(outdir,/directory) eq 0 then file_mkdir,outdir
  nfs=1
  uptheramp=0
  nocr=1
  cmjd=getcmjd(ims[0],mjd=mjd)
  getcal,mjd,dirs.calfile,dark=darkid,bpm=bpmid,det=detid

  ; read and process frames to 2D
  for ichip=0,2 do begin
    if darkid gt 0 then darkcorr=apogee_filename('Dark',num=darkid,chip=chip[ichip])
    if detid gt 0 then detcorr=apogee_filename('Detector',num=detid,chip=chip[ichip])
    for inum=0,n_elements(ims)-1 do begin
     num=ims[inum]
     ifile=apogee_filename('R',num=num,chip=chip[ichip])
     ofile=apogee_filename('2D',num=num,chip=chip[ichip],/base)
     ap3dproc,ifile,outdir+ofile,$
       detcorr=detcorr,$
       darkcorr=darkcorr,$
       nocr=nocr,uptheramp=uptheramp,nfowler=nfs,fitsdir=getlocaldir()
    endfor
  endfor
  
 ; sum up all of the individual flats
 ;  median nrep frames before summing if requested
 flats=fltarr(2048,2048,3,nframes)
 flatmasks=intarr(2048,2048,3)
 flatsum=fltarr(2048,2048,3)
 for ii=0,nframes-1,nrep do begin
  for irep=0,nrep-1 do begin
    i=ims[ii+irep]
    for ichip=0,2 do begin
      ofile=apogee_filename('2D',num=i,chip=chip[ichip],/base)
      f=mrdfits(outdir+ofile,0,head)
      flats[*,*,ichip,ii+irep]=mrdfits(outdir+ofile,1)
      flatmasks[*,*,ichip]=mrdfits(outdir+ofile,3)
    endfor
  endfor
  if ii eq 0 then head0=head
  if nrep gt 1 then flatsum+=median(flats[*,*,*,ii:ii+nrep-1],dimension=4) $
  else flatsum+=flats[*,*,*,ii]
 endfor
 
 ; create the superflat 
 for ichip=0,2 do begin
  flat=flatsum[*,*,ichip]
 
   ; normalize
  norm=median(flat[x1norm:x2norm,y1norm:y2norm])
  flat/=norm
  ; create mask
  sz=size(flat)
  mask=bytarr(sz[1],sz[2])

  ; mask from reductions, using last frame read
  bad=where((flatmasks[*,*,ichip] and badmask()) gt 0,nbad)
  if nbad gt 0 then flat[bad]=!values.f_nan
  if nbad gt 0 then mask[bad]=mask[bad] or 1

  ; set pixels to bad when below some fraction 
  localflat=flat/zap(flat,[100,10])
  ; relative to neighbors
  low=where(localflat lt perclow,nlow)
  if nlow gt 0 then mask[low]=mask[low] or 2
  ; absolute
  low=where(flat lt 0.1,nlow)
  if nlow gt 0 then mask[low]=mask[low] or 2
  ; high pixels
  hi=where(localflat gt perchi,nhi)
  if nhi gt 0 then mask[hi]=mask[hi] or 2

  ; set neighboring pixels to bad at slightly lower threshold, iteratively
  for iter=0,10  do begin
    low=where(mask gt 0)
    n=[-1,1,-2049,-2048,-2047,2047,2048,2049]
    for in=0,n_elements(n)-1 do begin
      neigh=low+n[in]
      off=where(neigh lt 0 or neigh gt 2048L*2048,noff)
      if noff gt 0 then neigh[off]=0
      lowneigh=where(localflat[neigh] lt nperclow,nlow)
      if nlow gt 0 then mask[neigh[lowneigh]]=mask[neigh[lowneigh]] or 4
      hineigh=where(localflat[neigh] gt nperchi,nhi)
      if nhi gt 0 then mask[neigh[hineigh]]=mask[neigh[hineigh]] or 4
    endfor
  endfor

  ; mask any zero values
  bad=where(flat eq 0.,nbad)
  if nbad gt 0 then mask[bad]=mask[bad] or 8
  if nbad gt 0 then flat[bad]=!values.f_nan
 
  if keyword_set(dithered) then begin
    ; get the large scale structure from smoothing, avoiding bad pixels (NaNs)
    sm=smooth(flat,100,/nan,/edge_truncate)
    rows=intarr(2048)+1
    smrows=(total(sm,1,/nan)/2048)##rows
    smcols=rows##(total(sm,2,/nan)/2048)
   
    ; median filter the median flat with a rectangular filter and 
   ;   divide flat by this to remove horizontal structure
   sflat=zap(flat,filter)
   flat/=sflat
  
   ; now put the large scale structure in apart from
   ;  structure that is horizontal (fibers) or vertical (SED)
   flat*=sm/smrows/smcols
  
  ; kludge to set unilluminated pixels to 1
;  if keyword_set(kludge) then begin
    for i=0,13 do  begin
      dark=where(flat[*,i] lt -99, ndark)
      if ndark gt 0 then flat[dark,i]=1.
      if ndark gt 0 then mask[dark,i]=0
      dark=where(flat[*,2047-i] lt -99, ndark)
      if ndark gt 0 then flat[dark,2047-i]=1.
      if ndark gt 0 then mask[dark,2047-i]=0
    endfor
;  endif
 endif else begin
 ; if not dithered, still take out spectral signature
   rows=intarr(2048)+1
   cols=fltarr(2048)

   ; spectral signature from median of each column
   for icol=0,2047 do cols[icol]=median(flat[icol,*])
   ; medfilt doesn't do much if intensity is varying across cols
   smrows=rows##medfilt1d(cols,100)
   sflat=smrows

   ; Dec 2018: don't take out median spectral signature, this leaves 
   ;  structure in spectra that is hard to normalize out
   ; instead, take out a low order polynomial fit to estimate spectral signature
   x=indgen(2048)
   gd=where(finite(cols))
   coef=robust_poly_fit(x[gd],cols[gd],2)
   smrows=rows##poly(x,coef)
   sflat=smrows

   ; divide out estimate of spectral signature
   flat/=smrows

   ;  do, however, adjust for small differences between quadrants
   ;  that are visible in high signal images and don't seem to 
   ;  be removed by reference pixel reduction....
   ;flat/=smrows

   ; WAIT !!! Doesn't this account for small gain differences?
   ;for i=1,3 do begin
   ;  left=median(flat[i*512-1,*])
   ;  right=median(flat[i*512,*])
   ;  flat[i*512:2047,*] *= left/right
   ;endfor
 endelse
 ; set bad values to -100 before writing to avoid NaNs in output file
 bad=where(finite(flat) eq 0,nbad)
 if nbad gt 0 then flat[bad]=0.

 ; write it out!
 file=apogee_filename('Flat',num=i1,chip=chip[ichip])
 leadstr = 'APMKFLAT: '
 sxaddhist,leadstr+systime(0),head0
 info = GET_LOGIN_INFO()
 sxaddhist,leadstr+info.user_name+' on '+info.machine_name,head0
 sxaddhist,leadstr+'IDL '+!version.release+' '+!version.os+' '+!version.arch,head0
 sxaddhist,leadstr+' APOGEE Reduction Pipeline Version: '+getvers(),head0
 mwrfits,0,file,head0,/create
 mwrfits,flat,file
 mwrfits,sflat,file
 mwrfits,mask,file

 ; make a jpg of the flat
 if not file_test(caldir+'flatcorr/plots',/dir) then file_mkdir,caldir+'flatcorr/plots'
 flatplot,flat,caldir+'flatcorr/plots/'+file_basename(file,'.fits')

 flatlog[ichip].name=file
 flatlog[ichip].num=i1
 flatlog[ichip].nframes=nframes

; flat=0
; time=systime(/seconds)
; print,'done '+chip[ichip],time-time0
endfor

; write out flat summary information
file=dirs.prefix+string(format='("Flat-",i8.8)',i1) 
mwrfits,flatlog,caldir+'flatcorr/'+file+'.tab',/create

; remove lock file
file_delete,flatfile+'.lock',/allow_non

; compile summary web page
flathtml,caldir

end
