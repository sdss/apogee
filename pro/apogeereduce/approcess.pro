; alternate  master reduction routine, not generally used!
; Holtz.

; approcess reduces a sequence of images, all 3 chips, and writes out
function approcess,nums,cmjd=cmjd,step=step,clobber=clobber,onedclobber=onedclobber,$
       detid=detid,darkid=darkid,flatid=flatid,traceid=traceid,psfid=psfid,fluxid=fluxid,$
       waveid=waveid,littrowid=littrowid,persistid=persistid,nocr=nocr,stp=stp,jchip=jchip,nfs=nfs,nofs=nofs,$
       doproc=doproc,doap3dproc=doap3dproc,doap2dproc=doap2dproc,logfile=logfile,outdir=outdir,maxread=maxread,$
       skywave=skywave

common savedepsf, savedepsffiles, epsfchip
savedepsffiles=[' ',' ',' ']
epsfchip=0

if not keyword_set(step) then step=1
if not keyword_set(stp) then stp=0
if not keyword_set(nocr) then nocr=0
if not keyword_set(nofs) then nofs=0
if not keyword_set(nfs) then nfs=0
if not keyword_set(detid) then detid=0
if not keyword_set(darkid) then darkid=0
if not keyword_set(flatid) then flatid=0
if not keyword_set(littrowid) then littrowid=0
if not keyword_set(persistid) then persistid=0
if not keyword_set(traceid) then traceid=0
if not keyword_set(psfid) then psfid=0
if not keyword_set(fluxid) then fluxid=0
if not keyword_set(waveid) then waveid=0
if keyword_set(clobber) then onedclobber=1
if not keyword_set(maxread) then maxread=[0,0,0]

dirs=getdir(apodir,caldir,spectrodir,vers,datadir=datadir)
leadstr = 'APPROCESS: '

print,'Using detector: ',detid,' dark: ', darkid,' flat: ', flatid, ' trace: ', traceid, ' psf: ',psfid
if keyword_set(jchip) then begin
  j1=jchip & j2=jchip
endif else begin
  j1=0 & j2=2
endelse
chip=['a','b','c']

if keyword_set(doproc) or keyword_set(doap3dproc) then begin 
 ; use approcess as front end for ap3dproc and ap2dproc only 
 for ichip=j1,j2 do begin
  ; set up calibration file names
  if darkid gt 0 then $
    bpmcorr=caldir+'bpm/'+dirs.prefix+string(format='("BPM-",a,"-",i8.8,".fits")',chip[ichip],darkid)
  if darkid gt 0 then $
    darkcorr=caldir+'darkcorr/'+dirs.prefix+string(format='("Dark-",a,"-",i8.8,".fits")',chip[ichip],darkid)
  if flatid gt 0 then $
    flatcorr=caldir+'flatcorr/'+dirs.prefix+string(format='("Flat-",a,"-",i8.8,".fits")',chip[ichip],flatid)
  apgundef,littrowcorr
  if littrowid gt 0 and ichip eq 1 then $
    littrowcorr=caldir+'littrow/'+dirs.prefix+string(format='("Littrow-",a,"-",i8.8,".fits")',chip[ichip],littrowid)
  if persistid gt 0 then $
    persistcorr=caldir+'persist/'+dirs.prefix+string(format='("Persist-",a,"-",i8.8,".fits")',chip[ichip],persistid)
  for inum=0,n_elements(nums)-1 do begin
   num=nums[inum]
   if not keyword_set(cmjd) then cmjd=getcmjd(num)
   if not keyword_set(outdir) then outdir=dirs.expdir+cmjd+'/'
   ifile=dirs.prefix+string(format='("R-",a,"-",i8.8)',chip[ichip],num) 
   ofile=dirs.prefix+string(format='("2D-",a,"-",i8.8)',chip[ichip],num) 
   if file_test(outdir,/directory) eq 0 then file_mkdir,outdir
   if nfs eq 0 then uptheramp=1 else uptheramp=0
print,'calling ap3dproc...'
   ap3dproc,datadir+cmjd+'/'+ifile+'.apz',outdir+ofile+'.fits',$
     flatcorr=flatcorr,darkcorr=darkcorr,bpmcorr=bpmcorr,littrowcorr=littrowcorr,$
     persistcorr=persistcorr,$
     nocr=nocr,uptheramp=uptheramp,nfowler=nfs,fitsdir=getlocaldir(),clobber=clobber,maxread=maxread[ichip]
  endfor
 endfor
 ; ap2dproc does all 3 chips together
 if keyword_set(doproc) then begin
  for inum=0,n_elements(nums)-1 do begin
   num=nums[inum]
   tracefile=caldir+'psf/'+string(format='(i8.8)',psfid)
   wavefile=0
   if keyword_set(waveid) then wavefile=caldir+'wave/'+string(format='(i8.8)',waveid)
   ap2dproc,dirs.expdir+cmjd+'/'+string(format='(i8.8)',num),$
            tracefile,4,outdir=outdir,wavefile=wavefile,clobber=clobber,skywave=skywave

   chiptag=['a','b','c']
   files = apogee_filename('2D',num=num,chip=chiptag)
   modfiles = apogee_filename('2Dmodel',num=num,chip=chiptag)
   for jj=0,n_elements(files)-1 do begin
      if file_test(files[jj]) then begin
        file_delete,files[jj]+'.fz',/allow_nonexistent
        ;SPAWN,['fpack','-D','-Y',files[jj]],/noshell
      endif
      if file_test(modfiles[jj]) then begin
        file_delete,modfiles[jj]+'.fz',/allow_nonexistent
        SPAWN,['fpack','-D','-Y',modfiles[jj]],/noshell
      endif
   endfor

  endfor
 endif
 return,0
endif

print,'rogue reduction....!?'

;cmjd=string(format='(i5.5)',mjd) 
dark=0 & flat=0 & out=0
for ichip=j1,j2 do begin
 time0=systime(/seconds)

 ; Prepare to read in calibration frames for this chip and set calmask
 calmask=bytarr(2048,2048)
 readdet=1 & readdark=1 & readflat=1 & readpsf=1 & readflux=1 & readwave=1
 print,'time: ',systime(/seconds)-time0

 ; loop over each requested input frame
 for inum=0,n_elements(nums)-1 do begin
  num=nums[inum]
  if not keyword_set(cmjd) then cmjd=getcmjd(num)
  outdir=spectrodir+'red/'+cmjd+'/'
  if file_test(outdir,/directory) eq 0 then file_mkdir,outdir

  ; 3D->2D (or just 3D to processed)
  if keyword_set(nofs) then file=dirs.prefix+string(format='("3D-",a,"-",i8.8)',chip[ichip],num) $
  else file=dirs.prefix+string(format='("2D-",a,"-",i8.8)',chip[ichip],num) 
  while file_test(outdir+'/'+file+'.lock') do apwait,outdir+'/'+file+'.lock',10
  if (keyword_set(clobber) or (not file_test(outdir+'/'+file+'.fits'))) then begin
    time1=systime(/seconds)
    openw,lock,/get_lun,outdir+'/'+file+'.lock'
    free_lun,lock
    if readdet and detid gt 0 then begin
      print,'reading detector'
      readdet=0
      detcorr=string(format='("apDetector-",a,"-",i8.8)',chip[ichip],detid)
      gain=mrdfits(detcorr+'.fits',1)
      rn=mrdfits(detcorr+'.fits',2)
      lincorr=mrdfits(detcorr+'.fits',3)
    endif else begin
      gain=[2.,2.,2.,2.]
      rn=[12.,12.,12.,12.]
      lincorr=fltarr(4,3)
      for iquad=0,3 do lincorr[iquad,*]=[0.,1.,0.]
    endelse
    if readdark and darkid gt 0 then begin
      print,'reading dark'
      readdark=0
      if step gt 1 then $
        darkcorr=string(format='("apDark",i1,"-",a,"-",i8.8)',step,chip[ichip],darkid) $
      else $
        darkcorr=string(format='("apDark-",a,"-",i8.8)',chip[ichip],darkid)
      bpmcorr=string(format='("apBPM-",a,"-",i8.8)',chip[ichip],darkid)
      dark=mrdfits(caldir+'darkcorr/'+darkcorr+'.fits',1)
      darkmask=mrdfits(caldir+'darkcorr/'+darkcorr+'.fits',3)
      bad=where(darkmask gt 0)
      calmask[bad] = calmask[bad] or maskval('BADDARK')
    endif
    if readflat and flatid gt 0 then begin
      print,'reading 2D flat'
      readflat=0
      flatcorr=string(format='("apFlat-",a,"-",i8.8)',chip[ichip],flatid)
      flat=mrdfits(caldir+'flatcorr/'+flatcorr+'.fits',1)
      flatmask=mrdfits(caldir+'flatcorr/'+flatcorr+'.fits',3)
      bad=where(flatmask gt 0)
      calmask[bad] = calmask[bad] or maskval('BADFLAT')
    endif
    time2=systime(/seconds)
    mask=calmask
    d=process(cmjd,num,chip[ichip],head,red,dark,flat,err,gain,rn,mask,$
      step=step,stp=stp,$
      nocr=nocr,nofs=nofs,nfs=nfs,ncr=ncr)
    sz=size(d)
    if sz[0] eq 0 then goto, nextim
    print,'time: ',systime(/seconds)-time0
    if keyword_set(nofs) then begin
      if ichip eq j1 then out=fltarr(sz[1],sz[2],sz[3],3)
      out[*,*,*,ichip]=d
    endif else begin
      ; add header information for master header
      sxaddhist,leadstr+'Reduction version '+vers,head
      sxaddhist,leadstr+systime(0),head
      sxaddhist,leadstr+'IDL '+!version.release+' '+!version.os+' '+$
         !version.arch,head
      sxaddhist,leadstr+'Output File:',head
      sxaddhist,leadstr+' HDU1 - image (ADU)',head
      sxaddhist,leadstr+' HDU2 - error (ADU)',head
      sxaddhist,leadstr+' HDU3 - flag mask',head
      sxaddhist,leadstr+'        1 - bad pixels',head
      sxaddhist,leadstr+'        2 - cosmic ray',head
      sxaddhist,leadstr+'        4 - saturated',head
      sxaddhist,leadstr+'        8 - bad dark',head
      sxaddhist,leadstr+'       16 - bad flat',head
      sxaddhist,leadstr+'       32 - bad err',head
      if detid gt 0 then sxaddhist,leadstr+'DETECTOR file="'+detcorr+'"',head
      sxaddpar,head,'DETID',detid
      if darkid gt 0 then sxaddhist,leadstr+'DARK file="'+darkcorr+'"',head
      sxaddpar,head,'DARKID',darkid
      if flatid gt 0 then sxaddhist,leadstr+'FLAT file="'+flatcorr+'"',head
      sxaddpar,head,'FLATID',flatid
      if not keyword_set(nocr) then sxaddhist,leadstr+'CR rejection ON',head
      sxaddhist,leadstr+' Fowler sampling: '+string(nfs)+' (0 = UTR)',head
      mwrfits,0,outdir+'/'+file+'.fits',head,/create
      ; set bad pixels to 0. for consistency with DLN
      dout=d & errout=err
      bad=where(mask and $
        (maskval('BADPIX') or maskval('SATPIX') or  $
         maskval('BADFLAT') or maskval('BADDARK')))
      if bad[0] ge 0 then begin
        d[bad]=!values.f_nan & dout[bad]=0. & errout[bad]=-1.
      endif
      bad=where(finite(d) eq 0) & if bad[0] ge 0 then dout[bad]=0.
      bad=where(finite(err) eq 0 or err eq 0.) 
      if bad[0] ge 0 then begin
        errout[bad]=-1.
        mask[bad] = mask[bad] or maskval('BADERR')
      endif
      ; header information for individual images
      mkhdr,head1,d,/image
      sxaddpar,head1,'CTYPE1','Pixel' & sxaddpar,head1,'CTYPE2','Pixel'
      sxaddpar,head1,'BUNIT','Flux (ADU)'
      mwrfits,dout,outdir+'/'+file+'.fits',head1,/silent
      sxaddpar,head1,'BUNIT','Error (ADU)'
      mwrfits,errout,outdir+'/'+file+'.fits',head1,/silent
      mkhdr,head1,mask,/image
      sxaddpar,head1,'CTYPE1','Pixel' & sxaddpar,head1,'CTYPE2','Pixel'
      sxaddpar,head1,'BUNIT','Flag mask (bitwise)'
      mwrfits,mask,outdir+'/'+file+'.fits',head1,/silent
      ; make some jpg output files
      impost,outdir+'/'+file+'.eps',d,z=-20,l=250
      if not file_test(outdir+'/plots',/dir) then file_mkdir,outdir+'/plots/'
      spawn,'convert '+outdir+'/'+file+'.eps '+outdir+'/plots/'+file+'a.jpg'
      impost,outdir+'/'+file+'.eps',d,z=-20,l=2500
      spawn,'convert '+outdir+'/'+file+'.eps '+outdir+'/plots/'+file+'b.jpg'
      file_delete,outdir+'/'+file+'.eps',/quiet
      if psfid eq 0 and ichip eq j1 then out=fltarr(2048,2048,3) 
      if keyword_set(doap3dproc) then begin
        ifile=string(format='("apR-",a,"-",i8.8)',chip[ichip],num) 
        ofile=string(format='("ap2D-",a,"-",i8.8)',chip[ichip],num) 
        if file_test(outdir+'dln',/directory) eq 0 then file_mkdir,outdir+'dln'
        ap3dproc,datadir+cmjd+'/'+ifile+'.fits',outdir+'dln/'+ofile+'.fits',flatcorr=caldir+'flatcorr/'+flatcorr+'.fits',darkcorr=caldir+'darkcorr/'+darkcorr+'.fits',bpmcorr=caldir+'bpm/'+bpmcorr+'.fits',/uptheramp,/clobber,maxread=maxread[ichip]
      endif
    endelse
    if keyword_set(logfile) then writelog,logfile,file+string(format='(2f8.2,i6)',systime(/seconds)-time1,time2-time1,ncr)
    file_delete,outdir+'/'+file+'.lock',/quiet
  endif else begin
    ; if ap2D exists, read it in for extraction
    print,'2D output file already exists! '+outdir+'/'+file+'.fits'
    if psfid eq 0 and ichip eq j1 then out=fltarr(2048,2048,3) 
    junk=mrdfits(outdir+'/'+file+'.fits',0,head)
    d=mrdfits(outdir+'/'+file+'.fits',1)
    err=mrdfits(outdir+'/'+file+'.fits',2)
    mask=mrdfits(outdir+'/'+file+'.fits',3)
    bad=where(mask and $
     (maskval('BADPIX') or maskval('SATPIX') or maskval('BADFLAT') or maskval('BADDARK') or maskval('BADFLAT')))
    if bad[0] ge 0 then d[bad]=!values.f_nan
    if bad[0] ge 0 then err[bad]=!values.f_nan
  endelse

  ; 2D->1D
  if not keyword_set(nofs) and psfid gt 0 then begin
    file=string(format='("ap1D-",a,"-",i8.8)',chip[ichip],num)
    while file_test(outdir+'/'+file+'.lock') do wait,10
    if (keyword_set(onedclobber) or (not file_test(outdir+'/'+file+'.fits'))) then begin
      time1=systime(/seconds)
      openw,lock,/get_lun,outdir+'/'+file+'.lock'
      free_lun,lock
      if readpsf and psfid gt 0 then begin
        print,'reading epsf '
        readpsf=0
        psfcorr=string(format='("apEPSF-",a,"-",i8.8)',chip[ichip],psfid)
        tmp=mrdfits(caldir+'/psf/'+psfcorr+'.fits',0,phead)
        ntrace=sxpar(phead,'NTRACE')
        img=ptrarr(ntrace,/allocate_heap)
        for i=0,ntrace-1 do begin
          ptmp=mrdfits(caldir+'/psf/'+psfcorr+'.fits',i+1,/silent)
          *img[i] = ptmp.img
          p ={lo: ptmp.lo, hi: ptmp.hi, img: img[i]}
          if i eq 0 then psf=replicate(p,ntrace)
          psf[i] = p
        endfor
        sxaddhist,leadstr+'EPSF file="'+psfcorr+'"',head
        sxaddpar,head,'EPSFID',psfid
      endif
      leadstr = 'APPROCESS: '
      if readflux and fluxid gt 0 then begin
        readflux=0
        fluxcorr=string(format='("apFlux-",a,"-",i8.8)',chip[ichip],fluxid)
        flux=mrdfits(caldir+'flux/'+fluxcorr+'.fits')
        sxaddhist,leadstr+'FLUX file="'+fluxcorr+'"',head
        sxaddpar,head,'FLUXID',fluxid
      endif
      if readwave and waveid gt 0 then begin
        readwave=0
        wavecorr=string(format='("apWave-",a,"-",i8.8)',chip[ichip],waveid)
        wave=mrdfits(caldir+'wave/'+wavecorr+'.fits',2)
        sxaddhist,leadstr+'WAVE file="'+wavecorr+'"',head
        sxaddpar,head,'WAVEID',waveid
      endif
      time2=systime(/seconds)
      print,'extract'
      if ichip eq j1 then out=fltarr(2048,ntrace,3)
      var=err^2
      model=d*0.
      out[*,*,ichip]=extract(d,ntrace,psf,var,model=model)
      if fluxid gt 0 then begin
        out[*,*,ichip]*=flux
        var*=flux
      endif
      mwrfits,0,outdir+file+'.fits',head,/create
      mkhdr,head1,d,/image
      mwrfits,out[*,*,ichip],outdir+file+'.fits',head1
      mwrfits,var,outdir+file+'.fits',head1,/silent
      mask=bytarr(2048,ntrace)  ; not yet populated!
      mwrfits,mask,outdir+file+'.fits',head1,/silent
      if waveid gt 0 then mwrfits,wave,outdir+file+'.fits',head1,/silent
      ; make some jpg output files
      if not file_test(outdir+'/plots',/dir) then file_mkdir,outdir+'/plots/'
      impost,outdir+'/'+file+'.eps',out[*,*,ichip],z=-20,l=500
      spawn,'convert '+outdir+'/'+file+'.eps '+outdir+'/plots/'+file+'a.jpg'
      impost,outdir+'/'+file+'.eps',out[*,*,ichip],z=-20,l=5000
      spawn,'convert '+outdir+'/'+file+'.eps '+outdir+'/plots/'+file+'b.jpg'
      file_delete,outdir+'/'+file+'.eps',/quiet
      ; output the 2Dmodel
      modfile=string(format='("ap2Dmodel-",a,"-",i8.8)',chip[ichip],num)
      mwrfits,model,outdir+modfile+'.fits',/create

      if keyword_set(doap2dproc) and ichip eq j2 then begin
        ; ap2dproc does all three chips
        if file_test(outdir+'dln',/directory) eq 0 then file_mkdir,outdir+'dln'
        tracefile=caldir+'trace/'+string(format='(i8.8)',psfid)
        ap2dproc,spectrodir+'red/'+cmjd+'/'+string(format='(i8.8)',num),$
                 tracefile,outdir+'dln/',3,/clobber
      endif
      if keyword_set(logfile) then writelog,logfile,file+string(format='(2f8.2)',systime(/seconds)-time1,time2-time1)
      file_delete,outdir+'/'+file+'.lock',/quiet
    endif else begin
      print,'1D output file already exists! '+outdir+'/'+file+'.fits'
      d=mrdfits(outdir+'/'+file+'.fits',1)
      sz=size(d)
      if ichip eq j1 then out=fltarr(sz[1],sz[2],3)
      out[*,*,ichip] = d
    endelse
  endif
  nextim:
  print,'time: ',systime(/seconds)-time0
 endfor
endfor
return,out

end


