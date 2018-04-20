;======================================================================
pro mkdet,detid,linid

 dirs=getdir()
 caldir=dirs.caldir
 file=dirs.prefix+string(format='("Detector-c-",i8.8)',detid)
 ;if another process is alreadying make this file, wait!
 while file_test(caldir+'detector/'+file+'.lock') do apwait,file,10
 ; does product already exist?
 print,'testing detector file: ',file+'.fits'
 if file_test(caldir+'/detector/'+file+'.fits') and not keyword_set(clobber) then begin
   print,' Detector file: ', file+'.fits', ' already made'
   return
 endif

 print,'making Detector: ', detid
 ; open .lock file
 openw,lock,/get_lun,caldir+'detector/'+file+'.lock'
 free_lun,lock

 lincorr=fltarr(4,3)
 for iquad=0,3 do lincorr[iquad,*]=[1.,0.,0.]
 if n_elements(linid) gt 0 then if linid gt 0 then begin
   linpar=mklinearity(linid)
   for iquad=0,3 do lincorr[iquad,*]=linpar
 endif

 if dirs.instrument eq 'apogee-n' then begin
   g=1.9
   ; these are supposed to be CDS DN!
   r=[20.,11.,16.] /sqrt(2.)
   ; 10/17 analysis get 12, 9, 10 single read in DN
   r=[12,8,8] ; actually achieved CDS
   r=[13,11,10] ; equivalent single read in UTR analysis (which gives lower rn overall)
 endif else begin
   g=3.0
   ; these are supposed to be CDS DN!
   r=[15.,15.,15.] /sqrt(2.)
   ; 10/17, get 6, 8, 4.5 single read in DN
   r=[4,5,3]   ; actually achieved CDS
   r=[7,8,4]   ; equivalent single read in UTR analysis (which gives lower rn overall)
   ; JCW 2/28/17 email
   ;Our current measurements are (blue, green, red):
   ;
   ;gain (e-/ADU): 3.0, 3.3, 2.7
   ;read noise (e-): 9.3, 15.2, 8.6
   ;dark current (e-/sec): 0.011, 0.014, 0.008
 endelse

 chips=['a','b','c']
 for ichip=0,2 do begin
   gain=[g,g,g,g]
   rn=[r[ichip],r[ichip],r[ichip],r[ichip]]*gain
   file=dirs.prefix+string(format='("Detector-",a,"-",i8.8)',chips[ichip],detid)
   mkhdr,head,0
   mwrfits,0,caldir+'/detector/'+file+'.fits',head,/create
   mwrfits,rn,caldir+'/detector/'+file+'.fits'
   mwrfits,gain,caldir+'/detector/'+file+'.fits'
   mwrfits,lincorr,caldir+'/detector/'+file+'.fits'
 endfor

 file_delete,caldir+'detector/'+file+'.lock'
end
