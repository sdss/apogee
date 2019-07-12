;+
; mkplan : makes plan files given input image numbers, MJD, psfid, fluxid
;          includes options for dark frames, calibration frames, sky frames,
;          ASDAF frames. This is called from the manually prepared MJD5.pro 
;          procedures
;-

pro mkplan,ims,plate,mjd,psfid,fluxid,$
   stars=stars,names=names,Hmags=Hmags,vers=vers,cal=cal,sky=sky,dark=dark,$
   plugid=plugid,ignore=ignore,test=test,onem=onem,fixfiberid=fixfiberid,suffix=suffix,mapper_data=mapper_data

print,'Making plan for MJD: ',mjd

; set up directories, plate, and MJD variables
dirs=getdir(apogee_dir,cal_dir,spectro_dir,vers,lib_dir,prefix,apred_vers=apred_vers,datadir=datadir,onem=onem)
if n_elements(mapper_data) eq 0 then mapper_data=dirs.mapperdir
telescope=dirs.telescope
calfile=dirs.calfile

if size(plate,/type) eq 7 then cplate=plate else cplate=strtrim(string(format='(i6.4)',plate),2)
cmjd=string(format='(i5.5)',mjd)

; planfile name and directory
if keyword_set(cal) then planfile=apogee_filename('CalPlan',mjd=mjd,instrument=dirs.instrument) else $
if keyword_set(dark) then planfile=apogee_filename('DarkPlan',mjd=mjd,instrument=dirs.instrument) else $
if keyword_set(onem) then begin
  planfile=apogee_filename('Plan',plate=plate,reduction=names[0],mjd=cmjd) 
  if suffix ne '' then planfile=file_dirname(planfile)+'/'+file_basename(planfile,'.par')+suffix+'.par'
endif else $
  planfile=apogee_filename('Plan',plate=plate,mjd=mjd)
outdir=file_dirname(planfile)+'/'
if not file_test(outdir,/dir) then file_mkdir,outdir

; get calibration files for this date
if n_elements(fixfiberid) gt 0 then fix0=fixfiberid else undefine,fix0
getcal,mjd,calfile,darkid=darkid,flatid=flatid,bpmid=bpmid,waveid=waveid,multiwaveid=multiwaveid,$
     responseid=responseid,lsfid=lsfid,detid=detid,sparseid=sparseid,fiberid=fiberid,badfiberid=badfiberid,$
     fixfiberid=fixfiberid,littrowid=littrowid,persistid=persistid,persistmodelid=persistmodelid
if n_elements(fix0) gt 0 then fixfiberid=fix0

; outplan plan file name
if keyword_set(stars) and not keyword_set(onem) then planfile=file_dirname(planfile)+'/'+file_basename(planfile,'.par')+'star.par' $
else if keyword_set(sky) then planfile=file_dirname(planfile)+'/'+file_basename(planfile,'.par')+'sky.par' 

if keyword_set(sky) then apdailycals,lsfs=ims,psf=psfid
print,planfile

; open plan file and write header
file_delete,planfile,/allow_nonexistent
openw,fplan,/get_lun,planfile
printf,fplan,'apogee_ver  ',getenv('APOGEE_VER')
;printf,fplan,'plateid ',plate
printf,fplan,'telescope '''+telescope+''''
printf,fplan,'instrument '''+dirs.instrument+''''
printf,fplan,'plateid '''+cplate+''''
printf,fplan,'mjd ',mjd
printf,fplan,'planfile '+file_basename(planfile)
printf,fplan,'logfile ''apDiag-'+cplate+'-'+cmjd+'.log'''
printf,fplan,'plotfile ''apDiag-'+cplate+'-'+cmjd+'.ps'''

; apred_vers keyword will override strict versioning using the plan file!
printf,fplan,'apred_vers '''+apred_vers+''''

if keyword_set(onem) then begin
    printf,fplan,'data_dir  '+datadir+'/'
    printf,fplan,'raw_dir  '+datadir+cmjd+'/'
    printf,fplan,'plate_dir '+outdir
    printf,fplan,'star_dir '+spectro_dir+'/fields/apo1m/'
    printf,fplan,'survey  apo1m'
    printf,fplan,'name '''+strtrim(string(names[0]),2)+'''' 
;    printf,fplan,'fiber '''+strtrim(string(stars[0]),2)+''''
;    if keyword_set(hmags) then printf,fplan,'hmag '''+strtrim(string(Hmags[0]),2)+''''    
    printf,fplan,'fiber ',stars[0]
    if keyword_set(hmags) then printf,fplan,'hmag ',Hmags[0]
    printf,fplan,'telliter      1'
    if suffix ne '' then printf,fplan,'mjdfrac     1'
endif
; platetype
if keyword_set(stars) then $
    printf,fplan,'platetype ''single''' $
else if keyword_set(cal) then $
  printf,fplan,'platetype ''cal''' $
else if keyword_set(sky) then $
  printf,fplan,'platetype ''sky''' $
else if keyword_set(dark) then $
  printf,fplan,'platetype ''dark''' $
else if keyword_set(test) then $
  printf,fplan,'platetype ''test''' $
else $
  printf,fplan,'platetype ''normal''' 

; note that q3fix is now done in ap3d.pro, not here!!
if mjd gt 56930L and mjd lt 57600L then printf,fplan,'q3fix    1'

file=apogee_filename('R',chip='a',num=ims[0])
if not file_test(file) then stop,'cant find file ', file
head=headfits(file,exten=1) 
plateid=sxpar(head,'PLATEID')
if not keyword_set(ignore) then if plate ne 0 and plateid ne plate then begin
  print,'plateid in header does not match plate!'
  stop
endif

; plugmap
print,keyword_set(plugid)
if not keyword_set(plugid) then begin
  file=apogee_filename('R',chip='a',num=ims[0])
  if file_test(file) eq 1 then begin
    head=headfits(file,exten=1)
    plugid=sxpar(head,'NAME')
    if size(plugid,/type) ne 7 then plugid='header'
  endif else plugid='header'
endif
print,ims[0]
print,plugid
;if not keyword_set(cal) and not keyword_set(dark) and not file_test(datadir+cmjd+'/plPlugMapA-'+plugid+'.par') then begin
if not keyword_set(cal) and not keyword_set(dark) and not keyword_set(onem) then begin
  tmp=strsplit(plugid,'-',/extract)
  if not file_test(mapper_data+'/'+tmp[1]+'/plPlugMapM-'+plugid+'.par') then begin
    print,'Cannot find plugmap file ', plugid
    spawn,'"ls" '+mapper_data+'/'+tmp[1]+'/plPlugMapA*'
    if not keyword_set(ignore) then stop
  endif
  if ~keyword_set(sky) then begin
    plug=getplatedata(cplate,cmjd,plugid=plugid,/noobj,mapper_data=mapper_data)
    cloc=strtrim(string(format='(i)',plug.locationid),2)
    file_mkdir,spectro_dir+'fields/'+telescope+'/'+cloc
    field=apogee_field(plug.locationid,plate,survey)
    printf,fplan,'survey   '+survey
    openw,file,spectro_dir+'fields/'+telescope+'/'+cloc+'/plan-'+cloc+'.lis',/get_lun,/append
    printf,file,telescope+'/'+cplate+'/'+cmjd+'/'+file_basename(planfile)
    free_lun,file
  endif
endif
printf,fplan,'plugmap '''+plugid+''''

; calibration frames to use
printf,fplan,'detid ', detid
printf,fplan,'bpmid ', bpmid
printf,fplan,'littrowid ', littrowid
printf,fplan,'persistid ', persistid
printf,fplan,'persistmodelid ', persistmodelid
printf,fplan,'darkid ', darkid
printf,fplan,'flatid ', flatid
printf,fplan,'sparseid ', sparseid
printf,fplan,'fiberid ', fiberid
printf,fplan,'badfiberid '+ string(format='(300i4)',badfiberid)
printf,fplan,'fixfiberid ',fixfiberid
printf,fplan,'psfid ', psfid
printf,fplan,'fluxid ', fluxid
printf,fplan,'responseid ', responseid
printf,fplan,'waveid ', multiwaveid
printf,fplan,'lsfid ', lsfid

; define plan structure
printf,fplan,'typedef struct {'
printf,fplan,' char plateid[20];'
printf,fplan,' int mjd;'
printf,fplan,' char flavor[8];'
printf,fplan,' char name[8];'
printf,fplan,' int single;'
printf,fplan,' char singlename[20];'
printf,fplan,'} APEXP;'
star=-1 & name='none'

; object frames
for i=0,n_elements(ims)-1 do begin
 if ims[i] gt 0 then begin
  if keyword_set(stars) then begin
    star=stars[i] & name=names[i]
  endif else begin
    star=-1 & name='none'
  endelse
  cid=string(format='(i8.8)',ims[i])
  printf,fplan,'APEXP '+cplate+' '+cmjd+' '+' object '+ $
    ' '+cid+string(format='(i6)',star)+' '+name
 endif
endfor
 
free_lun,fplan
file_chmod,planfile,'664'o

end
