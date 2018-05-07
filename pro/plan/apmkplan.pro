;+
; apmkplan:  automatically make apPlan files for a given MJD
;-

pro write_plateplan,out,plate,cmjd,exp,sky,dome,planfiles

 ; procedure to write out the information for a given plate, once the
 ;   variables have been loaded
 cplate=strtrim(string(format='(i6.4)',plate),2)
 printf,out,'plate='+cplate
 printf,out,'psfid='+dome
 printf,out,'fluxid='+dome
 printf,out,'ims=['+exp+']'
 printf,out,'mkplan,ims,plate,mjd,psfid,fluxid,vers=vers'
 printf,out,''
 planfile=apogee_filename('Plan',plate=cplate,mjd=cmjd,/base)
 if n_elements(planfiles) eq 0 then planfiles=planfile else planfiles=[planfiles,planfile]
 if n_elements(sky) gt 0 then begin
   printf,out,'ims=['+sky+']'
   printf,out,'mkplan,ims,plate,mjd,psfid,fluxid,vers=vers,/sky'
   planfile=file_basename(planfile,'.par')+'sky.par'
   if n_elements(planfiles) eq 0 then planfiles=planfile else planfiles=[planfiles,planfile]
 endif

end

pro apmkplan,mjd,planfiles=planfiles,apogees=apogees,vers=vers

if ~keyword_set(vers) then stop,'need to set vers'

if keyword_set(apogees) then begin
  prefix='as'
  telescope='lco25m'
  instrument='apogee-s'
  apogee_data=getenv('APOGEE_DATA_2S')
endif else begin
  prefix='ap'
  telescope='apo25m'
  instrument='apogee-n'
  apogee_data=getenv('APOGEE_DATA')
endelse

; main procedure to scan through apz files and collect exposure data

; get character MJD and search for files
cmjd=string(mjd,format='(i5.5)')
openw,out,getenv('APOGEEREDUCEPLAN_DIR')+'/pro/'+telescope+'/'+telescope+'_'+cmjd+'auto.pro',/get_lun

printf,out,"apsetver,telescope='"+telescope+"'"
;printf,out,"vers='"+vers+"'"
printf,out,'mjd='+cmjd

files=file_search(apogee_data+'/'+cmjd+'/*-c*.apz')

; loop over all files, accumulating the IDs of the different types
;vers='junk'
ndark=0
ncal=0
nexp=0
nsky=0
oldplate=0
for i=0,n_elements(files)-1 do begin
 h=headfits(files[i],ext=1)
 ; extract image number from file name
 s=strsplit(files[i],'-',/ext)
 num=0L
 ss=strsplit(s[2],'.',/ext)
 reads,ss[0],num
 ; get image exptype, plateid, and nreads from header
 exptype=strtrim(sxpar(h,'exptype'),2)
 plate=strtrim(sxpar(h,'plateid'),2)
 nread=sxpar(h,'nread')
 
 print,files[i],plate,oldplate

 ; if there has been a plate change, write the plan file from previous plate
 if plate ne oldplate and nexp gt 0 and n_elements(dome) gt 0 then begin
   write_plateplan,out,oldplate,cmjd,exp,sky,dome,planfiles
   nexp=0
   undefine,sky
 endif

 ; load image number in variable according to exptype and nreads
 ; discard images with nread<=3
 if exptype eq 'DARK' and sxpar(h,'NREAD') ge 3 then  begin
   if ndark eq 0 then dark=string(format='(i8.8)',num) else $
                     dark=dark+','+string(format='(i8.8)',num)
   ndark+=1
 endif
 if exptype eq 'INTERNALFLAT' and sxpar(h,'NREAD') ge 3 then  begin
   ; internal flats are reduced only to 2D, hence treated like darks
   if ndark eq 0 then dark=string(format='(i8.8)',num) else $
                     dark=dark+','+string(format='(i8.8)',num)
   ndark+=1
 endif
 if exptype eq 'QUARTZFLAT' and sxpar(h,'NREAD') ge 3 then  begin
   if ncal eq 0 then cal=string(format='(i8.8)',num) else $
                     cal=cal+','+string(format='(i8.8)',num)
   calpsfid=string(format='(i8.8)',num)
   ncal+=1
 endif
 if exptype eq 'ARCLAMP' and sxpar(h,'NREAD') gt 3 then  begin
   if ncal eq 0 then cal=string(format='(i8.8)',num) else $
                     cal=cal+','+string(format='(i8.8)',num)
   ncal+=1
 endif
 if exptype eq 'OBJECT' and (nread lt 15  and nread gt 10)  then  begin
   ; identify sky frames as object frames with 10<nread<15
   if nsky eq 0 then sky=string(format='(i8.8)',num) else $
                     sky=sky+','+string(format='(i8.8)',num)
   nsky+=1
 endif
 if exptype eq 'OBJECT' and sxpar(h,'NREAD') gt 15 then  begin
   if nexp eq 0 then exp=string(format='(i8.8)',num) else $
                     exp=exp+','+string(format='(i8.8)',num)
   nexp+=1
 endif
 if exptype eq 'DOMEFLAT' and sxpar(h,'NREAD') gt 3 then  begin
   dome=string(format='(i8.8)',num)
 endif

 oldplate=plate

endfor

; write out the last plate if we haven't already
if n_elements(dome) eq 0 then dome=0
if nexp gt 0 then write_plateplan,out,plate,cmjd,exp,sky,dome,planfiles

; write out the dark/calibration frame information
cplate='0000'
if n_elements(dark) gt 0 then begin
  printf,out,'plate=0'
  printf,out,'psfid=0'
  printf,out,'fluxid=0'
  printf,out,'ims=['+dark+']'
  printf,out,'mkplan,ims,plate,mjd,psfid,fluxid,vers=vers,/dark'
  planfile='apPlan-'+cplate+'-'+cmjd+'dark.par'
  planfile=apogee_filename('DarkPlan',mjd=cmjd,instrument=instrument,/base)
  if n_elements(planfiles) eq 0 then planfiles=planfile else planfiles=[planfiles,planfile]
endif
if n_elements(cal) gt 0 and n_elements(calpsfid) gt 0 then begin
  printf,out,''
  printf,out,'psfid='+calpsfid
  printf,out,'fluxid='+calpsfid
  printf,out,'ims=['+cal+']'
  printf,out,'mkplan,ims,plate,mjd,psfid,fluxid,vers=vers,/cal'
  planfile='apPlan-'+cplate+'-'+cmjd+'.par'
  planfile=apogee_filename('CalPlan',mjd=cmjd,instrument=instrument,/base)
  if n_elements(planfiles) eq 0 then planfiles=planfile else planfiles=[planfiles,planfile]
endif
free_lun,out
end
