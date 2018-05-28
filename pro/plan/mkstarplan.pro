pro mkstarplan,fields,apred_vers=apred_vers,apstar_vers=apstar_vers,mjdstart=mjdstart,mjdend=mjdend

 ; make plan file for apstar

 apsetver,vers=apred_vers
 dirs=getdir(a,c,s)
 outdir=a+'/'+apred_vers+'/'+apstar_vers+'/plan'
 file_mkdir,outdir

 for i=0,n_elements(fields)-1 do begin

 print,fields[i]
 field=file_basename(fields[i])
 telescope=file_basename(file_dirname(fields[i]))
 if telescope eq 'apo1m' then begin
   survey='apo1m'
   cfield=field 
 endif else begin
;   loc=0L
;   reads,field,loc
;   cfield=apogee_field(loc,0,survey)
   cfield=field
   apgundef,survey
   locid=apogee_locationid(field,survey)
 endelse
print,field
print,telescope
print,survey
 apsetver,telescope=telescope
 dirs=getdir()
 ;if size(field,/type) ne 7 then cfield=strtrim(string(format='(i)',field),2) else cfield=field
 openw,newplan,outdir+'/'+dirs.prefix+'Star-'+cfield+'.par',/get_lun
 printf,newplan,'apogee_ver  ',getenv('APOGEE_VER')
 printf,newplan,'apred_vers '''+apred_vers+''''
 printf,newplan,'apstar_vers '''+apstar_vers+''''
 if keyword_set(mjdstart) then printf,newplan,'mjdstart ',mjdstart
 if keyword_set(mjdend) then printf,newplan,'mjdend ',mjdend
 printf,newplan,'telescope  '+telescope
 printf,newplan,'survey  '+survey
 printf,newplan,'typedef struct {'
 printf,newplan,' char field[24];'
 printf,newplan,'} APFIELD'
 printf,newplan,'APFIELD '+cfield
 free_lun,newplan
 endfor

end
