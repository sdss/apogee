;+
;
; APRED
;
; This program runs various APOGEE pipeline modules
; Modules that are run depend on input parameters
;
; INPUTS:
;  planfile   The name of the plan file to reduce
;  flag       Character bitwise flag giving steps to call
;             for reduction option:
;                    1  ap3d
;                   10  ap2d
;                  100  ap1dvisit
;                 1000  apqa
;             for cal options:
;                    1  dark
;                   10  flat
;                  100  bpm
;                 1000  wave
;                10000  lsf
;  clobber    Whether to clobber existing files or use
;              them if they exist
;
; OUTPUTS:
;  The raw data will be reduced and output to
;  respective directories.
;
; USAGE:
;  IDL>apred,flag,planfile,clobber
;
; By J.Holtzman  Feb 2012
;-

pro apred,planfile,flag,sclobber

if n_elements(flag) eq 0 then flag = '1111111'
if n_elements(sclobber) eq 0 then sclobber = '1'

clobber=intarr(1)
reads,sclobber,clobber

bytes=byte(flag)
dim=n_elements(bytes)
charflag=strarr(dim)
for i=0,dim-1 do charflag[i] = string(bytes[i])
bin=Total( (Byte(Reverse(charflag)) EQ 49) * 2L^Indgen(dim) )
bin=long(bin)

aploadplan,planfile,planstr
override=0
if getenv('APOGEE_OVERRIDE_VERSION') eq '1' then override=1
if getenv('APOGEE_OVERRIDE_VERSION') ne '' then vers = getenv('APOGEE_OVERRIDE_VERSION') else vers=getenv('APOGEE_VER')
if (not override) and tag_exist(planstr,'apogee_ver') then $
 if planstr.apogee_ver ne vers and planstr.apogee_ver ne 'test' then  stop,'APOGEEREDUCE version does not match planfile!'
if tag_exist(planstr,'telescope') then telescope=planstr.telescope else telescope='apo25m'
if tag_exist(planstr,'apred_vers') then apsetver,vers=planstr.apred_vers,telescope=telescope

print,'bin: ', bin
print,'flag: ', flag
print,'planfile: ', planfile
print,'clobber: ', clobber

if planfile eq 'cal.par' then begin
print,'cal.par'
  if (bin and 1) ne 0 then makecal,/dark
  if (bin and 2) ne 0 then makecal,/flat
  if (bin and 4) ne 0 then makecal,/bpm
  if (bin and 8) ne 0 then makecal,/wave
  if (bin and 16) ne 0 then makecal,/lsf

endif else if planfile eq 'dailycal.par' then begin
print,'dailycal.par'
  if (bin and 1) ne 0 then makecal,/dark,file=planfile
  if (bin and 2) ne 0 then makecal,/flat,file=planfile
  if (bin and 4) ne 0 then makecal,/bpm,file=planfile
  if (bin and 8) ne 0 then makecal,/wave,file=planfile
  if (bin and 16) ne 0 then makecal,/lsf,file=planfile

endif else if strpos(planfile,'.cal') ge 0 then begin
print,'.cal'
  cmjd=strsplit(planfile,'.',/extract)
  mjd=cmjd[0]
  if (bin and 1) ne 0 then makecal,/dark,file='dailycal.par',mjd=mjd
  if (bin and 2) ne 0 then makecal,/flat,file='dailycal.par',mjd=mjd
  if (bin and 4) ne 0 then makecal,/bpm,file='dailycal.par',mjd=mjd
  if (bin and 8) ne 0 then makecal,/wave,file='dailycal.par',mjd=mjd
  if (bin and 16) ne 0 then makecal,/lsf,file='dailycal.par',mjd=mjd

endif else if strpos(planfile,'apPlan-') ge 0 or strpos(planfile,'apCalPlan-') ge 0 or strpos(planfile,'apDarkPlan-') ge 0 or $
              strpos(planfile,'asPlan-') ge 0 or strpos(planfile,'asCalPlan-') ge 0 or strpos(planfile,'asDarkPlan-') ge 0 then  begin
print,'.par'
  if (bin and 1) ne 0 then ap3d,planfile,clobber=clobber
  if (bin and 2) ne 0 then ap2d,planfile,clobber=clobber
  if (bin and 4) ne 0 then ap1dvisit,planfile,clobber=clobber
  if (bin and 8) ne 0 then apqa,planfile

endif else if strpos(planfile,'apStar-') ge 0 or strpos(planfile,'asStar-') ge 0 then  begin
print,'apstar'
  apstar,planfile,/sinc,/log,clobber=clobber

endif else if strpos(planfile,'apMJD-') ge 0 or strpos(planfile,'asMJD-') ge 0 then  begin
print,'apMJD'
  ;aploadplan,planfile,planstr
  ;if tag_exist(planstr,'apred_vers') then apsetver,planstr.apred_vers
  mkhtml,planstr.mjd

endif else begin
  print,'Unrecognized plan option: ', planfile
endelse

print,'apred completed successfully'

end
