;+
;
; ASPCAPRUN
;
; This program runs various APOGEE pipeline modules
; Modules that are run depend on input parameters
;
; INPUTS:
;  planfile   The name of the plan file to reduce
;  flag       Character bitwise flag giving steps to call
;  clobber    Whether to clobber existing files or use
;              them if they exist
;
; OUTPUTS:
;  The raw data will be reduced and output to
;  respective directories.
;
; USAGE:
;  IDL>aspcaprun,flag,planfile,clobber
;
; By J.Holtzman  Feb 2012
;-

pro aspcaprun,planfile,flag,sclobber,noelem=noelem,noplot=noplot,elemplot=elemplot,cal=cal

clobber=intarr(1)
reads,sclobber,clobber

bytes=byte(flag)
dim=n_elements(bytes)
charflag=strarr(dim)
for i=0,dim-1 do charflag[i] = string(bytes[i])
bin=Total( (Byte(Reverse(charflag)) EQ 49) * 2L^Indgen(dim) )
bin=long(bin)

print,'bin: ', bin
print,'flag: ', flag
print,'planfile: ', planfile
print,'clobber: ', clobber

aploadplan,planfile,planstr
override=0
if getenv('APOGEE_OVERRIDE_VERSION') eq '1' then override=1
if (not override) and tag_exist(planstr,'apogee_ver') then $
 if planstr.apogee_ver ne getenv('APOGEE_VER') and getenv('APOGEE_OVERRIDE_VERSION') ne getenv('APOGEE_VER') then $
 stop,'APOGEEREDUCE version does not match planfile!'

if keyword_set(cal) then aspcap_calibrate,planfile,caldir=cal $
else doaspcap,planfile,clobber=clobber,noelem=noelem,noplot=noplot,doelemplot=elemplot,/no_version_check
;doaspcap,planfile
print,'doaspcap completed successfully'

end
