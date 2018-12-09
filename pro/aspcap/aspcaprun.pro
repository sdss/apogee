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

pro aspcaprun,planfile,flag,sclobber,noplot=noplot,elemplot=elemplot

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

doaspcap,planfile,clobber=clobber,noplot=noplot,doelemplot=elemplot
;doaspcap,planfile
print,'doaspcap completed successfully'

end
