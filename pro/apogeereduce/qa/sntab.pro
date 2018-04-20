pro sntab,exptab=exptab,all=all,tabs=tabs,badims=badims,extra=extra,outfile=outfile,altsn=altsn

dirs=getdir(apodir,caldir,spectrodir,apovers)
if keyword_set(outfile) then openw,1,outfile else $
  openw,1,apodir+'reduction_sn_new.dat'

if keyword_set(badims) then begin
  for i=0,n_elements(badims)-1 do begin
    ;printf,1,format='(i,"	",f,"	",a,"   ",a)',badims[i],-1.,'manual' ,badims[i]
    printf,1,format='(i,"	",f,"	",a)',badims[i],-1.,'manual' 
  endfor
  close,1
  return
endif

if keyword_set(all) then begin
  spawn,'"ls" apo25m/*/*/apPlateSum*.fits',tabs
endif else if not keyword_set(tabs) then begin
  print,'you must either specify TABS=tabfiles or /ALL'
  stop
endif

if keyword_set(exptab) then begin
  readcol,'/net/apogee/exposure.dat',mjd,plate,snr,exposure_no,exposure_time,pk,format='(l,i,f,l,f,l)',delim='|'
  openw,2,'/net/apogee/noexp_new.dat'
  openw,3,'/net/apogee/multiexp_new.dat'
endif

x=0
y=0
vers=getenv('APOGEE_VERS')
for i=0,n_elements(tabs)-1 do begin

  a=mrdfits(tabs[i],1,/silent)
  for j=0,n_elements(a)-1 do begin
    pluginfo=strsplit(a[j].plugid,'-',/extract)
    if keyword_set(altsn) then sn=a[j].altsn[1] else sn=a[j].sn[1]
    if keyword_set(extra) then $
      printf,1,format='(i,"	",f,"	",a,"   ",a)',a[j].im,sn,vers,tabs[i] $
    else printf,1,format='(i,"	",f10.2,"	",a,"	",i4,"	",i8,"	",i6,"	",i8,"	",f24.6,"	",f10.2,"	",a)',$
      a[j].im,sn,vers,pluginfo[2],pluginfo[1],pluginfo[0],a[j].mjd,date_conv(a[j].dateobs,'MODIFIED')*86400.d0,a[j].nreads*10.647,'Object'

    if keyword_set(exptab) then begin
      k=where(exposure_no eq a[j].im)
      if k[0] lt 0 then printf,2,a[j].im
      if n_elements(k) gt 1 then printf,3,a[j].im
      if n_elements(k) eq 1 and k[0] ge 0 then begin
        if keyword_set(extra) then $
          printf,1,format='(i,i,i,i,"	",f,"	",a)',a[j].im,plate[k],mjd[k],pk[k],sn,vers $
        else print,format='(i,i,f,f)',pk[k],a[j].im,sn,snr[k]
        x=[x,sn]
        y=[y,snr[k]]
      endif
    endif
  endfor
endfor
close,1
if keyword_set(exptab) then begin
  close,2
  close,3
endif

end

