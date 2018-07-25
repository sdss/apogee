; procedure to derive linearity correction from an internal flat field frame

function mklinearity,frameid,clobber=clobber,chip=chip, stp=stp,lindata=lindata,nread=nread,minread=minread,norder=norder,nskip=nskip,inter=inter

if ~keyword_set(minread) then minread=2
if ~keyword_set(norder) then norder=2
if ~keyword_set(nskip) then nskip=4
cref=3000.

; get directories
dirs=getdir(apodir,caldir,specdir,apovers,libdir,datadir=datadir)

; character frameid
cframeid=string(format='(i8.8)',frameid)

; get calibratation file names for this MJD
cmjd=getcmjd(frameid,mjd=mjd)
getcal,mjd,libdir+'/cal/'+dirs.instrument+'.par',dark=darkid,bpm=bpmid,det=detid

; chip= keyword specifies, single chip, else use all 3 chips
chips=['a','b','c']
if n_elements(chip) gt 0 then begin
  ichip1=chip
  ichip2=chip
endif else begin
  ichip1=0
  ichip2=2
endelse

; get name of file for output linearity data
lindir=caldir+'/detector/'
file_mkdir,lindir+'/plots/'
if keyword_set(lindata) then $
linfile=lindir+dirs.prefix+'LinearityTest-'+cframeid+'.dat' else $
linfile=lindir+dirs.prefix+'Linearity-'+cframeid+'.dat'

; make sure file construction isn't already in process
while file_test(linfile+'.lock') do apwait,linfile,10

; does file already exist? 
if file_test(linfile) and ~keyword_set(clobber) then goto, havefile

; loop over the chips
openw,lun,linfile,/get_lun
for ichip=ichip1,ichip2 do begin

 ; uncompress data cube
 file=dirs.prefix+'R-'+chips[ichip]+'-'+cframeid
 info=apfileinfo(datadir+cmjd+'/'+file+'.apz')
 if keyword_set(nread) then info.nreads=nread

 if size(getlocaldir(),/type) eq 7 then $
   outdir=getlocaldir() else outdir='./'
 if not file_test(outdir+'/'+file+'.fits') then apunzip,datadir+cmjd+'/'+file+'.apz',fitsdir=getlocaldir()

 ; read the cube
 for i=1,info.nreads do begin
  FITS_READ,getlocaldir()+file+'.fits',im,head,exten_no=i,message=message,/no_abort   ; UINT
 
  if i eq 1 then cube=long(im) else cube=[[[cube]],[[im]]] 
  help,cube
 endfor

 ; do reference correction (is this right?)
 tmp=aprefcorr(cube,head,mask,indiv=0,/cds)
 cube=tmp

 ; if we have input linearity data, we will use it to test that things are working!
 if keyword_set(lindata) then begin
   oldcube=cube
   for iy=0,2047 do begin
     if iy mod 10 eq 0 then print,'linearity...',iy
     slice = float(reform(cube[*,iy,*]))
;if iy eq 1000 then aplincorr,slice,lindata,slice_out,/stp else $
     aplincorr,slice,lindata,slice_out
;if iy eq 1000 then stop
     cube[*,iy,*]=reform(slice_out,2048,1,info.nreads)
   endfor
 endif

  ; loop over different sections on chip
 for ix=0,39,5 do begin
  ;ix1=4+ix*50
  ;ix2=4+(ix+1)*50
  ix1=24+ix*50
  ix2=ix1+10
  for iy=0,39,5 do begin
   ; counts in section
   ;iy1=4+iy*50
   ;iy2=4+(iy+1)*50
   iy1=24+iy*50
   iy2=iy1+10

   ; get median in region
   cts=fltarr(info.nreads-2)
   rate=fltarr(info.nreads-2)
   instrate=fltarr(info.nreads-2)
   for i=2,info.nreads-1,nskip do begin
     cts[i-2] = median(float(cube[ix1:ix2,iy1:iy2,i])-cube[ix1:ix2,iy1:iy2,1])
     rate[i-2]=cts[i-2]/(i-1.)
     ; correct to "zero" read
     cts[i-2]*=(i+1.)/(i-1.)
     instrate[i-2]=median(float(cube[ix1:ix2,iy1:iy2,i])-cube[ix1:ix2,iy1:iy2,i-1])
   endfor

   ; normalize to rate at cref DN
   j=where(cts gt cref-2000 and cts lt cref+2000,nj)
   if min(cts[j]) gt cref or max(cts[j]) lt cref then nj=0
   if nj gt 2 then begin
     par=poly_fit(cts[j],rate[j],2)
     ref=par[0]+par[1]*cref+par[2]*cref^2
     for i=2,info.nreads-1,nskip do begin
       printf,lun,i,ichip,ix,iy,cts[i-2],rate[i-2]/ref,instrate[i-2]/ref,format='(i3,i3,i5,i5,3f12.4)'
       print,i,ichip,ix,iy,cts[i-2],rate[i-2]/ref,instrate[i-2]/ref,format='(i3,i3,i5,i5,3f12.4)'
     endfor
   endif

  endfor
 endfor
endfor
free_lun,lun

havefile:
; read the linearity data
readcol,linfile,read,chip,ix,iy,cts,rate,rate2,format='(i,i,i,f,f,f)'

if ~keyword_set(inter) then set_plot,'PS' else set_plot,'X'
!p.multi=[0,1,2]
smcolor

for ichip=ichip1,ichip2 do begin

gd=where(chip eq ichip)
ymax=50
if ichip eq 2 and dirs.instrument eq 'apogee-n' then ymax=18 
if ichip eq 2 and dirs.instrument eq 'apogee-n' then ymax=0
print,ichip,ymax

; make some plots
if ~keyword_set(inter) then begin
  if keyword_set(lindata) then $
  file=lindir+'plots/'+dirs.prefix+'LinearityTest-'+cframeid+'_'+string(format='(i1.1)',ichip)+'.eps' else $
  file=lindir+'plots/'+dirs.prefix+'Linearity-'+cframeid+'_'+string(format='(i1.1)',ichip)+'.eps'
  device,file=file,/encap,/color,xsize=12,ysize=12,/in
endif

; plot of instantaneous rate vs counts
plot,cts[gd],rate2[gd],psym=6,yr=[0.8,1.2],xr=[0,max(cts)],thick=3,charthick=3,charsize=1.5,xtit='DN',ytit='Relative count rate'
ii=1
for i=0,35,5 do begin
 j=where(ix[gd] eq i and iy[gd] lt ymax,nj)
 if nj gt 1 then oplot,cts[gd[j]],rate2[gd[j]],psym=6,color=(ii mod 6)+1
 ii+=1
endfor

; plot of accumulated rate, normalized by final rate, vs counts
; since illumination isn't uniform, this normalized rate is not the same for
;  all regions. Do a fit to get rate at cref DN
plot,cts[gd],rate[gd],psym=6,yr=[0.8,1.2],xr=[0,max(cts)],thick=3,charthick=3,charsize=1.5,xtit='DN',ytit='Relative count rate'
ii=1
for i=0,35,5 do begin
 j=where(ix[gd] eq i and iy[gd] lt ymax,nj)
 if nj gt 1 then oplot,cts[gd[j]],rate[gd[j]],psym=6,color=(ii mod 6)+1
 ii+=1
endfor

if ~keyword_set(inter) then begin
  device,/close
  ps2jpg,file,/eps
endif else stop
endfor

; now do the final linearity fit using all regions in chip a and non-persistence;  region of chip c
ii=0
for ichip=0,2 do begin
 ymax=50
 if ichip eq 1 and dirs.instrument eq 'apogee-n' then ymax=0 
 if ichip eq 2 and dirs.instrument eq 'apogee-n' then ymax=18
 gd=where((chip eq ichip) and (iy lt ymax) and (read ge minread) and (cts lt 50000),ngd)
 if ngd gt 0 then begin
  if ii eq 0 then begin
    x=cts[gd]
    y=rate[gd]
  endif else begin
    x=[x,cts[gd]]
    y=[y,rate[gd]]
  endelse
  ii+=1
    
 endif
endfor

; final fit and plot
if keyword_set(lindata) then $
file=lindir+'plots/'+dirs.prefix+'LinearityTest-'+cframeid+'.eps' else $
file=lindir+'plots/'+dirs.prefix+'Linearity-'+cframeid+'.eps' 
set_plot,'PS'
device,file=file,/encap,/color,xsize=12,ysize=8,/in
!p.multi=[0,0,0]

plot,x,y,psym=6,yr=[0.9,1.1],xtit='DN',ytit='Relative count rate', thick=2,charthick=2
par=poly_fit(x,y,norder)
xx=indgen(5000)*10.
yy=par[0]
tmp=xx
for iorder=1,norder do begin
  yy+=par[iorder]*tmp
  tmp*=xx
endfor
oplot,xx,yy,color=3,thick=3
device,/close
ps2jpg,file,/eps
set_plot,'X'

if keyword_set(stp) then stop

return,par

end

apsetver,vers='current',telescope='apo25m'
n=mklinearity(23900003,/clobber,nskip=2)
ntest=mklinearity(23900003,/clobber,nskip=2,lindata=n)
;apsetver,vers='current',telescope='lco25m'
;s=mklinearity(22820002,/clobber,nskip=1)
;stest=mklinearity(22820002,/clobber,nskip=1,lindata=s)
end
