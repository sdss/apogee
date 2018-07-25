
pro getrn,ims,cmjd=cmjd,noref=noref,indiv=indiv

; attempts to determine readout noise from several methods
;
;  rn1:  readout noise from single exposure (does not account for any bias structure)
;  rn2:  readout noise from difference of two images
;  rn3:  readout noise from difference of successive reads in single image
;
;  In all cases, noise is determined from various Fowler samplings as well as UTR sampling
;  The values are reduced to single-read readout noise based on the Rauscher et al formula

dirs=getdir(apogeedir,caldir,spectrodir,vers)

if n_elements(ims) ne 2 then begin
  print, 'ims should have 2 frames only'
  stop
endif

; output structure
rn={ im1: ims[0], im2: ims[1], n: intarr(6), m: intarr(6), rn1: fltarr(4,6), rn1corr: fltarr(4,6), rn2: fltarr(4,6), rn2corr: fltarr(4,6), rn3: fltarr(4,6), rn4: fltarr(4,6)}
rnlog=REPLICATE(rn,3)

; loop over chips
chips=['a','b','c']
for ichip=0,2 do begin
 if not keyword_set(cmjd) then cmjd=getcmjd(ims[0])
 d=process(cmjd,ims[0],chips[ichip],head,r,step=step,/nofs,/nofix,/nocr,noref=noref,indiv=indiv)
 if not keyword_set(cmjd) then cmjd=getcmjd(ims[1])
 d2=process(cmjd,ims[1],chips[ichip],head,r,step=step,/nofs,/nofix,/nocr,noref=noref,indiv=indiv)
 ; also read cube with any processing
 ;file=apogee_filename('R',num=ims[0],chip=chips[ichip])
 ;apunzip,file,fitsdir='./'
 ;cube=mrdfits(file_basename(file,'.apz')+'.fits',0)

 ; cut at 47 reads to get equivalent noise for typical exposure
 sz=size(d,/dim)
 if sz[2] gt 47 then begin
   d=d[*,*,0:46]
   d2=d2[*,*,0:46]
 endif

 ; loop over Fowler sampling options (0 = UTR)
 for i=0,5 do begin 
   sz=size(d,/dim)
   nreads=sz[2]
   if i eq 0 then m=1 else m=i
   if i eq 0 then n=nreads else n=2
   f=fsamp(d[*,*,0:46],i)
   f2=fsamp(d2[*,*,0:46],i)
   for j=0,3 do begin
     print,'FS: ', i,' Quadrant: ',j

     ; quadrant
     x1=j*512+5 & x2=(j+1)*512-5
     rnlog[ichip].n[i] = n
     rnlog[ichip].m[i] = m
     ; rn1 from single exposure gives single read readout noise in DN
     rnlog[ichip].rn1[j,i] = robust_sigma(f2[x1:x2,10:2040])
     rnlog[ichip].rn1corr[j,i] = rnlog[ichip].rn1[j,i] * sqrt(m * n * (n+1.) / 12. / (n-1.))

     ; rn2 from difference of 2 images, divide by sqrt(2) to account for difference 
     diff=f[x1:x2,10:2040]-f2[x1:x2,10:2040]
     rnlog[ichip].rn2[j,i] = robust_sigma(diff)/sqrt(2.)
     ; reduce to single read rn given fowler sampling per Rauscher et al.
     ; m as in Rauscher et al, with n=2, for Fowler sampling
     rnlog[ichip].rn2corr[j,i] = rnlog[ichip].rn2[j,i] * sqrt(m * n * (n+1.) / 12. / (n-1.))
     gd=where(finite(diff) eq 1)
     plothist,diff[gd],bin=1,xr=[-100,100]
     print,rnlog[ichip].rn2[j,i]

     ; successive reads readout noise, fitting histogram
     if i eq 1 then begin
       data=float(d[x1:x2,10:2040,*])
       del=data[*,*,1:*]-data[*,*,0:nreads-1]
       hist=histogram(del,bin=1,min=-100,max=100,locations=xhist)
       xhist+=0.5*1
       yfit=mpfitpeak(xhist,hist,par)
       print,par[2],par[2]*sqrt(2),par[2]/sqrt(2)*1.9,robust_sigma(del)
       ; rn3 from difference between successive reads, so divide by sqrt(2) to account for difference 
       rnlog[ichip].rn3[j,i]=robust_sigma(del)/sqrt(2.)
     endif
  
     ; successive reads readout noise, 100 columns from cube w/o process
     ;data=float(cube.data[x1:x2,10:2040,*])
     ;del=data[*,*,1:*]-data[*,*,0:nreads-1]
     ;hist=histogram(del,bin=1,min=-100,max=100,locations=xhist)
     ;xhist+=0.5*1
     ;yfit=mpfitpeak(xhist,hist,par)
     ;print,par[2],par[2]*sqrt(2),par[2]/sqrt(2)*1.9,robust_sigma(del)
     ;rnlog[ichip].rn4[j,i]=robust_sigma(del)

   endfor
 endfor
; print,'successive readout read noise: chip,quadrant, read, sigman'
; sz=size(d,/dim)
; for iread=2,sz[2]-1 do begin
;   for j=0,3 do begin
;     x1=j*512+5 & x2=(j+1)*512-5
;     print,ichip,j,iread,robust_sigma(d[x1:x2,10:2040,iread]-d[x1:x2,10:2040,iread-1])
;   endfor
; endfor
endfor
 

; output and plots
file=string(format='("apRN-",i8.8)',ims[0])
if keyword_set(noref) then file=file+'noref'
if keyword_set(indiv) then file=file+'indiv'+string(format='(i1)',indiv)
chip=['a','b','c']
hard=1
if hard eq 1 then begin
  set_plot,'ps'
  device,/color,file=caldir+'/detector/'+file+'.eps',/encap,xsize=30,ysize=20
endif
!p.multi=[0,2,3,0,0]
!p.charsize=2
smcolor
for ichip=0,2 do begin
  x=indgen(5)+1
  plot,x,rnlog[ichip].rn1[0,1:5],psym=1,xrange=[-1,6],yrange=[0,20],xtitle='Nfowler',$
    ytitle=textoidl('\sigma_{eff RN}')+' (DN)',thick=2,charthick=2,xthick=2,ythick=2,xstyle=1
  xyouts,5.5,17,'chip '+chip[ichip],align=1.,charsize=1.
  xyouts,-1,17,' Effective readout noise (DN)',charsize=1.
  for j=0,3 do begin
    oplot,x,rnlog[ichip].rn1[j,1:5],psym=1,color=j+2,thick=2
    oplot,[0],[rnlog[ichip].rn1[j,0]],psym=1,color=j+2,thick=2
  endfor
  ;oplot,x,rnlog[ichip].rn2[0,1:5],psym=6,xrange=[-1,6],yrange=[0,20],xtitle='Fowler',$
  ;  ytitle='Readout noise (DN)',thick=2,charthick=2,xthick=2,ythick=2
  for j=0,3 do begin
    oplot,x,rnlog[ichip].rn2[j,1:5],psym=6,color=j+2,thick=2
    oplot,[0],[rnlog[ichip].rn2[j,0]],psym=6,color=j+2,thick=2
  endfor
  for j=0,3 do begin
    oplot,x,rnlog[ichip].rn3[j,1:5],psym=5,color=j+2,thick=2
    oplot,[0],[rnlog[ichip].rn3[j,0]],psym=5,color=j+2,thick=2
  endfor
  ; plot with corrected noise, i.e. backed out to single read
  plot,x,rnlog[ichip].rn1corr[0,1:5],psym=1,xrange=[-1,6],yrange=[0,20],xtitle='Nfowler',$
    ytitle=textoidl('\sigma_{read}')+' (DN)',thick=2,charthick=2,xthick=2,ythick=2,xstyle=1
  xyouts,5.5,17,'chip '+chip[ichip],align=1.,charsize=1.
  xyouts,-1,17,' Derived single-read readout noise (DN)',charsize=1.
  for j=0,3 do begin
    oplot,x,rnlog[ichip].rn1corr[j,1:5],psym=1,color=j+2,thick=2
    oplot,[0],[rnlog[ichip].rn1corr[j,0]],psym=1,color=j+2,thick=2
  endfor
  for j=0,3 do begin
    oplot,x,rnlog[ichip].rn2corr[j,1:5],psym=6,color=j+2,thick=2
    oplot,[0],[rnlog[ichip].rn2corr[j,0]],psym=6,color=j+2,thick=2
  endfor
  ;oplot,x,rnlog[ichip].rn3[0,1:5],psym=5,xrange=[-1,6],yrange=[0,20],xtitle='Fowler',$
  ;  ytitle='Readout noise (DN)',thick=2,charthick=2,xthick=2,ythick=2
  for j=0,3 do begin
    oplot,x,rnlog[ichip].rn3[j,1:5],psym=5,color=j+2,thick=2
    oplot,[0],[rnlog[ichip].rn3[j,0]],psym=5,color=j+2,thick=2
  endfor
endfor
if hard eq 1 then begin
  device,/close
  if not file_test(caldir+'/detector/plots',/directory) then file_mkdir,caldir+'/detector/plots'
  spawn,'convert '+caldir+'/detector/'+file+'.eps '+caldir+'/detector/plots/'+file+'.jpg'
  ;file_delete,caldir+'/detector/'+file+'.eps'
  set_plot,'x'
endif
mwrfits,rnlog,caldir+'/detector/'+file+'.tab',/create

; single red chip plot
if hard eq 1 then begin
  set_plot,'ps'
  device,/color,file=caldir+'/detector/'+file+'_0.eps',/encap,xsize=30,ysize=20
endif
!p.multi=[0,2,3,0,0]
!p.charsize=2
smcolor
for ichip=0,2 do begin
  x=indgen(5)+1
  plot,x,rnlog[ichip].rn1[0,1:5],psym=1,xrange=[-1,6],yrange=[0,20],xtitle='Nfowler',$
    ytitle=textoidl('\sigma_{eff RN}')+' (DN)',thick=2,charthick=2,xthick=2,ythick=2,xstyle=1
  xyouts,5.5,17,'chip '+chip[ichip],align=1.,charsize=1.
  xyouts,-1,17,' Effective readout noise (DN)',charsize=1.
  ;for j=0,3 do begin
  ;  oplot,x,rnlog[ichip].rn1[j,1:5],psym=1,color=j+2,thick=2
  ;  oplot,[0],[rnlog[ichip].rn1[j,0]],psym=1,color=j+2,thick=2
  ;endfor
  ;oplot,x,rnlog[ichip].rn2[0,1:5],psym=6,xrange=[-1,6],yrange=[0,20],xtitle='Fowler',$
  ;  ytitle='Readout noise (DN)',thick=2,charthick=2,xthick=2,ythick=2
  for j=0,3 do begin
    oplot,x,rnlog[ichip].rn2[j,1:5],psym=6,color=j+2,thick=2
    oplot,[0],[rnlog[ichip].rn2[j,0]],psym=6,color=j+2,thick=2
  endfor
  ;for j=0,3 do begin
  ;  oplot,x,rnlog[ichip].rn3[j,1:5],psym=5,color=j+2,thick=2
  ;  oplot,[0],[rnlog[ichip].rn3[j,0]],psym=5,color=j+2,thick=2
  ;endfor
  ; plot with corrected noise, i.e. backed out to single read
  plot,x,rnlog[ichip].rn1corr[0,1:5],psym=1,xrange=[-1,6],yrange=[0,20],xtitle='Nfowler',$
    ytitle=textoidl('\sigma_{read}')+' (DN)',thick=2,charthick=2,xthick=2,ythick=2,xstyle=1
  xyouts,5.5,17,'chip '+chip[ichip],align=1.,charsize=1.
  xyouts,-1,17,' Derived single-read readout noise (DN)',charsize=1.
  ;for j=0,3 do begin
  ;  oplot,x,rnlog[ichip].rn1corr[j,1:5],psym=1,color=j+2,thick=2
  ;  oplot,[0],[rnlog[ichip].rn1corr[j,0]],psym=1,color=j+2,thick=2
  ;endfor
  for j=0,3 do begin
    oplot,x,rnlog[ichip].rn2corr[j,1:5],psym=6,color=j+2,thick=2
    oplot,[0],[rnlog[ichip].rn2corr[j,0]],psym=6,color=j+2,thick=2
  endfor
  ;oplot,x,rnlog[ichip].rn3[0,1:5],psym=5,xrange=[-1,6],yrange=[0,20],xtitle='Fowler',$
  ;  ytitle='Readout noise (DN)',thick=2,charthick=2,xthick=2,ythick=2
  ;for j=0,3 do begin
  ;  oplot,x,rnlog[ichip].rn3[j,1:5],psym=5,color=j+2,thick=2
  ;  oplot,[0],[rnlog[ichip].rn3[j,0]],psym=5,color=j+2,thick=2
  ;endfor
endfor
if hard eq 1 then begin
  device,/close
  if not file_test(caldir+'/detector/plots',/directory) then file_mkdir,caldir+'/detector/plots'
  spawn,'convert '+caldir+'/detector/'+file+'_0.eps '+caldir+'/detector/plots/'+file+'_0.jpg'
  ;file_delete,caldir+'/detector/'+file+'.eps'
  set_plot,'x'
endif

rnhtml,caldir

end

pro rnhtml,caldir

rnfiles=file_search(caldir+'/detector/apRN*.tab')
if not file_test(caldir+'/detector/html',/directory) then file_mkdir,caldir+'/detector/html'
openw,lun,/get_lun,caldir+'/detector/html/rn.html'

printf,lun,'<html><body>'
printf,lun,'Readout noise is given IN UNITS OF DN (as measured)'
printf,lun,'<table border=2>'
printf,lun,'<TR><TD>Image 1<TD>Image 2<TD>Chip'
for ifs=0,5 do begin
  for iquad=0,3 do begin
    printf,lun,'<TD>FS: ',ifs,' QUAD: ',iquad
  endfor
endfor
for i=0,n_elements(rnfiles)-1 do begin
  rnlog=mrdfits(rnfiles[i],1)
  for ichip=0,2 do begin
    if ichip eq 0 then printf,lun,'<TR><TD>',file_basename(rnfiles[i]),'<TD>',rnlog[ichip].im1,'<TD>',rnlog[ichip].im2 else $
      printf,lun,'<TR><TD><TD><TD>'
    printf,lun,'<TD>',ichip
    for ifs=0,5 do begin
     for iquad=0,3 do begin
       printf,lun,'<TD>',string(format='(4f8.2)',rnlog[ichip].rn2[iquad,ifs]),rnlog[ichip].rn3[iquad,ifs],rnlog[ichip].rn4[iquad,ifs]
     endfor
    endfor
  endfor
endfor
printf,lun,'</table></body></html>'
free_lun,lun

;openw,lun,/get_lun,'rn.html'
;printf,lun,'<HTML><BODY>'
;
;printf,lun,'Readout noise measurements, frames: ',im1,im2,'<br>'
;printf,lun,'Left panels are single frame measurements, i.e. rms of a quadrant on a single dark subtracted dark frame, using robust_sigma<br>'
;printf,lun,'Right panels are difference frame measurements, i.e. rms of a quadrant in the difference between two dark frames<br>'
;printf,lun,'Squares are measured noise, triangles are "corrected" to DCS (i.e. multiplied by sqrt(nfowler)) <br>'
;printf,lun,'Different colors are different quadrants<br>'
;printf,lun,'<b> UNITS ARE DN </b><br>'
;printf,lun,'<i> Comments: <br>'
;printf,lun,'difference frame readout noise is a bit smaller than single frame, suggesting some remaining bias structure (that repeats in these two frames)<br>'
;printf,lun,'increasing fowler sampling/UTR decreases noise almost, but not quite as much, as expected<br></i>'
;printf,lun,'<img src=rn.gif>'
;printf,lun,'<p>Readout noise measurements from frame ',im1
;printf,lun,' <b> in DN </b>; robust sigma after superdark subtraction'
;printf,lun,'<TABLE BORDER=2>'
;printf,lun,'<TR><TD>chip<TD>Fowler<TD>quad1<TD>quad2<TD>quad3<TD>quad4'
;for ichip=0,2 do begin
;;  for i=0,5 do begin
;   printf,lun,'<TR><TD>',ichip,'<TD>',i
;   for j=0,3 do printf,lun,format='(a,f8.2)','<TD>',rn[j,i,ichip]
;  endfor
;endfor
;printf,lun,'</table>'
;printf,lun,'<p>Readout noise measurements from difference frame ',im1, im2
;printf,lun,' <b> in DN </b>; robust sigma after superdark subtraction'
;printf,lun,'<TABLE BORDER=2>'
;printf,lun,'<TR><TD>chip<TD>Fowler<TD>quad1<TD>quad2<TD>quad3<TD>quad4'
;for ichip=0,2 do begin
;  for i=0,5 do begin
;   printf,lun,'<TR><TD>',ichip,'<TD>',i
;   for j=0,3 do printf,lun,format='(a,f8.2)','<TD>',rn2[j,i,ichip]/sqrt(2)
;  endfor
;endfor
;printf,lun,'</table><p>'
;
;free_lun,lun
;

;free_lun,lun
;  if keyword_set(dolog) then begin
;   set_plot,'ps'
;   device,file='bias.eps',/encap,/color
;   !p.multi=[0,2,2,0,0]
;   smcolor
;   m=fltarr(4,3)
;   s=fltarr(4,3)
;   for i=0,3 do begin
;    for fs=3,1,-1 do begin
;      d=fsamp(red,fs)
;      sect=d[i*512+5:(i+1)*512-5,5:2040]
;      good=where(abs(sect) lt 100 and sect ne 0)
;      m[i,fs-1]=mean(sect[good])
;      s[i,fs-1]=stddev(sect[good])
;      if fs eq 3 then plothist,sect[good],xrange=[-100,100],color=fs+1 $
;      else plothist,sect[good],/over,color=fs+1
;    endfor
;   endfor
;   device,/close
;   file=string(format='("bias-",a,"-",i8.8)',chip,num)
;   spawn,'convert bias.eps '+file+'.gif'
;   for fs=1,3 do print,num,chip[ichip],fs,m[*,fs-1],s[*,fs-1]
;   printf,lun,'<TR><TD>',num,'<TD>',chip[ichip],'<TD><TABLE border=1><TR><TD>NF<TD>Q1<TD>Q2<TD>Q3<TD>Q4'
;   for fs=1,3 do begin
;     line=string(format='(a,i2,a,f8.1,a,f8.1,a,f8.1,a,f8.1)','<TR><TD>',fs,'<TD>',m[0,fs-1],'<TD>',m[1,fs-1],'<TD>',m[2,fs-1],'<TD>',m[3,fs-1])
;     printf,lun,line
;   endfor
;   printf,lun,'</TABLE><TD><TABLE border=1><TR><TD>NF<TD>Q1<TD>Q2<TD>Q3<TD>Q4'
;   for fs=1,3 do begin
;     line=string(format='(a,i2,a,f8.1,a,f8.1,a,f8.1,a,f8.1)','<TR><TD>',fs,'<TD>',s[0,fs-1],'<TD>',s[1,fs-1],'<TD>',s[2,fs-1],'<TD>',s[3,fs-1])
;     printf,lun,line
;   endfor
;   printf,lun,'</TABLE><TD><IMG SRC='+file+'.gif>'
;  endif
end

apsetver,vers='current',telescope='apo25m'
;getrn,[15640004,15640005],indiv=1
;getrn,[15640004,15640005],indiv=3
;getrn,[12910009,12910010],indiv=1
getrn,[12910009,12910010],indiv=3
apsetver,telescope='lco25m'
;getrn,[22620002,22620003],indiv=1
;getrn,[22620002,22620003],indiv=3
end
