; makes some diagnostic plots for dark frames
pro darkplot,a,mask,darkfile,hard=hard

 sz=size(a)
 nr= sz[3]
 if keyword_set(hard) then begin
   set_plot,'ps'
   device,file=darkfile+'.eps',/color,/encap,xsize=20,ysize=20
 endif else set_plot,'X'
 smcolor
 !p.multi=[0,2,2,0,0]
 ; histogram of dark rates for each slice
 i1=2 & i2=nr-1 & del=1
 i1=(nr-1)/2
 i2=nr-1
 del=(nr-1)/2
 for i=i1,i2,del do begin
   if i eq i1 then plothist,a[*,*,i]-a[*,*,1],xh,yh,/nan,xrange=[0,100],xtitle='total dark',ytitle='Npixels' $
   else plothist,a[*,*,i]-a[*,*,1],xh,yh,/nan,/over,color=2+(i mod 6)
 endfor
 ; histogram of ratio of dark rate from full divided by dark rrate from first half
 d1=(a[*,*,nr-1]-a[*,*,1])
 d2=(a[*,*,(nr-1)/2]-a[*,*,1])
 d=where(d2 gt 100 and (mask and badmask()) eq 0)
 plothist,(d1[d]/d2[d]),/nan,bin=0.1,xrange=[0,4]

  ; logarithmic plot of dark rate
 for i=i2,i1,-del do begin
  print,i
   if i eq i2 then plothist,a[*,*,i]-a[*,*,1],xh,yh,/nan,/xlog,/ylog,xtitle='dark rate',ytitle='Npixels' $
   else plothist,a[*,*,i]-a[*,*,1],xh,yh,/nan,/over,color=2+(i mod 6)
 endfor

 ; cumulative histogram
 for i=i2,i1,-del do begin
   plothist,a[*,*,i]-a[*,*,1],xh,yh,/nan,xrange=[0,100],/noplot
   for j=1L,n_elements(yh)-1 do yh[j]+=yh[j-1]
   if i eq i2 then plot,xh,yh/2048./2048.,xrange=[0,50],yrange=[0.5,1.1] $
   else oplot,xh,yh/2048./2048.,color=2+(i mod 6)
 endfor
 if keyword_set(hard) then begin
   device,/close
   spawn,'convert '+darkfile+'.eps '+darkfile+'.gif'

   device,file=darkfile+'2.eps',/encap,xsize=20,ysize=20
   d=a[*,*,nr-1]-a[*,*,1]
   tvscl,d>(-20)<100
   device,/close
   spawn,'convert '+darkfile+'2.eps '+darkfile+'2.jpg'
   set_plot,'X'
 endif else stop
 ;file_delete,darkfile+'.eps'
end

