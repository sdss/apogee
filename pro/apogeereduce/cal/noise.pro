; see runs at bottom

pro noise,ims,dark=dark,bpmid=bpmid,hard=hard,rn=rn,fidgain=fidgain

; estimate noise properties by looking at mean and variance in stack of images
; input : ims is 2D array giving [[start,n],[start,n],...] for multiple intflat data sets

dirs=getdir()
if keyword_set(dark) then begin
  xr=[1,100]
  yr=[1,25]
endif else begin
  xr=[100,25000]
  yr=[1,10000]
endelse
sz=size(ims,/dim)
if n_elements(sz) eq 1 then nset=1 else nset=sz[1]

chips=['a','b','c']
lab=['red','middle','blue']

if keyword_set(hard) then begin
 set_plot,'PS'
 device,file=hard+'.eps',/encap,/color
 smcolor,/ps
endif else begin
 set_plot,'X'
 smcolor
endelse

!p.multi=[0,2,3,0,0]
!p.charsize=2
for ichip=0,2 do begin
  chip=chips[ichip]
  if n_elements(bpmid) gt 0 then begin
    file=apogee_filename('BPM',chip=chips[ichip],num=bpmid)
    bpm=mrdfits(file)
    bd=where(bpm gt 0,nbd)
  endif else nbd=0

  first=1
  for indiv=1,1,2 do begin
    ; loop over 3D-2D options, note that these have implications for gain!
    for nfs=0,0 do begin

      icolor=1+nfs*2

      x=[]
      y=[]
      ymed=[]
      ydiff=[]
      npix=[]
      ;plot,[0],[0],/nodata,xr=[10,100],yr=[0,0.1]

      ; load all the images into a 3D stack
      for iset=0,nset-1 do begin
        stack=[]
        errstack=[]
        start=ims[0,iset]
        nims=ims[1,iset]
        maxread=ims[2,iset]
        ;if nset eq 1 then nims=size(ims,/dim) else nims=size(ims[*,iset],/dim)
        for i=0,nims-1 do begin
          ;if nset eq 1 then im=ims[i] else im=ims[i,iset]
          im=start+i
          cmjd=getcmjd(im)
          nreads=0
          d=process(cmjd,im,chip,head,red,dark,flat,err,gain,norn,mask,nfs=nfs,/nocr,indiv=indiv,horz=horz,vert=vert,maxread=maxread,nread=nreads)
          ; set bad pixels to -1
          if nbd gt 0 then begin
             d[bd]=-1
             err[bd]=-1
          endif
          stack=[[[stack]],[[d]]]
          errstack=[[[errstack]],[[err]]]
        endfor

        if nfs eq 0 then m=1 else m=nfs
        if nfs eq 0 then n=nreads else n=2
        rneff = rn[ichip]*sqrt(12.*(n-1)/ ( m*n*(n+1) ) )
        gfactor = 6. * (n^2+1)/ ( 5. * n * (n+1) ) - 2. * (2*m-1)*(n-1) / (m*n*(n+1)) * (m-1) / (m * (n-1))

        ; get the mean and standard deviation, and difference image from first two
        mn=mean(stack,dim=3)
        ; note that the stddev is biased, the variance is not
        var=stddev(stack,dim=3)^2
        e=median(errstack,dim=3)
        diff=stack[*,*,0]-stack[*,*,1]

        ; find all pixels at various intensity levels (logaritmically spaced) and
        ;   calculate average noise in stack and standard deviation over difference image
        dl=0.01
        for lbin=2.,4.4,dl do begin
          b1=10^lbin
          b2=10^(lbin+dl)
          bin=(b1+b2)/2.
          j=where(mn gt b1 and mn lt b2,nj)
          if nj gt 30000 then begin
            x=[x,bin]
            npix=[npix,nj]
            ; for stack, fit parabola to peak of histogram
            ;yy=median(var[j])
            ;yhist=histogram(var[j],binsize=2,loc=xhist)
            ;ymax=max(yhist,imax)
            ;fit=poly_fit(xhist[imax/2:imax*3/2],yhist[imax/2:imax*3/2],2)
            ;yy=-fit[1]/2/fit[2]
            ; median variance
            ymed=[ymed,median(var[j])]
            ; use mean variance
            y=[y,mean(var[j])]
            ; for difference image, take robust_sigma
            ydiff=[ydiff,robust_sigma(diff[j])]
            print,lbin,b1,b2,bin,nj,bin/mean(var[j]),2*bin/robust_sigma(diff[j])^2
            ;plothist,var[j],xr=[0,2*mean(var[j])]
            ;for i=0,9 do begin
            ;  s=stddev(stack[*,*,0:i+2],dim=3)^2
            ;  plothist,s[j],/over,color=(i mod 6) +2
            ;  stop
            ;endfor
            ;stop
          endif
        endfor
      endfor

      ; plot the results up
      if ~first then begin
        !p=p1
        !x=x1
        !y=y1
      endif
      if first then plot,x,y,xr=xr,yr=yr,psym=5,color=icolor,xstyle=1,xtitle='Mean counts',ytitle='Variance' else oplot,x,y,psym=5,color=icolor,xtickformat='(i5)'
      oplot,x,ymed,psym=1,color=icolor
      oplot,x,ydiff^2/2.,psym=6,color=icolor
      ;xyouts,30000.,1000.,'chip '+chips[ichip],align=1
      xyouts,30000.,8000.,lab[ichip],align=1
      p1=!p & x1 = !x & y1 = !y
      if ~first then begin
        !p =p2
        !x=x2
        !y=y2
      endif
      if first then plot,x,x/y*gfactor,/xlog,xr=xr,yr=[1,4],psym=5,color=icolor,xstyle=1,xtitle='Mean counts',ytitle='Mean/variance' else oplot,x,x/y*gfactor,psym=5,color=icolor
      oplot,x,x/ymed*gfactor,psym=1,color=icolor
      oplot,x,x/(ydiff/sqrt(2.))^2*gfactor,psym=6,color=icolor
      ;xyouts,20000.,1.5,'chip '+chips[ichip],align=1
      xyouts,20000.,3.5,lab[ichip],align=1
      if keyword_set(rn) then begin
        oplot,x,x/(y-rneff^2)*gfactor,psym=5,color=icolor+1
        oplot,x,x/(ymed-rneff^2)*gfactor,psym=1,color=icolor+1
        oplot,x,2.*x/(ydiff^2-rneff^2)*gfactor,psym=6,color=icolor+1
      endif
      if keyword_set(fidgain) then oplot,[xr[0],xr[1]],[fidgain,fidgain]

      p2=!p & x2 = !x & y2 = !y
      first=0
      icolor+=1
    endfor
  endfor
endfor

if keyword_set(hard) then begin
  device,/close
  ps2jpg,hard+'.eps',/eps
endif else stop

end

; rn is equivalent single read in DN
;apogee-n
apsetver,vers='current',telescope='apo25m'
n=10
start=10
noise,[[13360054+start,n,0],[13360054+start,n,5],[13360054+start,n,10],[13360054+start,n,15],[13360054+start,n,30]],bpmid=15640003,rn=[12,8,8],hard='plots/gain_apogee-n',fidgain=1.9
;apogee-s
apsetver,vers='current',telescope='lco25m'
n=10
;noise,[[22800309+start,n,3],[22800309+start,n,0],[22810006+start,n,0],[22810037+start,n,0],[22820002+start,n,0],[22660002+start,n,0]],bpmid=22620001,rn=[4,5,3],hard='plots/gain_apogee-s',fidgain=3
end
