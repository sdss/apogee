pro aspcap_plotobssyn,filename,wave,frd,syn,err,spm,sn,objid,vvrad,vhel,vraderr,Jmag,Hmag,Kmag,$
xrange=xrange,yrange=yrange,ps=ps,obspath=obspath,label=label,chi2=chi2

loadct,39,/silent

!x.thick=5
!y.thick=5
!z.thick=5
!p.charsize=1.0

if n_elements(obspath) eq 0 then obspath='data'
angstrom=string(229b)

; Define the labels for the plot
if n_elements(ps) gt 0 then begin
     sschi2='!9' + String("143B) + '!X!U2!N!lr!N'
     ssmicro='!9' + String("170B) + '!X!lt!N'
     alpha='!9' + String("141B) + '!X'
  endif else begin 
     alpha='!4' + String("141B) + '!X'
     sschi2='!4' + String("166B) + '!X!U2!N!lr!N'
     ssmicro='!4' + String("156B) + '!X!lt!N
endelse
if n_elements(label) eq 0 then begin
 label='[Fe/H]/[C/Fe]/[N/Fe]/['+alpha+'/Fe]/'+ssmicro+'/T!leff!N/logg/'   ;+sschi2
endif

sspm=aspcap_stripar(spm,label)

; Define the range for the plots
if n_elements(yrange) eq 0 then yrange=[-0.2,1.2]

;if n_elements(ps) gt 0 then begin
;      set_plot,'ps'
;      device,filename='ASPCAP_obsvsyn.ps',/landscape,/color,/isolatin
;endif

!p.color=0

;  Loop over the detectors
gd=where(frd gt 0)
for i=0,2 do begin
        file=obspath+'/'+filename+'.fits'

        psave=!p
        if i eq 0 then begin
             plot,wave[gd]/10.,frd[gd],xtitle='Wavelength [nm]',ytitle='Flux',title=filename+'('+objid+')',yrange=yrange,charsize=2,$
             BACKGROUND = 255, COLOR = 0,xrange=[xrange[2*i],xrange[2*i+1]],xstyle=1,ystyle=9,psym=10
        endif else begin
             plot,wave[gd]/10.,frd[gd],xtitle='Wavelength [nm]',ytitle='Flux',yrange=yrange,charsize=2,$
             BACKGROUND = 255, COLOR = 0,xrange=[xrange[2*i],xrange[2*i+1]],xstyle=1,ystyle=9,psym=10
        endelse
        oplot,wave/10.,syn,color=250
        resid=(frd-syn)/syn
        oplot,wave[gd]/10.,resid[gd],lin=1

        ; calculate direct chi^2
        mychi2=TOTAL((frd[gd]-syn[gd])^2/err[gd]^2)/n_elements(gd)

        ; add information in labels
        if i eq 0 then begin
           !p.color=50
           al_legend,'SN='+string(sn,format='(F5.1)')+' H='+string(Hmag,format='(F5.2)')+' J-Ks='+strcompress(string(Jmag-Kmag,format='(F5.2)'),/remove_all)+' chi^2'+string(format='(2f8.2)',chi2,mychi2),$
           color=50,box=0,posit=[1520,0.3],/data,charsize=1;/top,/center
           !p.color=0
           al_legend,['obs','synth','rel. resid'],line=[0,0,1],color=[0,250,0],/center,/right,box=0,charsize=1.0
        endif

        if i eq 1 then begin
           !p.color=50
           al_legend,'Vrel= '+strcompress(string(vvrad,format='(F7.2)'),/remove_all)+' Vhelio='$
           +strcompress(string(vhel,format='(F7.2)'),/remove_all)+' +/-'$
           +strcompress(string(vraderr,format='(F7.2)'),/remove_all),color=50,box=0,posit=[1604,0.3],/data,charsize=1
           !p.color=0

        endif

        if i eq 2 then begin
           !p.color=50
           al_legend,label,position=[0.45,0.20],/normal,box=0,charsize=1.1,textcolor=50
           al_legend,sspm,position=[0.45,0.16],box=0,/normal,charsize=1.1,textcolor=50
           !p.color=0
        endif

        ; add chi^2 overplot
        !p=psave
        gd=where(frd gt 0)
        plot,wave[gd]/10.,(frd[gd]-syn[gd])^2/err[gd]^2,xrange=!x.crange,yrange=[0,100],/noerase,xstyle=5,ystyle=4,color=128,xtitle='',ytitle='',charsize=2
        axis,YAXIS=1, YRANGE=!y.crange, ytitle='Chi^2',color=128
        if !p.multi[0] eq 0 then !p.multi[0]=2 else if !p.multi[0] eq 2 then !p.multi[0]=1 else !p.multi[0]=0

        if n_elements(ps) eq 0 then wait,1
endfor

;  if n_elements(ps) gt 0 then begin
;    device,/close & set_plot,'x' 
;  endif

end
