pro aspcap_plateresults,obsID,class,res,name,path=path,star=star,$
png=png,vradcut=vradcut,ext=ext,interord=interord,version=version,ps=ps,single=single,$
sncut=sncut,plot=plot,append=append,versdat=versdat,obs_path=obs_path,$
out_path=out_path,apvisit=apvisit,zoom=zoom,libr_path=libr_path,queue=queue,ppath=ppath,oobspath=oobspath



; INPUT:
;
;      obsID           -string       the name specifying the
;                                    observation
;                                    (apVisit-obsID-fiberID.fits)
;
;
; KEYWORD:
;
;      apVisit         -            if set, the files will be open in
;                                   subdirectories specified by the
;                                   plate
;                                    and date





print,'out',out_path
if not keyword_set(ppath) then stop,'must set ppath'

if n_elements(vradcut) eq 0 then vradcut=3d10
if n_elements(sncut) eq 0 then sncut=-1
if n_elements(ncols) eq 0 then ncols=7
if n_elements(append) gt 0 then begin
  filename='all_data'
  if n_elements(out_path) gt 0 then filename=out_path+'/'+filename
  if n_elements(ext) gt 0 then filename=filename+'-'+ext
  openw,numf,filename+'.txt',/get_lun,/append
endif


; Set plots properties

!x.thick=3
!y.thick=3
!z.thick=3

xxrange=[1.585*1d4,1.585*1d4,1.645*1d4,1.645d4]


print,ppath+'/'+class+'-'+obsID+'.spm'

; Get the stellar parameters

value=file_test(ppath+'/'+class+'-'+obsID+'.spm')
print,'File',ppath+'/'+class+'-'+obsID+'.spm'
if value ne 0 then begin

    sschi2='log'+'!4' + String("166B) + '!X!U2!N!lr!N'
    ssmicro='!4' + String("156B) + '!X!lt!N'

    if n_elements(ps) gt 0 then begin
       sschi2='log'+'!9' + String("143B) + '!X!U2!N!lr!N'
       ssmicro='!9' + String("170B) + '!X!lt!N'
       alpha='!9' + String("141B) + '!X'
       
    endif

   ; output text file and HTML file (note with doaspcap, this is redone with aspcap_writehtml!)
    file_mkdir,out_path+'/plots'
    ;openw,num,out_path+'/'+obsID+'_interpol.txt',/get_lun
    openw,html,out_path+'/plots/'+obsID+'.html',/get_lun
    printf,html,'<HTML><HEAD><script type=text/javascript src=../../../../../../sorttable.js></script></head>'
    printf,html,'<BODY>'
    printf,html,'Click on column labels to sort'
    printf,html,'<TABLE BORDER=2 CLASS=sortable>'
    printf,html,'<TR><TD>Object<TD>S/N<TD>Chi^2<TD>P1<TD>P2<TD>P3<TD>P4<TD>P5<TD>P6<TD>P7'
    pfilename=out_path+'/plots/'+class+'-'+obsID
      
    ; get results 
    filename=ppath+'/'+class+'-'+obsID
    aspcap_strresults,filename,results
    aspcap_labels,filename,labels,slabel,npars,ps=ps,path=libr_path
    snpars=string(npars,format='(I1)')
    iindex=where(labels eq ssmicro,ncts)
    if ncts gt 0 then begin
       results.spm[iindex+1,*]=10.^(results.spm[iindex+1,*])
    endif
   
    if n_elements(zoom) gt 0 then aspcap_plotobssyn_zoom,results,ipar=ipar,$
      label=slabel,npars=npars,path=ppath,/ps
    invar=1d/(double(results.err)^2d)
    nfibers=n_elements(results.spm[0,*])
    files=strcompress(results.spm[0,*],/remove_all)
  
  ; Get the coordinates and radial velocities

    ii=-1
    index=-9999 & indexvrad=index
    vrad=dblarr(nfibers) & longi=vrad & lati=vrad & objtype=strarr(nfibers) 
    res_i={files:strarr(nfibers),vrad:vrad[0],vraderr:vrad[0],longi:longi[0],lati:lati[0],objtype:objtype[0],vrspm:replicate(0d,9),sn:-9999.,nfibers:0L}
    res_i=replicate(res_i,nfibers)

    ; Information for plots
    if n_elements(ps) gt 0  then begin
       set_plot,'ps'
       if keyword_set(single) then begin
         openw,cfile,/get_lun,out_path+'/plots/convert.csh'
         printf,cfile,'#!/bin/bash'
         printf,cfile,'#PBS -l select=1:ncpus=1:mem=100MB'
         printf,cfile,'#PBS -l walltime=99:00:00'
         printf,cfile,'#PBS -W group_list=apogee'
         printf,cfile,'#PBS -q apogee'
         printf,cfile,'#PBS -W umask=002'
         printf,cfile,'#PBS -N convert'
         printf,cfile,'#PBS -V'
         printf,cfile,'#PBS -z'
         printf,cfile,'#PBS -j oe'
       endif else device,filename=pfilename+'.ps',/landscape,/color,/isolatin      
       loadct,39,/silent
    endif
    ;printf,num,format='(A32,1X,A20,1X,2(A8,1X),A9,1X,A8,1X,5(A8,1X),'+snpars+'((A11,1X)),36X,A8,36X,1X,A12)',$
    ;         'File name','Star name','Gal. long.','Gal. latit','Vhel','SN','Jmag','Hmag','Kmag','JK0','H0mag',labels,'sigmas','log chi2'


;'[Fe/H]','[C/Fe]','[N/Fe]',$             '[alpha/Fe]','log micro','Teff','logg','sigmas','log chi2'


   for i=0,nfibers-1 do begin
      
; Get coordinates, vrad and photometry

      jpos=strpos(files[i],'_v')
      if jpos gt 0 then fitsname = strmid(files[i],0,jpos) else fitsname=files[i]
      aspcap_getinf4plt,oobspath+'/'+fitsname+'.fits',objtyp_i,objid_i,vvrad,vhel,vraderr,llongi,llati,Jmag,Hmag,Kmag

      posobj=strpos(objid_i,"'") & objid_i=strmid(objid_i,0,posobj)
      H0mag=-99 & JK0=-99
;      indexH0=where(strcompress(all.objid,/remove_all) eq objid_i,ncts)
;      if ncts gt 0 then begin
;         H0mag=all[indexH0].H0 
;         JK0=all[indexH0].jk0
;      endif

      if vvrad lt vradcut then begin             
         match,files[i],strcompress(results.ipf[0,*],/remove_all),index2,index1
         ; Estimate S/N
         sn=aspcap_snmedian(alog10(results.wav[*,index1]), results.frd[*,index1], invar[*,index1])
         if sn[0] gt sncut then begin 
 
            if keyword_set (single) then begin
              psfile=out_path+'/plots/'+files[i]+'.eps'
              jpgfile=out_path+'/plots/'+files[i]+'.jpg'
              device,filename=psfile,/color,xsize=48,ysize=16,/encap,/port
            endif

            ; text file output
            ;printf,num,format='(A32,1X,A20,1X,2(F8.3,1X),F9.3,1X,F8.1,1X,5(F8.3,1X),2('+snpars+'(F11.5,1X)),1X,F8.2)',$
            ;files[i],objid_i,llongi,llati,vhel,sn[0],Jmag,Hmag,Kmag,JK0,H0mag,results.spm[1:2*npars,i],results.spm[npars*2+3,i]
            if n_elements(append) gt 0 then printf,numf,format='(A32,1X,A20,1X,2(F8.3,1X),F9.3,1X,F8.1,1X,5(F8.3,1X),2('+snpars+'(F11.5,1X)),1X,F8.2)',$
            files[i],objid_i,llongi,llati,vhel,sn[0],Jmag,Hmag,Kmag,JK0,H0mag,results.spm[1:2*npars,i],results.spm[npars*2+3,i] 

            ; make three-panel plots
            !p.multi=[0,1,3]
            yrange=[-0.2,1.25]
            xrange=[min(results.wav),xxrange,max(results.wav)]/10.
            !p.color=0
            ;print,'plot',n_elements(plot)
            if n_elements(plot) gt 0 then aspcap_plotobssyn,files[i],results.wav[*,index1],results.frd[*,index1],results.syn[*,i],results.err[*,index1],results.spm[1:npars,i],$
            sn[0],objid_i,vvrad,vhel,vraderr,Jmag,Hmag,Kmag,xrange=xrange,yrange=yrange,ps=ps,label=slabel,obspath=oobspath,chi2=10.^results.spm[npars*2+3,i]

            if n_elements(ps) eq 0 then wait,1
         
            ii=ii+1L
            res_i[ii].vrad=vvrad
            res_i[ii].vraderr=vraderr
            res_i[ii].objtype=objtyp_i
            res_i[ii].longi=llongi
            res_i[ii].lati=llati
            res_i[ii].sn=sn[0]
            index=[index,i]
            if keyword_set (single) then begin
              device,/close
              printf,cfile,'echo '+psfile
              printf,cfile,'gs -sDEVICE=jpeg -sOutputFile='+jpgfile+' -r150 -dBATCH -dNOPAUSE -DDEVICEWIDTHPOINTS=1360 -dDEVICEHEIGHTPOINTS=453 '+psfile+'  >& /dev/null'
              printf,cfile,'''rm'' '+psfile
              printf,cfile,'chmod 664 '+jpgfile
            endif

            ; make single panel plot
            psfile=out_path+'/plots/'+files[i]+'_1.eps'
            jpgfile=out_path+'/plots/'+files[i]+'_1.jpg'
            device,filename=psfile,/color,xsize=60,ysize=8,/encap,/port
            gd=where(results.frd[*,index1[0]] gt 0)
            !p.multi=[0,0,0]
            plot,results.wav[gd,index1[0]]*1.e-4,results.frd[gd,index1[0]],yrange=[0,1.05],ystyle=9,xrange=[1.5100,1.7000],xstyle=1,xtitle='Wavelength',ytitle='Flux'
            oplot,results.wav[*,index1[0]]*1.e-4,results.syn[*,i],color=250
            plot,results.wav[gd,index1[0]]*1.e-4,(results.frd[gd,index1[0]]-results.syn[gd,i])^2/results.err[gd,index1[0]]^2,xrange=!x.crange,yrange=[0,50],/noerase,xstyle=5,ystyle=4,color=128,xtitle='',ytitle=''
            axis,YAXIS=1, YRANGE=!y.crange, ytitle='Chi^2',color=128
            if keyword_set(single) then begin
              device,/close
              printf,cfile,'echo '+psfile
              printf,cfile,'convert '+psfile+' '+jpgfile
              printf,cfile,'''rm'' '+psfile
              printf,cfile,'chmod 664 '+jpgfile
            endif

            ;html file output (note with doaspcap, this is redone with aspcap_writehtml!)
            printf,html,'<TR><TD><a href=apStar-'+strtrim(objid_i,2)+'.jpg>'+objid_i,'</a><TD>'+string(format='(f8.3)',sn[0])
            printf,html,'<TD>',10.^(results.spm[npars*2+3,i]),format='(a,f8.2)'
            for ipar=0,npars-1 do printf,html,'<TD>',results.spm[ipar+1,i],format='(a,6f8.2)'
            printf,html,'<TD><IMG src='+file_basename(jpgfile)+'>'

         endif
      endif
   endfor

   ;close,num & free_lun,num
   free_lun,html

   if n_elements(ps) gt 0 then begin
      cd,out_path+'/plots',current=cdir
      if keyword_set(single) then begin
        free_lun,cfile
	print,'converting plots to jpg...'
        if keyword_set(queue) then spawn,'qsub convert.csh' else spawn,'csh convert.csh'
        file_delete,'convert.csh',/allow_nonexistent
      endif else begin
        device,/close
        pstopdf,pfilename+'.ps'
      endelse
      cd,cdir
   endif

   if ii gt -1 then begin

         index=index[1:*]
         spm=[results.spm[2*npars+3,index],transpose([res_i[0:ii].sn]),results.spm[1:npars,index]]
         errspm=[transpose(replicate(0.,n_elements(index))),transpose(replicate(0.,n_elements(index))),double(results.spm[npars+1:npars*2,index])] 

  ; Add the radial velocity
         vrspm=[transpose([res_i[0:ii].vrad]),double(spm)]    
         errvrspm=[transpose([res_i[0:ii].vraderr]),errspm]
         res=res_i[index]

; Plots with parameters
      if n_elements(parplot) gt 0 then begin
         if n_elements(ps) gt 0 then begin
            set_plot,'ps'
            loadct,39,/silent
            device,filename=pfilename+'_hist.ps',/landscape,/color,/isolatin
         endif

         aspcap_parplots,vrspm,errvrspm,'Obs '+obsID,ps=ps,png=png,labelref=['Vrad [km/s]',sschi2,'SN',labels]
         aspcap_parhistog,vrspm,png=png,ps=ps,longitude=string(mean(res_i[ii].longi),format='(F6.2)'),$
         latitude=string(mean(res_i[ii].lati),format='(F6.2)'),labelref=['Vrad [km/s]',sschi2,'SN',labels]

         if n_elements(ps) gt 0 then begin
            device,/close
            set_plot,'x'
            loadct,39,/silent
            cd,out_path+'/plots',current=cdir
            pstopdf,pfilename+'_hist.ps',/deleteps
            cd,cdir
        endif
       
      endif
    endif else begin
           res={files:'NONE'}
    endelse

endif else begin 
   print,ppath+'/'+class+'-'+obsID+'.spm not found'
   res={files:'NONE'}
endelse 

if n_elements(append) gt 0 then begin
  close,numf & free_lun,numf
endif
end
