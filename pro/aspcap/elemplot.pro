;+
; NAME:
;   elemplot
; PURPOSE:
;   Make plots of spectra/fits in windows for each element
; CALLING SEQUENCE:
;   
; INPUTS:
;   
; OPTIONAL INPUTS:
;   
; OPTIONAL KEYWORDS:
;   
; COMMENTS:
;   
; REVISION HISTORY:
;   
;-

pro elemplot,str,istar,xrange=xrange,hard=hard,script=script,maskdir=maskdir,altmaskdir=altmaskdir

if keyword_set(hard) then begin
  set_plot,'PS'
  smcolor,/ps
  ;openw,convert,hard+str.param[istar].apogee_id+'.csh',/get_lun
  ;printf,convert,'cd '+hard
endif else begin
  set_plot,'X'
  smcolor
endelse
wave=10^str.lib.wave
spec=str.spec[istar].spec
;model=str.spec[istar].elem_bestfit
elem_mask=str.spec[istar].elem_mask
s=strsplit(hard,'/',/extract)
n=n_elements(s)
s[n-3]=s[n-2]
n-=2
froot=strjoin(s[0:n-1],'/')
field=file_basename(froot)

; need to set first pixel off or else ranges are screwed up
elem_mask[0]=0
plot,wave,str.spec[istar].spec,xrange=xrange

; loop over each element
if tag_exist(str.param,'file') then name=strtrim(str.param[istar].file,2) else name=strtrim(str.param[istar].apogee_id,2)
for i=0,n_elements(str.lib.elem_symbol)-1 do begin
 el=strtrim(str.lib.elem_symbol[i],2)
 if file_test('/'+froot+'/elem_'+el+'/'+el+'-'+strtrim(str.param[istar].class,2)+'-'+field+'.spm') then begin
  load,'/'+froot+'/elem_'+el+'/'+el+'-'+strtrim(str.param[istar].class,2)+'-'+field+'.spm',fspm
  load,'/'+froot+'/elem_'+el+'/'+el+'-'+strtrim(str.param[istar].class,2)+'-'+field+'.mdl',fmodel
  j=where(strtrim(fspm[0,*],2) eq name)
  model=fmodel[*,j]
  if keyword_set(hard) then device,file=hard+name+'.'+el+'.eps',/encap,/color,xsize=16,ysize=10,/in
  if keyword_set(maskdir) then readcol,maskdir+'/'+el+'.mask',mask,/silent
  if keyword_set(altmaskdir) then readcol,altmaskdir+'/'+el+'.filt',altmask,/silent
  imask=elem_mask and 2L^i
  istart=where(shift(imask,-1)-imask gt 0,nj)
  iend=where(shift(imask,1)-imask gt 0)
;print,el,nj,istart,iend
;stop
;if el eq 'Nd' then stop
  erase
  if nj le 6 then !p.multi=[0,3,2] else if nj le 12 then !p.multi=[0,4,3] else if  nj le 24 then !p.multi=[0,6,4] else if nj le 32 then !p.multi=[0,4,8] else begin
   for j=0,2 do begin
     junk=min(abs(str.lib.wavemin[j]-str.lib.wave),ipix) & istart[j]=ipix
     junk=min(abs(str.lib.wavemax[j]-str.lib.wave),ipix) & iend[j]=ipix
     nj=3 
     !p.multi=[0,1,3] 
   endfor
  endelse
  ;erase & multiplot,[4,3]
  if nj gt 0 then begin
   for j=0,nj-1 do begin
     is=max([0,istart[j]])
     ie=min([n_elements(spec)-1,iend[j]])
     plot,wave,spec,xtickformat='(i5)',charsize=1.75,yr=[0.5,1.25],thick=2,ystyle=1,xr=[wave[is]-2,wave[ie]+2],xstyle=1
     xyouts,wave[is]-1,0.5,el,charsize=2
     oplot,wave,model,color=2,thick=2
     if keyword_set(altmaskdir) then oplot,wave,altmask*0.5+0.5,color=4,thick=2
     if keyword_set(maskdir) then oplot,wave,mask*0.5+0.5,color=3,thick=2
     ;multiplot
   endfor
   if keyword_set(hard) then begin
     device,/close
     if keyword_set(script) then begin
     printf,script,'convert '+name+'.'+str.lib.elem_symbol[i]+'.eps '+name+'.'+el+'.jpg'
     printf,script,'rm -f '+name+'.'+el+'.eps'
     endif else begin
       ps2jpg,hard+name+'.'+el+'.eps',/eps,/delete
     endelse
   endif
  endif
  ;multiplot,/reset
 
;  j=where(((mask and 2^i) gt 0) or ((shift(mask,2) and 2^i) gt 0) or ((shift(mask,-2) and 2^i) gt 0),nj)
;  if nj gt 0 then begin
;    oplot,wave[j],spec[j],color=i+2
;    first=0
;  endif
 endif
endfor

;if keyword_set(hard) then begin
;  free_lun,convert
;  spawn,'csh '+hard+str.param[istar].apogee_id+'.csh'
;  file_delete,hard+str.param[istar].apogee_id+'.csh'
;endif
end
