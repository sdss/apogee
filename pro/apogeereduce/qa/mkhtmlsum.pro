pro mkhtmlsum,apred_vers=apred_vers,all=all,nocheck=nocheck,apstar_vers=apstar_vers,aspcap_vers=aspcap_vers,results_vers=results_vers,mjdmin=mjdmin,mjdmax=mjdmax,fieldfile=fieldfile,mjdfile=mjdfile,nomjd=nomjd,suffix=suffix

; make summary HTML pages both by MJD and by field

apsetver,vers=apred_vers
dirs=getdir(apodir,caldir,spectrodir,vers)
if not keyword_set(apred_vers) then stop,'Need to specify apred_vers'
if not keyword_set(nocheck) then nocheck=0
if not keyword_set(apstar_vers) then apstar_vers='stars'
if not keyword_set(aspcap_vers) then aspcap_vers='ASPCAP_VERS'
if not keyword_set(results_vers) then results_vers='RESULTS_VERS'
if not keyword_set(mjdmin) then mjdmin=0L
if not keyword_set(mjdmax) then mjdmax=99999L

if keyword_set(nomjd) then goto,dofield
; don't collide with another process....
while file_test(spectrodir+'index.html.lock') do $
   apwait,spectrodir+'index.html.lock',10

openw,lock,/get_lun,spectrodir+'index.html.lock'
free_lun,lock

; find all .log.html files, get all MJDs with data
logs=[file_search(getenv('APOGEE_DATA')+'/*/*.log.html'),$
      file_search(getenv('APOGEE_DATA_2S')+'/*/*.log.html'),$
      file_search(getenv('APOGEE_DATA_1M')+'/*/*.log.html')]
mjd=lonarr(n_elements(logs))
for i=0,n_elements(logs)-1 do begin
  dir=file_basename(logs[i])
  comp=strsplit(dir,'.',/extract)
  name=comp[0]
  m=0L
  reads,name,m
  mjd[i]=m
  if keyword_set(all) and m ge mjdmin  and m le mjdmax then mkhtml,m,vers=apred_vers,nocheck=nocheck
endfor
s=sort(mjd)
mjd=mjd[s]
logs=logs[s]

; MJD web page
if ~keyword_set(mjdfile) then mjdfile='mjd'
openw,html,/get_lun,spectrodir+'/'+mjdfile+'.html.tmp'
printf,html,'<HTML><BODY>'
printf,html,'<HEAD><script type=text/javascript src=html/sorttable.js></script></head>'
printf,html,'<p><A HREF=fields.html> FIELDS view </a><p>'
;printf,html,'<a href=http://sdss3.apo.nmsu.edu/sdssProcedures/> 2.5m Observers Procedures</a><a href=https://sdss3.org/mailman/private/sdss3-25mlog/all/date.html> (25mlog email archive) </A><br>'
;printf,html,'<a href=http://apogee-schedule.apo.nmsu.edu/AutoSched/cgi/summary.pl> Autoscheduler completion summary </A>'
printf,html, 'Blue: APO 2.5m, Green: LCO 2.5m, Red: APO 1m'
printf,html,'<br>Click on column headings to sort'
; create web page with entry for each MJD
printf,html,'<TABLE BORDER=2 CLASS=sortable>'
printf,html,'<TR><TD>Logs (data)<TD>Exposures<TD>Night QA<TD>Plate QA plots (indiv) <TD>Plate QA (summed plates) <TD>Indiv stars<TD>Dome flats'
; following perhaps more efficient, combined with lines below?
;print,'loading apo1m'
;plates_apo1m=file_search(spectrodir+'/visit/apo1m/*/*/html/*sum.html') 
;print,'loading apo25m'
;plates_apo25m=file_search(spectrodir+'/visit/apo25m/*/*/*/html/*sum.html') 
;print,'loading lco25m'
;plates_loc25m=file_search(spectrodir+'/visit/lco25m/*/*/*/html/*sum.html') 
for i=n_elements(mjd)-1,0,-1 do begin
 if mjd[i] ge mjdmin and mjd[i] le mjdmax then begin
  cmjd=string(format='(i5.5)',mjd[i])
  if strpos(logs[i],'data2s') ge 0 then apsetver,telescope='lco25m' $
  else if strpos(logs[i],'data1m') ge 0 then apsetver,telescope='apo1m' $
  else apsetver,telescope='apo25m'
  dirs=getdir(apodir,caldir,spectrodir,vers)
  if dirs.telescope eq 'apo25m' then color='b3b3ff' else if dirs.telescope eq 'apo1m' then color='ffb3b3' else color='b3ffb3'
  print,mjd[i],dirs.instrument,spectrodir,dirs.telescope
  printf,html,'<TR bgcolor='+color+'><TD><A HREF=../../'+file_basename(dirs.datadir)+'/'+cmjd+'/'+cmjd+'.log.html>',cmjd,'</A>'
  printf,html,'<A HREF=../../'+file_basename(dirs.datadir)+'/'+cmjd+'/>(raw)</A>'
  printf,html,'<TD><center><A HREF=exposures/'+dirs.instrument+'/'+cmjd+'/html/'+cmjd+'exp.html>',cmjd,'</A></center>'
  if file_test(spectrodir+'/exposures/'+dirs.instrument+'/'+cmjd+'/html/'+cmjd+'.html') then $
  printf,html,'<TD><center><A HREF=exposures/'+dirs.instrument+'/'+cmjd+'/html/'+cmjd+'.html> '+cmjd+' QA </a></center>' $
  else $
  printf,html,'<TD><center><FONT COLOR=red> '+cmjd+' QA </font></center>' 
  ; find plates reduced for this night
  if dirs.telescope eq 'apo1m' then $
  plates=file_search(spectrodir+'/visit/'+dirs.telescope+'/*/'+cmjd+'/html/'+cmjd+'*sum.html') else $
  plates=file_search(spectrodir+'/visit/'+dirs.telescope+'/*/*/'+cmjd+'/html/[1-9]*-'+cmjd+'.html')
;  if dirs.telescope eq 'apo1m' then allplates=plates_apo1m
;  if dirs.telescope eq 'apo25m' then allplates=plates_apo25m
;  if dirs.telescope eq 'lco25m' then allplates=plates_lco25m
;  j=where(strpos(allplates,'/'+cmjd+'/') ge 0)
;  plates=allplates[j]
  printf,html,'<TD>'
  for j=0,n_elements(plates)-1 do begin
   if plates[j] ne '' then begin
    plate=file_basename(plates[j])
    comp=strsplit(plate,'-',/extract)
    name=strsplit(comp[0],'.',/extract)
    field=apogee_field(0,fix(name[0]))
    p=0
    reads,name[0],p
    if dirs.telescope eq 'apo1m' then begin
      v=strpos(plates[j],'/visit/')
      tmp=strsplit(strmid(plates[j],v+1),'/',/extract)
      prog=tmp[2]
      v1=strpos(tmp[5],'-')
      v2=strpos(tmp[5],'sum')
      obj=strmid(tmp[5],v1+1,v2-v1-1)
      printf,html,'<A HREF='+strmid(plates[j],v+1)+'>',prog+'/'+obj,'</A>'
    endif else printf,html,'<A HREF=visit/'+dirs.telescope+'/'+field+'/'+name+'/'+cmjd+'/html/'+dirs.prefix+'QA-'+plate+'>',p,'</A>'
   endif
  endfor 

  ;find combined files for this night
  plates=file_search(spectrodir+'visit/'+dirs.telescope+'/*/*/'+cmjd+'/html/sum*-'+cmjd+'.html')
  printf,html,'<TD>'
  for j=0,n_elements(plates)-1 do begin
   if plates[j] ne '' then begin
    plate=file_basename(plates[j])
    comp=strsplit(plate,'-',/extract)
    name=strsplit(comp[0],'.',/extract)
    p=0
    cplate=strmid(name[0],3)
    reads,cplate,p
    printf,html,'<A HREF=visit/'+dirs.telescope+'/'+cplate+'/'+cmjd+'/html/'+plate+'>',p,'</A>'
   endif
  endfor 

  ; find single stars observed for this night
  stars=file_search(spectrodir+'visit/'+dirs.telescope+'/*/'+cmjd+'/html/[0-9]*-'+cmjd+'fiber.html')
  printf,html,'<TD>'
  for j=0,n_elements(stars)-1 do begin
   if stars[j] ne '' then begin
    star=file_basename(stars[j])
    comp=strsplit(star,'-',/extract)
    name=strsplit(comp[0],'.',/extract)
    p=0
    reads,name[0],p
    printf,html,'<A HREF=visit/'+dirs.telescope+'/'+name+'/'+cmjd+'/html/'+star+'>',p,'</A>'
   endif
  endfor 

  ; find dome flats observed for this night
  flats=file_search(spectrodir+'visit/'+dirs.telescope+'/*/'+cmjd+'/html/[0-9]?-'+cmjd+'flat.html')
  printf,html,'<TD>'
  for j=0,n_elements(flats)-1 do begin
   if flats[j] ne '' then begin
    flat=file_basename(flats[j])
    comp=strsplit(flat,'-',/extract)
    name=strsplit(comp[0],'.',/extract)
    p=0
    reads,name[0],p
    printf,html,'<A HREF=visit/'+dirs.telescope+'/'+name+'/'+cmjd+'/html/'+flat+'>',p,'</A>'
   endif
  endfor 
 endif
endfor
printf,html,'</table>'

; summary calibration data
caldir='cal/'
printf,html,'<P> calibration data '
printf,html,'<UL>'
printf,html,'<LI> <A HREF='+caldir+'/darkcorr/html/darks.html> Darks </A>'
printf,html,'<LI> <A HREF='+caldir+'/flatcorr/html/flats.html> Flats </A>'
printf,html,'<LI> <A HREF='+caldir+'/flux/html/flux.html> Fiber fluxes from petal flats </A>'
printf,html,'<LI> <A HREF='+caldir+'/trace/html/trace.html> Traces </A>'
printf,html,'<LI> <A HREF='+caldir+'/detector/html/rn.html> Readout noise </A>'
printf,html,'<LI> <A HREF='+caldir+'/detector/html/gain.html> Gain </A>'
printf,html,'<LI> <A HREF='+caldir+'/wave/html/wave.html> Wave cals </A>'
printf,html,'</UL>'

printf,html,'</body></html>'
free_lun,html
file_move,spectrodir+'/'+mjdfile+'.html.tmp',spectrodir+'/'+mjdfile+'.html',/overwrite

dofield:

if ~keyword_set(fieldfile) then fieldfile='fields'
plates=file_search(spectrodir+'/visit/*/*/*/*/'+'*PlateSum*.fits')

; fields view
openw,html,/get_lun,spectrodir+'/'+fieldfile+'.html.tmp'
printf,html,'<HTML><BODY>'
printf,html,'<HEAD><script type=text/javascript src=html/sorttable.js></script></head>'
printf,html,'<A HREF=mjd.html> MJD view</A>'
; should really get this next stuff direct from database!
plans= yanny_readone(getenv('PLATELIST_DIR')+'/platePlans.par')
plate_id=plans.plateid
location_id=plans.locationid
center_ra=plans.racen
center_dec=plans.deccen
name=strarr(n_elements(plans))
for i=0,n_elements(plans)-1 do name[i]=apogee_field(location_id[i],plate_id[i])
;readcol,spectrodir+'/plate.dat',delim='|',name,plate_id,location_id,center_ra,center_dec,format='(a,i,i,f,f)' 
aitoff,-1*(center_ra-180),center_dec,x,y
set_plot,'ps'
cleanplot,/silent
device,file='sky.eps',/encap,/color,xsize=11,ysize=8,/inches
smcolor
A = FINDGEN(17) * (!PI*2/16.)
usersym,cos(A),sin(A),/fill
plot,[0,1],[0,1],xstyle=4,ystyle=4,xrange=[-180,180],yrange=[-90,90]
aitoff_grid
iplate=intarr(n_elements(plates))
for i=0,n_elements(plates)-1 do begin
print,plates[i]
 platetab=mrdfits(plates[i],3,status=status)
 if status eq 0 then begin
  spl=strsplit(plates[i],'/',/extract)
  reads,spl[-3],plate
  mjd=0L
  reads,spl[-2],mjd
  iplate[i]=where(plate_id eq plate)
  color=1
  if plans[iplate[i]].survey eq 'apogee' then color=2  ; red
  if plans[iplate[i]].survey eq 'apogee2' then color=3 ; green
  if plans[iplate[i]].survey eq 'apogee2-manga' then color=3  ; green
  if plans[iplate[i]].survey eq 'manga-apogee2' then color=6  ; cyan
  if plans[iplate[i]].survey eq 'apogee2s' then color=5 ; magenta
  if mjd lt 55799L then color=6
  if iplate[i] ge 0 then oplot,[x[iplate[i]]],[y[iplate[i]]],psym=8,color=color
 endif
endfor
device,/close
ps2gif,'sky.eps',chmod='664'o,/delete,/eps
; galactic
glactc,center_ra,center_dec,2000.,l,b,1,/degree
jj=where(l gt 180)
l[jj]-=360
aitoff,-l,b,x,y
cleanplot,/silent
device,file='galactic.eps',/encap,/color,xsize=11,ysize=8,/inches
plot,[0,1],[0,1],xstyle=4,ystyle=4,xrange=[-180,180],yrange=[-90,90]
aitoff_grid
iplate=intarr(n_elements(plates))
for i=0,n_elements(plates)-1 do begin
 platetab=mrdfits(plates[i],3,status=status)
 platetab=mrdfits(plates[i],3,status=status)
 if status eq 0 then begin
  spl=strsplit(plates[i],'/',/extract)
  reads,spl[-3],plate
  iplate[i]=where(plate_id eq plate)
  color=1
  if plans[iplate[i]].survey eq 'apogee' then color=2
  if plans[iplate[i]].survey eq 'apogee2' then color=3
  if plans[iplate[i]].survey eq 'apogee2-manga' then color=3
  if plans[iplate[i]].survey eq 'manga-apogee2' then color=6
  if plans[iplate[i]].survey eq 'apogee2s' then color=5
  if iplate[i] ge 0 then oplot,[x[iplate[i]]],[y[iplate[i]]],psym=8,color=color
print,plates[i],x[iplate[i]],y[iplate[i]],color
 endif
endfor
device,/close
ps2gif,'galactic.eps',chmod='664'o,/delete,/eps

printf,html,'<p>APOGEE sky coverage: red=APOGEE1 (yellow: commissioning), green=APOGEE2, magenta=APOGEE2S, cyan=MaNGA-APOGEE2<p>'
printf,html,'<img src=sky.gif width=45%>'
printf,html,'<img src=galactic.gif width=45%>'
printf, html,'<p>Summary files:'

if ~keyword_set(suffix) then suffix='-'+apred_vers+'-'+aspcap_vers+'.fits'
printf,html,'<a href=../../aspcap/'+apred_vers+'/'+aspcap_vers+'/allStar'+suffix+'> allStar'+suffix+' file </a> '
printf,html,' and <a href=../../aspcap/'+apred_vers+'/'+aspcap_vers+'/allVisit'+suffix+'> allVisit'+suffix+' file </a>'

printf,html,'<br>Links on field name are to combined spectra plots and info'
printf,html,'<br>Links on plate name are to visit spectra plots and info'
printf,html,'<br>Links on MJD are to QA and summary plots for the visit'
printf,html,'<br>Click on column headings to sort'

printf,html,'<TABLE BORDER=2 CLASS=sortable>'
printf,html,'<TR><TD>FIELD<TD>Program<TD>ASPCAP<br>'+apred_vers+'/'+aspcap_vers+'<TD>PLATE<TD>MJD<TD>LOCATION<TD>RA<TD>DEC<TD>S/N(red)<TD>S/N(green)<TD>S/N(blue)'
;printf,html,'<TD>RV distrib'
oldname=''
locsort=sort(location_id[iplate])
locra=sort(center_ra[iplate])
;openw,2,'bad.dat'
; loop through plates, sorted by RA
for j=0,n_elements(plates)-1 do begin
 i=locsort[j]
 i=locra[j]
print,'doing plate: ',plates[i]
 platetab=mrdfits(plates[i],3,status=status)
 if iplate[i] ge 0 and status eq 0 then begin
  plate=0 & mjd=0L
  cloc=strtrim(string(format='(i)',location_id[iplate[i]]),2)
  cloc=name[iplate[i]]
  spl=strsplit(plates[i],'/',/extract)
  reads,spl[-3],plate
  reads,spl[-2],mjd
  if strpos(plans[iplate[i]].survey,'apogee2s') ge 0 then tele='lco25m' else tele='apo25m'
  ; is this a new location
  if name[iplate[i]] eq oldname then  begin
    outname='&nbsp;&nbsp;'+name[iplate[i]]
    outresults='' 
  endif else begin
    outname=name[iplate[i]]
    outresults='<A HREF=../../aspcap/'+apred_vers+'/'+aspcap_vers+'/'+tele+'/'+cloc+'>aspcapStar</A><br>'
    stardir=apstar_vers+'/'+strtrim(string(format='(i)',location_id[iplate[i]]),2)+'/'
    rvfile=stardir+'apFieldVisits-'+cloc+'.fits'
    print,rvfile
    rv=mrdfits(rvfile,1,status=status) 
    rvfile=stardir+'apFieldVisitsC-'+cloc+'.fits'
    print,rvfile
    rvc=mrdfits(rvfile,1,status=statusc) 
    if statusc eq 0 then begin
      if status eq 0 then rv=[rv,rvc] else rv=rvc
      if status ne 0 then status=statusc
    endif
  endelse

  if mjd lt 55799 then color='#b3ffff' else color='#ffb3b3'
  if plans[iplate[i]].survey eq 'apogee2' then color='#b3ffb3' 
  if plans[iplate[i]].survey eq 'apogee2-manga' then color='#b3ffb3' 
  if plans[iplate[i]].survey eq 'apogee2s' then color='#ffb3ff' 
  if plans[iplate[i]].survey eq 'manga-apogee2' then color='#b3ecff'
  apsetver,telescope=tele
  dirs=getdir()
  cplate=strtrim(string(format='(i6.4)',plate),2)
  cmjd=string(format='(i5.5)',mjd)
  platedir=spectrodir+'/plates/'+cplate+'/'+cmjd+'/'
  ; make RV histogram if we don't have one
  if not file_test(platedir+'/plots/rv.gif') then begin
;  if not tag_exist(rv,'objid') then begin
;    printf,2,'rm plates/',cplate,'/',cmjd,'/apPlate-?-',cplate,'-',cmjd,'.fits'  
;    printf,2,'rm log/apPlan-',cplate,'-',cmjd,'.par.done'  
;  endif
   platedir='/plates/'+cplate+'/'+cmjd+'/'
   if status eq 0 then begin
    if size(rv,/type) eq 8 then begin
     radvel=rv.vrel
     good=where(radvel lt 1e4 and rv.plate eq plate and rv.mjd eq mjd,ngood)
     if ngood gt 1 then begin
       set_plot,'PS'
       cleanplot,/silent
       device,file=platedir+'/plots/rv.eps',/encap,ysize=8
       plothist,radvel[good],bin=5,thick=2,xrange=[-600,600]
       device,/close
;       ps2gif,platedir+'/plots/rv.eps',chmod='664'o,/delete,/eps
     endif
    endif ; else printf,2,'apPlan-'+cplate+'-'+cmjd+'.par.done'
   endif
  endif

  platedir='plates/'+cplate+'/'+cmjd+'/'
  cid=strtrim(string(format='(i)',location_id[iplate[i]]),2)
  if name[iplate[i]] eq oldname then  $
  printf,html,'<TR bgcolor='+color+'><TD>',outname $
  else printf,html,'<TR bgcolor='+color+'><TD><A HREF='+apstar_vers+'/'+tele+'/'+outname+'/html/'+outname+'.html>',outname,'</A>'
  printf,html,'<TD>',plans[iplate[i]].programname
  printf,html,'<TD>',outresults
  ;printf,html,'<A HREF='+aspcap_vers+platedir+'>apVisit</A>'
  printf,html,'<TD><A href='+platedir+'/html/'+dirs.prefix+'QA-'+cplate+'-'+cmjd+'.html>',plate,'</a>'
 ; printf,html,'<TD><A href='+platedir+'/html/'+dirs.prefix+'QA-'+cplate+'-'+cmjd+'.html>',mjd,'</a>'
  printf,html,'<TD><center><A HREF=exposures/'+dirs.instrument+'/'+cmjd+'/html/'+cmjd+'.html> '+cmjd+' </a></center>' 
  printf,html,'<TD>',location_id[iplate[i]],'<TD>',center_ra[iplate[i]],'<TD>',center_dec[iplate[i]],$
    '<TD>',platetab.sn[0],platetab.altsn[0],'<TD>',platetab.sn[1],platetab.altsn[1],'<TD>',platetab.sn[2],platetab.altsn[2],$
    format='(a,i,a,f12.6,a,f12.6,a,f7.1,f7.1,a,f7.1,f7.1,a,f7.1,f7.1)'
;  printf,html,'<TD><A HREF='+platedir+'/plots/rv.gif><IMG src='+platedir+'/plots/rv.gif height=50></a>'
  oldname=name[iplate[i]]
 endif
endfor
printf,html,'</TABLE>'
printf,html,'</BODY></HTML>'
free_lun,html
file_move,spectrodir+'/'+fieldfile+'.html.tmp',spectrodir+'/'+fieldfile+'.html',/overwrite
;close,2

file_delete,spectrodir+'index.html.lock',/allow_nonexistent

end
