PRO counts1m,date=date,MJD=MJD,all=all,Mstart=Mstart,Mstop=Mstop
;CAN ONLY RUN IN redux/apver/apo1m

getdir,a,c,s
print,s

IF keyword_set(date) THEN BEGIN
; translate YYMMDD string to MJD via calendar date and date_conv
   year=2000+fix(strmid(date,0,2))
   imonth=fix(strmid(date,2,2))
   months=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
   month=months[imonth-1]
   day=fix(strmid(date,4,2))
   datestring=string(format='(i2.2,"-",a,"-",i4.4)',day,month,year)
   print,'datestring: ',datestring
   mjd = long(date_conv(datestring,'M'))
ENDIF

If keyword_set(all) THEN dir=file_search(s+'apo1m/*/?????') $
ELSE dir=file_search(s+'apo1m/*/'+MJD)
;dir=file_search(GETENV('APOGEE_REDUX')++'/'+apver+'/apo1m/*/'+MJD)
IF keyword_set(all) THEN files=file_search(s+'apo1m/*/*/apCframe-b-*.fits') $
ELSE files=file_search(s+'apo1m/*/'+MJD+'/apCframe-b-*.fits')
program=strarr(n_elements(dir))
MJDtmp=strarr(n_elements(dir))
FOR i=0,n_elements(dir)-1 DO BEGIN
   tmp=strsplit(dir[i],'/',/extract,/regex)
   ii=where(strtrim(tmp,2) eq 'apo1m')
   program[i]=tmp[ii+1]
   MJDtmp[i]=strmid(dir[i],strlen(dir[i])-5)
ENDFOR
print,program

IF keyword_set(Mstart) and keyword_set(Mstop) THEN $
   iuse=where(FIX(MJDtmp) ge FIX(Mstart) and FIX(MJDtmp) le FIX(Mstop))
IF keyword_set(Mstart) and not keyword_set(Mstop) THEN $
   iuse=where(FIX(MJDtmp) ge FIX(Mstart))
IF not keyword_set(Mstart) and keyword_set(Mstop) THEN $
   iuse=where(FIX(MJDtmp) le FIX(Mstop))

IF keyword_set(Mstart) or keyword_set(Mstop) THEN BEGIN
   dir=dir[iuse]
   
ENDIF

counts=replicate({obj:'a',hmag:-99.99,im:10000000,MJD:99999,program:'a',counts:-99.99,reads:1,sky74:-99.99,sky72:-99.99,sky79:-99.99,sky81:-99.99},n_elements(files))
IF strlen(files[0]) le 0 THEN print,'No apCframe files found on this night'

FOR i=0,n_elements(dir)-1 DO BEGIN
   plan=file_search(dir[i]+'/apPlan-*.par')
   FOR j=0,n_elements(plan)-1 DO BEGIN
      IF strlen(plan[j]) gt 0 THEN BEGIN
         close,1
         openr,1,plan[j]
         nlines=file_lines(plan[j])
         ss=strarr(nlines)
         readf,1,ss

         mag=where(strpos(ss,'hmag') ge 0)
         obj=where(strpos(ss,'object') ge 0)
         words=strsplit(ss[mag],/ext)
         hmag=words[1]
         ;IF n_elements(obj) eq 4 THEN BEGIN
         im=strarr(n_elements(obj))
         FOR k=0,n_elements(obj)-1 DO BEGIN
            words=strsplit(ss[obj[k]],/ext)
            im[k]=words[4]
         ENDFOR
         name=words[6]
         close,1
         
         xx=file_search(dir[i]+'/apCframe-b-'+im[0]+'.fits')
         IF strlen(xx) gt 0 THEN BEGIN
            FOR iim=0,n_elements(im)-1 DO BEGIN
               ;print,dir[i]+'/apCframe-b-'+im[iim]+'.fits'
               MJDtmp=strmid(dir[i],strlen(dir[i])-5)
               iind=where(files eq dir[i]+'/apCframe-b-'+im[iim]+'.fits')
               a=mrdfits(dir[i]+'/apCframe-b-'+im[iim]+'.fits',0,hdr,/SILENT)
               b=mrdfits(dir[i]+'/apCframe-b-'+im[iim]+'.fits',1,/SILENT)
               count=MEAN(b[*,77])
               sky81=MEDIAN(b[*,81])
               sky79=MEDIAN(b[*,79])
               sky72=MEDIAN(b[*,72])
               sky74=MEDIAN(b[*,74])
               rd=SXPAR(hdr,'NFRAMES')
               counts[iind].obj=name
               counts[iind].hmag=hmag
               counts[iind].im=im[iim]
               counts[iind].MJD=MJDtmp
               counts[iind].program=program[i]
               counts[iind].counts=count
               counts[iind].reads=rd
               counts[iind].sky72=sky72
               counts[iind].sky74=sky74
               counts[iind].sky79=sky79
               counts[iind].sky81=sky81
            ENDFOR      
         ENDIF
      ENDIF
   ENDFOR
ENDFOR

counts=counts[where(counts.Hmag gt -50.)]
counts=counts[sort(counts.im)]

getdir,a,c,s
print,'plotting'
print,s

IF keyword_set(all) THEN plotdir = s+'apo1m/' $
ELSE plotdir = s+'red/'+MJD+'/plots/'
IF keyword_set(Mstart) and keyword_set(Mstop) THEN $
   psfile1=plotdir+'zptim-'+Mstart+'-'+Mstop
IF not keyword_set(Mstart) and keyword_set(Mstop) THEN $
   psfile1=plotdir+'zptim-00000-'+Mstop
IF keyword_set(Mstart) and not keyword_set(Mstop) THEN $
   psfile1=plotdir+'zptim-'+Mstart+'-00000'
IF not keyword_set(Mstart) and not keyword_set(Mstop) THEN $
   psfile1=plotdir+'zptim-'+MJD

ps_close
ps_open,psfile1,thick=3,/color,/encap
loadct,39
IF keyword_set(all) THEN BEGIN
   imm = findgen(n_elements(counts.im))
   plotc,imm,counts.hmag+2.5*ALOG10(counts.counts)-2.5*ALOG10(counts.reads-2),counts.MJD,psym=symcat(16),yrange=[5,15],ytitle='Hmag + 2.5LOG(counts) - 2.5LOG(reads-2)',xtitle='Image #',title='MJD',/NAN
ENDIF ELSE $
   plotc,counts.im,counts.hmag+2.5*ALOG10(counts.counts)-2.5*ALOG10(counts.reads-2),counts.hmag,psym=symcat(16),yrange=[5,15],ytitle='Hmag + 2.5LOG(counts) - 2.5LOG(reads-2)',xtitle='Image #',title='H mag',/NAN
   oplot,[0,1E9],[12.,12.],linestyle=1
   oplot,[0,1E9],[12.,12.],linestyle=1
ps_close

SPAWN,'ls '+psfile1+'*'

IF keyword_set(Mstart) and keyword_set(Mstop) THEN $
   psfile2=plotdir+'zpthist-'+Mstart+'-'+Mstop
IF not keyword_set(Mstart) and keyword_set(Mstop) THEN $
   psfile2=plotdir+'zpthist-00000-'+Mstop
IF keyword_set(Mstart) and not keyword_set(Mstop) THEN $
   psfile2=plotdir+'zpthist-'+Mstart+'-00000'
IF not keyword_set(Mstart) and not keyword_set(Mstop) THEN $
   psfile2=plotdir+'zpthist-'+MJD

ps_open,psfile2,thick=3,/color,/encap
loadct,39
plothist,counts.hmag+2.5*ALOG10(counts.counts)-2.5*ALOG10(counts.reads-2),bin=0.2,/NAN,xtitle='Hmag + 2.5LOG(counts) - 2.5LOG(reads-2)',title='Objects'
ps_close


psfile=plotdir+'zptsky72-'+MJD
ps_open,psfile,thick=3,/color,/encap
loadct,39
plotc,counts.im,2.5*ALOG10(counts.sky72)-2.5*ALOG10(counts.reads-2),counts.hmag,psym=symcat(16),ytitle='2.5LOG(counts) - 2.5LOG(reads-2)',xtitle='Image #',title='Sky Fiber 72'
ps_close

psfile=plotdir+'zptsky74-'+MJD
ps_open,psfile,thick=3,/color,/encap
loadct,39
plotc,counts.im,2.5*ALOG10(counts.sky74)-2.5*ALOG10(counts.reads-2),counts.hmag,psym=symcat(16),ytitle='2.5LOG(counts) - 2.5LOG(reads-2)',xtitle='Image #',title='Sky Fiber 74'
ps_close

psfile=plotdir+'zptsky79-'+MJD
ps_open,psfile,thick=3,/color,/encap
loadct,39
plotc,counts.im,2.5*ALOG10(counts.sky79)-2.5*ALOG10(counts.reads-2),counts.hmag,psym=symcat(16),ytitle='2.5LOG(counts) - 2.5LOG(reads-2)',xtitle='Image #',title='Sky Fiber 79'
ps_close

psfile=plotdir+'zptsky81-'+MJD
ps_open,psfile,thick=3,/color,/encap
loadct,39
plotc,counts.im,2.5*ALOG10(counts.sky81)-2.5*ALOG10(counts.reads-2),counts.hmag,psym=symcat(16),ytitle='2.5LOG(counts) - 2.5LOG(reads-2)',xtitle='Image #',title='Sky Fiber 81'
ps_close


pltfiles=file_search(plotdir+'zpt*.eps')
FOR i=0,n_elements(pltfiles)-1 DO BEGIN
   ff=pltfiles[i]
   str=strsplit(ff,'/',/EXTRACT,/REGEX)
   spawn,'convert '+ff+' '+plotdir+strmid(str[n_elements(str)-1],0,strlen(str[n_elements(str)-1])-4)+'.jpg'
ENDFOR
;spawn,'convert '+plotdir+'zptim-'+MJD+'.eps '+plotdir+'zptim-'+MJD+'.jpg'
;spawn,'convert '+plotdir+'zpthist-'+MJD+'.eps '+plotdir+'zpthist-'+MJD+'.jpg'
;spawn,'convert '+plotdir+'zptsky72-'+MJD+'.eps '+plotdir+'zptsky72-'+MJD+'.jpg'
;spawn,'convert '+plotdir+'zptsky74-'+MJD+'.eps '+plotdir+'zptsky74-'+MJD+'.jpg'
;spawn,'convert '+plotdir+'zptsky79-'+MJD+'.eps '+plotdir+'zptsky79-'+MJD+'.jpg'
;spawn,'convert '+plotdir+'zptsky81-'+MJD+'.eps '+plotdir+'zptsky81-'+MJD+'.jpg'

plts=file_search(plotdir+'zpt*.eps')
file_delete,plts,/allow


END
