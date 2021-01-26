pro mkplan1m,date,write=write
;
; Makes plan files for 1m observations, in two stages
;  first stage, with /write, write the summary observation information into APOGEEREDUCEPLAN_DIR
;     this is run at NMSU, where there is access to the 1m observations files
;  second stage, without /write, makes the apPlan files with mkplan


; get directory with 1m reference files
if keyword_set(write) then  $
path=getenv('APOGEE_1M')+'/'+date+'/' else $
path=getenv('APOGEEREDUCEPLAN_DIR')+'/data/1m/'+date+'/' 

;Add ability to comment out stars in apogee_program.inp.dat

targdir=getenv('APOGEEREDUCEPLAN_DIR')+'/data/1m/'
;outdir=getenv('APOGEEREDUCE_DIR')+'/data/1m/'+date
if keyword_set(write) then begin
  outdir=write+date
  file_mkdir,outdir
endif

; translate YYMMDD string to MJD via calendar date and date_conv
year=2000+fix(strmid(date,0,2))
imonth=fix(strmid(date,2,2))
months=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
month=months[imonth-1]
day=fix(strmid(date,4,2))
datestring=string(format='(i2.2,"-",a,"-",i4.4)',day,month,year)
print,'datestring: ',datestring
; with updated astrolib date_conv, need time string or else MJD is off by 1!
mjd = long(date_conv(datestring+' 00:00:00','M'))

; get APOGEE dayno
day=getnum(mjd)*10000
print,mjd,day

if mjd lt 58047 then fiber=223 else fiber=221
if mjd ge 58190 then begin
  fixfiberid=1
  fiber=229
endif

;READ LIST OF FLATS, make sure to use LONG
print,path
readcol,path+'apogeeflat.dat',flat,tmp,date,F='L,L,L',/SILENT
if keyword_set(write) then file_copy,path+'apogeeflat.dat',outdir,/overwrite

;READ LISTS OF PROGRAMS
files=file_search(path+'apogee_*.dat')

; loop over programs
FOR ifile=0,n_elements(files)-1 DO BEGIN

  if keyword_set(write) then file_copy,files[ifile],outdir,/overwrite
  program=strmid(file_basename(files[ifile],'.inp.dat'),7)

  if (strmid(program,0,3) eq 'hip' and strlen(program) eq 5) then program='hip'
  if (strpos(program,'bright') ge 0) $
  then program=strmid(program,0,strpos(program,'bright')-1)
  if (strmid(program,0,11) eq 'calibration') then program='calibration'
  if (strpos(program,'fts') ge 0) then program='calibration'
  if (strpos(program,'standards') ge 0) then program='misc'
  if (strpos(program,'Bestars') ge 0) then program='Bestars'

  if file_test(targdir+program+'.fits') $
  then targ=mrdfits(targdir+program+'.fits',1) else $
  if file_test(targdir+'1mstars.fits') $
  then targ=mrdfits(targdir+'1mstars.fits',1) else targ=0

  ; get all stars for this program
  readcol,files[ifile],stars,F='A',/SILENT,COMMENT='#',count=ngood
  if ngood gt 0 then begin

    ;SORT FOR UNIQUE STAR NAMES
    ids = UNIQ(stars)
    star = STRING(FLTARR(n_elements(ids)))
    FOR x=0,n_elements(ids)-1 DO star[x] = stars[ids[x]]
    ;READ EACH STAR FILE, CONSTRUCT LIST OF IMAGES
    FOR i=0,n_elements(star)-1 DO BEGIN
     ; read image numbers, using LONG type
     ;readcol,path+star[i]+'.dat',imstart,imend,date,guide1,guide2, $
     ;  F='L,L,L,L,L',/SILENT
     if file_test(path+star[i]+'.dat') then begin
      readcol,path+star[i]+'.dat',imstart,imend,date, $
        F='L,L,L',/SILENT
      if keyword_set(write) then file_copy,path+star[i]+'.dat',outdir,/overwrite

      ; concatenate image numbers if they are >0
      nline=0
      for j=0,n_elements(imstart)-1 do begin
        if imstart[j] gt 0 and imend[j] gt 0 then begin
          if nline eq 0 then ims=imstart[j]+indgen(imend[j]-imstart[j]+1) else $
             ims=[ims,imstart[j]+indgen(imend[j]-imstart[j]+1)]
          nline+=1
          imflag=0
       endif else begin
          print,'!! Image number is <=0: '+star[i]
          imflag=1
       endelse
       endfor

       ; try to get an H mag
       nind=0
       if size(targ,/type) eq 8 then ind=where(strtrim(targ.name,2) eq strtrim(star[i],2),nind)
       if nind gt 0 then hmag=targ[ind].h else hmag=99.0
       print,program,star[i],hmag
;       if nind eq 0 then stop

       nset=1
       imstart=0
       imend=n_elements(ims)-1
       if (strmid(program,0,11) eq 'rrlyr') then begin
         ;separate exposure blocks to be reduced separately if cadence observations are desired
         del=ims-shift(ims,1)
         imstart=where(del gt 1,nset)    
         nset+=1
         imend=imstart-1
         if nset gt 1 then imstart=[0,imstart] else imstart=0
         if nset gt 1 then imend=[imend,n_elements(ims)-1] else imend=n_elements(ims)-1
       endif 
       for iset=0,nset-1 do begin

         ; construct arrays of fiber and starname
         nims=n_elements(ims[imstart[iset]:imend[iset]])
         fibers=replicate(fiber,nims)
         stars=replicate(star[i],nims)
         hmags=replicate(hmag,nims)

         ; find best flat
         medim=median(ims[imstart[iset]:imend[iset]])
         diff=ABS(medim-DOUBLE(flat))
         a=min(diff,imin)
         refflat=flat[imin] 
	
         ; write the plan file
         plate = 0001
         psfid = refflat
         fluxid = refflat
         print,refflat
         
         if nset eq 1 then suffix='' else suffix='_'+string(format='(i2.2)',iset+1)
         if not keyword_set(write) AND imflag eq 0 then mkplan,ims[imstart[iset]:imend[iset]],program,mjd,psfid,fluxid,$
              vers=vers,stars=fibers,fixfiberid=fixfiberid, $ 
              names=stars,/onem,/ignore,hmags=hmags,suffix=suffix
       endfor
     endif
    ENDFOR
  endif
ENDFOR
    
;FREE_LUN,outpln

END
