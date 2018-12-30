PRO mkwiki,program,apver,skipred=skipred,wrtall=wrtall,redonly=redonly

words = STRSPLIT(apver,'/',/EXTRACT)
apredver = words[0]
IF ~keyword_set(redonly) THEN BEGIN
apstarver = words[1]
aspver = words[2]
aspversion = words[3]
ENDIF

;FIND STARS THAT SUCCESSFULLY REDUCED
sumpath = getenv('APOGEE_REDUX')+'/'+apredver+'/fields/apo1m/'+program+'/'
apSum = file_search(sumpath+'apVisitSum-?????-*.fits')
mjdlist = STRMID(apSum,STRLEN(sumpath)+11,5)
tmp = STRSPLIT(STRMID(apSum,STRLEN(sumpath)+17),'.fits',/EXTRACT,/REGEX)
strlist = REPLICATE('a',n_elements(tmp))
FOR i=0,n_elements(strlist)-1 DO strlist[i] = tmp[i]

;IF program EQ 'calibration' THEN BEGIN
;   input_high = mrdfits(getenv('APOGEEREDUCEPLAN_DIR')+'/data/1m/calibration_high.fits',1)
;   input_low = mrdfits(getenv('APOGEEREDUCEPLAN_DIR')+'/data/1m/calibration_low.fits',1)
;ENDIF
input = mrdfits(getenv('APOGEEREDUCEPLAN_DIR')+'/data/1m/'+program+'.fits',1)

IF ~keyword_set(redonly) THEN BEGIN
;READ ASPCAP FILE
   IF FILE_SEARCH(getenv('APOGEE_REDUX')+'/'+apredver+'/'+apstarver+'/'+aspver+'/'+aspversion+'/'+program+'/aspcapField-'+program+'.fits') $
      NE '' THEN BEGIN
      print,getenv('APOGEE_REDUX')+'/'+apredver+'/'+apstarver+'/'+aspver+'/'+aspversion+'/'+program+'/aspcapField-'+program+'.fits'
      asp = mrdfits(getenv('APOGEE_REDUX')+'/'+apredver+'/'+apstarver+'/'+aspver+'/'+aspversion+'/'+program+'/aspcapField-'+program+'.fits',1)
      ASPCAP = replicate({name:'A',H:-99.0,NVISITS:99,MJD:'50000',SNR:-99.0,TYPE:'NONE',COMMENT:'NONE'},n_elements(asp.APOGEE_ID))
      
      FOR i=0,n_elements(asp.APOGEE_ID)-1 DO BEGIN
;FILL INFO FROM ASPCAP FILE
         ASPCAP[i].name = asp[i].APOGEE_ID
         ASPCAP[i].H = asp[i].H
         ASPCAP[i].NVISITS = asp[i].NVISITS
         ASPCAP[i].SNR = asp[i].SNR
         VISFL = STRSPLIT(asp[i].visitfiles,' ',/EXTRACT)
         MJD = REPLICATE('a',n_elements(VISFL))
         FOR j=0,n_elements(VISFL)-1 DO BEGIN
            tmp = STRSPLIT(VISFL[j],'-',/EXTRACT)
            MJD[j] = tmp[2]
         ENDFOR
         ASPCAP[i].MJD = STRJOIN(MJD,',')
                                ;FIND INFO FROM OBSERVING FILE
         IF program eq 'calibration' THEN BEGIN
         ;   index1 = WHERE(STRTRIM(input_high.name,2) EQ STRTRIM(asp[i].APOGEE_ID,2))
         ;   index2 = WHERE(STRTRIM(input_low.name,2) EQ STRTRIM(asp[i].APOGEE_ID,2))
            index = WHERE(STRTRIM(input.name,2) EQ STRTRIM(asp[i].APOGEE_ID,2))
            IF index[0] NE -1 THEN BEGIN
               ASPCAP[i].type = input[index[0]].comment
               ASPCAP[i].comment = input[index[0]].comment 
            ENDIF
        ;   ENDIF ELSE IF index2[0] NE -1 THEN BEGIN
         ;      ASPCAP[i].comment = input_low[index2[0]].comment 
         ;   ENDIF ELSE ASPCAP[i].comment = 'UNKNOWN'
         ENDIF ELSE ASPCAP[i].type = program
         
      ENDFOR
      
   ENDIF ELSE BEGIN 
      print,'NO ASPCAP FILE'
      ASPCAP = {name:'A',H:-99.0,NVISITS:99,MJD:'50000',SNR:-99.0,TYPE:'NONE',COMMENT:'NONE'}
   ENDELSE
;READ APSTAR FILES
   apField = mrdfits(getenv('APOGEE_REDUX')+'/'+apredver+'/'+apstarver+'/apo1m/'+program+'/apField-'+program+'.fits',1)
   APSTAR = replicate({name:'A',H:-99.0,NVISITS:99,MJD:'50000',SNR:-99.0,TYPE:'NONE',COMMENT:'NONE'},n_elements(apField.APOGEE_ID))
   
   j=0
   FOR i=0,n_elements(apField.APOGEE_ID)-1 DO BEGIN
      IF WHERE(STRTRIM(ASPCAP.NAME,2) EQ STRTRIM(apField[i].APOGEE_ID,2)) LT 0 THEN BEGIN
                                ;FILL INFO FROM APSTAR FILE
         APSTAR[j].name = apField[i].APOGEE_ID
         APSTAR[j].H = apField[i].H
                                ;NOT ALWAYS ACCURATE, REDUCED BUT NO APSTAR?
      ;APSTAR[j].NVISITS = apField[i].NVISITS
         APSTAR[j].SNR = apField[i].SNR
      ;GET MJDs FROM apVisitSum
         mjdind = where(STRTRIM(strlist,2) eq STRTRIM(apField[i].APOGEE_ID,2))
         MJD = mjdlist[mjdind]
         APSTAR[j].NVISITS = n_elements(MJD)
         APSTAR[j].MJD = STRJOIN(MJD,',')
         
   ;FIND INFO FROM OBSERVING FILE
         IF program eq 'calibration' THEN BEGIN
            index = WHERE(STRTRIM(input.name,2) EQ STRTRIM(APSTAR[j].name,2))
            IF index[0] NE -1 THEN BEGIN 
               APSTAR[j].comment = input[index[0]].comment 
                  APSTAR[j].type = input[index[0]].type 
            ENDIF 
         ENDIF ELSE APSTAR[j].type = program
         j+=1
      ENDIF
   ENDFOR
   APSTAR = APSTAR[0:j-1]
ENDIF ELSE BEGIN
   ASPCAP = {name:'A',H:-99.0,NVISITS:99,MJD:'50000',SNR:-99.0,TYPE:'NONE',COMMENT:'NONE'}
   APSTAR = {name:'A',H:-99.0,NVISITS:99,MJD:'50000',SNR:-99.0,TYPE:'NONE',COMMENT:'NONE'}
ENDELSE

;FIND STARS THAT ARE REDUCED BUT NO APSTAR OR ASPCAP
;THESE WILL BE CONSIDERED OBSERVED
IF (~keyword_set(skipred)) THEN BEGIN
OBSERVED = replicate({name:'A',H:-99.0,MJD:'50000',SNR:-99.0,TYPE:'NONE',COMMENT:'NONE'},n_elements(strlist))
j=0
FOR i=0,n_elements(strlist)-1 DO BEGIN
   IF WHERE(STRTRIM(ASPCAP.NAME,2) EQ STRTRIM(strlist[i],2)) LT 0 THEN BEGIN
      IF WHERE(STRTRIM(APSTAR.NAME,2) EQ STRTRIM(strlist[i],2)) LT 0 THEN BEGIN
         OBSERVED[j].name = strlist[i]
         OBSERVED[j].MJD = mjdlist[i]
         
         sum = mrdfits(apSum[i],1)
         OBSERVED[j].SNR = sum.SNR
         OBSERVED[j].H = sum.H
         
         ;FIND INFO FROM OBSERVING FILE
         IF program eq 'calibration' THEN BEGIN
            index = WHERE(STRTRIM(input.name,2) EQ STRTRIM(OBSERVED[j].name,2))
            IF index[0] NE -1 THEN BEGIN 
               OBSERVED[j].comment = input[index[0]].comment
               OBSERVED[j].type = input[index[0]].type 
            ENDIF 
         ENDIF ELSE OBSERVED[j].type = program
         j+=1
      ENDIF
   ENDIF
ENDFOR
OBSERVED = OBSERVED[0:j-1]
ENDIF

;FIND UNOBSERVED STARS

;IF program NE 'calibration' THEN input = mrdfits(getenv('APOGEEREDUCEPLAN_DIR')+'/data/1m/'+program+'.fits',1) $
;ELSE BEGIN
;   input = replicate({name:'A',H:-99.0,V:-99.0,RA:-99.0,DEC:-99.0,COMMENT:'NONE',PRIORITY:99},n_elements(input_high.name)+n_elements(input_low.name))
;   name = [input_high.name,input_low.name]
;   ra = [input_high.ra,input_low.ra]
;   dec = [input_high.dec,input_low.dec]
;   H = [input_high.H,input_low.H]
;   V = [input_high.V,input_low.V]
;   comment = [input_high.comment,input_low.comment]
;   priority = [replicate(1,n_elements(input_high.name)),replicate(2,n_elements(input_low.name))]
;   input.name = name
;   input.ra = ra
;   input.dec = dec
;   input.H = H
;   input.V = V
;   input.comment = comment
;   input.priority = priority
;ENDELSE

UNOBSERVED = replicate({name:'A',H:-99.0,PRIORITY:'NONE',TYPE:'NONE',COMMENT:'NONE'},n_elements(input.name))
j=0
FOR i=0,n_elements(input)-1 DO BEGIN
   strind = where(STRTRIM(strlist,2) eq STRTRIM(input[i].name,2))
   IF strind[0] lt 0 THEN BEGIN
      UNOBSERVED[j].name = input[i].name
      UNOBSERVED[j].H = input[i].H
      IF program eq 'calibration' THEN UNOBSERVED[j].PRIORITY = input[i].priority $
      ELSE UNOBSERVED[j].PRIORITY = 1
      IF program eq 'calibration' THEN UNOBSERVED[j].COMMENT = input[i].comment $
      ELSE UNOBSERVED[j].comment = 'none'
      IF program eq 'calibration' THEN UNOBSERVED[j].TYPE = input[i].type $
      ELSE UNOBSERVED[j].type = program
      j+=1
   ENDIF
ENDFOR
UNOBSERVED = UNOBSERVED[0:j-1]

IF keyword_set(wrtall) THEN BEGIN
   
   openw,file,getenv('APOGEEREDUCEPLAN_DIR')+'/data/1m/'+program+'.observed.txt',/GET_LUN
;WRITE OUT FILES
   printf,file,"        Star           H          NVISITS          MJD          SNR          TYPE          COMMENT      "
IF ~keyword_set(redonly) THEN BEGIN
   for i=0,n_elements(ASPCAP)-1 DO printf,file,ASPCAP[i].name,ASPCAP[i].H,ASPCAP[i].NVISITS,ASPCAP[i].MJD,ASPCAP[i].SNR,ASPCAP[i].TYPE,ASPCAP[i].COMMENT,$
                                          f='(A22,F12.6,I6,A30,F12.6,A20,A80)' 
   for i=0,n_elements(APSTAR)-1 DO printf,file,APSTAR[i].name,APSTAR[i].H,APSTAR[i].NVISITS,APSTAR[i].MJD,APSTAR[i].SNR,APSTAR[i].TYPE,APSTAR[i].COMMENT,$
                                          f='(A22,F12.6,I6,A30,F12.6,A20,A80)'
ENDIF
   IF (~keyword_set(skipred)) THEN BEGIN 
      for i=0,n_elements(OBSERVED)-1 DO printf,file,OBSERVED[i].name,OBSERVED[i].H,' 1 ',OBSERVED[i].MJD,OBSERVED[i].SNR,OBSERVED[i].TYPE,OBSERVED[i].COMMENT,$
                                               f='(A22,F12.6,A5,A30,F12.6,A20,A80)' 
   ENDIF
   
   FREE_LUN,file
   
   openw,file,getenv('APOGEEREDUCEPLAN_DIR')+'/data/1m/'+program+'.unobserved.txt',/GET_LUN
   ;WRITE OUT FILES
   printf,file,"        Star           H          Priority          type          comment     "
   for i=0,n_elements(UNOBSERVED)-1 DO printf,file,UNOBSERVED[i].name,UNOBSERVED[i].H,UNOBSERVED[i].PRIORITY,UNOBSERVED[i].TYPE,UNOBSERVED[i].COMMENT,$
      f='(A22,F12.6,I6,A20,A80)'
   FREE_LUN,file
   
   
ENDIF

openw,wikifile,getenv('APOGEEREDUCEPLAN_DIR')+'/data/1m/'+program+'.observed.wiki',/GET_LUN
printf,wikifile,"STARS WITH ASPCAP PARAMETERS"
printf,wikifile,""
printf,wikifile,"  || '''Star ''' || '''H''' || '''NVISITS''' || '''MJD''' || '''SNR''' || '''type''' || '''comment''' ||"
IF ~keyword_set(redonly) THEN BEGIN
for i=0,n_elements(ASPCAP)-1 DO printf,wikifile,'||',ASPCAP[i].name,'||',ASPCAP[i].H,'||',ASPCAP[i].NVISITS,'||',ASPCAP[i].MJD,'||',ASPCAP[i].SNR,'||',ASPCAP[i].TYPE,'||',$
                                       ASPCAP[i].COMMENT,'||',f='(A4,A22,A4,F12.6,A4,I6,A4,A30,A4,F12.6,A4,A20,A4,A80,A4)' 

printf,wikifile,""
printf,wikifile,"STARS RUN THROUGH APSTAR BUT NO ASPCAP"
printf,wikifile,""
printf,wikifile,"  || '''Star ''' || '''H''' || '''NVISITS''' || '''MJD''' || '''SNR''' || '''type''' || '''comment''' ||"

for i=0,n_elements(APSTAR)-1 DO printf,wikifile,'||',APSTAR[i].name,'||',APSTAR[i].H,'||',APSTAR[i].NVISITS,'||',APSTAR[i].MJD,'||',APSTAR[i].SNR,'||',APSTAR[i].TYPE,'||',$
                                       APSTAR[i].COMMENT,'||',f='(A4,A22,A4,F12.6,A4,I6,A4,A30,A4,F12.6,A4,A20,A4,A80,A4)' 
ENDIF
IF (~keyword_set(skipred)) THEN BEGIN
printf,wikifile,""
printf,wikifile,"STARS REDUCED BUT NO APSTAR OR ASPCAP"
printf,wikifile,""
printf,wikifile,"  || '''Star ''' || '''H''' || '''MJD''' || '''SNR''' || '''type''' || '''comment''' ||"

for i=0,n_elements(OBSERVED)-1 DO printf,wikifile,'||',OBSERVED[i].name,'||',OBSERVED[i].H,'||',OBSERVED[i].MJD,'||',OBSERVED[i].SNR,'||',OBSERVED[i].TYPE,'||',$
                                         OBSERVED[i].COMMENT,'||',f='(A4,A22,A4,F12.6,A4,A30,A4,F12.6,A4,A20,A4,A80,A4)' 
ENDIF

printf,wikifile,""
printf,wikifile,"STARS NOT YET OBSERVED"
printf,wikifile,""
printf,wikifile,"  || '''Star ''' || '''H mag''' || '''Priority''' || '''type''' || '''comment''' ||"

for i=0,n_elements(UNOBSERVED)-1 DO printf,wikifile,'||',UNOBSERVED[i].name,'||',UNOBSERVED[i].H,'||',UNOBSERVED[i].PRIORITY,'||',UNOBSERVED[i].TYPE,'||',$
   UNOBSERVED[i].COMMENT,'||',f='(A4,A22,A4,F12.6,A4,I5,A4,A20,A4,A80,A4)' 

FREE_LUN,wikifile



STOP

END
