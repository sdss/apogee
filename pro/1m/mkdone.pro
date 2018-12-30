PRO MKDONE,program,apver

;This procedure will check the reduced stars and write a file with
;the stars that have been successfully observed, or done. It will 
;check the old done file and update this list with newly observed 
;stars.
;This is hardcoded for 1m data.

;Grab path from reduction version
words = STRSPLIT(apver,'/',/EXTRACT)
apredver = words[0]
apstarver = words[1]

;Set SNR limit for 'done'
SNRlimit = 100.

;Read in all stars in programs that have been observed 
;and reduce to apStar status
apField = mrdfits(GETENV('APOGEE_REDUX')+'/'+apredver+'/'+apstarver+'/apo1m/'+program+'/apField-'+program+'.fits',1)

;Check for misobservations
badst=mrdfits(GETENV('APOGEEREDUCEPLAN_DIR')+'/data/1m/wrong_obs.fits',1)

;Read in existing -done file if it exists
IF file_test(GETENV('APOGEEREDUCEPLAN_DIR')+'/data/1m/'+program+'-done.fits') eq 1 THEN $
  olddone=mrdfits(GETENV('APOGEEREDUCEPLAN_DIR')+'/data/1m/'+program+'-done.fits',1) $
  ELSE olddone=!NULL

;Find where reduced stars have SNR above the SNR limit and check
;against misobserved stars
indext = where(apField.SNR gt SNRlimit)
FOR i=0,n_elements(indext)-1 DO $
   IF where(strtrim(badst.name,2) eq strtrim(apField[indext[i]].APOGEE_ID,2)) ge 0 THEN indext[i]=-1
index=indext[where(indext ge 0)]
help,index

;Fill in done structure, including previously done stars 
IF olddone eq !NULL THEN BEGIN 
;create structure for done stars
   done = replicate({name:'A',ra:-99.0,dec:-99.0,H:-99.0,SNR:-99.0},n_elements(index))
   done.name = apField[index].APOGEE_ID
   done.ra = apField[index].ra
   done.dec = apField[index].dec
   done.H = apField[index].H
   done.SNR = apField[index].SNR
ENDIF ELSE BEGIN
    ;Check for duplicates between apField and olddone
    apind = replicate(-1,n_elements(olddone.name))
    FOR i=0,n_elements(apind)-1 DO apind[i]=where(strtrim(apField[index].APOGEE_ID,2) eq strtrim(olddone[i].name,2))
    ndbl = where(apind lt 0)
;create structure for done stars
    done = replicate({name:'A',ra:-99.0,dec:-99.0,H:-99.0,SNR:-99.0},n_elements(index)+n_elements(ndbl))
    done.name = [apField[index].APOGEE_ID,olddone[ndbl].name]
    done.ra = [apField[index].ra,olddone[ndbl].ra]
    done.dec = [apField[index].dec,olddone[ndbl].dec]
    done.H = [apField[index].H,olddone[ndbl].H]
    done.SNR = [apField[index].SNR,olddone[ndbl].SNR]
ENDELSE

;Write out done structure to a fits file
mwrfits,done,'$APOGEEREDUCEPLAN_DIR/data/1m/'+program+'-done.fits',/create


END
