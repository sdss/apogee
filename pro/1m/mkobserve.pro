pro mkobserve,file,stars=stars,rdfile=rdfile,split=split,hours=hours,MJDlim=MJDlim,repall=repall,nobs=nobs
  
  if not keyword_set(split) then split=1

  ;1m commands
  ;If bright, then tell the 1m to check with the observer
  ;for exact target positioning, ie 101
  ;airmass limits, hour angle limits, 
  IF not keyword_set(MJDlim) THEN MJDlim=99999.00
  IF not keyword_set(nobs) THEN nobs=1
  p='2000.0  1.03 1.80 2007  0  0 12 31  24.0  0.00  24.00 -12.00  12.00   1.00  360.00    '+strtrim(nobs,2)+'   '+strtrim(MJDlim,2)+' 100   0   3  0  0'
  pb='2000.0  1.03 1.80 2007  0  0 12 31  24.0  0.00  24.00 -12.00  12.00   1.00  360.00    '+strtrim(nobs,2)+'   '+strtrim(MJDlim,2)+' 101   0   3  0  0'
  

  ;Read in program catalog file
  str=mrdfits(getenv('APOGEEREDUCEPLAN_DIR')+'/data/1m/'+file+'.fits',1)
  print,file
  help,str

  ;Check if select stars were provided
  ;either from an array or read in from
  ;a file. If not just use the whole catalog
  ;sort stars by name to check for duplicates
  if keyword_set(rdfile) then readcol,rdfile,stars,F='A',/SILENT
  if n_elements(stars) le 0 then stars=str[sort(str.name)].name $
  else stars=stars[sort(stars)]

  ;Find the catalog index for stars
  ;If select star is not in the program, do not use it
  ;Check for duplicate stars
  catind=fltarr(n_elements(stars))
  FOR i=0,n_elements(catind)-1 DO BEGIN
     catind[i]=where(strtrim(str.name,2) eq strtrim(stars[i],2))
     if catind[i] lt 0 then print,stars[i]+' not in program '+file
     if i gt 0 AND strtrim(stars[i],2) eq strtrim(stars[i-1],2) then $
        catind[i]=-1
  ENDFOR
  stars=stars[where(catind ge 0)]
  catind=catind[where(catind ge 0)]
  isort=sort(str[catind].ra)  ; sort by ra

  ;Check for stars that have been successfully observed
  if strpos(file,'calibration') ne -1 then fdone=getenv('APOGEEREDUCEPLAN_DIR')+'/data/1m/calibration-done.fits' else $
        fdone=getenv('APOGEEREDUCEPLAN_DIR')+'/data/1m/'+file+'-done.fits'
  if file_test(fdone) eq 1 then done=mrdfits(fdone,1) else done=!NULL
  ;option for repeat observations of all targets, even done ones
  if keyword_set(repall) then done=!NULL
  help,done

  nset=1
  hbright=2.  ;Bright stars require some manual observations
  hfaint=12.  ;Stars fainter than H of 10 should not get observed. Changed to 12 to allow low SN Be stars
  ;Open output files and print 2 blank lines
  openw,bright,'apogee_'+file+'_bright.inp',/get_lun
  printf,bright,' '
  printf,bright,' '
  ;split stars into multiple observing files if there are too many
  for istart=0,split-1 do begin
     if split eq 1 then  outfile='apogee_'+file else $
        outfile='apogee_'+file+string(format='(i2.2)',istart)
     openw,out,outfile+'.inp',/get_lun
     printf,out,' '
     printf,out,' '
     for j=istart,n_elements(stars)-1,split do begin
        i=catind[isort[j]]   ;index by sorted ra
        if done eq !NULL then idone=-1 $    ;check for done stars
        else idone=where(STRTRIM(done.name,2) eq STRTRIM(str[i].name,2))
        if idone ge 0 then print,str[i].name,'done'
        if idone eq -1 then begin
           ;Assign needed information from program catalog
           name=STRJOIN(STRSPLIT(str[i].name, /EXTRACT), '_')
           ;RA should be in degrees, but if not, then set keyword hours
           ra=str[i].ra/15.
           if keyword_set(hours) then ra*=15.
           dec=str[i].dec
           h=str[i].h
           rah=fix(ra)
           ram=fix((ra-rah)*60)
           ras=(ra-rah-ram/60.)*3600
           adec=abs(dec)
           decd=fix(adec)
           decm=fix((adec-decd)*60)
           decs=(adec-decd-decm/60.)*3600
           if dec/adec ge 0 then sign='+' else sign='-'
           ;Stars fainter than H of hfaint will not get observed.
           if h gt hbright and h lt hfaint then begin
              printf,out,name,rah,ram,ras,sign,decd,decm,decs,p, $
                     format='(a,3x,i2.2,":",i2.2,":",f05.2,3x,a,i2.2,":",i2.2,":",f04.1,3x,a)'
           endif else if h le hbright then $
              printf,bright,name,rah,ram,ras,sign,decd,decm,decs,pb, $
                     format='(a,3x,i2.2,":",i2.2,":",f05.2,3x,a,i2.2,":",i2.2,":",f04.1,3x,a)' 
           if h ge hfaint then print,name,' is too faint to observe with the 1m+APOGEE'
           ;SET NUMBER OF READS
           ; 10th mag needs 3*8*47 reads 
           ; 5th mag needs 3*8*.47 reads=4*3 reads
           ; nread=max([fix(10^(0.4*(h[i]-5))*5),5])
           ; Need to increase exposure times to reach SNR goal
           ; stars below 4.5 mag are reaching SNR 100, so just
           ; need to increase exposure times above 4.5 by ~1.5
           nread=max([fix(10^(0.4*(h-4.5))*7),5])
           nread=min([nread,47])
           if h le 0.8 then nread=3

           ;Print exposure information to observing file
           if h gt hbright and h lt hfaint then begin
              printf,out,nset,1,nread,2,nread,2,nread,1,nread,0,0,0,0,0,0,0,0,0,0,0,0, $
                     format='(i5,20i6)' 
           endif else if h le hbright then $
              printf,bright,nset,1,nread,2,nread,2,nread,1,nread,0,0,0,0,0,0,0,0,0,0,0,0, $
                     format='(i5,20i6)'
           
        endif
     endfor
     ;Close output files
     free_lun,out
  endfor
  free_lun,bright

end


