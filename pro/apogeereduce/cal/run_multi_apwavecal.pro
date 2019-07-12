pro run_multi_apwavecal

;; Create multi-night wavelength solutions

apsetver,vers=‘current’,telescope=‘apo25m’
dirs = getdir(apodir,caldir,spectrodir,vers)
wave_dir = dirs.caldir+'/wave/' 


;; Get list of wavelength solutions
wfiles = file_search(wave_dir+'apWave-????????.dat',count=nwfiles)
wbase = strmid(file_basename(wfiles,'.dat'),7)
;  remove 0000 ones
bd = where(strmid(wbase,4,4) eq '0000',nbd)
if nbd gt 0 then begin
  remove,bd,wbase,wfiles
  nwfiles = n_elements(wfiles)
endif
cmjd = getcmjd(long(wbase),mjd=mjd)

;; Load the wave.par file
readcol,caldir+'wave.par',caltype,mjd0,mjd1,name,frames,psfid,format='A,A,A,A,A,A',delim=' '
wavestr = replicate({mjd0:'',mjd1:'',name:'',frames:'',psfid:''},n_elements(mjd1))
wavestr.mjd0 = mjd0
wavestr.mjd1 = mjd1
wavestr.name = name
wavestr.frames = frames
wavestr.psfid = psfid

; The three periods of the cryostat being closed
; find the median wavelength array for each
mjd_periods = replicate({mjd0:0L,mjd1:0L,wave:dblarr(3*2048)},3)
mjd_periods[0].mjd0 = 55800
mjd_periods[0].mjd1 = 56900
mjd_periods[1].mjd0 = 56901
mjd_periods[1].mjd1 = 57951 ;   57935
mjd_periods[2].mjd0 = 57981 ;   58045
mjd_periods[2].mjd1 = 70000

; Summer shutdown period
shutmjd = [55800, 56130, 56512, 56876, 57230, 57600, 57966, 58335]

;; Loop over years
undefine,bad
undefine,wlines
for i=0,n_elements(shutmjd)-2 do begin
  ind = where(mjd ge shutmjd[i] and mjd le shutmjd[i+1],nind)
  print,'Year ',strtrim(i+1,2),' - ',strtrim(nind,2),' wavelength solutions'
  ;; Loop over multi-night wavelength solutions
  endflag = 0
  j = 0
  cnt = 0L
  while (endflag eq 0) do begin
    ;; 10 solutions up to 15 days apart
    wind = where(mjd[ind] ge mjd[ind[j]] and mjd[ind] le mjd[ind[j]]+15 and ind le min(ind)+j+9,nwind)
    frames = wbase[ind[wind]]
    print,frames     

    ;; Get the PSFID
    MATCH,wavestr.name,frames,ind1,ind2,/sort,count=nmatch
    if nmatch gt 0 then begin
      psfid = wavestr[ind1[0]].psfid
      print,psfid
      outfile = wave_dir+'apWave-'+strmid(frames[0],0,4)+'0000.dat'
      ;if file_test(outfile) eq 0 then push,bad,outfile
      ;if file_test(outfile) eq 0 then APWAVECAL,frames,psfid=psfid 
      ;APWAVECAL,frames,psfid=psfid
      ;wave       99999 99999 02450050 02450050,02450051  02450049
      cmjd1 = getcmjd(long(frames),mjd=mjd1)
      name = strmid(strtrim(frames[0],2),0,4)+'0000'
      push,wlines,'wave   '+strtrim(min(mjd1),2)+' '+strtrim(max(mjd1),2)+' '+name+' '+strjoin(strtrim(frames,2),',')+' '+strtrim(psfid,2)
    endif else print,'NO PSFID FOUND'

    ;; maybe we should always overlap the last night to make
    ;; sure that all nights will be covered?    

    j = max(wind)+1
    if j ge nind then endflag=1
    cnt++
  endwhile  

endfor

;; Output the wavelength file
outfile = 'mwave.par'
print,'Writing wavelength lines to ',outfile
writeline,outfile,wlines

stop

end
