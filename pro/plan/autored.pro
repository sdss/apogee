;+
; AUTORED
;      automated reduction procedure for APOGEE data
;      creates MJD5auto.pro files
;      runs runapred if plate not already done
;      runs mkhtml,mkhtmlsum,mkmonitor when all plates for a given MJD are done
;-
pro autored,mjds,vers=vers,norun=norun,apogees=apogees,override=override

; setup version and directories
if keyword_set(vers) then apsetver,vers=vers else stop,'need to set vers'

if keyword_set(apogees) then begin
  prefix='as' 
  telescope='lco25m'
  instrument='apogee-s'
  apogee_data=getenv('APOGEE_DATA_2S')
  mapper_data=getenv('MAPPER_DATA_2S')
endif else begin
  prefix='ap'
  telescope='apo25m'
  instrument='apogee-n'
  apogee_data=getenv('APOGEE_DATA')
  mapper_data=getenv('MAPPER_DATA')
endelse
apsetver,vers=vers,telescope=telescope

dirs=getdir(a,c,s)
dir=s+'/autored/'+telescope+'/'
file_mkdir,dir

; loop over each MJD
dosum=0
for i=0,n_elements(mjds)-1 do begin
  mjd=mjds[i]
  cmjd=string(format='(i5.5)',mjd)
  print,mjd,cmjd
  ; if this MJD is done, we are done
  if ~file_test(dir+cmjd+'.done') then begin
    ; if not done has it been started?
    if file_test(dir+cmjd+'.plans') then begin
      ; if it has been started, are all the jobs done? If so, run the MJD summary
      readcol,dir+cmjd+'.plans',format='(a)',plans
      done=1
      for j=0,n_elements(plans)-1 do begin
        junk=strsplit(file_basename(plans[j],'.par'),'-',/ext)
        cplate=junk[1]
        if strpos(plans[j],'Dark') ge 0 or strpos(plans[j],'Cal') ge 0 then $
           qafile=apogee_filename('QAcal',plate=cplate,mjd=cmjd) else $
           qafile=apogee_filename('QA',plate=cplate,mjd=cmjd)
        ;if cplate eq '0000' then qafile=prefix+'QAcal-'+cmjd+'.fits' else $
        ;                         qafile='html/'+prefix+'QA-'+cplate+'-'+junk[2]+'.html'
        ;if ~file_test(s+'/'+telescope+'/'+cplate+'/'+junk[2]+'/'+qafile) and $
        if ~file_test(qafile) and $
           strpos(file_basename(plans[j]),'sky') lt 0 and $
           strpos(file_basename(plans[j]),'dark') lt 0 $
           then done=0
        print,plans[j]
        print,cplate
        print,qafile
        print,file_test(qafile), done
      endfor
      if done then begin
        print,'all reductions complete'
        if ~file_test(s+'/exposures/'+instrument+'/'+cmjd+'/html/'+cmjd+'.html') then begin
          file_mkdir,s+'exposures/'+instrument+'/'+cmjd+'/plan'
          openw,plan,s+'exposures/'+instrument+'/'+cmjd+'/plan/'+prefix+'MJD-'+cmjd+'.par',/get_lun
          printf,plan,'apred_vers  '+vers
          printf,plan,'telescope  '+telescope
          printf,plan,'mjd  '+cmjd
          free_lun,plan
          openw,out,dir+cmjd+'.csh',/get_lun
          printf,out,'#!/bin/csh'
          printf,out,'cd $APOGEE_REDUX/'+vers
          printf,out,'runapred '+vers+' exposures/'+instrument+'/'+cmjd+'/plan/'+prefix+'MJD*.par'
          free_lun,out
          print,'running apMJD...'
          if ~keyword_set(norun) then spawn,'csh '+dir+cmjd+'.csh >&'+dir+cmjd+'.log'
          dosum=1
        endif

        ;file_delete,dir+cmjd+'.csh',/allow
        ;file_delete,dir+cmjd+'.plans',/allow
        openw,out,dir+cmjd+'.done',/get_lun
        free_lun,out
      endif else begin
        print,'reductions still running'
      endelse
    endif else begin
      ; if not started, make the plan files and start the reductions if transfer is complete
      if file_test(apogee_data+'/'+cmjd+'/'+cmjd+'.log') then begin
        readcol,apogee_data+'/'+cmjd+'/'+cmjd+'.log',n,f,skip=3,format='(i,a)'
        complete=1
        for j=0,n_elements(f)-1 do begin
          if ~file_test(apogee_data+'/'+cmjd+'/'+f[j]) then complete=0 else begin
            ; check to see if plugmap is available
            h=headfits(apogee_data+'/'+cmjd+'/'+f[j],ext=1)
            exptype=strtrim(sxpar(h,'exptype'),2)
            if exptype eq 'OBJECT' then begin
              plugid=strtrim(sxpar(h,'name'),2)
              tmp=strsplit(plugid,'-',/extract)
              if ~file_test(mapper_data+'/'+tmp[1]+'/plPlugMapM-'+plugid+'.par') then begin
                print,'No plugmap found for: ', f[j],' ',plugid
                if ~keyword_set(override) then complete=0
              endif
            endif
          endelse
        endfor
        if complete then begin
          print,cmjd+' not done and transferred, creating plan files and running'
          undefine,planfiles
          ; make the automatic reduction file and copy to MJD5.pro if it doesn't exist
          apmkplan,mjd,planfiles=planfiles,apogees=apogees,vers=vers
          if ~file_test(getenv('APOGEEREDUCEPLAN_DIR')+'/pro/'+telescope+'/'+telescope+'_'+cmjd+'.pro') then begin
            file_copy,getenv('APOGEEREDUCEPLAN_DIR')+'/pro/'+telescope+'/'+telescope+'_'+cmjd+'auto.pro',$
                      getenv('APOGEEREDUCEPLAN_DIR')+'/pro/'+telescope+'/'+telescope+'_'+cmjd+'.pro' 
          endif
          openw,out,dir+cmjd+'.plans' ,/get_lun
          for j=0,n_elements(planfiles)-1 do printf,out,planfiles[j]
          free_lun,out
          openw,out,dir+cmjd+'.csh',/get_lun
          printf,out,'#!/bin/csh'
          printf,out,'idl << endidl'
          printf,out," vers='"+vers+"'"
          printf,out,' apsetver,vers=vers'
          printf,out,' @'+telescope+'_'+cmjd+'.pro'
          printf,out,'endidl'
          printf,out,'cd $APOGEE_REDUX/'+vers
          if ~keyword_set(norun) then printf,out,'runapred '+vers+' visit/'+telescope+'/*/*/'+cmjd+'/'+prefix+'Plan*.par'+' cal/'+instrument+'/'+cmjd+'/*Plan*.par'
          free_lun,out
          spawn,'csh '+dir+cmjd+'.csh >&'+dir+cmjd+'.log &'
        endif else print,cmjd+' not done, but still transferring'
      endif else print,'no data log file yet...'
    endelse
  endif else print,' already done'
endfor

help,dosum

if dosum then begin
;  openw,out,dir+'sum.csh',/get_lun
;  printf,out,'#!/bin/csh'
;  printf,out,'cd $APOGEE_REDUX/'+vers
;  printf,out,'idl << endidl'
;  printf,out," mkhtmlsum,/nocheck,apred='"+vers+"',apstar='stars',aspcap='a',results='v'"
;  printf,out,' mkmonitor'
;  printf,out,'endidl'
;  free_lun,out
;  print,'running summary...'
;  if ~keyword_set(norun) then spawn,'csh '+dir+'sum.csh >&'+dir+'sum.log'
mkmonitor
mkhtmlsum,/nocheck,apred=vers,apstar='stars',aspcap='a',results='v'
endif

end
