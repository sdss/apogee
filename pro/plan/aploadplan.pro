pro aploadplan,planfile,planstr,verbose=verbose,silent=silent,struct=struct,$
               expand=expand,plugmapm=plugmapm,stp=stp,error=error,newlog=newlog

;+
;
; APLOADPLAN
;
; This program loads an APOGEE plan file
;
; INPUTS:
;  planfile  The absolute path of the plan file
;  /expand   Expand the calibration ID paths and add in
;              directories.
;
; OUTPUTS:
;  planstr   The plan structure with all the
;              relevant information
;  /verbose  Print a lot of information to the screen
;  /silent   Don't print anything to the screen
;  /stp      Stop at the end of the program
;  =error    The error message if one occurred
;
; USAGE;
;  IDL>aploadplan,planfile,planstr
;
; By D.Nidever  May 2010
;-

apgundef,error,planstr

; Not enough inputs
nplanfile = n_elements(planfile)
if nplanfile eq 0 then begin
  if not keyword_set(silent) then begin
    print,'Syntax - aploadplan,planfile,planstr,verbose=verbose,silent=silent,'
    print,'                    expand=expand,stp=stp,error=error'
  endif
  error = 'Not enough inputs'
  return
endif

; More than one file
if nplanfile gt 1 then begin
  error = 'More than one plan file. Only ONE allowed'
  if not keyword_set(silent) then print,error
  return
endif


; Check that the plan file exists
if file_test(planfile) eq 0 then begin
  error = planfile+' NOT FOUND'
  if not keyword_set(silent) then print,error
  return
endif

;----------
; Strip path from plan file name, and change to that directory
thisplan = file_basename(planfile)
thispath = file_dirname(planfile)
CD, thispath, current=origdir

;----------
; Find the desired structure
if not keyword_set(struct) then struct='APEXP'
YANNY_READ, thisplan, pdata, hdr=hdr, errcode=errcode ;,/anonymous
if errcode ne 0 then begin
  error = 'ERROR LOADING '+planfile
  if not keyword_set(silent) then print,error
  return
endif
if size(pdata,/type) eq 10 then begin
  for i=0, N_elements(pdata)-1 do begin
   if (tag_names(*pdata[i], /structure_name) EQ struct) then $
    allseq = *pdata[i]
  endfor
endif else allseq=0
YANNY_FREE, pdata

if (N_elements(allseq) EQ 0) then begin
  ;splog, 'ABORT: No APEXP structures in plan file ' + thisplan
  error = 'ABORT: No APEXP structure in plan file ' + planfile
  if not keyword_set(silent) then print,error
  CD, origdir
;  return
endif
CD,origdir

; Make the PLANSTR structure
;----------------------------
;planstr = {hdr:hdr,apexp:allseq,plateid:'',mjd:'',planfile:'',logfile:'',plotfile:'',$
;           plugmap:'',detid:'',bpmid:'',darkid:'',flatid:'',psfid:'',waveid:'',lsfid:'',fluxid:'',$
;           extract_type:0,raw_dir:'',red_dir:'',plate_dir:'',star_dir:'',apogeereduceversion:''}
planstr = {hdr:hdr}
if n_elements(allseq) ne 0 then add_tag,planstr,struct,allseq,planstr
tags = TAG_NAMES(planstr)
apgundef,hdr,allseq

; Load all keywords like aploadplugmap.pro

;;----------
;; Find keywords from the header
;keywords = ['plateid','mjd','planfile','logfile','plotfile','plugmap',$
;            'detid','bpmid','darkid','flatid','psfid','waveid','lsfid',$
;            'fluxid','extract_type','raw_dir','red_dir','plate_dir','star_dir','apogeereduceVersion']
;nkeywords = n_elements(keywords)
;for j=0,nkeywords-1 do begin
;  ; NOTE, the keywords are CASE SENSITIVE
;  value = YANNY_PAR(planstr.hdr,keywords[j],count=count)
;  if count eq 0 then begin
;    error = 'KEYWORD '+keywords[j]+' NOT FOUND IN '+planfile
;    if not keyword_set(silent) then print,error
;    return
;  endif
;  tagind = where(tags eq strupcase(keywords[j]))
;  type = size(planstr.(tagind),/type)
;  planstr.(tagind) = fix(value,type=type)
;  ; Print out the keyword/value pair
;  if keyword_set(verbose) then print,keywords[j],' = ',fix(value,type=type)
;end
 

; Get the header keywords
gdline = where(strtrim(strmid(planstr.hdr,0,1),2) ne '#' and $
               strtrim(planstr.hdr,2) ne '',ngdline)
keywords = strarr(ngdline)
keytype = lonarr(ngdline)
for i=0,ngdline-1 do begin
  line = strtrim(planstr.hdr[gdline[i]],2)
  dum = strsplit(line,' ',/extract)
  keywords[i] = idl_validname(dum[0],/convert_all)
endfor

; Check that we can make a valid name
valid_name = IDL_VALIDNAME(keywords)
bdkeywords = where(valid_name eq '',nbdkeywords,ncomp=ngdkeywords)
if nbdkeywords gt 0 then begin
  print,'CANNOT create IDL variable for: ',strjoin(keywords[bdkeywords],', ')
  if ngdkeywords eq 0 then begin
    apgundef,keywords
    nkeywords = 0
  endif else begin
    REMOVE,bdkeywords,keywords
    nkeywords = ngdkeywords
  endelse
endif

;----------
; Find keywords from the header
;keywords = ['reddeningMed','tileId','raCen','decCen','plateId','locationId','temp','haMin','haMax']
;keytype = [4,7,5,5,7,7,4,4,4]
nkeywords = n_elements(keywords)
for j=0,nkeywords-1 do begin

  ; NOTE, the keywords are CASE SENSITIVE
  value = YANNY_PAR(planstr.hdr,keywords[j],count=count)

  ; NO match, check without case sensitivity
  if count eq 0 then begin
    ; must be at the beginning and have a space afterwards
    ind = where(stregex(strtrim(planstr.hdr,2),'^'+keywords[j]+' ',/boolean,/fold_case) eq 1,nind)
    if nind gt 0 then begin
      line = strtrim(planstr.hdr[ind[0]],2)
      len = strlen(keywords[j])
      value = strtrim(strmid(line,len),2)
      count = 1
    end
  endif

  ; No match
  if count eq 0 then begin
    error = 'KEYWORD '+keywords[j]+' NOT FOUND IN '+planfile
    if not keyword_set(silent) then print,error
    return
  endif

  ; Get type
  len = strlen(strtrim(value,2))
  type = 7
  if min(valid_num(value)) eq 1 and max(len) lt 7 then type=4
  if min(valid_num(value)) eq 1 and max(len) ge 7 then type=5
  if min(valid_num(value,/integer)) eq 1 then type=3   ; use min in case value is an array

  ; Use STRING for calIDs
  ;calkeys = ['DETID','BPMID','LITTROWID','PERSISTID','PERSISTMODELID','HISTID','DARKID','FLATID','SPARSEID','FIBERID','PSFID','FLUXID','RESPONSEID','WAVEID','LSFID','PLATEID']
  ;dum = where(calkeys eq strupcase(keywords[j]),ncalid)
  ;if ncalid gt 0 then type=7

  ;type = keytype[j]
  planstr = CREATE_STRUCT(planstr,keywords[j],fix(value,type=type))
  tags = tag_names(planstr)
  ;tagind = where(tags eq strupcase(keywords[j]))
  ;type = size(planstr.(tagind),/type)
  ;planstr.(tagind) = fix(value,type=type)
  ; Print out the keyword/value pair
  if keyword_set(verbose) then print,keywords[j],' = ',fix(value,type=type)
end

if tag_exist(planstr,'apred_vers') then begin
  apsetver,vers=planstr.apred_vers
  if tag_exist(planstr,'telescope') then apsetver,vers=planstr.apred_vers,telescope=planstr.telescope
endif

; Add paths to the calibration IDs and add directories
if keyword_set(expand) then begin

  ; Get data directories/version from environment
  if tag_exist(planstr,'data_dir') then apsetver,datadir=planstr.data_dir
  dirs=getdir(apogee_dir,cal_dir,spectro_dir,apred_vers,datadir=datadir)

  calkeys = ['DETID','BPMID','LITTROWID','PERSISTID','PERSISTMODELID','DARKID','FLATID','PSFID','FLUXID','RESPONSEID','WAVEID','LSFID']
  caltype = ['Detector','BPM','Littrow','Persist','PersistModel','Dark','Flat','PSF','Flux','Response','Wave','LSF']
  caldirs = cal_dir+['detector/','bpm/','littrow/','persist/','persist/','darkcorr/','flatcorr/','psf/','flux/','flux/','wave/','lsf/']
  ncalkeys = n_elements(calkeys)

  tags = tag_names(planstr)
  ; Loop through Calibration IDs
  for i=0,ncalkeys-1 do begin
    tagind = where(tags eq calkeys[i],ntagind)
    ; Tag exists for this calid
    if ntagind gt 0 then begin
      ; Check that it has no path information
      val = planstr.(tagind[0])
      ;if val ne '' and strpos(val,'/') eq -1 then planstr.(tagind[0]) = caldirs[i]+val
      if size(val,/type) ne 7 then begin
        if val ne 0 then planstr.(tagind[0]) = apogee_filename(caltype[i],num=val) else planstr.(tagind[0]) = ''
      endif
    endif
  endfor

  ; Add directores
  if tag_exist(planstr,'telescope') then telescope=planstr.telescope else telescope='apo25m'
  dum = where(tags eq 'PLATE_DIR',nplate_dir)
  if nplate_dir eq 0 then planstr = CREATE_STRUCT(planstr,'PLATE_DIR',spectro_dir+'/visit/'+telescope+'/'+$
                           strtrim(planstr.plateid,2)+'/'+strtrim(long(planstr.mjd),2)+'/')
  dum = where(tags eq 'STAR_DIR',nstar_dir)
  if nstar_dir eq 0 then planstr = CREATE_STRUCT(planstr,'STAR_DIR',spectro_dir+'fields/apo25m/')

  if tag_exist(planstr,'instrument') then instrument=planstr.instrument else instrument='apogee-n'
  dum = where(tags eq 'MJD',nmjd)
  if nmjd gt 0 then begin
    dum = where(tags eq 'RED_DIR',nred_dir)
    ;if nred_dir eq 0 then planstr = CREATE_STRUCT(planstr,'RED_DIR',spectro_dir+'red/'+strtrim(long(planstr.mjd),2)+'/')
    if nred_dir eq 0 then planstr = CREATE_STRUCT(planstr,'RED_DIR',dirs.expdir+'/'+strtrim(long(planstr.mjd),2)+'/')
  endif else begin
    print,'NO MJD.  CANNOT make RED_DIR or RAW_DIR'
    return
  endelse

  dum = where(tags eq 'RAW_DIR',nraw_dir)
  if nraw_dir eq 0 then planstr = CREATE_STRUCT(planstr,'RAW_DIR',dirs.datadir+strtrim(long(planstr.mjd),2)+'/')

  ; Expand plugmap
  tagind = where(tags eq 'PLUGMAP',nplugmap)
  if nplugmap gt 0 then begin
    ; Check that it has no path information
    val = planstr.(tagind[0])
    if strpos(val,'/') eq -1 and strmid(val,2) ne '' then begin
      ; Add directory
      plugfile = datadir+strtrim(long(planstr.mjd),2)+'/'
      ; Add plPlugMap prefix if necessary
      if strmid(val,0,2) ne 'pl' then begin
        if not keyword_set(plugmapm) then plugfile+='plPlugMapA-' else plugfile+='plPlugMapM-'
      endif
      plugfile += val
      ; Add .par ending if necessary
      len = strlen(val)
      if strmid(val,len-4,4) ne '.par' then plugfile+='.par'
      planstr.(tagind[0]) = plugfile
    endif
  endif

endif ; expand path/directories

; create log file if one is specified so we can append to it
if keyword_set(newlog) then begin
  logfile=apogee_filename('Diag',plate=planstr.plateid,mjd=planstr.mjd)
  tags = tag_names(planstr)
  dum = where(tags eq 'LOGFILE',nlog)
  if nlog gt 0 then begin
    openw,log,/get_lun,logfile
    free_lun,log
  endif
endif

; Make sure MJD is there
tags = tag_names(planstr)
dum = where(tags eq 'MJD',nmjd)
if nmjd eq 0 then planstr=CREATE_STRUCT(planstr,'MJD','')

if keyword_set(stp) then stop

end
