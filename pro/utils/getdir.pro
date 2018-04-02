;======================================================================
; 
; getdir : returns directories given current environment
;
;  USAGE:  getdir,apogeedir,spectrodir,caldir,vers,libdir
;
;    INPUT:   none
;    OUTPUT:  returns top level directories and version
;

function getdir,apogeedir,caldir,spectrodir,vers,libdir,prefix,v=v,apred_vers=apred_vers,datadir=datadir,onem=onem,apogees=apogees

  ;common apver,ver,oldvers,ddir,telescop,instrume
  common apver,ver,telescop,instrume

  ; override
  if n_elements(ver) gt 0 then apred_vers=ver else begin
    ; strict versioning
    vers=apogeereduce_version()
    apred_vers=vers
  endelse
  vers=apred_vers

  apogeedir=getenv('APOGEE_REDUX')
  speclib=getenv('APOGEE_SPECLIB')
  aspcap=getenv('APOGEE_ASPCAP')
  pipedir = getenv('APOGEEREDUCE_DIR')
  libdir=pipedir+'/lib/'
  if telescop eq 'apo1m' then begin
    datadir=getenv('APOGEE_DATA_1M')+'/'
  endif else if telescop eq 'lco25m' then begin
    datadir=getenv('APOGEE_DATA_2S')+'/' 
  endif else begin
    datadir=getenv('APOGEE_DATA')+'/' 
  endelse
  if instrume eq 'apogee-n' then begin
    prefix='ap'
    mapperdir=getenv('MAPPER_DATA')+'/'
    calfile=libdir+'cal/apogee-n.par'
  endif else begin
    prefix='as'
    mapperdir=getenv('MAPPER_DATA_2S')+'/'
    calfile=libdir+'cal/apogee-s.par'
  endelse
  spectrodir=apogeedir+'/'+apred_vers+'/'
  caldir=spectrodir+'cal/'
  expdir=spectrodir+'/exposures/'+instrume+'/'

  ; if keyword_set(onem) then datadir=getenv('APOGEE_DATA')+'1m/' else $
  ;if keyword_set(apogees) then datadir=getenv('APOGEE_DATA_2S')+'/' else $
  ;datadir=getenv('APOGEE_DATA')+'/'
  ;if n_elements(ddir) gt 0 then datadir=ddir
  ;if keyword_set(apogees) then prefix='as' else prefix='ap'

  ; get version, uses apogeereduce_version script except for trunk/local version,
  ;   where environment variable is used
;  if not keyword_set(v) then vers=apogeereduce_version() else vers=v
;  comp=strsplit(vers,'.',/extract)
;  apred_vers=comp[0]
;  if apred_vers eq 'v0' then apred_vers=vers


; override using apsetver command (don't use in pipeline!)
;  if n_elements(ver) eq 0 then begin
;    print,'You must set the version with apsetver,ver before proceeding!'
;  endif
;  vers=ver
;  apred_vers=vers

;  if oldvers then spectrodir=apogeedir+'spectro/'+apred_vers+'/'

  return,{datadir: datadir, apogeedir: apogeedir, expdir: expdir, caldir: caldir, spectrodir: spectrodir, libdir: libdir, prefix: prefix, calfile: calfile, mapperdir: mapperdir, apred: apred_vers, telescope: telescop, instrument: instrume, redux: apogeedir, speclib: speclib, aspcap: aspcap}

end

function getspectrodir
  getdir,a,c,s,v
  return,s
end

