;+
; NAME:
;   aspcap_allstar
; PURPOSE:
;   Make allStar and allVisit files, roll-ups of APOGEE reduction results
; CALLING SEQUENCE:
;   aspcap_allstar, aspcap_version=, apred_version=, apstar_version=, $
;     results_version= [, outdir=, [/unsafe_synth ]
; INPUTS:
;   apred_version - version of ASPCAP ("rN")
;   apstar_version - version of ASPCAP ("sN")
;   aspcap_version - version of ASPCAP ("aN")
;   results_version - version of ASPCAP ("vQQQ")
; OPTIONAL INPUTS:
;   outdir - output directory, instead of $APOGEE_REDUX/rN
; OPTIONAL KEYWORDS:
;   /unsafe_synth - work-around for missing SYNTHRVEL values
; COMMENTS:
;   Creates the files:
;      $APOGEE_REDUX/rN/allFields-vQQQ.fits
;      $APOGEE_REDUX/rN/allPlates-vQQQ.fits
;      $APOGEE_REDUX/rN/allStar-vQQQ.fits
;      $APOGEE_REDUX/rN/allVisit-vQQQ.fits
;   From the files:
;      $APOGEE_REDUX/rN/sN/LLLL/apFieldVisits-LLLL.fits
;      $APOGEE_REDUX/rN/sN/LLLL/apField-LLLL.fits
;      $APOGEE_REDUX/rN/sN/LLLL/apFieldVisitsC-LLLL.fits
;      $APOGEE_REDUX/rN/sN/LLLL/apFieldC-LLLL.fits
;      $APOGEE_REDUX/rN/sN/aN/l6/LLLL/aspcapField-LLLL.fits
; DEPENDS ON:
;   getdir.pro [in apogeereduce product]
;   apgundef.pro [in apogeereduce product]
;   field.pro [in apogeereduce product]
;   rdlibhead.pro [in apogeereduce product]
; REVISION HISTORY:
;   J. Holtzman, NMSU, May 2012
;   Re-organized, M. Blanton, NYU, February 2013
;-
;
function catalog_info_blank

  cat0=create_struct('reduction_id',' ',$
                     'src_h', ' ', $
                     'wash_m', 0.,$
                     'wash_m_err', 0.,$
                     'wash_t2', 0., $
                     'wash_t2_err', 0., $
                     'ddo51', 0., $
                     'ddo51_err', 0., $
                     'irac_3_6', 0., $
                     'irac_3_6_err', 0., $
                     'irac_4_5', 0., $
                     'irac_4_5_err', 0., $
                     'irac_5_8', 0., $
                     'irac_5_8_err', 0., $
                     'irac_8_0', 0., $
                     'irac_8_0_err', 0., $
                     'wise_4_5', 0., $
                     'wise_4_5_err', 0., $
                     'targ_4_5', 0., $
                     'targ_4_5_err', 0., $
                     'ak_targ', -9999.99, $
                     'ak_targ_method', '', $
                     'ak_wise', -9999.99, $
                     'sfd_ebv', -9999.99, $
                     'wash_ddo51_giant_flag', 0, $
                     'wash_ddo51_star_flag', 0, $
                     'pmra', 0., $
                     'pmdec', 0., $
                     'pm_src', ' ')
  return, cat0
end
; Use this function to replace existing data with catalog data
pro catalog_info_replace,src,dest,missing=missing
  if keyword_set(missing) then begin
    if abs(dest.ra-src.ra) gt 2./3600. or abs(dest.dec-src.dec) gt 2./3600. and src.ra gt 0. then $
      printf,missing,'coords not matching: '+src.apogee_id+ string(format='(4f12.8)', src.ra, src.dec, dest.ra,dest.dec)
    dest.ra = src.ra
    dest.dec = src.dec
  endif
  dest.j = src.j
  dest.j_err = src.j_err
  dest.h = src.h
  dest.h_err = src.h_err
  dest.k = src.k
  dest.k_err = src.k_err
  dest.ak_targ = src.ak_targ
  dest.ak_targ_method = src.ak_targ_method
  dest.ak_wise = src.ak_wise
  dest.sfd_ebv = src.sfd_ebv
end
;
function allstar_blank, apfield0, maxvisit=maxvisit, nparam=nparam, nelem=nelem, nclass=nclass, nwind=nwind

  if(NOT keyword_set(maxvisit)) then $
     maxvisit=50L
  if(NOT keyword_set(nparam)) then $
     nparam=7L
  if(NOT keyword_set(nelem)) then $
     nelem=15L
  if(NOT keyword_set(nwind)) then $
     nwind=1L

  struct_assign, {junk:0}, apfield0
  
  allstar0=create_struct('apstar_id', ' ', $
                         'target_id', ' ', $
                         'aspcap_id', ' ', $
                         apfield0, $
  ;                       'apogee2_target1', 0L, $
  ;                       'apogee2_target2', 0L, $
  ;                       'apogee2_target3', 0L, $
                         'meanfib', -1. , $
                         'sigfib', -1., $
                         'snrev', -1., $
                         'apstar_version', ' ', $
                         'aspcap_version', ' ', $
                         'results_version', ' ', $
                         'extratarg', -1, $
                         'min_h', -9999.99, $
                         'max_h', 9999.99, $
                         'min_jk', -9999.99, $
                         'max_jk', 9999.99, $
                         'param', fltarr(nparam)-9999.99, $
                         'fparam', fltarr(nparam)-999999., $
                         'param_cov', fltarr(nparam, nparam)-9999.99, $
                         'fparam_cov', fltarr(nparam, nparam)-9999.99,$
                         ;'elem', fltarr(nelem)-9999.99, $
                         ;'elem_err', fltarr(nelem)-9999.99, $
                         'teff', -9999.99, $
                         'teff_err', -9999.99, $
                         'logg', -9999.99, $
                         'logg_err', -9999.99, $
                         ;'param_teff', -1., $
                         ;'logvmicro', -9999.99, $
                         ;'vmacro', -9999.99, $
                         ;'lgvsini', -9999.99, $    ; covers both rotation and vmacro
                         'vmicro', -9999.99, $    
                         'vmacro', -9999.99, $   
                         'vsini', -9999.99, $   
                         ;'param_logg', -9999.99, $
                         ;'param_logvmicro', -1., $
                         'm_h', -9999.99, $
                         'm_h_err', -9999.99, $
                         ;'param_c_m', -9999.99, $
                         ;'param_n_m', -9999.99, $
                         'alpha_m', -9999.99,$
                         'alpha_m_err', -9999.99, $
                         ;'param_teff_err', -1., $
                         ;'param_logg_err', -1., $
                         ;'param_logvmicro_err', -1., $
                         ;'param_c_m_err', -1., $
                         ;'param_n_m_err', -1., $
                         'aspcap_chi2', 0., $
                         'aspcap_class', ' ',$
                         'aspcapflag',  long(aspcapflagval('NO_ASPCAP_RESULT')),$
                         'aspcapflags', 'NO_ASPCAP_RESULT',$
                         'paramflag', lonarr(nparam),$
                         'felem', reform(fltarr(nelem,nwind))-9999.99, $
                         'felem_err', reform(fltarr(nelem,nwind))-9999.99, $
                         'X_H', fltarr(nelem)-9999.99, $
                         'X_H_ERR', fltarr(nelem)-9999.99, $
                         'X_M', fltarr(nelem)-9999.99, $
                         'X_M_ERR', fltarr(nelem)-9999.99, $
                         'C_Fe',-9999.99,$
                         'CI_Fe',-9999.99,$
                         'N_Fe',-9999.99,$
                         'O_Fe',-9999.99,$
                         'Na_Fe',-9999.99,$
                         'Mg_Fe',-9999.99,$
                         'Al_Fe',-9999.99,$
                         'Si_Fe',-9999.99,$
                         'P_Fe',-9999.99,$
                         'S_Fe',-9999.99,$
                         'K_Fe',-9999.99,$
                         'Ca_Fe',-9999.99,$
                         'Ti_Fe',-9999.99,$
                         'TiII_Fe',-9999.99,$
                         'V_Fe',-9999.99,$
                         'Cr_Fe',-9999.99,$
                         'Mn_Fe',-9999.99,$
                         'Fe_H',-9999.99,$
                         'Co_Fe',-9999.99,$
                         'Ni_Fe',-9999.99,$
                         'Cu_Fe',-9999.99,$
                         'Ge_Fe',-9999.99,$
                         ;'Ce_Fe',-9999.99,$
                         'Rb_Fe',-9999.99,$
                         'Y_Fe',-9999.99,$
                         'Nd_Fe',-9999.99,$
                         'C_Fe_err',-9999.99,$
                         'CI_Fe_err',-9999.99,$
                         'N_Fe_err',-9999.99,$
                         'O_Fe_err',-9999.99,$
                         'Na_Fe_err',-9999.99,$
                         'Mg_Fe_err',-9999.99,$
                         'Al_Fe_err',-9999.99,$
                         'Si_Fe_err',-9999.99,$
                         'P_Fe_err',-9999.99,$
                         'S_Fe_err',-9999.99,$
                         'K_Fe_err',-9999.99,$
                         'Ca_Fe_err',-9999.99,$
                         'Ti_Fe_err',-9999.99,$
                         'TiII_Fe_err',-9999.99,$
                         'V_Fe_err',-9999.99,$
                         'Cr_Fe_err',-9999.99,$
                         'Mn_Fe_err',-9999.99,$
                         'Fe_H_err',-9999.99,$
                         'Co_Fe_err',-9999.99,$
                         'Ni_Fe_err',-9999.99,$
                         'Cu_Fe_err',-9999.99,$
                         'Ge_Fe_err',-9999.99,$
                         ;'Ce_Fe_err',-9999.99,$
                         'Rb_Fe_err',-9999.99,$
                         'Y_Fe_err',-9999.99,$
                         'Nd_Fe_err',-9999.99,$
                         'C_Fe_flag',0L,$
                         'CI_Fe_flag',0L,$
                         'N_Fe_flag',0L,$
                         'O_Fe_flag',0L,$
                         'Na_Fe_flag',0L,$
                         'Mg_Fe_flag',0L,$
                         'Al_Fe_flag',0L,$
                         'Si_Fe_flag',0L,$
                         'P_Fe_flag',0L,$
                         'S_Fe_flag',0L,$
                         'K_Fe_flag',0L,$
                         'Ca_Fe_flag',0L,$
                         'Ti_Fe_flag',0L,$
                         'TiII_Fe_flag',0L,$
                         'V_Fe_flag',0L,$
                         'Cr_Fe_flag',0L,$
                         'Mn_Fe_flag',0L,$
                         'Fe_H_flag',0L,$
                         'Co_Fe_flag',0L,$
                         'Ni_Fe_flag',0L,$
                         'Cu_Fe_flag',0L,$
                         'Ge_Fe_flag',0L,$
                         ;'Ce_Fe_flag',0L,$
                         'Rb_Fe_flag',0L,$
                         'Y_Fe_flag',0L,$
                         'Nd_Fe_flag',0L,$
                         'elem_chi2', fltarr(nelem), $
                         'elemflag', lonarr(nelem),$
                         catalog_info_blank(), $
                         'all_visits', ' ', $
                         'visits', ' ', $
                         'all_visit_pk', lonarr(maxvisit)-1L, $
                         'visit_pk', lonarr(maxvisit)-1L )
  if n_elements(nclass) gt 0 then begin
    add_tag,allstar0,'fparam_class',fltarr(nparam,nclass),allstar0
    add_tag,allstar0,'chi2_class',fltarr(nclass),allstar0
  endif

  return, allstar0

end
;
function allvisit_blank, apfieldvisit0, unsafe_synth=unsafe_synth

  common com_allvisit_blank, allstar0_test

  struct_assign, {junk:0}, apfieldvisit0
  
  allvisit0=create_struct('visit_id', ' ', $
                         ;'target_id', ' ', $
                         'apred_version', ' ', $
                         apfieldvisit0, $
                      ;   'field', ' ', $
                         'apogee2_target1', 0L, $
                         'apogee2_target2', 0L, $
                         'apogee2_target3', 0L, $
                         'commiss', 0, $
                         'extratarg', -1, $
                         'min_h', -9999.99, $
                         'max_h', 9999.99, $
                         'min_jk', -9999.99, $
                         'max_jk', 9999.99, $
                         catalog_info_blank())
;                         'wash_m', 0., $
;                         'wash_m_err', 0., $
;                         'wash_t2', 0., $
;                         'wash_t2_err', 0., $
;                         'ddo51', 0., $
;                         'ddo51_err', 0., $
;                         'irac_3_6', 0., $
;                         'irac_3_6_err', 0., $
;                         'mag_4_5', 0., $
;                         'mag_4_5_err', 0., $
;                         'irac_5_8', 0., $
;                         'irac_5_8_err', 0., $
;                         'irac_8_0', 0., $
;                         'irac_8_0_err', 0., $ 
;                         'wise_4_5', 0., $
;                         'wise_4_5_err', 0., $
;                         'giant', 0, $
;                         'star', 0,$ 
;                         'pmra', 0., $
;                         'pmdec', 0., $
;                         'pmcat', ' ')

  if(keyword_set(unsafe_synth)) then $
     allvisit0=struct_trimtags(allvisit0, $
                              except=['SYNTHVREL', 'SYNTHVRELERR', $
                                      'SYNTHVHELIO'])

  if(n_tags(allvisit0_test) eq 0) then begin
     allvisit0_test=allvisit0
  endif else begin
     diff= compare_struct(allvisit0_test, allvisit0)
     if(diff.ndiff gt 0) then $
        message, 'Structures differ!'
  endelse

  return, allvisit0
end
;
function allfields_blank

allfields0= {name:' ', $
             location_id:0L, $
             racen:0.D, $
             deccen:0.D, $
             radius:0., $
             shared:0, $
             field_type:0}

return, allfields0
  
end
;
function allplate_blank

allplate0= {plate_visit_id:' ', $
            location_id:0L, $
            plate:0L, $
            mjd:0L, $
            apred_version:' ', $
            name:' ', $
            racen:0.D, $
            deccen:0.D, $
            radius:0., $
            shared:0, $
            field_type:0, $
            survey:' ', $
            programname:' ', $
            platerun:' ', $
            chunk:' ', $
            ha:fltarr(6), $
            designid:0L, $
            nstandard:0L, $
            nscience:0L, $
            nsky:0L, $
            platedesign_version:' '}

return, allplate0

end
;
pro aspcap_allstar, aspcap_version=aspcap_version, $
                    apred_version=apred_version, $
                    apstar_version=apstar_version, $
                    results_version=results_version, $
                    locationdirs=locationdirs, $
                    outdir=outdir, $
                    caldir=caldir, $
                    unsafe_synth=unsafe_synth, override=override, fixplate=fixplate, notelescope=notelescope, $
                    maxlocations=maxlocations, test=test, params=params, nelem=nelem, apogee2=apogee2, nwind=nwind

;; Set up parameters
maxvisit=50 ;; maximum number of visits to track per star
nparam=7L ;; number of parameters
if ~keyword_set(nelem) then nelem=15L ;; number of elements
telescope='apo25m' ;; northern telescope, placeholder!
if ~keyword_set(apstar_version) then apstar_version='stars'
; specify nclass if fparam_class is desired
nclass=23
if ~keyword_set(nwind) then nwind=1

; set up directories
apsetver, vers=apred_version
dirs=getdir(apogee_dir, cal_dir, spectro_dir, apogee_vers, lib_dir)
aspcap_dir=dirs.aspcap+'/'+apred_version+'/'+aspcap_version+'/'
if not keyword_set(outdir) then outdir=aspcap_dir
if n_elements(results_version) eq 0 then results_version=apred_version+'-'+aspcap_version

;; Set up file for tracking missing visits
apgundef,allvisit,allstar
file_mkdir,outdir+'/log'
openw,missing,outdir+'/log/missing-'+results_version+'.txt',/get_lun
openw,altname,outdir+'/log/altname-'+results_version+'.txt',/get_lun
openw,dup,outdir+'/log/dup-'+results_version+'.txt',/get_lun

get_date, date

;; define index order used by structure for named parameters
if n_elements(params) eq 0 then params=aspcap_params(paramtags)
nparam=n_elements(params)
elems=aspcap_elems(elemtags,elemtoh,nelem=nelem)
nelem=n_elements(elems)

;for itelescope = 0,2 do begin
; if itelescope eq 0 then telescope = 'apo25m'
; if itelescope eq 1 then telescope = 'lco25m'
; if itelescope eq 2 then telescope = 'apo1m'

targetdir=getenv('APOGEE_TARGET')
design1=mrdfits(targetdir+'/apogeeDesign.fits',1)
design2=mrdfits(targetdir+'/apogee2Design.fits',1)
min_jk=-9999.99
max_jk=9999.99

;; Get list of locations by looking at directories
if ~keyword_set(locationdirs) then begin
  ;if itelescope eq 0 then $
  ;locationdirs = file_search(visits_dir+'[0-9][0-9][0-9][0-9]', $
  ;                         /test_directory, count=nlocationdirs) $
  ;else if itelescope eq 1 then $
  locationdirs = file_search(apogee_dir+'/'+apred_version+'/'+apstar_version+'/*/*', $
                           /test_directory, count=nlocationdirs) 
  splog,strtrim(nlocationdirs,2),' Location directories found'
endif else locationdirs=apogee_dir+'/'+apred_version+'/'+apstar_version+'/'+locationdirs
nlocationdirs=n_elements(locationdirs)

; maxlocations keyword for testing limited number of fields
if keyword_set(maxlocations) then nlocationdirs=min([maxlocations,nlocationdirs])

;; Get list of plate plans
plans= yanny_readone(getenv('PLATELIST_DIR')+'/platePlans.par')

; Loop through the locationdirs
if file_test('.exclude') then readcol,'.exclude',excludedirs,format='(a)'
For i=0L,nlocationdirs-1L do begin 
   ;; location id 
   idir = locationdirs[i]+'/'
   ilocationid = file_basename(locationdirs[i])
   print,'ilocationid: ',ilocationid
   printf,missing,'ilocationid: ',ilocationid

   ; get telescope name from directory
   junk=strsplit(idir,'/',/extract)
   telescope=junk[-2]

   ;; Get APOGEE directories
   apsetver,telescope=telescope
   visits_dir=apogee_dir+'/'+apred_version+'/'+apstar_version+'/'+telescope+'/'
   stars_dir=apogee_dir+'/'+apred_version+'/'+apstar_version+'/'+telescope+'/'

   ; do we want to skip this one?
   nskip=0
   if n_elements(excludedirs) gt 0 then skip=where(strtrim(excludedirs,2) eq strtrim(ilocationid,2),nskip)
   if nskip gt 0 then goto,nextloc

   ;; Get all of the apFieldVisits files
   rvisitfiles = file_search(idir+'apFieldVisits-*.fits',count=nrvisitfiles)
   cvisitfiles = file_search(idir+'apFieldVisitsC-*.fits',count=ncvisitfiles)
   apgundef, visitfiles
   if(nrvisitfiles gt 0 and ncvisitfiles gt 0) then begin
      visitfiles= [rvisitfiles, cvisitfiles]
      visitcommiss= [bytarr(nrvisitfiles), bytarr(nrvisitfiles)+1]
   endif else begin
      if(nrvisitfiles gt 0) then begin
         visitfiles= rvisitfiles 
         visitcommiss= bytarr(nrvisitfiles)
      endif
      if(ncvisitfiles gt 0) then begin
         visitfiles= cvisitfiles 
         visitcommiss= bytarr(ncvisitfiles)+1
      endif
   endelse 
   nvisitfiles= n_elements(visitfiles)

   ;; Get all of the apField files 
   idir=stars_dir+'/'+ilocationid+'/'
   rstarfiles = file_search(idir+'apField-*.fits',count=nrstarfiles)
   cstarfiles = file_search(idir+'apFieldC-*.fits',count=ncstarfiles)
   apgundef,starfiles
   if(nrstarfiles gt 0 and ncstarfiles gt 0) then begin
      starfiles= [rstarfiles, cstarfiles]
   endif else begin
      if(nrstarfiles gt 0) then $
         starfiles= rstarfiles 
      if(ncstarfiles gt 0) then $
         starfiles= cstarfiles 
   endelse 
   nstarfiles= n_elements(starfiles)

   if(nstarfiles ne nvisitfiles) then begin
      splog, 'Inconsistency in which apField and apFieldVisits files '+ $
             'exist for location= '+ilocationid
stop
      printf,missing, ilocationid, 'nstarfiles: ', nstarfiles, $
             ' nvisitfiles: ',nvisitfiles
      ;continue
   endif
   if(nstarfiles eq 0 and nvisitfiles eq 0)  then begin
      splog, 'No apField or apFieldVisits files '+ $
             'exist for location= '+ilocationid
      continue
   endif

   splog, strtrim(string(i+1),2),'/',strtrim(nlocationdirs,2), $
          ' LocationID=',ilocationid,'  ',strtrim(nvisitfiles,2), $
          ' apFieldVisit files'
   apgundef,allvisitloc
   apgundef,allstarloc
   ;; Loop the apFieldVisits files
   nplates=0
   for j=0L, nvisitfiles-1L do begin

      splog,visitfiles[j]
      apgundef,str,status
      help, visitfiles[j]
      str = MRDFITS(visitfiles[j],1,/silent,status=status)

      ; remove any targetting tags, as we will repopulate from apogeeObject
      if keyword_set(notelescope) then $
      remove_tags,str,['ak_targ','ak_targ_method','ak_wise','sfd_ebv','telescope'],new else $
      remove_tags,str,['ak_targ','ak_targ_method','ak_wise','sfd_ebv'],new 
      str=new

      ; kludge for fixing old tag names
      if apred_version eq 'r3'  or apred_version eq 'r' then begin
        oldtags=['obj','fiber','locid','aktarg','akmethod','akwise','targ1','targ2']
        newtags=['apogee_id','fiberid','location_id','ak_targ','ak_targ_method','ak_wise','apogee_target1','apogee_target2']
        new=struct_rename_tags(str,oldtags,newtags)
        str=new
      endif
      ; kludge for fixing plate tags
      if keyword_set(fixplate) then begin
        oldtags=['plate']
        newtags=['cplate']
        new=struct_rename_tags(str,oldtags,newtags,addtag='plate',addval=0)
        new.plate=fix(new.cplate)
        remove_tags,new,'cplate',str
      endif
      junk= tag_exist(str,'apogee_id',index=objind) 
      junk= tag_exist(str,'location_id',index=locind) 
      if ~tag_exist(str,'obsvrel') then begin
        add_tag,str,'obsvrel',0.,str
        add_tag,str,'obsvrelerr',0.,str
        add_tag,str,'obsvhelio',0.,str
      endif
      ;; fix bad survey for apo1m
      fix = where(str.location_id eq 1, nfix)
      if nfix gt 0 then str[fix].survey = 'apo1m'
 
      if (status lt 0) or (size(str,/type) ne 8) then begin
         message,'Problems reading '+visitfiles[j]
      endif

      ;; Get character field name and other info
      ;if itelescope eq 0 then begin
      ;  field=strtrim(apogee_field(long(ilocationid),str[0].plate),2) 
      ;  if ilocationid ne str[0].(locind) then begin
      ;    print, 'location id does not match: ', ilocationid, str[0].(locind)
      ;    stop
      ;    for  jj=0,n_elements(str) do str[i].(locind) = ilocationid
      ;  endif
      ;endif else field=ilocationid
      field = ilocationid

      ;; Get uniq plate/MJD combinations for this location (not apo1m)
      if telescope ne 'apo1m' then begin
       plate_mjd = strtrim(string(str.plate),2)+'_'+strtrim(string(str.mjd),2)
       ind=uniq(plate_mjd,sort(plate_mjd))
       nplates+=n_elements(ind)
       for iplate=0,n_elements(ind)-1 do begin
        print,str[ind[iplate]].plate,str[ind[iplate]].mjd
        plate=str[ind[iplate]].plate
        mjd=str[ind[iplate]].mjd

        ;; make location structure
        iplan= where(plans.plateid eq plate, nplan)
        if(nplan ne 1) then  printf,missing, 'Missing or multiple plate in platePlans.par!',plate
        tmp_plates= allplate_blank()
        tmp_plates.name= field
        tmp_plates.plate_visit_id= $
           apogee_visit_id(plate=plate, $
                         mjd=mjd, $
                         apred_version=apred_version, $
                         commissioning=visitcommiss[j], $
                         telescope=telescope)
        ;tmp_plates.location_id= long(ilocationid)
        tmp_plates.location_id= plans[iplan].locationid
        tmp_plates.plate= plate
        tmp_plates.mjd= mjd
        tmp_plates.apred_version= apred_version
        if nplan eq 1 then begin
         tmp_plates.racen= plans[iplan].racen
         tmp_plates.deccen= plans[iplan].deccen
         tmp_plates.survey= plans[iplan].survey
         tmp_plates.programname= plans[iplan].programname
         tmp_plates.platerun= plans[iplan].platerun
         tmp_plates.chunk= plans[iplan].chunk
         tmp_plates.ha= plans[iplan].ha
         tmp_plates.designid= plans[iplan].designid
        endif
        p100= plate/100
        holefile=getenv('PLATELIST_DIR')+'/plates/'+ $
                 strtrim(string(p100, f='(i4.4)'),2)+'XX/'+ $
                 strtrim(string(plate, f='(i6.6)'),2)+'/'+ $
                 'plateHolesSorted-'+ $
                 strtrim(string(plate, f='(i6.6)'),2)+'.par'
        if(NOT file_test(holefile)) then $
         message, 'No such file: '+holefile
        yanny_read, holefile, hdr=hdr, /anon
        tmp_plates.nstandard= long(yanny_par(hdr,'napogee_standard'))
        tmp_plates.nscience= long(yanny_par(hdr,'napogee_science'))
        tmp_plates.nsky= long(yanny_par(hdr,'napogee_sky'))
        if(keyword_set(yanny_par(hdr,'tilerad'))) then $
           tmp_plates.radius= float(yanny_par(hdr,'tilerad')) $
        else $
           tmp_plates.radius= 1.49
        tmp_plates.platedesign_version= yanny_par(hdr,'platedesign_version')
        push,allplates,tmp_plates,count=count
        if count lt 0 then $
           message, 'Error in pushing'
       endfor

       idesign=where(design1.design_id eq plans[iplan].designid,nd)
       idesign=idesign[0]
       if nd eq 0 then begin
        ; APOGEE-2
        idesign=where(design2.design_id eq plans[iplan].designid,nd)
        idesign=idesign[0]
        if nd eq 0 then printf,missing,'no design found for ',field
        ;if nd eq 0 then stop,'no design found!'
        design=design2
        min_h=design[idesign].cohort_min_h
        max_h=design[idesign].cohort_max_h
       endif else begin
        ; APOGEE-1
        design=design1
        min_h=[design[idesign].short_cohort_min_h,design[idesign].medium_cohort_min_h,design[idesign].long_cohort_min_h]
        max_h=[design[idesign].short_cohort_max_h,design[idesign].medium_cohort_max_h,design[idesign].long_cohort_max_h]
        min_jk=design[idesign].DEREDDENED_MIN_J_KS_COLOR
        max_jk=9999.99
       endelse
      endif

      ;; Make structure of this visit for this location
      tmp_allvisitloc= replicate(allvisit_blank(str[0], unsafe_synth=unsafe_synth), $
                                 n_elements(str))
      struct_assign, str, tmp_allvisitloc, /nozero
      ; reset field name with latest conventions from platePlans
      tmp_allvisitloc.field= field
      ; add commissioniong tag info
      k=where(tmp_allvisitloc.mjd le 55761L,nk)
      if nk gt 0 then tmp_allvisitloc[k].commiss=1

      ; add extratarg tag info
      for k=0,n_elements(tmp_allvisitloc)-1 do begin
        extratarg=0
        if (tmp_allvisitloc[k].commiss eq 0 and  $
           ((tmp_allvisitloc[k].apogee_target1 and 2L^11) gt 0 or $
            (tmp_allvisitloc[k].apogee_target1 and 2L^12) gt 0 or $
            (tmp_allvisitloc[k].apogee_target1 and 2L^13) gt 0)) then begin
          extratarg=0
        endif else begin
          if tmp_allvisitloc[k].commiss eq 1 then extratarg = extratarg or 2
          if strpos(tmp_allvisitloc[k].targflags,'_TELLURIC') ge 0 then extratarg= extratarg or 4
          if strtrim(tmp_allvisitloc[k].telescope,2) eq 'apo1m' then extratarg= extratarg or 8
          if extratarg eq 0 then extratarg = extratarg or 1
        endelse
        tmp_allvisitloc[k].extratarg = extratarg   
        ; get cohort maximum magnitude
        icohort=-1
        if  (tmp_allvisitloc[k].apogee_target1 and 2L^11) gt 0  then icohort=0
        if  (tmp_allvisitloc[k].apogee_target1 and 2L^12) gt 0  then icohort=1
        if  (tmp_allvisitloc[k].apogee_target1 and 2L^13) gt 0  then icohort=2
        if icohort ge 0 then begin
          if strpos(tmp_allvisitloc[k].survey,'apogee2') ge 0 then begin
            if  (tmp_allvisitloc[k].apogee_target1 and 2L^0) gt 0  then min_jk=0.5
            if  (tmp_allvisitloc[k].apogee_target1 and 2L^1) gt 0  then min_jk=0.5
            if  (tmp_allvisitloc[k].apogee_target1 and 2L^2) gt 0  then min_jk=0.8
            if  (tmp_allvisitloc[k].apogee_target1 and 2L^16) gt 0  then min_jk=0.3
            if  (tmp_allvisitloc[k].apogee_target1 and 2L^1) gt 0  then max_jk=0.8 else max_jk=9999.99
          endif 
          tmp_allvisitloc[k].min_h = min_h[icohort]
          tmp_allvisitloc[k].max_h = max_h[icohort]
          tmp_allvisitloc[k].min_jk = min_jk
          tmp_allvisitloc[k].max_jk = max_jk
        endif
      endfor
      ; move targflags for apogee2
      fix = where(strpos(str.survey,'apogee2') ge 0, nfix)
      if nfix gt 0 then begin
        tmp_allvisitloc[fix].apogee2_target1 = tmp_allvisitloc[fix].apogee_target1
        tmp_allvisitloc[fix].apogee2_target2 = tmp_allvisitloc[fix].apogee_target2
        tmp_allvisitloc[fix].apogee2_target3 = tmp_allvisitloc[fix].apogee_target3
        tmp_allvisitloc[fix].apogee_target1 = 0
        tmp_allvisitloc[fix].apogee_target2 = 0
        tmp_allvisitloc[fix].apogee_target3 = 0
      endif
 

      push,allvisitloc,tmp_allvisitloc,count=count
      if count lt 0 then $
         message, 'Error in pushing'
   endfor

   splog, strtrim(i+1,2),'/',strtrim(nlocationdirs,2), $
          ' LocationID=',ilocationid,'  ',strtrim(nstarfiles,2), $
          ' apField files   ',n_elements(allstar)
   if (size(allvisitloc,/type) gt 0) then begin
     junk=tag_exist(allvisitloc,'apogee_id',index=objind) 
     if objind lt 0 then junk=tag_exist(allvisitloc,'obj',index=objind) 
     junk=tag_exist(allvisitloc,'location_id',index=locind) 
     if locind lt 0 then junk=tag_exist(allvisitloc,'locid',index=locind) 
   endif
   
   ;; Loop through the apField files
   for j=0L, nstarfiles-1L do begin
      apgundef,str,status
      str = MRDFITS(starfiles[j],1,/silent,status=status)
      if status lt 0 or size(str,/type) ne 8 then begin
         message, 'Problems reading '+starfiles[j]
      endif

      if keyword_set(notelescope) then $
      remove_tags,str,['ak_targ','ak_targ_method','ak_wise','sfd_ebv','telescope'],new else $
      remove_tags,str,['ak_targ','ak_targ_method','ak_wise','sfd_ebv'],new
      str=new
      ;; fix tag names and the SUSPECT_RV_COMBINATION flag for s3
      if apstar_version eq 's3' then begin
        oldtags=['obj','locid','aktarg','akmethod','akwise']
        newtags=['apogee_id','location_id','ak_targ','ak_targ_method','ak_wise']
        new=struct_rename_tags(str,oldtags,newtags,addtag='sfd_ebv',addval=-99.)
        str=new
        bd=where((str.starflag and starflagval('SUSPECT_RV_COMBINATION')) gt 0 and $
                  str.synthscatter lt 1,nbd)
        if nbd gt 0 then str[bd].starflag = str[bd].starflag and not starflagval('SUSPECT_RV_COMBINATION')
      endif
      ;; fix bad survey for location_id 5042
      fix = where(str.location_id eq 5042, nfix)
      if nfix gt 0 then str[fix].survey = 'apogee2'

      ; make location id compatible types -- should be done in reduction
      if size(str.location_id,/type) eq 2 then begin
;stop
        locid=long(str[0].location_id)
        remove_tags,str,['location_id'],new
        str=new
        add_tag,str,'location_id',locid,str
      endif
    
      ;; see if we can find information from an aspcapField files
      if strpos(starfiles[j],'apFieldC') ge 0 then $
        file = aspcap_dir+telescope+'/'+ilocationid+'C/aspcapFieldC-'+ilocationid+'C.fits' else $
        file = aspcap_dir+telescope+'/'+ilocationid+'/aspcapField-'+ilocationid+'.fits' 
      splog,'Looking for aspcap: ', file
      haveaspcap = 0
      if file_test(file) ne 0 then begin
         haveaspcap=1
         ;; read library parameters
         aspcap=mrdfits(file,2)
         ;libr_path= apogee_dir+'/speclib/'
         ;libfile= strtrim(libr_path+file_basename(aspcap[0].grid,'.dat'),2)
         ;rdlibhead,libfile,libhead0,libhead

         ;; read ASPCAP output
         aspcap=mrdfits(file,1)
         aspcaplabs=mrdfits(file,3)

         ; calibrate now (again?) if requested, but only for final summary files!
         if keyword_set(caldir) then aspcap_correct,aspcap,aspcaplabs.elem_symbol,aspcap_dir+'/'+caldir+'/'

         index=intarr(nparam)
         ;lindex=intarr(nparam)
         for ipar=0L, nparam-1L do begin
           index[ipar]= where(strtrim(aspcaplabs.param_symbol,2) eq params[ipar])
           ;lindex[ipar]=where(strtrim(libhead0.label,2) eq params[ipar])
         endfor
         eindex=intarr(nelem)
         if tag_exist(aspcaplabs,'elem_symbol') then begin
           for ipar=0L, nelem-1L do eindex[ipar]= where(strtrim(aspcaplabs.elem_symbol,2) eq elems[ipar])
         endif else haveaspcap=0
      endif
      if ~haveaspcap then printf,missing,'No ASPCAP file for ', field

      ;; loop over the stars and transfer to new structure
      tmp_str=allstar_blank(str[0], maxvisit=maxvisit, nparam=nparam, nelem=nelem, nclass=nclass, nwind=nwind)
      tmp_allstarloc= replicate(tmp_str,n_elements(str))
      struct_assign, str, tmp_allstarloc, /nozero

      ; reset field name with latest conventions from platePlans
      tmp_allstarloc.field = field
      ;if strtrim(tmp_allstarloc[0].telescope,2) ne 'apo1m' then begin
      ;  if tmp_allstarloc[0].location_id ne ilocationid then begin
      ;    print,'star location does not match: ', tmp_allstarloc[0].location_id,ilocationid
      ;    stop
      ;  endif
      ;endif
      
      ;; Add visits reference indices
      for k=0L, n_elements(tmp_allstarloc)-1L do begin
         if (size(allvisitloc,/type) gt 0) then begin
           ; fill all_visits and visits (latter only gives good visits) - NOW DONE BELOW
           pk=where(allvisitloc.(objind) eq tmp_allstarloc[k].apogee_id,nk)
           if nk gt maxvisit then begin
              message, 'Error: need to increase maxvisit! '+ maxvisit+ nk
           endif
           ;if nk gt 0 then visits = strmid(file_basename(allvisitloc[pk].file,'.fits'),8) else visits = ' '
           ;tmp_allstarloc[k].all_visits = strjoin(visits,',')
           ;pk=where(allvisitloc.(objind) eq tmp_allstarloc[k].apogee_id and $
           ;        (finite(allvisitloc.vrel) eq 1) and (allvisitloc.commiss eq tmp_allstarloc[k].commiss),nk)
           ;if nk gt 0 then visits = strmid(file_basename(allvisitloc[pk].file,'.fits'),8) else visits = ' '
           ;tmp_allstarloc[k].visits = strjoin(visits,',')

           ; fix NVISIT=1 RVs for s3
           if tmp_allstarloc[k].nvisits eq 1 and (apstar_version eq 's3' or apstar_version eq 's') then begin
             gpk=where(finite(allvisitloc[pk].vrel) eq 1,ngpk)
             if ngpk gt 0 then begin
              tmp_allstarloc[k].vhelio_avg = allvisitloc[pk[gpk[0]]].vhelio
              tmp_allstarloc[k].verr = allvisitloc[pk[gpk[0]]].synthvrelerr
              tmp_allstarloc[k].verr_med = allvisitloc[pk[gpk[0]]].synthvrelerr
              tmp_allstarloc[k].synthvhelio_avg = allvisitloc[pk[gpk[0]]].synthvhelio
              tmp_allstarloc[k].synthverr = allvisitloc[pk[gpk[0]]].synthvrelerr
              tmp_allstarloc[k].synthverr_med = allvisitloc[pk[gpk[0]]].synthvrelerr
             endif
           endif

         endif
          ; add extratarg tag info
         extratarg=0
         if (tmp_allstarloc[k].commiss eq 0 and  $
            ((tmp_allstarloc[k].apogee_target1 and 2L^11) gt 0 or $
             (tmp_allstarloc[k].apogee_target1 and 2L^12) gt 0 or $
             (tmp_allstarloc[k].apogee_target1 and 2L^13) gt 0)) then begin
           extratarg=0
         endif else begin
           if tmp_allstarloc[k].commiss eq 1 then extratarg = extratarg or 2
           if strpos(tmp_allstarloc[k].targflags,'_TELLURIC') ge 0 then extratarg= extratarg or 4
           if strtrim(tmp_allstarloc[k].telescope,2) eq 'apo1m' then extratarg= extratarg or 8
           if extratarg eq 0 then extratarg = extratarg or 1
         endelse
         tmp_allstarloc[k].extratarg = extratarg   
         tmp_allstarloc[k].min_h = min(allvisitloc[pk].min_h)
         tmp_allstarloc[k].max_h = max(allvisitloc[pk].max_h)
         tmp_allstarloc[k].min_jk = min(allvisitloc[pk].min_jk)
         tmp_allstarloc[k].max_jk = max(allvisitloc[pk].max_jk)
         if max(allvisitloc[pk].max_h) ne min(allvisitloc[pk].max_h) then begin
           print,'inconsistent max_h: ',tmp_allstarloc[k].apogee_id,allvisitloc[pk].max_h
           printf,dup,'inconsistent max_h: ',tmp_allstarloc[k].apogee_id,allvisitloc[pk].max_h
         endif

         ;; add ASPCAP parameters if we have them and this is not commissioning data
         ;if file_test(file) and tmp_allstarloc[k].commiss eq 0 then begin
         ;; add ASPCAP parameters if we have them
         if haveaspcap then begin
            kk=where(strtrim(aspcap.apogee_id,2) eq strtrim(tmp_allstarloc[k].apogee_id,2),nk)
            if nk gt 1 then  begin
               if not keyword_set(override) then $
                  message, 'Something went wrong: multiple cases of this object' $
               else begin
                  kk=kk[nk-1]
                  nk=1
               endelse
            endif
            if nk eq 1 then begin
              ; Assign values to the _named_ parameter tags
              tagnames=tag_names(tmp_allstarloc)
              for iparam=0,n_elements(params)-1 do begin
                itag=where(strtrim(tagnames,2) eq strtrim(paramtags[iparam],2),ntag)
                ierrtag=where(strtrim(tagnames,2) eq strtrim(paramtags[iparam],2)+'_ERR')
                if ntag gt 0 and index[iparam] ge 0 then begin
                  tmp_allstarloc[k].(itag)=aspcap[kk].param[index[iparam]]
                  if aspcap[kk].param_cov[index[iparam],index[iparam]] gt 0 then $
                      tmp_allstarloc[k].(ierrtag) = $
                  sqrt(aspcap[kk].param_cov[index[iparam],index[iparam]])
                endif
              endfor
              ; vmicro, vmacro/vsini different for dwarfs and giants
              tmp_allstarloc[k].vmicro=10.^aspcap[kk].param[2]
              if (strpos(aspcap[kk].class,'GKd') ge 0  or strpos(aspcap[kk].class,'Fd') ge 0 or strpos(aspcap[kk].class,'Md') ge 0) then begin
                  tmp_allstarloc[k].vmacro=10.^0.
                  tmp_allstarloc[k].vsini=10.^aspcap[kk].param[7]
              endif else begin
                  tmp_allstarloc[k].vmacro=10.^aspcap[kk].param[7]
              endelse

              ; Assign values to the _named_ element tags
              ife = where(strtrim(aspcaplabs.elem_symbol,2) eq 'Fe')
              for ielem=0,n_elements(elems)-1 do begin
                itag=where(strtrim(tagnames,2) eq strtrim(elemtags[ielem],2))
                ierrtag=where(strtrim(tagnames,2) eq strtrim(elemtags[ielem],2)+'_ERR')
                iflagtag=where(strtrim(tagnames,2) eq strtrim(elemtags[ielem],2)+'_FLAG')
                if eindex[ielem] ge 0 then begin
                  ;tmp_allstarloc[k].x_h[eindex[ielem]]=aspcap[kk].elem[eindex[ielem]]
                  ;tmp_allstarloc[k].x_h_err[eindex[ielem]]=aspcap[kk].x_h_err[eindex[ielem]]
                  ;; if this is C or N and we are in dwarf grid, then parameter is already [X/H]
                  ;if elems[ielem] eq 'C' or elems[ielem] eq 'CI' or elems[ielem] eq 'N' and $
                  ;  (strpos(aspcap[kk].class,'GKd') ge 0  or strpos(aspcap[kk].class,'Fd') ge 0 or strpos(aspcap[kk].class,'Md') ge 0) then begin
                  ;   tmp_allstarloc[k].x_h[eindex[ielem]]+=0.
                  ;endif else begin
                  ;  if ~elemtoh[ielem] and aspcap[kk].fparam[3] gt -90 and tmp_allstarloc[k].x_h[eindex[ielem]] gt -90 then begin
                  ;   tmp_allstarloc[k].x_h[eindex[ielem]]+=aspcap[kk].fparam[3]
                  ;  endif
                  ;endelse
                  ;tmp_allstarloc[k].x_m[eindex[ielem]]=tmp_allstarloc[k].x_h[eindex[ielem]]-aspcap[kk].fparam[3]
                  ;tmp_allstarloc[k].x_m_err[eindex[ielem]]=aspcap[kk].x_m_err[eindex[ielem]]
                  tmp_allstarloc[k].x_h[eindex[ielem]]=aspcap[kk].x_h[eindex[ielem]]
                  tmp_allstarloc[k].x_h_err[eindex[ielem]]=aspcap[kk].x_h_err[eindex[ielem]]
                  tmp_allstarloc[k].x_m[eindex[ielem]]=aspcap[kk].x_m[eindex[ielem]]
                  tmp_allstarloc[k].x_m_err[eindex[ielem]]=aspcap[kk].x_m_err[eindex[ielem]]

                  if itag ge 0 then begin
                    if ielem eq ife then begin
                      ; Fe is special, since we don't want [Fe/Fe]!
                      tmp_allstarloc[k].(itag)=aspcap[kk].x_h[eindex[ielem]]
                      tmp_allstarloc[k].(ierrtag) = aspcap[kk].x_h_err[eindex[ielem]]
                      tmp_allstarloc[k].(iflagtag) = aspcap[kk].elemflag[eindex[ielem]]
                    endif else begin
                      if aspcap[kk].x_h[eindex[ielem]] gt -9998. then begin  
                        tmp_allstarloc[k].(itag)=aspcap[kk].x_h[eindex[ielem]]-aspcap[kk].x_h[ife]
                        tmp_allstarloc[k].(ierrtag) = aspcap[kk].x_m_err[eindex[ielem]]
                      endif
                      tmp_allstarloc[k].(iflagtag) = aspcap[kk].elemflag[eindex[ielem]]
                    endelse

                    ; special handling for Rb in l31c.2 --> remove it!
                    if strpos(results_version,'l31c') ge 0  then begin
                      irb = where(strtrim(aspcaplabs.elem_symbol,2) eq 'Rb')
                      if ielem eq irb then begin
                        tmp_allstarloc[k].(itag) = -9999.99
                        tmp_allstarloc[k].(ierrtag) = -9999.99
                        tmp_allstarloc[k].x_h[eindex[ielem]] = -9999.99
                        tmp_allstarloc[k].x_h_err[eindex[ielem]] = -9999.99
                        tmp_allstarloc[k].x_m[eindex[ielem]] = -9999.99
                        tmp_allstarloc[k].x_m_err[eindex[ielem]] = -9999.99
                      endif
                      irb = where(strtrim(aspcaplabs.elem_symbol,2) eq 'Na')
                      if ielem eq irb and (tmp_allstarloc[k].x_m[eindex[ielem]] lt -1 or tmp_allstarloc[k].x_h[eindex[ielem]] lt -1) then begin
                        tmp_allstarloc[k].(itag) = -9999.99
                        tmp_allstarloc[k].(ierrtag) = -9999.99
                        tmp_allstarloc[k].x_h[eindex[ielem]] = -9999.99
                        tmp_allstarloc[k].x_h_err[eindex[ielem]] = -9999.99
                        tmp_allstarloc[k].x_m[eindex[ielem]] = -9999.99
                        tmp_allstarloc[k].x_m_err[eindex[ielem]] = -9999.99
                      endif
                    endif
                  endif
                endif
              endfor

              ; add the ASPCAP output tags
              if tag_exist(aspcap,'meanfib') then tmp_allstarloc[k].meanfib=aspcap[kk].meanfib
              if tag_exist(aspcap,'sigfib') then tmp_allstarloc[k].sigfib=aspcap[kk].sigfib
              if tag_exist(aspcap,'snrev') then tmp_allstarloc[k].snrev=aspcap[kk].snrev
              tmp_allstarloc[k].aspcapflag=aspcap[kk].aspcapflag
              tmp_allstarloc[k].aspcapflags=aspcapflag(aspcap[kk].aspcapflag,0)
              tmp_allstarloc[k].paramflag=aspcap[kk].paramflag
              tmp_allstarloc[k].aspcap_chi2=aspcap[kk].param_chi2
              tmp_allstarloc[k].aspcap_class=aspcap[kk].class
              tmp_allstarloc[k].param=aspcap[kk].param
              tmp_allstarloc[k].fparam=aspcap[kk].fparam
              if n_elements(nclass) gt 0 and tag_exist(aspcap,'fparam_class') then begin
                tmp_allstarloc[k].fparam_class=aspcap[kk].fparam_class
                tmp_allstarloc[k].chi2_class=aspcap[kk].chi2_class
              endif
              tmp_allstarloc[k].param_cov=aspcap[kk].param_cov
              tmp_allstarloc[k].fparam_cov=aspcap[kk].fparam_cov
              ;if tag_exist(aspcap,'elem') then tmp_allstarloc[k].elem=aspcap[kk].elem
              ;if tag_exist(aspcap,'elem_err') then tmp_allstarloc[k].elem_err=aspcap[kk].elem_err
              if tag_exist(aspcap,'felem') then tmp_allstarloc[k].felem=aspcap[kk].felem
              if tag_exist(aspcap,'felem_err') then tmp_allstarloc[k].felem_err=aspcap[kk].felem_err
              if tag_exist(aspcap,'elem_chi2') then tmp_allstarloc[k].elem_chi2=aspcap[kk].elem_chi2
              if tag_exist(aspcap,'elemflag') then tmp_allstarloc[k].elemflag=aspcap[kk].elemflag
            endif
         endif 

      endfor
      ; move targflags for apogee2, now already done in apstar!
      ;fix = where(strpos(str.survey,'apogee2') ge 0, nfix)
      ;if nfix gt 0 then begin
      ;  tmp_allstarloc[fix].apogee2_target1 = tmp_allstarloc[fix].apogee_target1
      ;  tmp_allstarloc[fix].apogee2_target2 = tmp_allstarloc[fix].apogee_target2
      ;  tmp_allstarloc[fix].apogee2_target3 = tmp_allstarloc[fix].apogee_target3
      ;  tmp_allstarloc[fix].apogee_target1 = 0
      ;  tmp_allstarloc[fix].apogee_target2 = 0
      ;  tmp_allstarloc[fix].apogee_target3 = 0
      ;endif

      PUSH,allstarloc,tmp_allstarloc,count=count
   endfor
   ;; populate the structures with catalogdb information

   ;; First visits
  if (size(allvisitloc,/type) gt 0) then begin
   splog,'Getting unique star names'
   stars=uniq(allvisitloc.(objind), sort(allvisitloc.(objind)))
   catmin=fltarr(n_elements(stars))+9999.99
   
   ; get apogeeObject catalog info for this field
;   if strpos(allvisitloc[0].survey,'apogee2') ge 0 then apogeeobject='apogee2Object' $
;     else if strpos(allvisitloc[0].survey,'apo1m') ge 0 then apogeeobject='apogee1mObject' else $
;     apogeeobject='apogeeObject'
;   objectfile=targetdir+'/'+apogeeobject+'/'+apogeeobject+'_'+strtrim(apogee_field(allvisitloc[0].(locind),allvisitloc[0].plate,/addloc),2)+'.fits' 

   ;if itelescope eq 0 then $
   ;objectfile=targetdir+'/'+apogeeobject+'/'+apogeeobject+'_'+strtrim(apogee_field(allvisitloc[0].(locind),allvisitloc[0].plate,/addloc),2)+'.fits' else $
   ;objectfile=targetdir+'/apogeeObject/apogeeObject_'+strtrim(field,2)+'.fits'
   ;if not file_test(objectfile) then begin
   ;  ;print,'cant find apogeeObject file: ', objectfile, allvisitloc[0].(locind), allvisitloc[0].plate
   ;  ;printf,missing,'cant find apogeeObject file: ', objectfile, allvisitloc[0].(locind), allvisitloc[0].plate
   ;  ; try apogeeObject
   ;  objectfile=targetdir+'/apogeeObject/apogeeObject_'+strtrim(apogee_field(allvisitloc[0].(locind),allvisitloc[0].plate),2)+'.fits' 
   ;endif 
;   if not file_test(objectfile) then begin

     ; find all matching apogeeObject files and loop through them looking for matches
     files=file_search(getenv('APOGEE_TARGET')+'/apogee*Object/*'+field+'*')
     if files[0] eq '' then begin
       print,'cant find apogeeObject file: '+field, allvisitloc[0].(locind), allvisitloc[0].plate
       printf,missing,'cant find apogeeObject file: '+ field, allvisitloc[0].(locind), allvisitloc[0].plate
     endif else begin
       ;objects=mrdfits(files[0],1)
       ;for ifile=1,n_elements(files)-1 do begin
       ;  tmp=mrdfits(files[ifile],1)
       ;  objects=[objects,tmp]
       ;endfor
      if n_elements(files) gt 1 then printf,missing,'using multiple apogeeObject files: '+ files

      for ifile=0,n_elements(files)-1 do begin
       objects=mrdfits(files[ifile],1)
;   endif else begin
;    objects=mrdfits(objectfile,1)
;    locid=allvisitloc[0].(locind)
;    if locid eq 2111 or locid eq  2119 or locid eq  2120 or locid eq  2121 or locid eq  2122 or locid eq  2382 then begin
;      files=file_search(getenv('APOGEE_TARGET')+'/apogee2Object/*'+string(format='(i4.4)',locid)+'*')
;      objects=mrdfits(files[0],1)
;      for ifile=1,n_elements(files)-1 do begin
;        tmp=mrdfits(files[ifile],1)
;        objects=[objects,tmp]
;      endfor
;    endif

    ; fix NaNs, etc.
    aspcap_fixobject,objects

    ; fix any negative RAs
    ra=allvisitloc.ra
    dec=allvisitloc.dec
    bd=where(ra lt 0, nbd)
    if nbd gt 0 then ra[bd] = 0.
    if nbd gt 0 then dec[bd] = 0.
    ;spherematch,objects.ra,objects.dec,ra,dec,10./3600.,match1,match2,dist,maxmatch=nplates
    spherematch,objects.ra,objects.dec,ra,dec,2./3600.,match1,match2,dist,maxmatch=nplates

    for istar=0L,n_elements(stars)-1L do begin
      if istar mod 20 eq 0 then $
         splog,istar,n_elements(stars)
      ; first try to match by position (in case names are screwed up)
      j=where(match2 eq stars[istar],nj)
      if allvisitloc[stars[istar]].ra gt 0 and nj gt 0 then begin
        ; if more than one match, take the closest
        if nj gt 1 then junk=min(dist(j),jj) else jj=0
        if dist[j[jj]] lt catmin[istar] then begin
          iobject=match1[j[jj]] 
          dmin=dist[j[jj]]
        endif else iobject=-1
      endif else begin
        ; if position match fails, try to match by name
        iobject=where((strtrim(objects.apogee_id,2) eq strtrim(allvisitloc[stars[istar]].apogee_id,2)) or $
                      (strtrim(objects.alt_id,2) eq strtrim(allvisitloc[stars[istar]].apogee_id,2)),nj) 
        dmin=-1
        if nj gt 1 then printf,missing,'more than one object found by name!'
        if nj eq 0 then iobject=-1
      endelse
      objname=strtrim(allvisitloc[stars[istar]].apogee_id,2)
      if iobject ge 0 then begin
         catmin[istar]=dmin
         if strtrim(allvisitloc[stars[istar]].apogee_id,2) ne strtrim(objects[iobject].apogee_id,2) then begin
           printf,altname, 'alt name: ',field,' ',$
                allvisitloc[stars[istar]].(locind),' ',$
                allvisitloc[stars[istar]].(objind),' ', $
                objects[iobject].apogee_id
           have_altname = 1 
         endif else have_altname=0
         cat=objects[iobject]
         if strlen(strtrim(cat.apogee_id,2)) eq 19 and strmid(cat.apogee_id,0,2) eq '2M' then begin
           cat.apogee_id = strmid(cat.apogee_id,0,18)
           printf,altname, 'alt name trim character: ',field,' ',$
                allvisitloc[stars[istar]].(locind),' ',$
                allvisitloc[stars[istar]].(objind),' ', $
                cat.apogee_id
           have_altname = 1 
         endif
         j=where(strtrim(allvisitloc.(objind),2) eq strtrim(allvisitloc[stars[istar]].(objind),2),nj)
         for jj=0,n_elements(j)-1 do begin
            jjj=j[jj]
            tmp_cat= catalog_info_blank()
            tmp_allvisitloc= allvisitloc[jjj]
            struct_assign, cat, tmp_cat
            struct_assign, tmp_cat, tmp_allvisitloc, /nozero
            catalog_info_replace,cat,tmp_allvisitloc,missing=missing
            allvisitloc[jjj]= tmp_allvisitloc
            allvisitloc[jjj].reduction_id = allvisitloc[jjj].apogee_id
            if have_altname then begin
              allvisitloc[jjj].target_id = allvisitloc[jjj].apogee_id
              allvisitloc[jjj].apogee_id = cat.apogee_id
            endif
         endfor
     
         if nstarfiles gt 0 then begin
            ;j=where((strtrim(allvisitloc.(objind),2) eq strtrim(allvisitloc[stars[istar]].(objind),2)) and $
            ;         (allvisitloc.vtype gt 0),nj)
            ;if nj gt 0 then visits = strmid(file_basename(allvisitloc[j].file,'.fits'),8)
            
            j=where(strtrim(allstarloc.apogee_id,2) eq objname,nj)
            if nj eq 0 then begin
               splog,'No apstar found for object ',allvisitloc[stars[istar]].(objind),allvisitloc[stars[istar]].mjd
               splog,'flag: ',starflag(allvisitloc[stars[istar]].starflag)
               if allvisitloc[stars[istar]].mjd gt 55761 and $
                  (allvisitloc[stars[istar]].starflag and $
                   badstarflag() eq 0) then begin
                      message, 'Uh oh! Expected apstar information here!'
                  printf,message,'No apstar found for object ',allvisitloc[stars[istar]].(objind),allvisitloc[stars[istar]].mjd
               endif
            endif
            for jj=0,nj-1 do begin
               jjj=j[jj]
               tmp_cat= catalog_info_blank()
               tmp_allstarloc= allstarloc[jjj]
               struct_assign, cat, tmp_cat
               struct_assign, tmp_cat, tmp_allstarloc, /nozero
               catalog_info_replace,cat,tmp_allstarloc,missing=missing
               allstarloc[jjj]= tmp_allstarloc
               allstarloc[jjj].reduction_id = allstarloc[jjj].apogee_id
               if have_altname then begin
                 allstarloc[jjj].apogee_id = cat.apogee_id
               endif
            endfor
         endif
      endif ;else begin
      ;   splog,'missing ',allvisitloc[stars[istar]].(objind),' ', $
      ;         allvisitloc[stars[istar]].(locind)
      ;   printf,missing, 'missing: ', $
      ;          field,' ',$
      ;          allvisitloc[stars[istar]].(locind),' ',$
      ;          allvisitloc[stars[istar]].(objind),' ',$
      ;          allvisitloc[stars[istar]].ra,' ',$
      ;          allvisitloc[stars[istar]].dec,' ',$
      ;          allvisitloc[stars[istar]].h,' ',$
      ;          allvisitloc[stars[istar]].snr
      ;endelse
     endfor  ; loop over stars
    endfor  ; loop over object files
    for istar=0,n_elements(stars)-1 do begin
      if catmin[istar] gt 9999 then begin
         splog,'missing ',allvisitloc[stars[istar]].(objind),' ', $
               allvisitloc[stars[istar]].(locind)
         printf,missing, 'missing: ', $
                field,' ',$
                allvisitloc[stars[istar]].(locind),' ',$
                allvisitloc[stars[istar]].(objind),' ',$
                allvisitloc[stars[istar]].ra,' ',$
                allvisitloc[stars[istar]].dec,' ',$
                allvisitloc[stars[istar]].h,' ',$
                allvisitloc[stars[istar]].snr
      endif
    endfor
   endelse  ; if have any object files
  endif

   PUSH,allvisit,allvisitloc,count=count
   if count lt 0 then message, 'Error in combining with push'
   if nstarfiles gt 0 then begin
      PUSH,allstar,allstarloc,count=count
      if count lt 0 then message, 'Error in combining with push'
   endif

   nextloc:   
endfor ; locationdirs loop

;endfor ; telescope loop

free_lun,missing
free_lun,altname

;; set IDs based on final names
for k=0L,n_elements(allstar)-1 do begin
   if allstar[k].location_id eq 1 then telescope = 'apo1m' else telescope = 'apo25m'
   allstar[k].apogee_id = strtrim(allstar[k].apogee_id,2)
   allstar[k].apstar_id= $
            apogee_apstar_id(locid=allstar[k].field, $;locid= allstar[k].location_id, $
                             star= allstar[k].apogee_id, $
                             apstar_version=apstar_version, $
                             commissioning=allstar[k].commiss, $
                             telescope=telescope)
   allstar[k].target_id= $
            apogee_target_id(locid=allstar[k].field, $;locid= allstar[k].location_id, $
                             star= allstar[k].apogee_id, $
                             field= allstar[k].field)
   allstar[k].aspcap_id= $
            apogee_aspcap_id(locid=allstar[k].field, $;locid= allstar[k].location_id, $
                             star= allstar[k].apogee_id, $
                             results_version=results_version, $
                             commissioning=allstar[k].commiss, $
                             telescope=telescope)
endfor
for k=0L,n_elements(allvisit)-1 do begin
   ;if allvisit[k].location_id eq 1 then telescope = 'apo1m' else telescope = 'apo25m'
   allvisit[k].apogee_id = strtrim(allvisit[k].apogee_id,2)
   allvisit[k].target_id= $
            apogee_target_id(locid=allvisit[k].field, $;locid= allvisit[k].location_id, $
                             star= allvisit[k].apogee_id, $
                             field= allvisit[k].field)
   allvisit[k].visit_id= $
         apogee_visit_id(plate=allvisit[k].plate, $
                         mjd=allvisit[k].mjd, $
                         file=allvisit[k].file, $
                         fiberid=allvisit[k].fiberid, $
                         apogeeid=allvisit[k].apogee_id, $
                         apred_version=apred_version, $
                         commissioning=allvisit[k].commiss, $
                         telescope=allvisit[k].telescope)
endfor

;; Set APOSPECOBJID and reduction version
red= replicate(long(strmid(apred_version, 1)), n_elements(allvisit))
allvisit.apred_version= apred_version

;; Set apstar, aspcap, results versions
red= long(strmid(apred_version, 1))
allstar.apstar_version= apstar_version
allstar.aspcap_version= aspcap_version
allstar.results_version= results_version

; Remove duplicates (that shouldn't be there!)
u=rem_dup(allstar.apstar_id,allstar.snr)
if n_elements(u) ne n_elements(allstar) then begin
  for i=0L,n_elements(allstar)-1 do begin
    junk=where(u eq i,nj)
    if nj eq 0 then begin
      j=where(allstar[u].apstar_id eq allstar[i].apstar_id)
      printf,dup,'dup: ',i,' ',allstar[i].apogee_id,' ',allstar[i].reduction_id,' ',allstar[i].location_id,' ',allstar[i].field
      printf,dup,'           ',allstar[u[j]].apogee_id,' ',allstar[u[j]].reduction_id,' ',allstar[u[j]].location_id,' ',allstar[u[j]].field
    endif
  endfor
  allstar=allstar[u]
endif

; Flag duplicate objects from different fields (but don't remove them)
u=rem_dup(allstar.apogee_id,allstar.snr)
ndup=0
if n_elements(u) ne n_elements(allstar) then begin
  for i=0L,n_elements(allstar)-1 do begin
    junk=where(u eq i,nj)
    if nj eq 0 then begin
      j=where(allstar[u].apogee_id eq allstar[i].apogee_id)
      print,'dup: ',i,' ',allstar[i].apogee_id,' ',allstar[i].reduction_id,' ',allstar[i].location_id,' ',allstar[i].field,allstar[i].snr
      print,'           ',allstar[u[j]].apogee_id,' ',allstar[u[j]].reduction_id,' ',allstar[u[j]].location_id,' ',allstar[u[j]].field,allstar[u[j]].snr
      ; if stars are from same field, are they just commissioning?
      if allstar[i].location_id ne allstar[u[j]].location_id then $
        allstar[i].extratarg = allstar[i].extratarg or 16 $
      else begin
        if allstar[i].commiss eq 0 and allstar[u[j]].commiss eq 0 then begin
          print,'WARNING: duplicate non-commissioning stars in same field' 
          printf,dup,'WARNING: duplicate non-commissioning stars in same field' 
        endif
      endelse
      ndup+=1
    endif
  endfor
endif
free_lun,dup
print,'flagged duplicate objects: ', ndup

; sort visit file by RA
s=sort(allvisit.ra)
allvisitsort=allvisit[s]
ind=lonarr(360)
for i=1L,359L do begin
   j=where(allvisitsort.ra ge double(i))
   ind[i]=j[0]
endfor

; Sort the allstar table, and put pointers into (now sorted) allvisit table
openw,missing,outdir+'/log/visits-'+results_version+'.txt',/get_lun
s=sort(allstar.ra)
allstarsort=allstar[s]
for i=0L,n_elements(allstarsort)-1 do begin
   ira=fix(allstarsort[i].ra)
   i1=ind[ira]
   if ira lt 359 then i2=ind[ira+1] else i2=n_elements(allvisitsort)-1
   for kkk=0L,n_elements(maxvisit)-1 do $
      allstarsort[i].visit_pk[kkk]=-1
   pk=where(allvisitsort[i1:i2].(objind) eq allstarsort[i].apogee_id,nk)
   if nk gt maxvisit then begin
     printf,missing,'Problem: too many visits? '+ allstarsort[i].field+' '+allstarsort[i].apogee_id
     nk=maxvisit
   endif
   pk+=i1
   if nk gt 0 then begin
     for kkk=0L,nk-1 do allstarsort[i].all_visit_pk[kkk]=pk[kkk]
     visits = strmid(file_basename(allvisitsort[pk].file,'.fits'),8) 
     allstarsort[i].all_visits = strjoin(visits,',')
   endif else begin
     allstarsort[i].visits = ' '
     printf,missing,'No good visits: '+allstarsort[i].field+' '+allstarsort[i].apogee_id
   endelse
   junk=where(allvisitsort[i1:i2].(objind) eq allstarsort[i].apogee_id and $
            (finite(allvisitsort[i1:i2].vrel) eq 1) ,nk2)
   if allstarsort[i].nvisits ne nk then printf,missing,'nvisits does not match all visits: ',allstarsort[i].nvisits,nk,nk2
   pk=where((allvisitsort[i1:i2].target_id eq allstarsort[i].target_id) and $
            (finite(allvisitsort[i1:i2].vrel) eq 1) and $
            (allvisitsort[i1:i2].vrel lt 99999) and $
            (allvisitsort[i1:i2].commiss eq allstarsort[i].commiss),nk)
   pk+=i1
   if nk gt maxvisit then begin
     printf,missing,'Problem: too many visits? '+ allstarsort[i].field+' '+allstarsort[i].apogee_id
     nk=maxvisit
   endif
   if nk gt 0 then begin
     for kkk=0L,nk-1 do allstarsort[i].visit_pk[kkk]=pk[kkk]
     visits = strmid(file_basename(allvisitsort[pk].file,'.fits'),8) 
     allstarsort[i].visits = strjoin(visits,',')
   endif else begin
     allstarsort[i].visits = ' '
     printf,missing,'No visits: '+allstarsort[i].field+' '+allstarsort[i].apogee_id
   endelse
   if allstarsort[i].nvisits ne nk then begin
     printf,missing,'nvisits does not match good visits: '+allstarsort[i].field+' '+string(allstarsort[i].nvisits)+' '+string(nk)
     allstarsort[i].nvisits=nk
   endif
 
endfor
free_lun,missing
; get the RA indices in the sorted array
starind=lonarr(360)
for i=1,359 do begin
   j=where(allstarsort.ra ge i)
   starind[i]=j[0]
endfor

; write out the files: allPlates, allVisit, allStar
if keyword_set(test) then begin
  prefix='test' 
  stop
endif else prefix = 'all'

mjdlast=strtrim(max(allvisitsort.mjd),2)
file_move,outdir+'/log/missing-'+results_version+'.txt',outdir+'/log/missing-'+results_version+'-'+mjdlast+'.txt',/over
file_move,outdir+'/log/altname-'+results_version+'.txt',outdir+'/log/altname-'+results_version+'-'+mjdlast+'.txt',/over
file_move,outdir+'/log/dup-'+results_version+'.txt',outdir+'/log/dup-'+results_version+'-'+mjdlast+'.txt',/over
file_move,outdir+'/log/visits-'+results_version+'.txt',outdir+'/log/visits-'+results_version+'-'+mjdlast+'.txt',/over
if n_elements(allplates) gt 0 then begin
  isort=sort(allplates.plate*100000LL+allplates.mjd)
  allplates=allplates[isort]
  outfile = outdir+'/'+prefix+'Plates-'+results_version+'-'+mjdlast+'.fits'
  mkhdr,header,0
  sxaddhist,'APRED VERSION: '+apred_version,header
  sxaddhist,'APSTAR VERSION: '+apstar_version,header
  sxaddhist,'ASPCAP VERSION: '+aspcap_version,header
  sxaddhist,'RESULTS VERSION: '+results_version,header
  sxaddhist,'APOGEE SOFTWARE VERSION: '+getvers(),header
  sxaddhist,'DATE: '+systime(),header
  MWRFITS,0,outfile,header,/create
  MWRFITS,allplates,outfile
endif

outfile = outdir+'/'+prefix+'Visit-'+results_version+'-'+mjdlast+'.fits'
mkhdr,header,0
sxaddhist,'APRED VERSION: '+apred_version,header
sxaddhist,'APSTAR VERSION: '+apstar_version,header
sxaddhist,'ASPCAP VERSION: '+aspcap_version,header
sxaddhist,'RESULTS VERSION: '+results_version,header
sxaddhist,'APOGEE SOFTWARE VERSION: '+getvers(),header
sxaddhist,'DATE: '+systime(),header
MWRFITS,0,outfile,header,/create
MWRFITS,allvisitsort,outfile
MWRFITS,ind,outfile

outfile = outdir+'/'+prefix+'Star-'+results_version+'-'+mjdlast+'.fits'
MWRFITS,0,outfile,header,/create
MWRFITS,allstarsort,outfile
MWRFITS,starind,outfile
if n_elements(aspcaplabs) gt 0 then begin
  labs={param_symbol: aspcaplabs.param_symbol, elem_symbol: aspcaplabs.elem_symbol, elem_value: aspcaplabs.elem_value, elemtoh: aspcaplabs.elemtoh, classes: aspcaplabs.classes}
  MWRFITS,labs,outfile
endif

; make summary plots
if keyword_set(test) then stop
aspcap_plot,allstarsort,aspcap_dir+'/html/all',/hard,results_vers=results_version,/fit,elems=aspcaplabs.elem_symbol

; make ASCII output files
;aspcap_maketxt,allvisitsort,allstarsort,date=date,outdir=outdir


end
