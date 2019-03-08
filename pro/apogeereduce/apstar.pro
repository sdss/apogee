pro apstar,planfile,clobber=clobber,sinc=sinc,log=log,noplot=noplot,obj=inpobj,rvrefine=rvrefine,snmin=snmin,snsig=snsig,snmax=snmax,nstars=nstars,nres=nres,bccomb=bccomb,stars_dir=stars_dir, trimgrid = trimgrid

;+
;
; APSTAR
;
; This program combines different visits into a single spectrum, given
;   input spectra and measured RVs, with option to refine RVs based on combination
;   Also does binary check from final RVs
;
; INPUT:
;  locations    input list of location IDs
;  /clobber     recombine/overwrite existing files even if no new data exists
;  sinc=        use sinc interpolation (default, set to 0 to turn off)
;  log=         output on log-lambda scale (default, set to 0 to turn off)
;  /noplot      suppresses plot generations
;  obj=obj      input single object name to process (must be in location directory)
;  rvrefine=    refine RVs based on cross correlation between input spectra
;               with final cross correlation of combined spectrum against grid (default, set to 0 to turn off)
;               If set to 2 - refine RVs with the ability to choose between synthetic grid templates and 
;               observed spectra templates on each iteration.
;  snmin=       SNR ratio below which visits are rejected (regardless of any sigma clipping) Set to 0 to never reject a visit. Default: 3 (1 if snsig set)
;  snsig =      If set, reject visits with snr below a certain sigma (calculated from mad.pro) from the median). Default = 0 (No sigma clipping)
;  snmax =      Specify maximum snr for which visits can be rejected. Used only if snsig set. Default = 10. 
;  trimgrid =   If set, this keyword will activiate the aptrimgrid routine in aprvrefine. 
;               This limits the grid templates that the code can select based on the J-K color (for Teff) and Washington DDO photometery (for logg)
;
; OUTPUT:
;   The combined spectra for each star
;   HDU0: flux: combined (row 0), +individual visits (rows 1-nvisits)
;   HDU1: error/variance: combined (row 0), +individual visits (rows 1-nvisits)
;   HDU2: sky: combined (row 0), +individual visits (rows 1-nvisits)
;   HDU3: telluric: combined (row 0), +individual visits (rows 1-nvisits)
;   HDU4: wavelength: combined (row 0), +individual visits (rows 1-nvisits)
;   HDU5: flags: combined (row 0), +individual visits (rows 1-nvisits)
;   HDU6: LSF fit coefficients 
;  The names are apStar-STAR8.fits
;
; USAGE:
;  IDL>ap1dobject,planfiles
;  or  IDL>ap1dobject,platefiles=platfiles
;
; Written by J.Holtzman October 2011
; New keywords added by N.Troup June 2016
;-

; get root directory
aploadplan,planfile,planstr,struct='APFIELD'
if tag_exist(planstr,'telescope') then telescope=planstr.telescope else telescope = 'apo25m'
if tag_exist(planstr,'apred_vers') then apsetver,vers=planstr.apred_vers,telescope=telescope
dirs=getdir(apodir,caldir,spectrodir,apovers,apred_vers=apred_vers)
if tag_exist(planstr,'apstar_vers') then apstar_vers=planstr.apstar_vers else apstar_vers='stars_'+apovers
locations=planstr.apfield.field
if tag_exist(planstr,'apred_vers') then in_vers=planstr.apred_vers else in_vers=apovers
top_dir=apodir+'/'+in_vers+'/'
if tag_exist(planstr,'mjdstart') then mjdstart=planstr.mjdstart else mjdstart=0
if tag_exist(planstr,'mjdend') then mjdend=planstr.mjdend else mjdend=99999
if tag_exist(planstr,'indir') then in_dir=planstr.indir else in_dir=top_dir+'/visit/'
;;if tag_exist(planstr,'starsdir') then stars_dir=planstr.starsdir else stars_dir=top_dir+'/'+apstar_vers+'/'
;if tag_exist(planstr,'starsdir') then stars_dir=planstr.starsdir else stars_dir=top_dir+'/stars/'
in_dir = in_dir + '/' +telescope+'/'
if tag_exist(planstr,'survey') then survey=planstr.survey else survey='apogee'

; set parameter defaults: keyword has highest priority, then parameter file, finally default value
if n_elements(sinc) eq 0 then sinc=apsetpar(planstr,'sinc',1)
if n_elements(log) eq 0 then log=apsetpar(planstr,'log',1)
if n_elements(rvrefine) eq 0 then rvrefine=apsetpar(planstr,'rvrefine',2)  ; 1 for DR13 
if n_elements(snmax) eq 0 then snmax = apsetpar(planstr,'snmax',10) ; 0 for DR13
if n_elements(snsig) eq 0 then snsig = apsetpar(planstr,'snsig',3)  ; 0 for DR13
   if keyword_set(snsig) then defsnmin = 1 else defsnmin = 5
if n_elements(snmin) eq 0 then snmin=apsetpar(planstr,'snmin',defsnmin)
if n_elements(trimgrid) eq 0 then trimgrid=apsetpar(planstr,'trimgrid',0)
if n_elements(nstars) eq 0 then nstars=apsetpar(planstr,'nstars',0)
if n_elements(clobber) eq 0 then clobber=apsetpar(planstr,'clobber',0)
if n_elements(nres) eq 0 then nres=apsetpar(planstr,'nres',0)
if n_elements(bccomb) eq 0 then bccomb=apsetpar(planstr,'bccomb',0)
if n_elements(stars_dir) eq 0 then stars_dir=apsetpar(planstr,'stars_dir',top_dir+'/stars/')
stars_dir = stars_dir + '/' +telescope+'/'

; combine spectra for each of the specified  locations
for iloc=0,n_elements(locations)-1 do begin
  ;cloc=string(format='(i4.4)',locations[iloc])
  clocation=strtrim(locations[iloc],2)
  ncommiss=0
  nsurvey=0
  plot_dir=stars_dir+clocation+'/plots/'
  file_mkdir,plot_dir
  html_dir=stars_dir+clocation+'/html/'
  file_mkdir,html_dir
  openw,html,/get_lun,html_dir+clocation+'.html.tmp'
  if telescope eq 'apo1m' then locid=1L else locid=apogee_locationid(clocation,survey)

  ; find the apVisitSum files for this location, read in visit structure, and append them all
  files=file_search(in_dir+clocation+'/'+'apVisitSum*.fits')
  ifirst=1
  for i=0,n_elements(files)-1 do begin
    print,file_basename(files[i])
    v=mrdfits(files[i],1,/silent)
    if v[0].mjd ge mjdstart and v[0].mjd le mjdend then begin
      if ifirst eq 1 then begin
        allvisits=v 
        ifirst=0
      endif else allvisits=struct_append(allvisits,v)
    endif
  endfor
  if keyword_set(rvrefine) then begin  ; put visit RV estimate in separate structure tags
    if not tag_exist(allvisits,'ESTVREL') then begin
      add_tag,allvisits,'ESTVTYPE',0,allvisits
      add_tag,allvisits,'ESTVREL',0.0,allvisits
      add_tag,allvisits,'ESTVRELERR',0.0,allvisits
      add_tag,allvisits,'ESTVHELIO',0.0,allvisits
      allvisits.estvtype = allvisits.vtype
      allvisits.estvrel = allvisits.vrel
      allvisits.estvrelerr = allvisits.vrelerr
      allvisits.estvhelio = allvisits.vhelio
      allvisits.vtype = 0
      allvisits.vrel = !values.f_nan
      allvisits.vrelerr = !values.f_nan
      allvisits.vhelio = !values.f_nan
    endif
    if not tag_exist(allvisits,'SYNTHVREL') then begin
      ADD_TAG,allvisits,'SYNTHVREL',0.0,allvisits
      ADD_TAG,allvisits,'SYNTHVRELERR',0.0,allvisits
      ADD_TAG,allvisits,'SYNTHVHELIO',0.0,allvisits
    endif
    if (rvrefine EQ 2) and (not tag_exist(allvisits,'OBSVREL')) then begin
      ADD_TAG,allvisits,'OBSVREL',0.0,allvisits
      ADD_TAG,allvisits,'OBSVRELERR',0.0,allvisits
      ADD_TAG,allvisits,'OBSVHELIO',0.0,allvisits
    endif
  endif

  ; HTML file
  printf,html,'<HTML><BODY>'
  if telescope eq 'apo1m' then f=allvisits[0].plate else begin
    lplate=0L
    reads,allvisits[0].plate,lplate
    f=apogee_field(locid,lplate)
  endelse
  printf,html,'<H2> Field: ',f,' Location ID: ', clocation,'</H2><p>'
  printf,html,'<TABLE BORDER=2>'

  ; find the uniq object names (or ra if we don't have them)
  ;good=where(allvisits.plate gt 0 and strpos(allvisits.file,'Visit') ge 0 and $
  ;good=where(strpos(allvisits.file,'Visit') ge 0 and $
  ;           (allvisits.starflag and badstarflag()) eq 0)
  good=where(strpos(allvisits.file,'Visit') ge 0,ngood)
  if ngood eq 0 then stop,'halt: no good visits?'
  objuniq=uniq(strtrim(allvisits[good].apogee_id,2),sort(strtrim(allvisits[good].apogee_id,2)))
  ; for each unique object, do the combination 
  openw,csh,stars_dir+clocation+'/a.csh',/get_lun

  if keyword_set(inpobj) then begin
    i=where(strtrim(allvisits[good[objuniq]].apogee_id,2) eq inpobj)
    i1=i[0]
    i2=i1
  endif else begin
    i1=0
    i2=n_elements(objuniq)-1
    if keyword_set(nstars) then i2=nstars-1
  endelse

  for i=i1,i2 do begin
   objname=strtrim(allvisits[good[objuniq[i]]].apogee_id,2)

   ; allobj has allvisit indices of all visits of this object without regard to S/N and commissioning
   objgood=where(strtrim(allvisits[good].apogee_id,2) eq objname ,ngood)
   allobj=good[objgood]
   print,'OBJECT: ',objname,'   HMAG: ',allvisits[allobj[0]].h

   ; recalculate barycentric correction with more accurate calculation
   if telescope eq 'lco25m' then obs='LCO' else obs='APO'
   openw,fp,stars_dir+clocation+'/'+objname+'.bcin',/get_lun
   for j=0,n_elements(allobj)-1 do $
     printf,fp,string(format='(2f14.8,f18.8,1x,a)',allvisits[allobj[j]].ra,allvisits[allobj[j]].dec,allvisits[allobj[j]].jd,obs)
   free_lun,fp
   cmd=['bc',stars_dir+clocation+'/'+objname+'.bcin','--out',stars_dir+clocation+'/'+objname+'.bc']
   print,cmd
   spawn,cmd,/noshell
   readcol,stars_dir+clocation+'/'+objname+'.bc',bc
   for j=0,n_elements(allobj)-1 do begin
     oldbc=allvisits[allobj[j]].bc
     if oldbc ne 0. then begin
       allvisits[allobj[j]].bc=bc[j]
       allvisits[allobj[j]].vhelio+=(bc[j]-oldbc)
       if abs(bc[j]-oldbc) gt 0.1 then stop,'halt: recalculated BC off by more than 100 m/s'
     endif
   endfor
   file_delete,stars_dir+clocation+'/'+objname+'.bcin',/allow_non
   file_delete,stars_dir+clocation+'/'+objname+'.bc',/allow_non

   for isurvey=0,1 do begin
    ; find all of the matching objects, separately for commissioning
    if isurvey eq 0 then begin
      objgood=where(allvisits[good].mjd lt 55800 and $
                    strtrim(allvisits[good].apogee_id,2) eq objname  and allvisits[good].snr gt snmin,ngood)
      root='apStarC-'+apred_vers+'-'
      commiss=1
      if keyword_set(sinc) then if n_elements(nres) ne 3 then sinc=[4.,2.,2.] else sinc=nres
    endif else begin
      objgood=where(allvisits[good].mjd gt 55800 and $
                    strtrim(allvisits[good].apogee_id,2) eq objname and allvisits[good].snr gt snmin,ngood)
      root='apStar-'+apred_vers+'-'
      commiss=0
      if keyword_set(sinc) then if n_elements(nres) ne 3 then sinc=[5.,4.25,3.5] else sinc=nres
    endelse
    if ngood eq 0 then goto, bomb

    ; obj has allvisit indices of all visits of this object that exceed snmin and are/are not commissioning
    obj=good[objgood]
    
    IF keyword_set(snsig) and ngood gt 1 THEN BEGIN
      newsnmin = median([allvisits[obj].snr],/even) - snsig*mad([allvisits[obj].snr])
      IF (newsnmin LT snmin) OR (n_elements(obj) LT 2) THEN newsnmin = snmin
      IF newsnmin GT snmax THEN newsnmin = snmax; IF lower sigma clip is above snmax, then set isnmin to that value
    ENDIF ELSE newsnmin = snmin
    snrgood = where(allvisits[obj].snr GE newsnmin,ngood)
    if ngood eq 0 then goto, bomb
    obj = obj[snrgood]
    
    ; check to see if we have already done this stars with all of the visits
    ;starfile=stars_dir+clocation+'/'+root+objname+'.fits'
    starfile=apogee_filename('Star',field=clocation,obj=objname)
    done=0
    if file_test(starfile) then begin
      done=1
      star=mrdfits(starfile,0,header,/silent)
      nvisits=sxpar(header,'NVISITS')
      if nvisits ne n_elements(obj) then done=0 else begin
       for ivisit=1,nvisits do begin
        card='SFILE'+strtrim(ivisit,2)
        sfile=sxpar(header,card)
        if sfile ne file_basename(allvisits[obj[ivisit-1]].file) then done=0
       endfor
      endelse
    endif

    ; if we haven't done it, do it!
    if keyword_set(clobber) or not done then begin
      ; load up the spectra into str structure  
      nvisits=n_elements(obj)
      for ivisit=0,nvisits-1 do begin
        ;dir=platedir(allvisits[obj[ivisit]].plate,allvisits[obj[ivisit]].mjd,vers=in_vers,platetype=platetype)
        ;aploadvisit,dir+allvisits[obj[ivisit]].file,str
        v=allvisits[obj[ivisit]]
        file=apogee_filename('Visit',plate=v.plate,mjd=v.mjd,fiber=v.fiberid,reduction=v.apogee_id)
        aploadvisit,file,str
        if ivisit eq 0 then allstr=ptr_new(str) else allstr=[allstr,ptr_new(str)]
      endfor
      ; do the resamling and visit combination using the apVisit relative velocities
      ;if keyword_set(sinc) then sinc=nres

      ; do relative velocities by cross-correlating spectra against each other
      if keyword_set(rvrefine) then begin
        if not tag_exist(allvisits,'SYNTHVREL') then begin
          ADD_TAG,allvisits,'SYNTHVREL',-999999.,allvisits
          ADD_TAG,allvisits,'SYNTHVRELERR',-999999.,allvisits
          ADD_TAG,allvisits,'SYNTHVHELIO',-999999.,allvisits
        endif
        if not tag_exist(allvisits,'OBSVREL') then begin
          ADD_TAG,allvisits,'OBSVREL',-999999.,allvisits
          ADD_TAG,allvisits,'OBSVRELERR',-999999.,allvisits
          ADD_TAG,allvisits,'OBSVHELIO',-999999.,allvisits
        endif
        ; IF rvrefine keyword is set to 2, then insert synthetic templates into the iterations with the observed templates. 
        IF rvrefine EQ 2 THEN BEGIN
           ;Right now this is a package deal, and requires no new keywords in the plan files, but these could be seperate keywords
           synthiter = 1
           refinefixmasked = 2 ; Set whether or not to use fixmasked keyword in aprvprep in aprvrefine. Not having it set worked well with the synth option. 
           allsmooth = 0
        ENDIF ;If rvrefine is set but not specifically to 2, then run aprvrefine using default values from DR13.
        vtemp = allvisits[obj]
        aprvrefine,allstr,vtemp,starstr,sinc=sinc,log=log,plotdir=plot_dir,bccomb=bccomb, synth = synthiter, fixmasked = refinefixmasked, trimgrid=trimgrid, allsmooth = allsmooth
        ; check for failure
        gdrv = where(vtemp.vtype gt 0,ngdrv,comp=bdrv,ncomp=nbdrv)
        if ngdrv eq 0 then goto, bomb
        allvisits[obj] = vtemp  ; update allvisits structure - MIGHT NEED TO UPDATE allvisits to have obsvhelio
      endif else begin
        apvisitcomb,allstr,allvisits[obj],starstr,sinc=sinc,log=log  ;/nolsffit
      endelse

      ; run binary check statistics 
      for is=0,1 do begin
        apg_binarycheck,allvisits[obj],binary,stablerv_chi2,stablerv_rchi2,$
               chi2_threshold,stablerv_chi2_prob,synth=is
        if n_elements(stablerv_chi2) eq 0 then stablerv_chi2 = -1.;
        if n_elements(stablerv_rchi2) eq 0 then stablerv_rchi2 = -1.;
        if n_elements(chi2_threshold) eq 0 then chi2_threshold = -1.;
        if n_elements(stablerv_chi2_prob) eq 0 then stablerv_chi2_prob = -1.;
        bstr = { binary: binary, stablerv_chi2: stablerv_chi2, stablerv_rchi2: stablerv_rchi2, $
                 chi2_threshold: chi2_threshold,  stablerv_chi2_prob: stablerv_chi2_prob}
        if is eq 0 then binstr=bstr else binstr=[binstr,bstr]
      endfor
      ; Output apStar file and make plots
      print,'Writing output apStar and apLSF files'
      ;apstar_output,starstr,allvisits[obj],binstr,stars_dir,commiss=commiss,localdir=getlocaldir(),apstar_vers=apstar_vers,locationdir=clocation  ;/nolsf
      apstar_output,starstr,allvisits[obj],binstr,stars_dir,commiss=commiss,localdir=getlocaldir(),apstar_vers=apred_vers,locationdir=clocation,survey=survey  ;/nolsf

      ; link and plot conversion script
      finaldir=stars_dir+clocation+'/' 
      if size(getlocaldir(),/type) eq 7 then $
        outdir=getlocaldir()+'/'+clocation+'/' else $
        outdir=stars_dir+clocation+'/' 
      
      MJD = max(allvisits[obj].jd - 2400000.5 )  ; last visit
      mjd5 = long(mjd)  ; clip the decimals
      vhead1 = reform(starstr.header[0,*])
      objid = sxpar(vhead1,'OBJID')
      printf,csh,'convert '+outdir+'plots/'+file_basename(starfile,'.fits')+'.eps '+finaldir+'plots/'+file_basename(starfile,'.fits')+'.jpg'
      printf,csh,'convert '+outdir+'plots/'+file_basename(starfile,'.fits')+'SN.eps '+finaldir+'plots/'+file_basename(starfile,'.fits')+'SN.jpg'
      printf,csh,'"rm" '+outdir+'plots/'+file_basename(starfile,'.fits')+'.eps'
      printf,csh,'"rm" '+outdir+'plots/'+file_basename(starfile,'.fits')+'SN.eps'
    endif

    ; load the apStar file  
    aploadstar,starfile,apstr
    ; replace the initial RVs in visit structure with refined ones (already done if we ran
    ;   this object, but not if we're just reading it back from a previous computation)
    
    for ivisit=0,apstr.nvisits-1 do begin
      num = strtrim(ivisit+1,2)
      ;vrad=sxpar(apstr.head,'VRAD'+num)
      ;vhelio=sxpar(apstr.head,'VHELIO'+num)
      ;verr=sxpar(apstr.head,'VERR'+num)
      ;teff=sxpar(apstr.head,'RVTEFF'+num)
      ;logg=sxpar(apstr.head,'RVLOGG'+num)
      ;feh=sxpar(apstr.head,'RVFEH'+num)
      ;alpha=sxpar(apstr.head,'RVALPH'+num)
      ;carbon=sxpar(apstr.head,'RVCARB'+num)
      vdiff=apstr.rv.vhelio[ivisit]-allvisits[obj[ivisit]].vhelio
      allvisits[obj[ivisit]].vhelio=apstr.rv.vhelio[ivisit]
      allvisits[obj[ivisit]].vrel=apstr.rv.vrel[ivisit]
      allvisits[obj[ivisit]].vrelerr=apstr.rv.vrelerr[ivisit]
      allvisits[obj[ivisit]].synthvrel=apstr.rv.synthvrel[ivisit]
      allvisits[obj[ivisit]].synthvrelerr=apstr.rv.synthvrelerr[ivisit]
      allvisits[obj[ivisit]].synthvhelio=apstr.rv.synthvhelio[ivisit]
      if tag_exist(allvisits,'OBSVREL') then begin
        allvisits[obj[ivisit]].obsvrel=apstr.rv.obsvrel[ivisit]
        allvisits[obj[ivisit]].obsvrelerr=apstr.rv.obsvrelerr[ivisit]
        allvisits[obj[ivisit]].obsvhelio=apstr.rv.obsvhelio[ivisit]
      endif   

      if tag_exist(apstr.rv,'CHISQ') then allvisits[obj[ivisit]].CHISQ = apstr.rv.CHISQ[ivisit]
      allvisits[obj[ivisit]].RV_TEFF = apstr.rv.TEFF
      allvisits[obj[ivisit]].RV_LOGG = apstr.rv.LOGG
      allvisits[obj[ivisit]].RV_FEH = apstr.rv.FEH
      allvisits[obj[ivisit]].RV_ALPHA = apstr.rv.ALPHA
      allvisits[obj[ivisit]].RV_CARB = apstr.rv.CARBON
      if tag_exist(apstr,'VLSR') then begin 
       vlsr = apstr.vlsr
       vgsr = apstr.vgsr
      endif else begin
       VCONV,apstr.rv.vhelio[ivisit],apstr.glon,apstr.glat,vlsr,vgsr
      endelse
      allvisits[obj[ivisit]].VLSR = vlsr
      allvisits[obj[ivisit]].VGSR = vgsr
    endfor
    ; load up information for an output table for apField
    junk=where(allvisits[obj].mjd lt 55800,n1)
    junk=where(allvisits[allobj].mjd gt 55800 and allvisits[allobj].mjd lt 56860,n2)
    junk=where(allvisits[allobj].mjd gt 56860,n3)

    ; does this star appear from  multiple surveys?
    isurveys = uniq(allvisits[obj].survey,sort(allvisits[obj].survey))
    nsurveys = n_elements(isurveys)
    survey=allvisits[obj[isurveys[0]]].survey
    for jsurvey=1,nsurveys-1 do survey=survey+','+allvisits[obj[isurveys[jsurvey]]].survey
    ; use location ID from first visit
    locid=allvisits[obj[0]].location_id
    ; separate targetflags into apogee and apogee2
    apogee_target1=0L
    apogee_target2=0L
    apogee_target3=0L
    apogee_targflags=''
    j=where(allvisits[obj].survey eq 'apogee',nj)
    if nj gt 0 then begin
      apogee_target1=allvisits[obj[j[0]]].apogee_target1
      apogee_target2=allvisits[obj[j[0]]].apogee_target2
      apogee_target3=allvisits[obj[j[0]]].apogee_target3
      apogee_targflags=targflag(apogee_target1,apogee_target2,apogee_target3,survey='apogee')
    endif 
    apogee2_target1=0L
    apogee2_target2=0L
    apogee2_target3=0L
    apogee2_targflags=''
    j=where(strpos(allvisits[obj].survey,'apogee2') ge 0,nj)
    if nj gt 0 then begin
      apogee2_target1=allvisits[obj[j[0]]].apogee_target1
      apogee2_target2=allvisits[obj[j[0]]].apogee_target2
      apogee2_target3=allvisits[obj[j[0]]].apogee_target3
      apogee2_targflags=targflag(apogee2_target1,apogee2_target2,apogee2_target3,survey='apogee2')
    endif
    ; combine targflags
    if apogee_targflags ne '' then begin
      if apogee2_targflags ne '' then targflags=apogee_targflags+','+apogee2_targflags else targflags=apogee_targflags
    endif else targflags=apogee2_targflags
    a={file: file_basename(starfile), apogee_id: objname, telescope: apstr.telescope, $
       location_id: long(locid), field: apstr.field, $
       j: apstr.j, j_err: apstr.j_err, h: apstr.h, h_err: apstr.h_err, k: apstr.k, k_err: apstr.k_err, $
       ra: apstr.ra, dec: apstr.dec, glon: apstr.glon, glat: apstr.glat, $
       ak_targ: float(apstr.ak_targ), ak_targ_method: apstr.ak_targ_method, $
       ak_wise: float(apstr.ak_wise), sfd_ebv: float(apstr.sfd_ebv),$
       apogee_target1: apogee_target1, apogee_target2: apogee_target2, apogee_target3: apogee_target3, $
       apogee2_target1: apogee2_target1, apogee2_target2: apogee2_target2, apogee2_target3: apogee2_target3, $
       targflags: targflags, survey: survey, programname: allvisits[obj[j[0]]].programname, $
       ninst: [n1,n2,n3],$
       nvisits: apstr.nvisits, combtype:apstr.combtype, commiss: commiss, $
       snr: float(apstr.snr), $
       starflag: apstr.starflag, starflags: starflag(apstr.starflag), $
       andflag: apstr.andflag, andflags: starflag(apstr.andflag), $
       vhelio_avg: float(apstr.vhelio), vscatter: float(apstr.vscatter),$
       verr: float(apstr.verr), verr_med: float(apstr.verr_med),$
       obsvhelio_avg: float(apstr.obsvhelio),  obsvscatter: float(apstr.obsvscatter),$
       obsverr: float(apstr.obsverr), obsverr_med: float(apstr.obsverr_med),$
       synthvhelio_avg: float(apstr.synthvhelio),  synthvscatter: float(apstr.synthvscatter),$
       synthverr: float(apstr.synthverr), synthverr_med: float(apstr.synthverr_med),$
       rv_teff: float(apstr.rv_teff), rv_logg: float(apstr.rv_logg), rv_feh: float(apstr.rv_feh),$
       rv_alpha: float(apstr.rv_alpha), rv_carb: float(apstr.rv_carb), $
       rv_ccfwhm: float(apstr.ccpfwhm), rv_autofwhm: float(apstr.autofwhm), synthscatter: float(apstr.synthscatter),$
       ;binary: apstr.rv.binary, $
       stablerv_chi2: apstr.rv.stablerv_chi2, $
       stablerv_rchi2: apstr.rv.stablerv_rchi2, chi2_threshold: apstr.rv.chi2_threshold, $
       stablerv_chi2_prob: apstr.rv.stablerv_chi2_prob}
    b={spec: apstr.spec[*,0], err: apstr.err[*,0], mask: apstr.mask[*,0]}
    c={wave: apstr.wavelength[*,0]}

    ; concatenate results for this star with previous stars
    comb=apstr.spec[*,0]
    if commiss eq 1 then begin
      if ncommiss eq 0 then begin
        allp=a
        alldatap=b 
      endif else begin
        allp=[allp,a]
        alldatap=[alldatap,b]
      endelse
      ncommiss+=1
    endif
    if isurvey eq 1 then begin
      if nsurvey eq 0 then begin
        all=a
        alldata=b 
      endif else begin
        all=[all,a]
        alldata=[alldata,b]
      endelse
      nsurvey+=1
    endif

    ; RV plot
    set_plot,'PS'
    mjd=allvisits[obj].mjd+0L
    vhelio=allvisits[obj].vhelio
    outfile=stars_dir+clocation+'/plots/apRV-'+objname
    if size(getlocaldir(),/type) eq 7 then $
    plotfile=getlocaldir()+'apRV-'+objname else $
    plotfile=stars_dir+clocation+'/plots/apRV-'+objname
    device,file=plotfile+'.eps',/encap,xsize=12,ysize=12
    ymin=min([min(vhelio)-1,apstr.vhelio-1])
    ymax=max([max(vhelio)+1,apstr.vhelio+1])
    xmin=min(mjd)-5
    xmax=max(mjd)+5
    ; plot both commissioning and survey data
    if nvisits gt 1 then plot,allvisits[allobj].mjd+0L,allvisits[allobj].vhelio,psym=6,xtitle='MJD',ytitle='V(helio)',yrange=[ymin,ymax],xrange=[xmin,xmax],xtickformat='(i5)' else $
      plot,[allvisits[allobj].mjd+0L],[allvisits[allobj].vhelio],psym=6,xtitle='MJD',ytitle='V(helio)',yrange=[ymin,ymax],xrange=[xmin,xmax],xtickformat='(i5)'
    oplot,!x.crange,[apstr.vhelio,apstr.vhelio],linestyle=1
    device,/close
    set_plot,'X'
    printf,csh,'convert '+plotfile+'.eps '+outfile+'.jpg'
    printf,csh,'"rm" '+plotfile+'.eps'

    ; HTML output
    if isurvey eq 0 then commiss='<FONT COLOR=RED> COMMISSIONING </FONT> ' else commiss=' '
    if is_bit_set(apstr.apogee_target2,9) then $
      printf,html,'<TR><TD bgcolor=lightblue>' else  printf,html,'<TR><TD>'
    printf,html,commiss,'<BR>'

    printf,html,'<A HREF=../'+file_basename(starfile)+'>',objname,'</a>'
    printf,html,'(<A HREF=../plots/'+objname+'_rvccf.gif>RV template + CCFs</a>)'
    rastring=stringize(apstr.ra,ndec=5)
    decstring=stringize(apstr.dec,ndec=5)
    printf,html,'(<A HREF="http://simbad.cfa.harvard.edu/simbad/sim-basic?Ident='+RASTRING+'+%09'+DECSTRING+'++&submit=SIMBAD+search"> SIMBAD </A>)'

    printf,html,'<br> H='+string(format='(f8.2)',apstr.h)
    printf,html,'<br>'+targflag(apstr.apogee_target1,apstr.apogee_target2,apstr.apogee_target3,survey=survey)+'<BR>'+starflag(apstr.starflag)
    printf,html,'<br> SNR='+string(format='(f8.2)',apstr.snr)
    printf,html,'<br> VSCATTER='+string(format='(f8.2)',apstr.vscatter)
    printf,html,'<br> RV_TEFF='+string(format='(f8.2)',apstr.rv_teff),$
                    ' RV_LOGG='+string(format='(f8.2)',apstr.rv_logg),$
                    ' RV_FEH='+string(format='(f8.2)',apstr.rv_feh)
    printf,html,' <TABLE BORDER=2>'
    printf,html,'<TR><TD>Plate<TD>MJD<TD>VHELIO<TD>S/N'
    for ivisit=0,n_elements(mjd)-1 do begin
      printf,html,'<TR><TD>',apstr.rv.plate[ivisit],'<TD>',apstr.rv.mjd[ivisit],'<TD>',apstr.rv.vhelio[ivisit],'<TD>',allvisits[obj[ivisit]].snr
    endfor
    printf,html,'</TABLE>'
    printf,html,'<TD> <IMG SRC=../plots/apRV-'+objname+'.jpg>'
    printf,html,'<TD><IMG SRC=../plots/'+file_basename(starfile,'.fits')+'.jpg>'
    printf,html,'<IMG SRC=../plots/'+file_basename(starfile,'.fits')+'SN.jpg>'
 
    bomb: 
   endfor
  endfor

  if keyword_set(inpobj) then return

  ; write out the apField file(s)
  mkhdr,hdr,0
  sxaddpar,hdr,'V_APSTAR',getvers()
  if nsurvey gt 0 then begin
    outfile=apogee_filename('Field',field=clocation)
    mwrfits,0,outfile,hdr,/create
    mwrfits,all,outfile
    mwrfits,alldata,outfile
    mwrfits,c,outfile
  endif
  if ncommiss gt 0 then begin
    outfile=stars_dir+clocation+'/apFieldC-'+clocation
    mwrfits,0,outfile+'.fits',hdr,/create
    mwrfits,allp,outfile+'.fits'
    mwrfits,alldatap,outfile+'.fits'
    mwrfits,c,outfile+'.fits'
  endif

  ; write out the apFieldVisits file, which has refined RVs
  if ncommiss gt 0 then begin
    j=where(allvisits.mjd lt 55800)
    for jj=0,n_elements(j)-1 do begin
      i=j[jj]
      ;dir=platedir(allvisits[i].plate,allvisits[i].mjd,platetype=platetype)
      ;aploadvisit,dir+allvisits[i].file,str
      v=allvisits[i]
      file=apogee_filename('Visit',plate=v.plate,mjd=v.mjd,fiber=v.fiberid,reduction=v.apogee_id)
      aploadvisit,file,str
      if jj eq 0 then begin
        tmp={spec: str.spec*0., err: str.err*0., mask: str.mask*0, wave: str.wave*0.}
        allvisitspec=replicate(tmp,n_elements(j))
      endif
      allvisitspec[jj].spec=str.spec
      allvisitspec[jj].err=str.err
      allvisitspec[jj].mask=str.mask
    endfor
    outfile=stars_dir+clocation+'/apFieldVisitsC-'+clocation
    mwrfits,0,outfile+'.fits',hdr,/create
    mwrfits,allvisits[j],outfile+'.fits'
    mwrfits,allvisitspec,outfile+'.fits'
  endif
  if nsurvey gt 0 then begin
    j=where(allvisits.mjd gt 55800)
    jjj=0
    for jj=0,n_elements(j)-1 do begin
      i=j[jj]
      ;dir=platedir(allvisits[i].plate,allvisits[i].mjd,platetype=platetype)
      ;aploadvisit,dir+allvisits[i].file,str
      v=allvisits[i]
      file=apogee_filename('Visit',plate=v.plate,mjd=v.mjd,fiber=v.fiberid,reduction=v.apogee_id)
      aploadvisit,file,str
      sz=size(str.spec)
      ; only take observations with successful dithers, i.e. 4096 pixels
      if sz[1] eq 4096 then begin
       if jjj eq 0 then begin
        tmp={spec: str.spec*0., err: str.err*0., mask: str.mask*0, wave: str.wave*0.}
        allvisitspec=replicate(tmp,n_elements(j))
        jjj+=1
       endif
       allvisitspec[jj].spec=str.spec
       allvisitspec[jj].err=str.err
       allvisitspec[jj].mask=str.mask
      endif else begin
       print,'not halted, but 2048 spectrum ignored: '+allvisits[i].file
     endelse
    endfor
    outfile=apogee_filename('FieldVisits',field=clocation)
    mwrfits,0,outfile,hdr,/create
    mwrfits,allvisits[j],outfile
    mwrfits,allvisitspec,outfile
  endif

  free_lun,csh
  if not keyword_set(noplot) then spawn,'csh '+stars_dir+clocation+'/a.csh'
  file_delete,stars_dir+clocation+'/a.csh',/allow_non

  printf,html,'</TABLE></BODY></HTML>'
  free_lun,html
  file_move,html_dir+clocation+'.html.tmp',$
          html_dir+clocation+'.html',/overwrite
endfor

end
