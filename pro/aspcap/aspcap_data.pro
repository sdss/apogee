pro aspcap_data,str,files,getobject=getobject

if keyword_set(getobject) then begin
  targetdir=getenv('APOGEE_TARGET')
  ; preload ra, dec, and obj arrays for use with spherematch
  for i=0,n_elements(files)-1 do begin
    a=mrdfits(files[i],0,hdr,/silent)
    if i eq 0 then begin
      ra=sxpar(hdr,'RA')
      dec=sxpar(hdr,'DEC')
      obj=sxpar(hdr,'OBJID')
    endif else begin
      ra=[ra,sxpar(hdr,'RA')]
      dec=[dec,sxpar(hdr,'DEC')]
      obj=[obj,sxpar(hdr,'OBJID')]
    endelse
  endfor
endif

; routine to grab information from header of apStar file and 
;   incorporate into output ASPCAP structure

oldfield=''
for i=0,n_elements(str.param)-1 do begin
  ; find and read the corresponding apStar file
  jpos=strpos(str.param[i].apogee_id,'_v')
  visit=0
  if jpos gt 0 then begin
    reads,strmid(str.param[i].apogee_id,jpos+2,jpos+4),visit
    fitsname = strmid(str.param[i].apogee_id,0,jpos) 
  endif else fitsname=str.param[i].apogee_id

  jj=where(strpos(files,fitsname) ge 0)
  j=jj[0]
  a=mrdfits(files[j],0,hdr,/silent)
  ; calculate revised SNR based on green wavelength range
  spec=mrdfits(files[j],1,whdr,/silent)
  err=mrdfits(files[j],2,/silent)
  wave=wave_fits(whdr,/log)
  gd=where(wave gt 15850 and wave lt 16440)
  snrev=median(spec[gd,0]/err[gd,0])

  ; add catalog info if we have it
  if keyword_set(getobject) then begin
    ;field=apogee_fixfield(sxpar(hdr,'FIELD'))
    locid=sxpar(hdr,'LOCID')
    survey=sxpar(hdr,'SURVEY')
    ;if strtrim(field,2) eq 'Unknown' and locid ne 0 then field=apogee_field(locid,0,/addloc)
    field=apogee_field(locid,0,/addloc)
    if strtrim(field,2) ne strtrim(oldfield,2) then begin
      ; if necessary, read new apogeeObject field and match 
      ;objectfile=targetdir+'/apogeeObject/apogeeObject_'+strtrim(field,2)+'.fits'
      if strpos(survey,'apogee2') ge 0 then apogeeobject='apogee2Object' $
      else if strpos(survey,'apo1m') ge 0 then apogeeobject='apogee1mObject' else $
      apogeeobject='apogeeObject'
      objectfile=targetdir+'/'+apogeeobject+'/'+apogeeobject+'_'+strtrim(field,2)+'.fits'

      if not file_test(objectfile) then begin
        print,'cant find apogeeObject file: ', objectfile
      endif else begin
        print,'reading apogeeObject file: ', objectfile
        objects=mrdfits(objectfile,1)
        bd=where(objects.ra lt 0,nbd)
        if nbd gt 0 then objects[bd].ra = 0.
        if nbd gt 0 then objects[bd].dec = 0.
        bd=where(ra lt 0,nbd)
        if nbd gt 0 then ra[bd] = 0.
        if nbd gt 0 then dec[bd] = 0.
        spherematch,objects.ra,objects.dec,ra,dec,10./3600.,match1,match2,dist,maxmatch=1
        oldfield=field
      endelse
    endif

    tmp=strsplit(str.param[i].apogee_id,'-')
    ntmp=n_elements(tmp)
    name=file_basename(strmid(str.param[i].apogee_id,tmp[ntmp-1]),'.fits')
    ;iobj=where(strtrim(objects.apogee_id,2) eq name,nobj)
    if file_test(objectfile) then begin
     jobj=where(match2 eq j,nobj)
     if nobj gt 0 then begin
      iobj=match1(jobj)
      sxaddpar,hdr,'J',objects[iobj].j
      sxaddpar,hdr,'H',objects[iobj].h
      sxaddpar,hdr,'K',objects[iobj].k
      sxaddpar,hdr,'AKTARG',objects[iobj].ak_targ
      sxaddpar,hdr,'AKMETHOD',objects[iobj].ak_targ_method
      sxaddpar,hdr,'AKWISE',objects[iobj].ak_wise
      sxaddpar,hdr,'SFD_EBV',objects[iobj].sfd_ebv
      if dist(jobj) gt 2./3600. then $
        print,'Large distance match for ',name,': ', dist(jobj)*3600.
     endif else begin
      print,'not halted'
      print,'Missing object information for ',name,' ', objectfile
     endelse
    endif
  endif

  field=apogee_fixfield(sxpar(hdr,'FIELD'))
  objid=sxpar(hdr,'OBJID')
  if size(objid,/type) eq 3 then objid=str.param[i].apogee_id
  loc=sxpar(hdr,'LOC')
  if loc eq 0 then loc=sxpar(hdr,'LOCATION')
  if loc eq 0 then loc=sxpar(hdr,'LOCID')
  jmag=sxpar(hdr,'J')
  if jmag eq 0 then jmag=sxpar(hdr,'JMAG')
  hmag=sxpar(hdr,'H')
  if hmag eq 0 then hmag=sxpar(hdr,'HMAG')
  kmag=sxpar(hdr,'K')
  if Kmag eq 0 then Kmag=sxpar(hdr,'KMAG')
  targ1=sxpar(hdr,'TARG1')
  if targ1 eq 0 then targ1=sxpar(hdr,'PRIMTARG')
  targ2=sxpar(hdr,'TARG2')
  if targ2 eq 0 then targ2=sxpar(hdr,'SECTARG')
  targ3=sxpar(hdr,'TARG3')
  survey=sxpar(hdr,'SURVEY',count=count)
  if count le 0 then survey='apogee'
  ak_targ_method= sxpar(hdr,'AKMETHOD')
  if size(ak_targ_method,/type) eq 3 then ak_targ_method='Unknown '

  nvisits=sxpar(hdr,'NVISITS')
  sfile=''
  if nvisits gt 0 then begin
    for ivisit=0,nvisits-1 do begin
      num=strtrim(ivisit+1,2)
      card=sxpar(hdr,'SFILE'+num,count=count)
      if count gt 0 then sfile=sfile+sxpar(hdr,'SFILE'+num)+' '
    endfor
  endif

  ; extract a reduction version (not optimal given multiple upstream products!)
  if strpos(files[j],'apStar') gt 0 then begin
    hist='HISTORY APSTAR:  APOGEE Reduction Pipeline Version: '
  endif else begin
    hist='HISTORY AP1DVISIT:  APOGEE Reduction Pipeline Version: '
  endelse
  ppos=strpos(hdr,hist) 
  pos=where(ppos eq 0, ncts)
  if ncts gt 0 then begin
    red_vers=hdr[pos] & pos=strpos(red_vers,'Version:')
    red_vers=strcompress(strmid(red_vers,pos+8))
  endif else red_vers = 'Unknown'

  if visit eq 0 then begin
    snr = float(sxpar(hdr,'SNR'))
    vrad = float(sxpar(hdr,'VRAD'))
    vhelio = float(sxpar(hdr,'VHELIO'))
    vscatter = float(sxpar(hdr,'VSCATTER'))
    verr =  float(sxpar(hdr,'VERR'))
    nv = nvisits
  endif else begin
    num=strtrim(visit,2)
    snr = sxpar(hdr,'SNRVIS'+num)
    vrad = sxpar(hdr,'VRAD'+num)
    vhelio = sxpar(hdr,'VHELIO'+num)
    vscatter = 0.
    verr =  sxpar(hdr,'VERR'+num)
    nv = 1
  endelse

  objparam={$
	    file: str.param[i].apogee_id,$
	    red_vers: strtrim(red_vers,2),$
            field: field,$
            location_id: loc,$
	    apogee_id: objid,$
            ra: double(sxpar(hdr,'RA')),$
            dec: double(sxpar(hdr,'DEC')),$
            glon: double(sxpar(hdr,'GLON')),$
            glat: double(sxpar(hdr,'GLAT')),$
            j: float(jmag),$
            h: float(hmag),$
            k: float(kmag),$
            ak_targ: float(sxpar(hdr,'AKTARG')),$
            ak_targ_method: ak_targ_method,$
            ak_wise: float(sxpar(hdr,'AKWISE')),$
            sfd_ebv: float(sxpar(hdr,'SFD_EBV')),$
            apogee_target1: long(targ1),$
            apogee_target2: long(targ2),$
            apogee_target3: long(targ3),$
            targflags: targflag(targ1,targ2,targ3,survey=survey),$
            survey: survey,$
            starflag: long(sxpar(hdr,'STARFLAG')),$
            andflag: long(sxpar(hdr,'ANDFLAG')),$
            snr: float(snr),$
            snrev: float(snrev),$
            vrad: float(vrad),$
            vhelio: float(vhelio),$
            vscatter: float(vscatter),$
            verr: float(verr),$
            ;vlsr: float(sxpar(hdr,'VLSR')),$
            ;vgsr: float(sxpar(hdr,'VGSR')),$
            nvisits: fix(nv),$
            visit: fix(visit),$
            rv_ccfwhm: float(sxpar(hdr,'CCPFWHM')),$
            rv_autofwhm: float(sxpar(hdr,'AUTOFWHM')),$
            meanfib: float(sxpar(hdr,'MEANFIB')),$
            sigfib: float(sxpar(hdr,'SIGFIB')),$
            visitfiles: sfile}

  ; add in the ASPCAP tags
  tags=tag_names(str.param)
  for itag=0,n_elements(tags)-1 do begin
    if tags[itag] ne 'OBJ' and tags[itag] ne 'APOGEE_ID' then $
    add_tag,objparam,tags[itag],str.param[i].(itag),objparam
  endfor

  ; create or add to the final parameter structure
  if i eq 0 then new=objparam else new=[new,objparam]
endfor
; create and return the modified final structure
str={param: new, spec: str.spec, lib: str.lib}
end
