PRO mk1mtarg,program,addtarg=addtarg

;what values are in all 1m target files?
;name  ;ra  ;dec  ;h   ;v   ;comments   ;type    ;priority
;plx   ;eplx   ;b   ;j   ;k   ;HIP   ;HD   ;HR   ;PLXFLAG
;APOGEE_ID    ;alt_name    ;period(for rrlyr)   ;PMRA    ;PMDEC

 a=mrdfits(program+'.fits',1,/SILENT)
 tn = tag_names(a)
;hh=mrdfits('/home/regulus/feuilldk/1msurvey/XHIP_Anderson12.fits',1,/SILENT)
hh=mrdfits('/regulus/feuilldk/1msurvey/new_hipparcos_cat.fits',1,/SILENT)

;------------------------------------------------------
;ADD NEW TARGETS AND FILL ARAYS WITH KNOWN PARAMETERS
;------------------------------------------------------
IF keyword_set(addtarg) THEN BEGIN
   print,'ADDING NEW TARGETS'
   IF strtrim(program,2) eq 'rrlyr' THEN $
      readcol,program+'.add',id2,ra2,dec2,h2,v2,type2,priority,period2,comment2,f='a,f,f,f,f,a,i,f,a' $
   ELSE readcol,program+'.add',id2,ra2,dec2,h2,v2,comment2,priority,f='a,f,f,f,f,a,i'
   targind = fltarr(n_elements(id2))
   FOR i=0,n_elements(targind)-1 DO targind[i]=where(strtrim(a.name,2) eq strtrim(id2[i],2))
   newind = where(targind lt 0)
   help,newind
   name=[a.name,id2[newind]]
   ra=[a.ra,ra2[newind]]
   dec=[a.dec,dec2[newind]]
   h=[a.h,h2[newind]]
   v=[a.v,v2[newind]]
   comment=[a.comment,comment2[newind]]
   IF strtrim(program,2) eq 'rrlyr' THEN period=[a.period,period2[newind]]
   IF where(tn eq 'TYPE') GT 0 THEN type=[a.type,replicate('none',n_elements(newind))] $
   ELSE type=[replicate(program,n_elements(a.name)),type2[newind]] 
   IF where(tn eq 'PRIORITY') GT 0 THEN priority=[a.priority,priority[newind]] $
   ELSE priority=[replicate(1,n_elements(a.name)),priority[newind]]

ENDIF ELSE BEGIN
   name=a.name
   ra=a.ra
   dec=a.dec
   h=a.h
   IF where(tn eq 'V') GT 0 THEN v=a.v $
   ELSE v=replicate(!VALUES.F_NAN,n_elements(name))
   IF where(tn eq 'COMMENT') GT 0 THEN comment=a.comment $
   ELSE comment=replicate(program,n_elements(name))
   IF where(tn eq 'PRIORITY') GT 0 THEN priority=a.priority $
   ELSE priority=replicate(1,n_elements(a.name))
   IF where(tn eq 'TYPE') GT 0 THEN type=a.type $
   ELSE type=replicate(program,n_elements(a.name))
   IF strtrim(program,2) eq 'rrlyr' THEN period=a.period
ENDELSE


;STOP
;------------------------------------------------------
;MATCH TARGETS TO 2MASS AND FILL VALUES
;------------------------------------------------------

print,'MATCHING TO 2MASS'
APOGEE_ID=replicate('a',n_elements(name))
j=dblarr(n_elements(name))
k=dblarr(n_elements(name))
h_err=dblarr(n_elements(name))
j_err=dblarr(n_elements(name))
k_err=dblarr(n_elements(name))
FOR i=0,n_elements(name)-1 DO BEGIN
   ;print,name[i],' ',h[i]
   coords=[ra[i],dec[i]]
   data=QUERYVIZIER('2MASS-PSC',coords,25/60.,/ALLCOL,/SILENT)
   ;sz=SIZE(data) & if sz[0] eq 0 then data=QUERYVIZIER('2MASS-PSC',coords,25/60.,/ALLCOL,/SILENT,/CANADA)
   ;sz=SIZE(data) & if sz[0] eq 0 then data=QUERYVIZIER('2MASS-PSC',coords,25/60.,/ALLCOL,/SILENT)
   sz=SIZE(data) & if sz[0] ne 0 then BEGIN
      nfound=N_ELEMENTS(data.raj2000) & if nfound gt 1 then gd=WHERE(data._R eq MIN(data._R),ngd) else gd=0
      IF ABS(data[gd].hmag-h[i]) lt 0.2 THEN BEGIN
         apogee_id[i]='2M'+data[gd]._2mass
         j[i]=data[gd].jmag
         h[i]=data[gd].hmag
         k[i]=data[gd].kmag
         j_err[i]=data[gd].e_jmag
         h_err[i]=data[gd].e_hmag
         k_err[i]=data[gd].e_kmag
      ENDIF ELSE BEGIN
         IF ABS(data[gd].raj2000-ra[i]) lt 0.0003 and ABS(data[gd].dej2000-dec[i]) lt 0.0003 THEN BEGIN
            print,'No H match, but coords match'
            data=data[gd]
            apogee_id[i]='2M'+data[gd]._2mass
            j[i]=data[gd].jmag
            h[i]=data[gd].hmag
            k[i]=data[gd].kmag
            j_err[i]=data[gd].e_jmag
            h_err[i]=data[gd].e_hmag
            k_err[i]=data[gd].e_kmag
         ENDIF ELSE BEGIN
            print,'Found in 2MASS but H mag and coords do not match: '+strtrim(name[i],2)
            print,'Manually making APOGEE_ID'
            rh = FIX(ra[i]/15.)
            rm = FIX((ra[i]/15.-rh)*60.)
            rs = ((((ra[i]/15.-rh)*60.)-rm)*60.)
            rastr = STRTRIM(STRING(rh,F='(I02)'),2)+STRTRIM(STRING(rm,F='(I02)'),2)+STRTRIM(STRING(FIX(rs*100.),F='(I04)'),2)
            dd = FIX(dec[i])
            dm = FIX(ABS((dec[i]-dd)*60.))
            ds = ((ABS((dec[i]-dd)*60.)-dm)*60.)
            decstr = STRTRIM(STRING(ABS(dd),F='(I02)'),2)+STRTRIM(STRING(dm,F='(I02)'),2)+STRTRIM(STRING(FIX(ds*10.),F='(I03)'),2)
            if dec[i] ge 0. then id = 'AP'+rastr+'+'+decstr $
            else id = 'AP'+rastr+'-'+decstr
            IF name[i] eq 'VESTA' THEN id='VESTA'
            print,id
            apogee_id[i] = id
            j[i]=!VALUES.F_NAN
            k[i]=!VALUES.F_NAN
            j_err[i]=!VALUES.F_NAN
            h_err[i]=!VALUES.F_NAN
            k_err[i]=!VALUES.F_NAN            
            
         ENDELSE
      ENDELSE
   ENDIF ELSE IF name[i] eq 'VESTA' THEN apogee_id='VESTA' ELSE print,'PROBELM!!!: '+name[i]
   
ENDFOR

tmpstr=replicate({name:'a',apogee_id:'a',ra:-99.,dec:-99.},n_elements(name))
tmpstr.name=name
tmpstr.apogee_id=apogee_id
tmpstr.ra=ra
tmpstr.dec=dec

mwrfits,tmpstr,program+'.tmp.fits',/create

;-------------------
;NEED TO FIND b,v mags and HIP ID, HR,HD
;USE SIMBAD??
;--------------------------------

print,'MATCHING TO SIMBAD'
;try calling python script to query simbad?
SPAWN,'python simbad_names.py "'+program+'"'
;file_delete,program+'.tmp.fits'

readcol,program+'_simbadnames.txt',z1,hdname,z2,hrname,z3,hipname,z4,twomname,F='a,a,a,a,a,a,a,a'

hd=hdname
hip=hipname
hr=hrname
tmass=replicate('none',n_elements(name))
FOR i=0,n_elements(tmass)-1 DO IF twomname[i] eq 'none' THEN tmass[i]='none' $
ELSE tmass[i]='2M'+STRMID(STRTRIM(twomname[i],2),1)

FOR i=0,n_elements(name)-1 DO IF strtrim(tmass[i],2) ne strtrim(apogee_id[i],2) $
   AND tmass[i] ne 'none' THEN $
      print,'2MASS not match SIMBAD: '+apogee_id[i]+' '+tmass[i]

;------------------------------------------------------
;MATCH TARGETS TO HIPPARCOS CATALOG AND FILL VALUES
;------------------------------------------------------
print,'MATCHING TO HIPPARCOS'


;!!! ADD INFO FROM OLD HIP: VT MAG
;   ADD PM

hii=dblarr(n_elements(name))
FOR i=0,n_elements(hii)-1 DO BEGIN
   tmpii=where(strtrim(hh.hip,2) eq strtrim(hip[i],2))
   IF tmpii lt 0 THEN BEGIN
      print,'no hip id match: '+name[i]
      tmpii=where(ABS(hh._RAJ2000-ra[i]) lt 0.002 AND ABS(hh._DEJ2000-dec[i]) lt 0.002,num)
      IF num gt 1 THEN print,'multiple ra,dec, hip matches: '+strtrim(name[i],2)+'  '+strtrim(tmpii,2)
      IF tmpii lt 0 THEN print,'no hip coord match' $
      ELSE print,'hip coord match: HIP'+STRTRIM(hh[tmpii].hip,2) 
   ENDIF
   IF dec[i] lt -20. THEN print,'****** dec is < -20. '+strtrim(name[i],2)+strtrim(dec[i],2)+' ******'
   hii[i]=tmpii
ENDFOR

plx=dblarr(n_elements(name))
e_plx=dblarr(n_elements(name))
plxflag=fix(fltarr(n_elements(name)))

FOR i=0,n_elements(plx)-1 DO BEGIN
   IF hii[i] ge 0 THEN BEGIN
      IF strtrim(hip[i],2) eq 'none' then begin
         hip[i]=strtrim(hh[hii[i]].hip,2)
         plxflagtmp=1
      ENDIF ELSE plxflagtmp=0
      plx[i]=hh[hii[i]].plx
      e_plx[i]=hh[hii[i]].e_plx
      plxflag[i]=plxflagtmp
   ENDIF ELSE BEGIN
      plx[i]=!VALUES.F_NAN
      e_plx[i]=!VALUES.F_NAN
      plxflag[i]=0
   ENDELSE
ENDFOR


;------------------------------------------------------
;CREATE AND FILL STRUCTURE FOR NEW TARGET FILE
;------------------------------------------------------
print,'FILLING OTHER VARIABLES (TYPE, PERIOD, ETC)'
IF program ne 'rrlyr' THEN period=replicate(!VALUES.F_NAN,n_elements(name))

IF where(tn eq 'TYPE') GT 0 THEN BEGIN
   IF keyword_set(addtarg) THEN print,'May need add type to new targs'
ENDIF ELSE type=replicate('none',n_elements(name))

IF where(tn eq 'B') GT 0 THEN BEGIN
   IF keyword_set(addtarg) THEN b=[a.b,replicate(!VALUES.F_NAN,n_elements(newind))] $
   ELSE b=a.b
ENDIF ELSE b=replicate(!VALUES.F_NAN,n_elements(name))

;IF where(tn eq 'V') GT 0 THEN BEGIN
;   IF keyword_set(addtarg) THEN print,'need add targ for v' $
;   ELSE v=a.v
;ENDIF ELSE b=replicate(!VALUES.F_NAN,n_elements(name))

alt_name=replicate('none',n_elements(name))

print,'FILLING FULL STRUCTURE'
str = replicate({NAME:'a',HIP:'a',HD:'a',HR:'a',APOGEE_ID:'a',ALT_NAME:'a',$
                 RA:-99.,DEC:-99.,B:-99.,V:-99.,J:-99.,H:-99.,K:-99.,$
                 PLX:-99.,E_PLX:-99.,PLXFLAG:0,PERIOD:-99.,TYPE:'a',$
                 COMMENT:'a',PRIORITY:0},n_elements(name))

rsort=SORT(ra)
   ;print,'CHECKING FOR DUPLICATES'
   ;FOR i=1,n_elements(ra)-1 DO IF ABS(ra[rsort[i]]-ra[rsort[i-1]]) lt 0.01 $
   ;   AND ABS(dec[rsort[i]]-dec[rsort[i-1]]) lt 0.01 THEN $
   ;      print,'Possible Duplicate: ',strtrim(apogee_id[rsort[i]],2),'  ',strtrim(apogee_id[rsort[i-1]],2)

print,'CHECKING FOR DUPLICATES'
FOR nind=0,n_elements(name)-1 DO BEGIN 
   IF nind gt 0 AND ABS(ra[rsort[nind]]-ra[rsort[nind-1]]) lt 0.01 $
      AND ABS(dec[rsort[nind]]-dec[rsort[nind-1]]) lt 0.01 THEN $
         print,'Possible Duplicate: ',strtrim(apogee_id[rsort[nind]],2),$
               '  ',strtrim(apogee_id[rsort[nind-1]],2)
   IF strtrim(apogee_id[rsort[nind]],2) eq strtrim(apogee_id[rsort[nind-1]],2) $
      THEN print,'Removing second incidence of star: '$
                 +strtrim(apogee_id[rsort[nind]],2) $
   ELSE BEGIN
      str[nind].name=name[rsort[nind]]
      str[nind].hip=hip[rsort[nind]]
      str[nind].hd=hd[rsort[nind]]
      str[nind].hr=hr[rsort[nind]]
      str[nind].APOGEE_ID=APOGEE_ID[rsort[nind]]
      str[nind].alt_name=alt_name[rsort[nind]]
      str[nind].ra=ra[rsort[nind]]
      str[nind].dec=dec[rsort[nind]]
      str[nind].b=b[rsort[nind]]
      str[nind].v=v[rsort[nind]]
      str[nind].j=j[rsort[nind]]
      str[nind].h=h[rsort[nind]]
      str[nind].k=k[rsort[nind]]
      str[nind].plx=plx[rsort[nind]]
      str[nind].e_plx=e_plx[rsort[nind]]
      str[nind].plxflag=plxflag[rsort[nind]]
      str[nind].period=period[rsort[nind]]
      str[nind].type=type[rsort[nind]]
      str[nind].comment=comment[rsort[nind]]
      str[nind].priority=priority[rsort[nind]]
   ENDELSE
ENDFOR

str=str[where(strtrim(str.name,2) ne 'a')]

mwrfits,str,program+'.t.fits',/create

tt=mrdfits(program+'.t.fits',1)

help,tt,/str

END
