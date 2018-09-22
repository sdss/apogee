;+
;
; APQUICKRED_DBINSERT 
;
; This routine is called by the quicklook program to load information
; in the database. The information is passed as a data structure
; or as a savefile with the data structure written to it (one or the other
; is required).
;
; The savefile is used for now when this routine is called through 
; IDL_BRIDGE->EXECUTE because of the limitations of passing structures 
; and pointers.
; 
;-
;
pro apquickred_dbinsert, instruct=dbstr, savefile=savefile, exp_pk=exp_pk

   t0 = systime(1)  

   if n_elements(dbstr) eq 0 and n_elements(savefile) eq 0 then begin
      print,'Usage:  apquickred_dbinsert, instruct=str, savefile=savestr, exp_pk=exp_pk'
      return
   endif

   print,'starting apquickred_dbinsert ...'
   if n_elements(savefile) gt 0 then begin
      ; restore the data
      res=file_search(savefile, count=count, /TEST_READ)
      if count eq 0 then begin
         print,'Error: file '+savefile+' not found'
         return
      endif else if count gt 1 then begin
         print,'Error: more than one found ('+savefile+')'
         return
      endif else begin
         sObj = OBJ_NEW('IDL_Savefile', savefile)
         sNames = sObj->Names()
         ; look for the structure named "str"
         p=strcmp(sNames, 'dbstr', /fold_case)
         pos=where(p eq 1)
         if pos lt 0 then begin
            print,'Error: expected structure named "str" is not in the savefile'
            return
         endif
         sObj->Restore, 'dbstr'
         ; close the file
         OBJ_DESTROY, sObj
      endelse
      ; we delete the save file since we won't need it anymore
      ;FILE_DELETE, savefile, /quiet
   endif    

   ; get all the tag_names for the structure
   tnames = tag_names(dbstr)

   ; get the latest quicklook pk to put in the table
   ; get the pk for the newly created row for this eposure
   print,'exp_pk: ', exp_pk
   get_sql_col,'select pk from apogeeqldb.quicklook where exposure_pk='+strtrim(string(exp_pk,format="(i0)"))+$
       'order by pk DESC limit 1',quicklook_pk,/long
   if n_elements(quicklook_pk) eq 0 then begin
      ; the row just entered was not found -> a problem occured during the previous insert
      print,'Error: the row for '+dbstr.frameid+' was not successfully INSERTED in quicklook'
      return
   endif
   quicklook_pk=quicklook_pk[0]
   qlpk=string(quicklook_pk,format='(i0)')

   print,'Inserting Quickred data'

   ; insert the data in the database using the idl-sql routines
   ; first all the scalars
   p=where(strpos(tnames,'SNR_STANDARD') ge 0,count)
   if count gt 0 then snr_standard=dbstr.snr_standard else snr_standard=0.0
   FORMAT='("INSERT INTO apogeeqldb.quickred (exposure_pk, last_quicklook_pk, bzero, bscale, '+ $
          'zscale1, zscale2, dither_pixpos, snr_standard) VALUES (",2(I0,","),5(F0.8,","),F0.8,")")'

   exec_string = string(exp_pk, quicklook_pk, dbstr.arraydisplay.bzero, dbstr.arraydisplay.bscale, $
      dbstr.arraydisplay.zscale[0], dbstr.arraydisplay.zscale[1], dbstr.dithpix, snr_standard,$
      FORMAT=format)

   exec_string = REPSTR(exec_string,"NaN","'NaN'")  ; put any NaNs in correct format
   exec_sql, exec_string

   get_sql_col,'select pk from apogeeqldb.quickred order by pk DESC limit 1',quickred_pk,/long
   if n_elements(quickred_pk) eq 0 then begin
      ; the row just entered was not found -> a problem occured during the previous insert
      print,'Error: the row for '+dbstr.frameid+' was not successfully INSERTED quickred'
      return
   endif
   qrpk=string(quickred_pk,format='(i0)')
   set_sql_colarray,'UPDATE apogeeqldb.quickred SET data=? WHERE pk='+qrpk, dbstr.arraydisplay.data
   ; may or may be in structure (depends if a plate is associated with exposure)
   p=where(strpos(tnames,'LOGSNR_HMAG_COEF') ge 0,count)
   if count gt 0 then logsnr=dbstr.logsnr_hmag_coef else logsnr=[0.0,0.0]
   set_sql_colarray,'UPDATE apogeeqldb.quickred SET logsnr_hmag_coef=? WHERE pk='+qrpk, logsnr
   ; add snr_goals_version
   ; may or may be in structure (depends if a plate is associated with exposure)
   p=where(strpos(tnames,'HMAG_STANDARD_VERSION') ge 0,count)
   if count gt 0 then snr_goals_version=dbstr.hmag_standard_version else snr_goals_version=0
   set_sql_colarray,'UPDATE apogeeqldb.quickred SET snr_goals_version=? WHERE pk='+qrpk, snr_goals_version


   ; load the QUICKRED_IMBINZOOM
   print,'Inserting Quickred IMBINZOOM'
   for zpos=0, n_elements(dbstr.arraydisplay_sub)-1 do begin
      FORMAT='("INSERT INTO apogeeqldb.quickred_imbinzoom (quickred_pk, ylo, yhi, bzero, bscale, ' + $
          'zscale1, zscale2) VALUES (",3(I0,","),3(F0.8,","),F0.8,")")'
      exec_string = string(quickred_pk, dbstr.arraydisplay_sub[zpos].yrange[0], dbstr.arraydisplay_sub[zpos].yrange[1], $
          dbstr.arraydisplay_sub[zpos].bzero, dbstr.arraydisplay_sub[zpos].bscale, $
          dbstr.arraydisplay_sub[zpos].zscale[0], dbstr.arraydisplay_sub[zpos].zscale[1], FORMAT=format)
      exec_sql, exec_string

      ; get the pk of previously inserted row
      get_sql_col,'select pk from apogeeqldb.quickred_imbinzoom order by pk DESC limit 1',quickred_imbinzoom_pk,/long
      if n_elements(quickred_imbinzoom_pk) eq 0 then begin
          ; the row just entered was not found -> a problem occured during the previous insert
          print,'Error: the row for '+dbstr.frameid+' was not successfully INSERTED quickred_imbinzoom'
          return
      endif
      qrzpk=string(quickred_imbinzoom_pk,format='(i0)')
      set_sql_colarray,'UPDATE apogeeqldb.quickred_imbinzoom SET data=? WHERE pk='+qrzpk, dbstr.arraydisplay_sub[zpos].data

   endfor

   ; load the QUICKRED_SPECTRUM
   p=where(strpos(tnames,'QR_SPECTRUM') ge 0,count)
   if count gt 0 then begin
       print,'Inserting Quickred SPECTRUM'

       for fid=0, n_elements(dbstr.qr_spectrum)-1 do begin
          FORMAT='("INSERT INTO apogeeqldb.quickred_spectrum (quickred_pk, fiberid, bzero, bscale, ' + $
              'medsnr) VALUES (",2(I0,","),2(F0.8,","),F0.8,")")'
          exec_string = string(quickred_pk, dbstr.qr_spectrum[fid].fiberid, $
              dbstr.qr_spectrum[fid].bzero, dbstr.qr_spectrum[fid].bscale, dbstr.qr_spectrum[fid].medsnr, FORMAT=format)
          exec_sql, exec_string

          ; get the pk of previously inserted row
          get_sql_col,'select pk from apogeeqldb.quickred_spectrum order by pk DESC limit 1',quickred_spectrum_pk,/long
          if n_elements(quickred_spectrum_pk) eq 0 then begin
              ; the row just entered was not found -> a problem occured during the previous insert
              print,'Error: the row for '+dbstr.frameid+' was not successfully INSERTED quickred_spectrum'
              return
          endif
          qrspk=string(quickred_spectrum_pk,format='(i0)')
          set_sql_colarray,'UPDATE apogeeqldb.quickred_spectrum SET spectrum=? WHERE pk='+qrspk, dbstr.qr_spectrum[fid].spectrum
       endfor
   endif

   ; load the EXPOSURE HEADER information
   ; exposure_header_keyword, exposure_header_value
   ; exposure_header_value table
   ;  pk, value, comment, index, exposure_pk, exposure_header_keyword_pk
   ; 
   ; -get a list of all the keywords that we have
   ; -get all of the rows in exposure_header_keyword and see if we have all of them
   ;   if not add the new ones that we need
   ; -get exposure_header_keyword_pk for each line
   p=where(strpos(tnames,'HEADER') ge 0,count)
   if count gt 0 then begin   

     print,'Inserting exposure header information'

     ; Only get good lines, don't want COMMENT, HISTORY, END, etc.
     ninechar = strmid(dbstr.header,8,1)
     gdhd = where(ninechar eq '=',ngdhd)
     if ngdhd eq 0 then begin
       print,'No good header lines'
       return
     endif

     ; Get the keyword/value/comments
     head = dbstr.header[gdhd]
     nhead = n_elements(head)
     headstr = replicate({index:-1,keyword:'',value:'',comment:'',keyword_pk:-1L},nhead)
     headstr.index = indgen(ngdhd)
     for i=0,nhead-1 do begin
       key = strtrim(strmid(head[i],0,8),2)
       headstr[i].keyword = key
       value = sxpar(dbstr.header,key,comment=comment,count=ncount)
       if ncount gt 0 then begin
         headstr[i].value = value
         headstr[i].comment = comment
       endif
       ;left = strmid(head[i],9)
       ;dum = strsplit(left,'/',/extract)
       ;headstr[i].value = strtrim(dum[0],2)
       ;if n_elements(dum) gt 1 then headstr[i].comment=strtrim(dum[1],2)
     end

     ; Get keywords already in database
     get_sql_col,'select pk,label from platedb.exposure_header_keyword order by pk',keyword_pk,keyword_label,/string
     if n_elements(keyword_pk) gt 0 then dbkeywords = {pk:keyword_pk, label:keyword_label}

     ; Get the keyword_pk
     if n_elements(keyword_pk) gt 0 then begin
       for i=0,nhead-1 do begin
         gdkey = where(dbkeywords.label eq headstr[i].keyword,ngdkey)
         if ngdkey gt 0 then headstr[i].keyword_pk = dbkeywords.pk[gdkey[0]]
       end
     end

     ; Some keywords missing
     nokey = where(headstr.keyword_pk eq -1,n_nokey)
     if n_nokey gt 0 then begin

       ; Add new keywords to database
       for i=0,n_nokey-1 do begin
         newkey = headstr[nokey[i]].keyword
         exec_string = "INSERT INTO platedb.exposure_header_keyword (label) VALUES ('"+newkey+"')"
         print,'Adding keyword ',newkey,' to platedb.exposure_header_keyword'
         exec_sql, exec_string
       endfor 

       ; RE-Get keywords already in database
       apgundef,keyword_pk,keyword_label
       get_sql_col,'select pk,label from platedb.exposure_header_keyword order by pk',keyword_pk,keyword_label,/string
       if n_elements(keyword_pk) gt 0 then dbkeywords = {pk:keyword_pk, label:keyword_label}

       ; Get the keyword_pk
       if n_elements(keyword_pk) gt 0 then begin
         for i=0,nhead-1 do begin
           gdkey = where(dbkeywords.label eq headstr[i].keyword,ngdkey)
           if ngdkey gt 0 then headstr[i].keyword_pk = dbkeywords.pk[gdkey[0]]
         end
       end

       ; Still some missing keywords
       bdkey = where(headstr.keyword_pk eq -1,nbdkey)
       if nbdkey gt 0 then begin
         print,'Still some missing keywords'
       endif

     end ; need to add keywords to database

     ; Add values to exposure_header_value
     for i=0,nhead-1 do begin
       if headstr[i].keyword_pk ne -1 then begin
         ; the fits header values already have '' if they are strings
         value = headstr[i].value
         if strmid(value,0,1) ne "'" then value="'"+value+"'"  ; add tick marks if necessary
         ; we need to clean up the comment field in case in contains quotes (') since will mess up the sql INSERT
         this_cmnt = headstr[i].comment
         ppos=0
         ; handle single quotes
         repeat begin
             pos = strpos(this_cmnt,"'")
             if pos ge 0 then this_cmnt = strmid(this_cmnt,ppos,pos)+strmid(this_cmnt,pos+1)
             ppos = (pos-1)>0
         endrep until (pos eq -1)
         ; now handle double quotes
         repeat begin
             pos = strpos(this_cmnt,'"')
             if pos ge 0 then this_cmnt = strmid(this_cmnt,ppos,pos)+strmid(this_cmnt,pos+1)
             ppos = (pos-1)>0
         endrep until (pos eq -1)
         exec_string = "INSERT INTO platedb.exposure_header_value (value,comment,index,exposure_pk,exposure_header_keyword_pk) "+$
                      " VALUES ("+value+",'"+this_cmnt+"',"+strtrim(headstr[i].index,2)+","+$
                       strtrim(exp_pk,2)+","+strtrim(headstr[i].keyword_pk,2)+")"
         ;print,'EXEC_STRING='+exec_string
         exec_sql, exec_string
       endif
     end

   endif ; header exists

   print,'apquickred_dbinsert completed'
   dt = systime(1)-t0
   print,'dt = ',strtrim(dt,2),' sec.'

   ; This takes ~200 sec for a full extracted exposure

   ;stop

END
