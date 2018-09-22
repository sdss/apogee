;+
;
; APQL_DBINSERT 
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
pro apql_dbinsert, instruct=str, savefile=savefile, exp_pk=exp_pk, $
   fitskw_version=fitskw_version,obs=obs

   ; setup SDSS database parameters
   SDSS_DB_PARAMS,obs=obs

   if n_elements(str) eq 0 and n_elements(savefile) eq 0 then begin
      print,'Usage:  apql_dbinsert, instruct=str, savefile=savestr, /UTR'
      return
   endif

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
         p=strcmp(sNames, 'str', /fold_case)
         pos=where(p eq 1)
         if pos lt 0 then begin
            print,'Error: expected structure named "str" is not in the savefile'
            OBJ_DESTROY, sObj
            return
         endif
         sObj->Restore, 'str'
         ; look for the structure named "fitskw_err"
         p=strcmp(sNames, 'fitskw_err', /fold_case)
         pos=where(p eq 1)
         if pos ge 0 then begin
            sObj->Restore, 'fitskw_err'
         endif
         ; look for the structure named "predict_str"
         p=strcmp(sNames, 'predict_str', /fold_case)
         pos=where(p eq 1)
         if pos ge 0 then begin
            sObj->Restore, 'predict_str'
         endif
         ; close the file
         OBJ_DESTROY, sObj
      endelse
      ; we delete the save file since we won't need it anymore
      FILE_DELETE, savefile, /quiet
   endif
   
   if n_elements(snrgoals_pk) eq 0 then begin
      ; look for the apogee_snr_goals_pk as the largest version in the table
      get_sql_col,'select pk from apogeeqldb.apogee_snr_goals order by version DESC limit 1',snrgoals_pk,/long
      if n_elements(snrgoals_pk) eq 0 then begin
         print, 'Error: no rows in the apogee_snr_goals table found'
         return
      endif
      snrgoals_pk=snrgoals_pk[0]
   endif

   
   ; start by inserting info relevant to all exptypes
   ; insert the data in the database using the idl-sql routines
   ; first all the scalars
   print,'Inserting basic information for  -> ',FILE_BASENAME(str.filename)
   FORMAT='("INSERT INTO apogeeqldb.quicklook (filename, exptype, readnum, exptime, jd, fitsheader_status,'+$
      'required_fitskeywords_version) VALUES (","''",A,"''",",",2("''",A,"'',"),2(F0.5,","),I0,",",I0,")")'

   exec_string = string(FILE_BASENAME(str.filename), str.exptype, str.readnum, str.exptime, str.jd, $
      str.fitsheader_status, fitskw_version, FORMAT=format)

   exec_sql, exec_string

   ; get the pk for the newly created row
   get_sql_col,'select pk from apogeeqldb.quicklook order by pk DESC limit 1',quicklook_pk,/long
   if n_elements(quicklook_pk) eq 0 then begin
      ; the row just entered was not found -> a problem occured during the previous insert
      print,'Error: the row for '+str.filename+' was not successfully INSERTED'
      return
   endif
   quicklook_pk=quicklook_pk[0]
   qlpk=string(quicklook_pk,format='(i0)')

   ;stop

   ; Update S/N information
   ;-----------------------
   print,'Inserting S/N information'
   if strupcase(str.exptype) eq 'OBJECT' then begin

      ; OBJECT: logsnr_hmag_coef, logsnr_hmag_coef_goal
      ; ALL: medsnr, snr_standard, snr2_time_coef, expected_total_readnum, delta_snr2_standard,
      ;         snr2_time_recent_coef, expected_total_readnum_recent


      FORMAT='("UPDATE apogeeqldb.quicklook SET snr_standard=",F0.3,", delta_snr2_standard=",F0.3,",'+$
               ' expected_total_readnum=",F0.3,", expected_total_readnum_recent=",F0.3," WHERE pk=",I0)'
      exec_string = string(str.snr_standard, str.delta_snr2_standard,str.expected_total_readnum, $
                           str.expected_total_readnum_recent, quicklook_pk, FORMAT=format)
      exec_sql, exec_string

      set_sql_colarray,'UPDATE apogeeqldb.quicklook SET medsnr=? WHERE pk='+qlpk, str.medsnr
      set_sql_colarray,'UPDATE apogeeqldb.quicklook SET logsnr_hmag_coef=? WHERE pk='+qlpk, str.logsnr_hmag_coef
      set_sql_colarray,'UPDATE apogeeqldb.quicklook SET logsnr_hmag_coef_goal=? WHERE pk='+qlpk, str.logsnr_hmag_coef_goal
      set_sql_colarray,'UPDATE apogeeqldb.quicklook SET snr2_time_coef=? WHERE pk='+qlpk,str.snr2_time_coef
      set_sql_colarray,'UPDATE apogeeqldb.quicklook SET snr2_time_coef_recent=? WHERE pk='+qlpk, str.snr2_time_recent_coef
   end
   if strupcase(str.exptype) eq 'QUARTZFLAT' or strupcase(str.exptype) eq 'DOMEFLAT' or $
      strupcase(str.exptype) eq 'ASDAF' then begin

      ; OBJECT: logsnr_hmag_coef, logsnr_hmag_coef_goal
      ; ALL: medsnr, snr_standard, snr2_time_coef, expected_total_readnum, delta_snr2_standard,
      ;         snr2_time_recent_coef, expected_total_readnum_recent

      FORMAT='("UPDATE apogeeqldb.quicklook SET snr_standard=",F0.3,", delta_snr2_standard=",F0.3,",'+$
               ' expected_total_readnum=",F0.3,", expected_total_readnum_recent=",F0.3," WHERE pk=",I0)'
      exec_string = string(str.snr_standard, str.delta_snr2_standard,str.expected_total_readnum, $
                           str.expected_total_readnum_recent, quicklook_pk, FORMAT=format)
      exec_sql, exec_string

      set_sql_colarray,'UPDATE apogeeqldb.quicklook SET medsnr=? WHERE pk='+qlpk, str.medsnr
      set_sql_colarray,'UPDATE apogeeqldb.quicklook SET snr2_time_coef=? WHERE pk='+qlpk,str.snr2_time_coef
      set_sql_colarray,'UPDATE apogeeqldb.quicklook SET snr2_time_coef_recent=? WHERE pk='+qlpk, str.snr2_time_recent_coef
   end

   ;stop

   ; Update wavelength information 
   ;------------------------------
   if strupcase(str.exptype) eq 'OBJECT' or strupcase(str.exptype) eq 'ASDAF' then begin
      print,'Inserting wavelength information'
      FORMAT='("UPDATE apogeeqldb.quicklook SET wavefit_rms=",F0.3,",wavelength_status=",I0," WHERE pk=",I0)'
      exec_string = string(str.wavefit_rms,str.wavelength_status, quicklook_pk, FORMAT=format)
      exec_sql, exec_string

      set_sql_colarray,'UPDATE apogeeqldb.quicklook SET wavefit_pars=? WHERE pk='+qlpk, str.wavefit_pars
      set_sql_colarray,'UPDATE apogeeqldb.quicklook SET waverange_exp=? WHERE pk='+qlpk, str.waverange_exp
      set_sql_colarray,'UPDATE apogeeqldb.quicklook SET waverange_meas=? WHERE pk='+qlpk, str.waverange_meas
      set_sql_colarray,'UPDATE apogeeqldb.quicklook SET waverange_diff=? WHERE pk='+qlpk, str.waverange_diff
   endif

   ;stop

   ; Update sky and dither information
   if strupcase(str.exptype) eq 'OBJECT' or strupcase(str.exptype) eq 'ASDAF' then begin
      print,'Inserting dither information'
      FORMAT='("UPDATE apogeeqldb.quicklook SET skyvar_meddev_perc=",F0.3," '+$
         ',skyvar_stddev_perc=",F0.3,",skyvar_contflux=",F0,",skyvar_contflux_rate=",F0," '+$
         ',skyvar_avglineflux=",F0,",skyvar_avglineflux_rate=",F0,",dither_prevexp_measured=",F0.2," '+ $
         ',dither_prevexp_header=",F0.2,",dither_relative=",F0.2,",sky_status=",A0," '+$
         ',dither_status=",A0," WHERE pk=",I0)'

      exec_string = string(str.skyvar_meddev_perc, str.skyvar_stddev_perc, str.skyvar_contflux, $
         str.skyvar_contflux_rate, str.skyvar_avglineflux, str.skyvar_avglineflux_rate, $
         str.dither_prevexp_measured, str.dither_prevexp_header, str.dither_relative, str.sky_status, $
         str.dither_status, quicklook_pk, FORMAT=format)

      ;print,'exec_string=',exec_string
      exec_sql, exec_string
   endif

   ;stop

   ; Finished Status
   if strupcase(str.exptype) eq 'OBJECT' then begin
      print,'Inserting finished status information'
      FORMAT='("UPDATE apogeeqldb.quicklook SET exp_finished_status=",I0,", visit_finished_status=",I0,'+$
         '" WHERE pk=",I0)'
      exec_string = string(str.exp_finished_status, str.visit_finished_status, quicklook_pk, FORMAT=format)
      exec_sql, exec_string
   endif

   ;stop

   ;if strupcase(str.exptype) eq 'OBJECT' or strupcase(str.exptype) eq 'ASDAF' then begin
   ;   FORMAT='("UPDATE apogeeqldb.quicklook SET wavefit_rms=",F0.3,",skyvar_meddev_perc=",F0.3," '+$
   ;      ',skyvar_stddev_perc=",F0.3,",skyvar_contflux=",F0,",skyvar_contflux_rate=",F0," '+$
   ;      ',skyvar_avglineflux=",F0,",skyvar_avglineflux_rate=",F0,",dither_prevexp_measured=",F0.2," '+ $
   ;      ',dither_prevexp_header=",F0.2,",dither_relative=",F0.2,",sky_status=",I0," '+$
   ;      ',wavelength_status=",I0,",dither_status=",I0,",snr_goals_version=",I0," WHERE pk=",I0)'
   ;
   ;   exec_string = string(str.wavefit_rms, str.skyvar_meddev_perc, $
   ;      str.skyvar_stddev_perc, str.skyvar_contflux, str.skyvar_contflux_rate, $
   ;      str.skyvar_avglineflux, str.skyvar_avglineflux_rate, str.dither_prevexp_measured, $
   ;      str.dither_prevexp_header, str.dither_relative, str.sky_status, $
   ;      str.wavelength_status, str.dither_status, str.snr_goals_version, quicklook_pk, FORMAT=format)
   ;
   ;   exec_sql, exec_string
   ;
   ;   set_sql_colarray,'UPDATE apogeeqldb.quicklook SET wavefit_pars=? WHERE pk='+qlpk, str.wavefit_pars
   ;   set_sql_colarray,'UPDATE apogeeqldb.quicklook SET waverange_exp=? WHERE pk='+qlpk, str.waverange_exp
   ;   set_sql_colarray,'UPDATE apogeeqldb.quicklook SET waverange_meas=? WHERE pk='+qlpk, str.waverange_meas
   ;   set_sql_colarray,'UPDATE apogeeqldb.quicklook SET waverange_diff=? WHERE pk='+qlpk, str.waverange_diff
   ;endif

   ;if strupcase(str.exptype) eq 'OBJECT' or strpos(strupcase(str.exptype),'FLAT') ge 0 then begin
   ;
   ;   FORMAT='("UPDATE apogeeqldb.quicklook SET snr_standard=",F0.3,", expected_total_readnum=",F0.3,",'+$
   ;      ' expected_total_readnum_recent=",F0.3,", exp_finished_status=",I0,", visit_finished_status=",I0,'+$
   ;      '" WHERE pk=",I0)'
   ;   exec_string = string(str.snr_standard, str.expected_total_readnum, str.expected_total_readnum_recent, $
   ;      str.exp_finished_status, str.visit_finished_status, quicklook_pk, FORMAT=format)
   ;   exec_sql, exec_string
   ;
   ;   set_sql_colarray,'UPDATE apogeeqldb.quicklook SET medsnr=? WHERE pk='+qlpk, str.medsnr
   ;   set_sql_colarray,'UPDATE apogeeqldb.quicklook SET logsnr_hmag_coef=? WHERE pk='+qlpk, str.logsnr_hmag_coef
   ;   set_sql_colarray,'UPDATE apogeeqldb.quicklook SET logsnr_hmag_coef_goal=? WHERE pk='+qlpk, str.logsnr_hmag_coef_goal
   ;   set_sql_colarray,'UPDATE apogeeqldb.quicklook SET snr2_time_coef=? WHERE pk='+qlpk,str.snr2_time_coef
   ;   set_sql_colarray,'UPDATE apogeeqldb.quicklook SET snr2_time_coef_recent=? WHERE pk='+qlpk, str.snr2_time_recent_coef
   ;
   ;endif


   if n_elements(exp_pk) gt 0 then begin
      print,'Inserting exposure_pk information ', exp_pk
      exec_sql,'UPDATE apogeeqldb.quicklook SET exposure_pk='+strtrim(string(exp_pk),2)+' WHERE pk='+qlpk
   endif

   ;stop

   ; Update the required_fitskeywords_error table if error were found
   if PTR_VALID(str.fitsheader_errors) then begin
      print,'Inserting required fits keywords error information'
      FORMAT='("INSERT INTO apogeeqldb.required_fitskeywords_error (quicklook_pk, name, errortype_pk) ' + $
             'VALUES (",I0,",''",A,"'',",I0,")")'
      for i=0, n_elements((*str.fitsheader_errors).name)-1 do begin
         p = where((*str.fitsheader_errors).name[i] eq fitskw_err.name,count)
         if count ge 1 then begin
            exec_string = string(quicklook_pk, (*str.fitsheader_errors).name[i], $
               fitskw_err.pk[p[0]], FORMAT=format)
            exec_sql, exec_string
         endif
      endfor
   endif

   ;stop

   ; Check if we need to insert data in the quicklook60 table 
   ;  binned images and representative spectra
   if PTR_VALID(str.arraydisplay_sub[0].data) then begin
      print,'Inserting quicklook60 information'

      FORMAT='("INSERT INTO apogeeqldb.quicklook60 (quicklook_pk, bzero, bscale, zscale1, zscale2) '+$
             'VALUES (",I0,",",3(F0.8,","),F0.8,")")'
      exec_string = string(quicklook_pk, str.arraydisplay.bzero, str.arraydisplay.bscale, str.arraydisplay.zscale[0], $
         str.arraydisplay.zscale[1], FORMAT=format)
      exec_sql, exec_string

      ; get the pk for the newly created row
      get_sql_col,'select pk from apogeeqldb.quicklook60 order by pk DESC limit 1',quicklook60_pk,/long

      ql60pk=string(quicklook60_pk,format='(i0)')

      ; finish loading quicklook60
      set_sql_colarray,'UPDATE apogeeqldb.quicklook60 SET data=? WHERE pk='+ql60pk, str.arraydisplay.data
      dim1=n_elements(str.medsky)
      dim2=n_elements(str.medsky[0].flux)
      medsky=fltarr(dim1*dim2)
      p1=0L
      p2=dim2-1L
      for i=0L,dim1-1L do begin
         medsky[p1:p2]=str.medsky[i].flux
         p1=p2+1L
         p2=(p1+dim2-1L)
      endfor
      set_sql_colarray,'UPDATE apogeeqldb.quicklook60 SET medsky=? WHERE pk='+ql60pk, medsky

      ; load the quicklook60_imbinzoom
      for nimzoom=0, n_elements(str.arraydisplay_sub)-1 do begin
         FORMAT='("INSERT INTO apogeeqldb.quicklook60_imbinzoom (quicklook60_pk, ylo, yhi, bzero, bscale, zscale1, zscale2) '+$
                'VALUES (",3(I0,","),3(F0.8,","),F0.8,")")'
         exec_string = string(quicklook60_pk, str.arraydisplay_sub[nimzoom].yrange[0], str.arraydisplay_sub[nimzoom].yrange[1], $
            str.arraydisplay_sub[nimzoom].bzero, str.arraydisplay_sub[nimzoom].bscale, str.arraydisplay_sub[nimzoom].zscale[0], $
            str.arraydisplay_sub[nimzoom].zscale[1], FORMAT=format)
         exec_sql, exec_string

         ; get the pk for the newly created row
         get_sql_col,'select pk from apogeeqldb.quicklook60_imbinzoom order by pk DESC limit 1',quicklook60_imbinzoom_pk,/string

         ; finish loading quicklook60_imbinzoom
         set_sql_colarray,'UPDATE apogeeqldb.quicklook60_imbinzoom SET data=? WHERE pk='+quicklook60_imbinzoom_pk[0], $
            *(str.arraydisplay_sub[nimzoom].data)
      endfor
 
      ; load the quicklook60_repspec
      if strupcase(str.exptype) eq 'OBJECT' or strupcase(str.exptype) eq 'QUARTZFLAT' or strupcase(str.exptype) eq 'DOMEFLAT' or $
         strupcase(str.exptype) eq 'ARCLAMP' or strupcase(str.exptype) eq 'ASDAF' then begin
         for nspec=0, n_elements(str.representative_spectra)-1 do begin
            FORMAT='("INSERT INTO apogeeqldb.quicklook60_repspec (quicklook60_pk, fiberid, bzero, bscale, medsnr) '+$
                   'VALUES (",2(I0,","),2(F0.8,","),F0.8,")")'
            exec_string = string(quicklook60_pk, str.representative_spectra[nspec].fiberid, $
               str.representative_spectra[nspec].bzero, str.representative_spectra[nspec].bscale, $
               str.representative_spectra[nspec].medsnr, FORMAT=format)
            ;   str.medsnr[str.representative_spectra[nspec].fiberid], FORMAT=format)
            exec_sql, exec_string

            ; get the pk for the newly created row
            get_sql_col,'select pk from apogeeqldb.quicklook60_repspec order by pk DESC limit 1',qlk60_repspec_pk,/string

            ; finish loading quicklook60_repspec
            if PTR_VALID(str.representative_spectra[nspec].data) then begin
                set_sql_colarray,'UPDATE apogeeqldb.quicklook60_repspec SET spectrum=? WHERE pk='+qlk60_repspec_pk[0], $
                   *(str.representative_spectra[nspec].data)
            endif
         endfor
      endif
 
   endif

   ;stop

   ; Check if we need to insert data in the quicklook_prediction table 
   if n_elements(predict_str) ne 0 then begin
      print,'Inserting prediction information'

      for predpos=0, n_elements(predict_str)-1 do begin
         FORMAT='("INSERT INTO apogeeqldb.quicklook_prediction (quicklook_pk, frameid, exptime, snr_standard, nread, '+$
            ' ditherpos, exp_stopcode, visit_stopcode) VALUES (",I0,",''",A,"'',",2(F0.4,","),3(I0,","),I0,")")'
         ; handle the case where snr = -NaN (which the database doesn't like)
         if FINITE(predict_str[predpos].snr_standard) eq 0 then predict_str[predpos].snr_standard=0.0
         exec_string = string(quicklook_pk, predict_str[predpos].frameid, predict_str[predpos].exptime, $
             predict_str[predpos].snr_standard, predict_str[predpos].nreads, predict_str[predpos].ditherexposure, $
             predict_str[predpos].exp_finished_status, predict_str[predpos].visit_finished_status, FORMAT=format)
         ;print,'exec_string='
         ;print,exec_string
         exec_sql, exec_string
      endfor

   endif

   ;stop

END
