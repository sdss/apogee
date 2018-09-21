;+
;
; APQL_WRITEOUTPUTFILE
;
; This routine is called by the quicklook program to print out information
; in to a output file. The information is passed as a data structure.
; 
;-
;
pro apql_writeoutputfile, outfile, allstr, all_predictstr, exp_pk=exp_pk, fitskw_version=fitskw_version

   if n_elements(allstr) eq 0 or n_elements(outfile) eq 0 then begin
      print,'Usage:  apql_writeoutputfile, outfile, allstr, all_predictstr'
      return
   endif

   openw,loglun,outfile,/append,/get_lun

   printf,loglun,'EXPOSURE: ',allstr[0].frameid

   ql_line = ['TABLE: QUICKLOOK']
   ql_line = [ql_line,'filename, exposure_pk, exptype, readnum, exptime, jd, fitsheader_status, required_fits, ' $
       +'keywords_version, snr_standard, delta_snr2_standard, expected_total_readnum, expected_total_readnum_recent, '$
       +'logsnr_hmag_coef, logsnr_hmag_coef_goal, snr2_time_coef, snr2_time_coef_recent, wavefit_rms, wavelength_status, '$
       +'skyvar_meddev_perc, skyvar_stddev_perc, skyvar_contflux, skyvar_contflux_rate, skyvar_avglineflux, skyvar_avglineflux_rate, '$
       +'sky_status, dither_prevexp_measured, dither_prevexp_header, dither_relative, dither_status, exp_finished_status, '$
       +'visit_finished_status']
   ql_format='(4(A,","),2(F0.5,","),I0,",",I0,",",4(F0.3,","),4("(",F0.3,",",F0.3,"),"),F0.3,",",I0,",",2(F0.2,","),4(F0,","),' $
       +'I0,","3(F0.2,","),I0,",",I0,",",I0)'

   ql60_line = ['TABLE: QUICKLOOK60']
   ql60_line = [ql60_line,'filename, bzero, bscale, zscale1, zscale2']
   ql60_format = '(A,",",3(F0.8,","),F0.8)'

   ql60i_line = ['TABLE: QUICKLOOK60_IMBINZOOM']
   ql60i_line = [ql60i_line,'filename, ylo, yhi, bzero, bscale, zscale1, zscale2']
   ql60i_format = '(A,","2(I0,","),3(F0.8,","),F0.8)'

   ql60r_line = ['TABLE: QUICKLOOK60_REPSPEC']
   ql60r_line = [ql60r_line,'filename, fiberid, bzero, bscale, medsnr']
   ql60r_format = '(A,",",I0,",",2(F0.8,","),F0.8)'
 
   qlp_line = ['TABLE: QUICKLOOK_PREDICTION']
   qlp_line = [qlp_line,'filename, frameid, exptime, snr_standard, nread, ditherpos, exp_stopcode, visit_stopcode']
   qlp_format = '(A,",",I0,",",2(F0.4,","),3(I0,","),I0)'


   exptype = allstr[0].exptype 
   nreads = n_elements(allstr)

   for i = 0, nreads -1 do begin
         str = allstr[i]
         filename =  file_basename(str.filename)
         print,'Writing QUICKLOOK information for  -> ',FILE_BASENAME(str.filename)

         ql_line = [ql_line,string(filename, strtrim(exp_pk,2), str.exptype, str.readnum, str.exptime, str.jd, $
                           str.fitsheader_status, fitskw_version, str.snr_standard, str.delta_snr2_standard,str.expected_total_readnum, $
                           str.expected_total_readnum_recent, str.logsnr_hmag_coef[0], str.logsnr_hmag_coef[1], $
                           str.logsnr_hmag_coef_goal[0], str.logsnr_hmag_coef_goal[1], str.snr2_time_coef[0], $
                           str.snr2_time_coef[1], str.snr2_time_recent_coef[0], str.snr2_time_recent_coef[1], $ 
                           str.wavefit_rms, str.wavelength_status, str.skyvar_meddev_perc, str.skyvar_stddev_perc, $
                           str.skyvar_contflux, str.skyvar_contflux_rate, str.skyvar_avglineflux, str.skyvar_avglineflux_rate, $
                           str.sky_status, str.dither_prevexp_measured, str.dither_prevexp_header, str.dither_relative,  $
                           str.dither_status, str.exp_finished_status, str.visit_finished_status, $
                           FORMAT=ql_format)]

         ; Check if we need to insert data in the quicklook60 table 
         ; binned images and representative spectra
         if str.arraydisplay_sub[0].bscale ne 0 then begin  
            print,'Inserting quicklook60 information'

            ql60_line = [ql60_line,string(filename, str.arraydisplay.bzero, str.arraydisplay.bscale, $
                           str.arraydisplay.zscale[0], str.arraydisplay.zscale[1], $
                           FORMAT=ql60_format)]


            ; load the quicklook60_imbinzoom
            for nimzoom=0, n_elements(str.arraydisplay_sub)-1 do begin
               ql60i_line = [ql60i_line,string(filename,str.arraydisplay_sub[nimzoom].yrange[0], $
                              str.arraydisplay_sub[nimzoom].yrange[0], str.arraydisplay_sub[nimzoom].bzero, $
                              str.arraydisplay_sub[nimzoom].bscale, str.arraydisplay_sub[nimzoom].zscale[0], $
                              str.arraydisplay_sub[nimzoom].zscale[1], $
                              FORMAT=ql60i_format)]
            endfor

            ; load the quicklook60_repspec
            if exptype eq 'OBJECT' or exptype eq 'QUARTZFLAT' or exptype eq 'DOMEFLAT' or $
               exptype eq 'ARCLAMP' or exptype eq 'ASDAF' then begin
                  for nspec=0, n_elements(str.representative_spectra)-1 do begin
                     ql60r_line = [ql60r_line, string(filename, str.representative_spectra[nspec].fiberid, $
                                 str.representative_spectra[nspec].bzero, str.representative_spectra[nspec].bscale, $
                                 str.representative_spectra[nspec].medsnr, $
                                 FORMAT=ql60r_format)]
                  endfor
            endif



         endif


      ; Check if we need to insert data in the quicklook_prediction table 
      if n_elements(all_predictstr) gt 0 then begin
         tmp = where(all_predictstr.filename eq filename,npredict)
         if npredict ne 0 then begin
            predict_str = all_predictstr[tmp]
            print,'Inserting prediction information'
            for predpos=0, npredict-1 do begin
               ; handle the case where snr = -NaN (which the database doesn't like)
               if FINITE(predict_str[predpos].snr_standard) eq 0 then predict_str[predpos].snr_standard=0.0
               qlp_line = [qlp_line, string(filename, predict_str[predpos].frameid, predict_str[predpos].exptime, $
                           predict_str[predpos].snr_standard, predict_str[predpos].nreads, predict_str[predpos].ditherexposure, $
                           predict_str[predpos].exp_finished_status, predict_str[predpos].visit_finished_status, $
                           FORMAT=qlp_format)]
             endfor

         endif
      endif  
   endfor

   printf,loglun, ql_line,format='(A)'
   printf,loglun, ''
   printf,loglun, ql60_line,format='(A)'
   printf,loglun, ''
   printf,loglun, ql60i_line,format='(A)'
   printf,loglun, ''
   printf,loglun, ql60r_line,format='(A)'
   printf,loglun, ''
   printf,loglun, qlp_line,format='(A)'
   printf,loglun, ''
   printf,loglun, ''
   close,loglun
   

end
