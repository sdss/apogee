;+
;
; APQUICKRED_WRITEOUPUTFILE 
;
; This routine is called by the quickred program to print out information
; in to a output file. The information is passed as a data structure.
; 
;-
;
pro apquickred_writeoutputfile, outfile, dbstr, exp_pk=exp_pk

   if n_elements(dbstr) eq 0 or n_elements(outfile) eq 0 then begin
      print,'Usage:  apquickred_writeoutputfile, outfile, dbstr, exp_pk=exp_pk'
      return
   endif

   openw,loglun,outfile,/append,/get_lun
 
   ; get all the tag_names for the structure
   tnames = tag_names(dbstr)

   ; find quicklook_pk for this exposure, if it has been inserted.
   if exp_pk eq 'None' then begin
      qlpk = 'None'
   endif else begin
     get_sql_col,'select pk from apogeeqldb.quicklook where exposure_pk='+strtrim(string(exp_pk,format="(i0)"))+$
         'order by pk DESC limit 1',quicklook_pk,/long
      if n_elements(quicklook_pk) eq 0 then begin
        ; the row just entered was not found -> a problem occured during the previous insert
        qlpk = 'None'
      endif else begin
        quicklook_pk=quicklook_pk[0]
        qlpk=string(quicklook_pk,format='(i0)')
      endelse
   endelse

   printf,loglun,'EXPOSURE: ',dbstr[0].frameid

   qr_line = ['TABLE: QUICKRED']
   qr_line = [qr_line,'exposure_pk, last_quicklook_pk, bzero, bscale, zscale1, zscale2, dither_pixpos, snr_standard, ' $
       +'logsnr_hmag_coef, snr_goals_version']
   qr_format = '(2(A,","),5(F0.8,","),F0.8,",[",F0.8,",",F0.8,"],",I0)'

   p=where(strpos(tnames,'SNR_STANDARD') ge 0,count)
   if count gt 0 then snr_standard=dbstr.snr_standard else snr_standard=0.0

   ; logsnr_hmag_coef may or may be in structure (depends if a plate is associated with exposure)
   p=where(strpos(tnames,'LOGSNR_HMAG_COEF') ge 0,count)
   if count gt 0 then logsnr=dbstr.logsnr_hmag_coef else logsnr=[0.0,0.0]

   ; snr_goals_version may or may be in structure (depends if a plate is associated with exposure)
   p=where(strpos(tnames,'HMAG_STANDARD_VERSION') ge 0,count)
   if count gt 0 then snr_goals_version=dbstr.hmag_standard_version else snr_goals_version=0

   qr_line = [qr_line, string(exp_pk, qlpk, dbstr.arraydisplay.bzero, dbstr.arraydisplay.bscale, $
      dbstr.arraydisplay.zscale[0], dbstr.arraydisplay.zscale[1], dbstr.dithpix, snr_standard, logsnr[0], $
      logsnr[1], snr_goals_version, $
      FORMAT=qr_format)]
   
   printf,loglun, qr_line,format='(A)'

   ; load the QUICKRED_IMBINZOOM
   qri_line = ['TABLE: QUICKRED_IMBINZOOM']
   qri_line = [qri_line,'ylo, yhi, bzero, bscale, zscale1, zscale2']
   qri_format = '(2(I0,","),3(F0.8,","),F0.8)'

   for zpos=0, n_elements(dbstr.arraydisplay_sub)-1 do begin
      qri_line = [qri_line, string(dbstr.arraydisplay_sub[zpos].yrange[0], dbstr.arraydisplay_sub[zpos].yrange[1], $
          dbstr.arraydisplay_sub[zpos].bzero, dbstr.arraydisplay_sub[zpos].bscale, $
          dbstr.arraydisplay_sub[zpos].zscale[0], dbstr.arraydisplay_sub[zpos].zscale[1], $
          FORMAT=qri_format)]
   endfor
   printf,loglun, ''
   printf,loglun, qri_line,format='(A)'

   ; load the QUICKRED_SPECTRUM
   p=where(strpos(tnames,'QR_SPECTRUM') ge 0,count)
   if count gt 0 then begin
     qrs_line = ['TABLE: QUICKRED_SPECTRUM']
     qrs_line = [qrs_line,'fiberid, bzero, bscale, medsnr']
     qrs_format = '(I0,",",2(F0.8,","),F0.8)'

     for fid=n_elements(dbstr.qr_spectrum)-1,0,-1 do begin
        qrs_line = [qrs_line,string(dbstr.qr_spectrum[fid].fiberid, dbstr.qr_spectrum[fid].bzero, $
            dbstr.qr_spectrum[fid].bscale, dbstr.qr_spectrum[fid].medsnr, $
            FORMAT=qrs_format)]
     endfor

     printf,loglun, ''
     printf,loglun, qrs_line,format='(A)'
   endif

   printf,loglun, ''
   printf,loglun, ''
   close,loglun

end
