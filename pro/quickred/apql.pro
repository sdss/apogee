;+
;
; APQL
;
; APOGEE Quicklook software. This program checks exposures in progress
;
; INPUTS:
;  filename          The Raw APOGEE image filename
;  allstr            All of the STR structures for this EXPOSURE
;  prevstr           Information from previous exposures
;  =snr_goals        A structure giving the S/N goals, normally
;                      obtained from the apogeeql database.
;  =fitskw_err       A structure giving the required FITS keyword
;                      error codes.  This is normally obtained from
;                      the apogeeql database.
;  =required_fitskw  A structure giving the required FITS keywords.
;                      This is normally obtained from the apogeeql database.
;  /silent           Don't print anything to the screen.
;  /stp              Stop at the end of the program.
;
; OUTPUTS:
;  str               The output structure
;  =predict_str      The prediction structure.  This is only
;                      output every 6th read.
;  =error            The error message if one occurred.
;
; USAGE:
;  IDL>apql,'/net/stream/apogee/data/raw/apRaw-00000118-001.fits',str,allstr,prevstr,predict_str=predict_str
;
; By D.Nidever  Aug 2010
;-
pro apql,filename,str,allstr,prevstr,obs=obs,predict_str=predict_str, snr_goals=snr_goals, $
   fitskw_err=fitskw_err, required_fitskw=required_fitskw, error=error, silent=silent,$
   stp=stp

   ; setup	SDSS database parameters
   SDSS_DB_PARAMS,obs=obs

   t0 = systime(1)

   ; Not enough inputs
   if n_elements(filename) eq 0 then begin
     error = 'Not enough inputs'
     if not keyword_set(silent) then begin
       print,'Syntax - apql,filename,str,allstr,prevstr,predict_str=predict_str, snr_goals=snr_goals,'
       print,'              fitskw_err=fitskw_err, required_fitskw=required_fitskw, error=error, silent=silent,'
       print,'              stp=stp'
     endif
     return
   endif

   ; Default S/N goal parameters
   if keyword_set(snr_goals) then begin
      hmag_standard = snr_goals.hmag
      snr_standard_goal = snr_goals.snr_goals
      hmag_standard_version= snr_goals.version
   endif else begin
      hmag_standard = 12.0       ; Hmag used for getting S/N
      snr_standard_goal = 30.0   ; S/N goal at standard Hmag
      hmag_standard_version = 0
   endelse

   ; Frameid and readnum
   filedir = FILE_DIRNAME(filename)
   filebase = FILE_BASENAME(filename,'.fits')
   dum = strsplit(filebase,'-',/extract)
   frameid = dum[1]
   readnum = dum[2]

   if not keyword_set(silent) then $
     print,'Processing ',filename

   ; Default parameters
   npix = 2048L
   nfibers = 300L
   nchips = 3L
   nbin = 8
   npixbin = npix/nbin


   ; Load the current read
   ;-----------------------
   if file_test(filename) eq 0 then begin
     error = 'Last read file '+filename+' NOT FOUND'
     if not keyword_set(silent) then print,error
     return
   endif
   FITS_READ,filename,im1,head1,message=message1,/no_abort
   if message1 ne '' then begin
     error = 'Error reading file '+filename
     if not keyword_set(silent) then print,error
     return
   endif
   im1 = long(im1)
   nhead1 = n_elements(head1)

   ; Initialize the structure that will hold all the information on this dither
   str = {filename:filename,frameid:strtrim(frameid,2),readnum:readnum,$
          exptime:0.0,header:strarr(1000),exptype:'',plate:'',date:'',$
          rawimage:PTR_NEW(),cds_image:PTR_NEW(),jd:0.0d0,mjd5:'',$
          representative_spectra:REPLICATE({data:PTR_NEW(),bscale:0.0,bzero:0.0,yrange:[0.0,0.0],$
          objtype:'',fiberid:0,mag:fltarr(3),medsnr:0.0},9),$
          arraydisplay_nbin:nbin,arraydisplay:{data:bytarr(npixbin,npixbin*nchips),bscale:0.0,bzero:0.0,zscale:[0.0,0.0]},$
          arraydisplay_sub:REPLICATE({data:PTR_NEW(),bscale:0.0,bzero:0.0,zscale:[0.0,0.0],yrange:[0L,0L]},3),$
          snr_goals_version:hmag_standard_version, hmag_standard:hmag_standard,snr_standard_goal:snr_standard_goal,$
          snr_standard:0.0,delta_snr2_standard:0.0,logsnr_hmag_coef:fltarr(2),logsnr_hmag_coef_goal:fltarr(2),$
          snr2_time_coef:fltarr(2),expected_total_readnum:0.0,snr2_time_recent_coef:fltarr(2),$
          expected_total_readnum_recent:0.0,wavefit_pars:dblarr(15),wavefit_rms:0.0,$
          waverange_exp:dblarr(nchips,2),waverange_meas:dblarr(nchips,2),waverange_diff:dblarr(nchips,2),$
          wavelength_status:-1,skyvar_meddev_perc:0.0,skyvar_stddev_perc:0.0,skyvar_contflux:0.0,$
          skyvar_contflux_rate:0.0,skyvar_avglineflux:0.0,skyvar_avglineflux_rate:0.0,sky_status:-1,$
          dither_prevexp_measured:99.0,dither_prevexp_header:99.0,dither_relative:0.0,$
          dither_status:-1,fitsheader_errors:PTR_NEW(), fitsheader_status:-1,$
          exp_finished_status:-1,visit_finished_status:-1,medflux:fltarr(nchips,nfibers),frame:PTR_NEW(),$
          medsnr:fltarr(nfibers),hmag:fltarr(nfibers),objtype:strarr(nfibers),$
          medsky:REPLICATE({flux:fltarr(npix),err:fltarr(npix),mask:lonarr(npix)},nchips),$
          medskylinestr:PTR_NEW(),skyfiberstr:PTR_NEW(),skylinestr:PTR_NEW(),$
          nmissingfibers:0L,missingfibers:PTR_NEW()}

   ; object, flat, dark, sky, calib (comps), localflat, superdark, superflat
   exptype = strtrim(strupcase(sxpar(head1,'EXPTYPE')),2)


   ; Check if this is an Any-Star-Down-Any-Fiber (ASDAF) exposure
   ;  ASDAF if it's object and ra/dec coordinates do match those in plugmap
   if exptype eq 'OBJECT' and PTR_VALID(!apql.plugmap.datastr) then begin
     rahd = sxpar(head1,'RA',count=nra)
     dechd = sxpar(head1,'DEC',count=ndec)
     rapl = (*!apql.plugmap.datastr).racen
     decpl = (*!apql.plugmap.datastr).deccen
     ; we have valid coordinates
     if nra gt 0 and ndec gt 0 and min(valid_num([rahd,dechd,rapl,decpl])) eq 1 then begin
       dist = sphdist(rahd,dechd,rapl,decpl,/deg)
       if dist gt 2 then exptype='ASDAF'
       ;print,'This is an ASDAF exposure'
     endif
   endif

   ; Check that there are APOGEE fibers in the plugmap file
   napogeefibers = 0
   if PTR_VALID(!apql.plugmap.datastr) then begin
     apogeefibers = where( (*!apql.plugmap.datastr).fiberdata.SPECTROGRAPHID eq 2 AND $
                           (*!apql.plugmap.datastr).fiberdata.holetype eq 'OBJECT' AND $
                           (*!apql.plugmap.datastr).fiberdata.fiberId ge 1,napogeefibers)
     if napogeefibers eq 0 then print,'NO APOGEE fibers in this plugmap.'
   endif


   ; print,'APQL  -> EXPTYPE =',exptype

   ; Don't extract if it's a dark
   doextract = 0
   if exptype eq 'OBJECT' or exptype eq 'QUARTZFLAT' or exptype eq 'DOMEFLAT' or $
      exptype eq 'ARCLAMP' or exptype eq 'ASDAF' then doextract=1

   ; Change S/N goals for different exposure types
   if exptype eq 'FLAT' or exptype eq 'SUPERFLAT' then str.snr_standard_goal = 140.
   if exptype eq 'LOCALFLAT' then str.snr_standard_goal = 15.


   ; Load the current read information into the structure
   str.rawimage = PTR_NEW(im1)
   str.header[0:nhead1-1] = head1
   ;str.exptime = float(sxpar(head1,'EXPTIME'))
   str.exptime = long(readnum)*!apql.timePerRead
   str.exptype = exptype
   str.plate = strtrim(sxpar(head1,'PLATEID'),2)
   str.date = sxpar(head1,'DATE-OBS')
   str.jd = date2jd(str.date)
   str.mjd5 = getmjd5(str.header,/sdss)
   if keyword_set(doextract) then $
     str.frame = PTR_NEW(REPLICATE({flux:lonarr(npix,nfibers),err:lonarr(npix,nfibers),mask:intarr(npix,nfibers)},nchips))

   ; Check that we have the first read data
   if not PTR_VALID(!apql.firstread.data) then begin
     error = 'APQL: No first read data'
     if not keyword_set(silent) then print,error
     return
   endif

   ; Do double-correlated sampling
   ;-------------------------------
   str.cds_image = PTR_NEW( long( im1 - (*!apql.firstread.data) ) )

   ; Extract the spectra
   ;---------------------
   if keyword_set(doextract) then begin

     ; Loop through the chips
     FOR i=0,nchips-1 DO BEGIN

       ; Chip portion of 2D array
       chipim = (*str.cds_image)[i*npix:i*npix+npix-1,*]
       chiperr = long( sqrt(chipim>1) )  ; should we add rdnoise?
       mask = fix(chipim*0)
       chipstr = {flux:chipim, err:chiperr, mask:mask}

       ; Get the trace structure
       tracestr = *(!apql.psf.data[i].tracestr)
       psfim = *(!apql.psf.data[i].psfim)

       ; We only need all three chips for the sky fibers.
       ; For the object fibers we just need the middle chip.

       ; OBJECT Exposure
       ;-----------------------
       if (exptype eq 'OBJECT') then begin

         ; "Regular" read
         ;----------------
         ;   extract all fiber with BOXCAR and sky fibers
         ;   with PSF
         if (long(readnum) mod 6 ne 0) then begin

           ; Boxcar extract all fibers
           if i eq 1 then begin
              print, 'APEXTRACT...'
              APEXTRACT,chipstr,tracestr,outstr
           endif else begin
             ntraces = n_elements(tracestr)
             outstr = {flux:lonarr(npix,ntraces),err:lonarr(npix,ntraces),mask:intarr(npix,ntraces)}
           endelse

           ; Extract SKY fibers
           ;  If we have many sky fibers then we should just use a subsample
           if PTR_VALID(!apql.plugmap.datastr) then begin
              ; only keep SKY for the APOGEE instrument with a valid fiberId
              skyind = where((*!apql.plugmap.datastr).fiberdata.objtype eq 'SKY' AND $
                             (*!apql.plugmap.datastr).fiberdata.SPECTROGRAPHID eq 2 AND $
                             (*!apql.plugmap.datastr).fiberdata.holetype eq 'OBJECT' AND $
                             (*!apql.plugmap.datastr).fiberdata.fiberId ge 1,nskyind)

              if nskyind gt 0 then begin
                 skyfiberind = (*!apql.plugmap.datastr).fiberdata[skyind].fiberid-1
                 ;APEXTRACTPSF,chipstr,tracestr,psfim,skyoutstr,fibers=skyfiberind
                 APEXTRACT,chipstr,tracestr,skyoutstr,fibers=skyfiberind
     
                 ; Stick the PSF extracted sky fibers into OUTSTR
                 outstr.flux[*,skyfiberind] = skyoutstr.flux
                 outstr.err[*,skyfiberind] = skyoutstr.err
                 outstr.mask[*,skyfiberind] = skyoutstr.mask
          
              endif else begin
                if not keyword_set(silent) then print,'APQL: NO SKY fibers'
              endelse
           endif else begin
             if not keyword_set(silent) then print,'APQL: No plugmap data'
           endelse

         ; 60sec, extract all fibers with PSF
         ;-----------------------------------
         endif else begin

            ; PSF extract all fibers
            ;print, 'APEXTRACTPSF...'
            ;APEXTRACTPSF,chipstr,tracestr,psfim,outstr

            ; Boxcar extract all fibers
            APEXTRACT,chipstr,tracestr,outstr

         endelse

       ; NON-OBJECT Exposures
       ;------------------------
       ; use boxcar for flat, comp and dome
       endif else begin
          print, 'APEXTRACT...'
          APEXTRACT,chipstr,tracestr,outstr
       endelse

       ; Get median fluxed at center of chip
       str.medflux[i,*] = median(outstr.flux[npix/2-50:npix/2+50,*],dim=1)

       ; Save the extracted spectra
       (*str.frame)[i] = outstr[0]

     ENDFOR  ; chip loop

   endif ; doextract=1


   ;---------------------------
   ; RUN ALL OF THE ROUTINES
   ;---------------------------

   ; Representative spectra
   ;-----------------------
   if keyword_set(doextract) and (long(readnum) mod 6 eq 0) then begin
      print,'APQL_PLOTSPEC ...'
      APQL_PLOTSPEC,str
   endif

   ; Array display
   ;--------------
   if (long(readnum) mod 6 eq 0) then begin
     print,'APQL_PLOTARRAYS ...'
     APQL_PLOTARRAYS,str
   endif

   ; S/N vs. Hmag
   ;-------------
   ; object and exposures with continuum
   if (exptype eq 'OBJECT' and napogeefibers gt 0) or (exptype eq 'ASDAF' and napogeefibers gt 0) or $
       exptype eq 'QUARTZFLAT' or exptype eq 'DOMEFLAT' then begin
     print,'APQL_SNRMAG ...'
     APQL_SNRMAG,str,allstr
   endif

   ; Check the wavelength solution, skyline positions
   ;-------------------------------------------------
   if (exptype eq 'OBJECT' or exptype eq 'ASDAF') and (napogeefibers gt 0) then begin
     print,'APQL_WAVESOL ...'
     APQL_WAVESOL,str,allstr
   endif

   ; Check skyline positions
   ;------------------------
   if (exptype eq 'OBJECT' or exptype eq 'ASDAF') and (napogeefibers gt 0) then begin
     print,'APQL_SKYCHECK ...'
     APQL_SKYCHECK,str,allstr
   endif

   ; Sky brightness vs. time
   ;------------------------
   if (exptype eq 'OBJECT' or exptype eq 'ASDAF') and (napogeefibers gt 0) then begin
      print,'APQL_SKYVAR ...'
      APQL_SKYVAR,str,allstr
   endif

   ; Check Dither position
   ;----------------------
   if (exptype eq 'OBJECT' or exptype eq 'ASDAF') and (napogeefibers gt 0) then begin
      print,'APQL_DITHER ...'
      APQL_DITHER,str,allstr,prevstr
   endif

   ; Check FITS header
   ;------------------
   print,'APQL_FITSHEADER ...'
   APQL_FITSHEADER,str, required_fitskw=required_fitskw, fitskw_err=fitskw_err

   ; Check completion status
   ;------------------------
   if (exptype eq 'OBJECT') and (napogeefibers gt 0) then begin
      print,'APQL_STATUSCHECK ...'
      APQL_STATUSCHECK,str,allstr,prevstr,status_str
   endif

   ; Make predictions
   ;-----------------
   ;   only run every 6th read, takes ~0.12 sec to run
   if (exptype eq 'OBJECT') and (long(readnum) mod 6 eq 0) and (napogeefibers gt 0) then begin
      print,'APQL_PREDICT ...'
      APQL_PREDICT,str,allstr,prevstr,predict_str,/summary ;/verbose
   endif

   ; Runtime
   dt = systime(1)-t0
   print,'dt = ',strtrim(dt,2),' sec'


   ; This now takes about 6-8 sec for a normal read and
   ; ~11.5 sec for every 6th read.
   ; The number of threads/cores/cpus used doesn't seem to really
   ; change this!!1

   if keyword_set(stp) then stop

end
