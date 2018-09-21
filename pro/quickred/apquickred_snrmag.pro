;+
;
; APQUICKRED_SNRMAG
;
; The calculates the S/N for the fiducial magnitude
; for the quick reduction spectra.
;
; What information to return:
;  -SNR and Hmag per object fiber
;  -fit to SNR vs. magnitude
;  -coefficients for fitted SNR at H=11.5
;
; INPUTS:
;  dbstr     The database structure created by apquickred.pro
;  =plugmap  The plugmap structure
;  /silent   Don't print anything to the screen.
;
; OUTPUTS:
;  The S/N values are updated in the STR structure.
;  =error   The error message if one occurred.
;
; USAGE:
;  IDL>apquickred_snrmag,str
;
; By D.Nidever  2010
;-
pro apquickred_snrmag,dbstr,plugmap=plugmap,silent=silent,error=error

   ; Error Handling
   ;------------------
   ; Establish error handler. When errors occur, the index of the
   ; error is returned in the variable Error_status:  
   ;CATCH, Error_status 

   ;This statement begins the error handler:  
   ;if (Error_status ne 0) then begin 
   ;   error = !ERROR_STATE.MSG  
   ;   if not keyword_set(silent) then print,error
   ;   CATCH, /CANCEL 
   ;   return
   ;endif

   ; Not enough inputs
   if n_elements(dbstr) eq 0 then begin
     error = 'Not enough inputs'
     if not keyword_set(silent) then $
       print,'Syntax - apquickred_snrmag,dbstr,plugmap=plugmap,debug=debug,silent=silent,error=error'
     return
   endif

   timestep = 10.6  ; seconds
   npix = 2048L
   nfibers = 300L
   ;nfibers = n_elements(plugmap.fiberdata)

   apogeefibers = where(plugmap.fiberdata.fiberid ge 0 AND $
                        plugmap.fiberdata.holetype eq 'OBJECT' AND $
                        plugmap.fiberdata.spectrographid eq 2,napogeefibers)
   
   ; Calculate the median S/N in the middle chip
   flux = (dbstr.frame).(1).flux
   err = (dbstr.frame).(1).err

   ; OBJECT Exposures
   if (dbstr.exptype eq 'OBJECT') then begin

     ; No apogee fibers
     if napogeefibers eq 0 then begin
       error = 'APQUICKRED_SNRMAG: No APOGEE fibers in this plugmap'
       if not keyword_set(silent) then print,error
       return
     endif

     ; Get the apogee fibers
     apogee_fibers = where(plugmap.fiberdata.spectrographid eq 2 and $
                    plugmap.fiberdata.holetype eq 'OBJECT' AND $
                    plugmap.fiberdata.fiberid ge 1,napogee_fibers)
     if napogee_fibers eq 0 then begin
       print,'No apogee fibers for this plugmap'
       dbstr.snr_standard = -1
       return
     endif

     fiberid = plugmap.fiberdata[apogee_fibers].fiberid ; 1-300
     ; fiberid=1 is at the top of the detector or index=299
     ; index = 300-fiberid
     fiberind = 300-fiberid
     hmag = plugmap.fiberdata[apogee_fibers].mag[1] ; J, H, Ks
     objtype = plugmap.fiberdata[apogee_fibers].objtype
     ;snr_match = snr[fiberind]   ; S/N matched to the fiber data
     dbstr.hmag[fiberind] = hmag
     dbstr.objtype[fiberind] = objtype  ; we want objtype for ALL apogee fibers

     ; Get non-sky fibers
     objind = where(dbstr.objtype ne '' and dbstr.objtype ne 'SKY',nobjind)
     if nobjind eq 0 then begin
       print,'No object fibers'
       dbstr.snr_standard = -1
       return
     endif

    ;Get sky fibers
     skyind = where(dbstr.objtype eq 'SKY',nskyind)
     obs=fltarr(300)
     for j=0,nfibers-1 do begin
        ;for each fiber, get an observed mag from a median value
        obs[j]=median(flux[*,j])
     endfor
     medsky=median(obs[skyind])
     ;perform sky subtraction whilst calculating signal-to-noise
     snr = MEDIAN((flux-medsky)/(err>1),dim=1) 
     dbstr.snr = snr

     obj_hmag = dbstr.hmag[objind]
     obj_snr = dbstr.snr[objind]

     snstars=where(obj_hmag gt 12 and obj_hmag lt 12.2,nsn)
     scale=1
     if nsn lt 3 then begin
        bright=where(obj_hmag lt 12)
        hmax=max(obj_hmag[bright])
        snstars=where(obj_hmag gt hmax-0.2 and obj_hmag le hmax,nsn)
        scale=sqrt(10^(0.4*(hmax-12.2)))
     endif
     achievedsn=median(obj_snr[snstars],dimension=1)*scale
     dbstr.snr_standard = achievedsn

;     No longer used
;     dbstr.logsnr_hmag_coef = reform(icoef)   


   ; NON-OBJECT Exposures
   ;----------------------
   ; flat, dome
   Endif else begin

     ; For non-object exposures SNR_STANDARD means SNR_MEDIAN
     dbstr.snr_standard = median([snr])

   Endelse


  ;stop

end
