;+
;
; APQL_PLOTSPEC
;
; This prepares some representative spectra for the APOGEE
; Quicklook software.
;
; INPUTS:
;  str       The quicklook structure.
;  /silent   Don't print anything to the screen
;
; OUTPUTS:
;  The representative spectra are put into the STR structure.
;  =error    The error message if one occurred.
;
; USAGE:
;  IDL>apql_plotspec,str
;
; By D. Nidever  2010
;-
pro apql_plotspec,str,silent=silent,error=error

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
   if n_elements(str) eq 0 then begin
     error = 'Not enough inputs'
     if not keyword_set(silent) then print,'Syntax - apql_plotspec,str,silent=silent,error=error'
     return
   endif

   ; Check that the extracted spectra are there
   if not PTR_VALID(str.frame) then begin
     error = 'APQL_PLOTSPEC: NO extracted spectra'
     if not keyword_set(silent) then print,error
     return
   endif

   ; Check that there are APOGEE fibers in the plugmap file
   napogeefibers = 0
   if PTR_VALID(!apql.plugmap.datastr) then $
     apogeefibers = where( (*!apql.plugmap.datastr).fiberdata.SPECTROGRAPHID eq 2 AND $
                           (*!apql.plugmap.datastr).fiberdata.fiberId ge 1,napogeefibers)   

   ; Default parameters
   npix = 2048L
   chipgap = 150 ;300
   ;npick = 3
   nspec = n_elements(str.representative_spectra)
   npick = nspec/3

   nfibers = n_elements( ((*str.frame)[0]).flux[0,*])


   ; OBJECT Exposure
   ;-----------------
   If str.exptype eq 'OBJECT' and napogeefibers gt 0 then begin

     if not PTR_VALID(!apql.plugmap.datastr) then begin
       ; skipping papql_plotspec since we have no fiber information available
       error = 'APQL_PLOTSPEC: No valid plugmap available'
       if not keyword_set(silent) then print,error
       return
     endif

     ; Object spectra
     ;-----------------
     ;objind = where(plugmap.fiberdata.objtype eq 'STAR',nobjind)
     ;objind = where((*!apql.plugmap.datastr).fiberdata.objtype eq 'STAR',nobjind)
     ;objind = where((*!apql.plugmap.datastr).fiberdata.objtype eq 'SCIENCE',nobjind)
     objind = where(((*!apql.plugmap.datastr).fiberdata.objtype eq 'STAR' or $
                    (*!apql.plugmap.datastr).fiberdata.objtype eq 'EXTOBJ') AND $
                    (*!apql.plugmap.datastr).fiberdata.holetype eq 'OBJECT' AND $
                    (*!apql.plugmap.datastr).fiberdata.fiberid ne -1 AND $
                    (*!apql.plugmap.datastr).fiberdata.spectrographid eq 2,nobjind)
     ; Pick three spectra at random
     if nobjind gt 0 then begin
       objrnd = randomu(seed,nobjind)
       siobjrnd = sort(objrnd)
       objindrnd = objind[siobjrnd[0:(npick-1)<(nobjind-1)]]
       push,allind,objindrnd
     endif else begin
       print,'No STAR or EXTOBJ found in this plugmap for exptype=',str.exptype
     endelse

     ; Telluric spectra
     ;telind = where(plugmap.fiberdata.objtype eq 'HOT_STD',ntelind)
     telind = where((*!apql.plugmap.datastr).fiberdata.objtype eq 'HOT_STD' AND $
        (*!apql.plugmap.datastr).fiberdata.fiberid ne -1 AND $
        (*!apql.plugmap.datastr).fiberdata.holetype eq 'OBJECT' AND $
        (*!apql.plugmap.datastr).fiberdata.spectrographid eq 2,ntelind)
     if ntelind gt 0 then begin
       telrnd = randomu(seed,ntelind)
       sitelrnd = sort(telrnd)
       telindrnd = telind[sitelrnd[0:(npick-1)<(ntelind-1)]]
       push,allind,telindrnd
     endif else begin
       print,'No HOT_STD found in this plugmap for exptype=',str.exptype
     endelse

     ; Sky spectra
     ;skyind = where(plugmap.fiberdata.objtype eq 'SKY',nskyind)
     skyind = where((*!apql.plugmap.datastr).fiberdata.objtype eq 'SKY' AND $
        (*!apql.plugmap.datastr).fiberdata.fiberid ne -1 AND $
        (*!apql.plugmap.datastr).fiberdata.holetype eq 'OBJECT' AND $
        (*!apql.plugmap.datastr).fiberdata.spectrographid eq 2,nskyind)
     if nskyind gt 0 then begin
       skyrnd = randomu(seed,nskyind)
       siskyrnd = sort(skyrnd)
       skyindrnd = skyind[siskyrnd[0:(npick-1)<(nskyind-1)]]
       push,allind,skyindrnd
     endif else begin
       print,'No SKY found in this plugmap for exptype=',str.exptype
     endelse

   ; NON-OBJECT Exposure
   ;------------------------
   ; flat, comp or dome
   Endif Else Begin

     ; Pick fibers evenly spaced
     sep = float(nfibers)/nspec
     allind = round( indgen(nspec)*sep+sep*0.5 )  ; these are the indices for FLUX

     ; if ASDAF pick the brightest one first
     if str.exptype eq 'ASDAF' then begin
       medflux = median((*str.frame)[1].flux,dim=1)
       maxind = where(medflux eq max(medflux))
       allind[0] = maxind
     endif

   Endelse


   ; Loop through the spectra
   for i=0,n_elements(allind)-1 do begin

     if str.exptype eq 'OBJECT' and napogeefibers gt 0 then begin

       ind = allind[i]
       fiberid = (*!apql.plugmap.datastr).fiberdata[ind].fiberid
       fiberind = 300-fiberid

       ;fiberid = (*!apql.plugmap.datastr).fiberdata[ind].fiberid
       ;fiberid = plugmap.fiberdata[ind].fiberid
       ;fiberind = fiberid-1
       ; fiberid=1 is at the top of the detector or index=299
       ;index = 300-fiberid
       ;fiberind = 300-fiberid

       flux = fltarr(3*npix)
       ;err = fltarr(3*npix)
       ;for j=0,2 do flux[j*npix:j*npix+npix-1] = frame.(j).flux[*,fiberind]
       ;for j=0,2 do err[j*npix:j*npix+npix-1] = frame.(j).err[*,fiberind]
       ;mag = plugmap.fiberdata[ind].mag  ; J, H, Ks
help,flux
help,npix
help,(*str.frame),/str
help,fiberind
       for j=0,2 do flux[j*npix:j*npix+npix-1] = (*str.frame)[j].flux[*,fiberind]
       ;for j=0,2 do err[j*npix:j*npix+npix-1] = str.frame[j].err[*,fiberind]
       mag = (*!apql.plugmap.datastr).fiberdata[ind].mag  ; J, H, Ks
       mag = mag[0:2]  ; in case there are extra values
       hmag = mag[1]
       jk = mag[0]-mag[2]

       ;objtype = plugmap.fiberdata[ind].objtype
       objtype = (*!apql.plugmap.datastr).fiberdata[ind].objtype
       case objtype of
        'STAR': spectype='Object'
        'EXTOBJ': spectype='Object'
        ;'SCIENCE': spectype='Object'
        'SKY': spectype='Sky'
        'HOT_STD': spectype='Telluric'
        else: spectype=objtype
       endcase

     ; Non-object exposure
     endif else begin

       fiberind = allind[i]
       ;fiberind = ind
       ;fiberid = 300-fiberind
       fiberid = 300-fiberind
       flux = fltarr(3*npix)
       for j=0,2 do flux[j*npix:j*npix+npix-1] = (*str.frame)[j].flux[*,fiberind]
       spectype = str.exptype
       mag = 99.99
       objtype = spectype

     endelse

     ; Calculate the median S/N in the middle chip
     flux1 = (*str.frame)[1].flux[*,fiberind] > 0.0
     err1 = (*str.frame)[1].err[*,fiberind] > 1.0
     medsnr = MEDIAN(flux1/err1,dim=1)

     ;midx = n_elements(flux)/2
     maxy = max(flux)
     if spectype eq 'Object' or spectype eq 'Telluric' then maxy=max(median(flux,11))*1.3
     ;yr = [0,max(flux)*1.05]
     yr = [0,maxy*1.15]

     ; Scale the data
     minflux = min(flux)
     maxflux = max(flux)
     ; real_values = byte_values * bscale + bzero
     bzero = minflux
     ; UINT goes from 0-65535
     ; integer 
     bscale = (maxflux-minflux)/65535. > 1  ; don't scale up
     uint_flux = ( flux - bzero )/bscale
     uint_flux = uint( round( uint_flux ) ) ; round and make uint type
     ; stuff into the structure
     if PTR_VALID(str.representative_spectra[i].data) then PTR_FREE,str.representative_spectra[i].data
     str.representative_spectra[i].data = PTR_NEW(uint_flux,/no_copy)
     str.representative_spectra[i].bscale = bscale
     str.representative_spectra[i].bzero = bzero
     str.representative_spectra[i].yrange = yr
     str.representative_spectra[i].objtype = objtype
     str.representative_spectra[i].fiberid = fiberid
     str.representative_spectra[i].mag = mag
     str.representative_spectra[i].medsnr = medsnr

   endfor

end
