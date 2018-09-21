;+
;
; APQL_DITHER
;
; This checks the dither position for the quicklook software
;
; INPUTS:
;  str      The quicklook structure.
;  allstr   The quicklook structure of all previous reads of this exposure.
;  prevstr  The structure of all previous exposures for this plate visit.
;  /silent  Don't print anything to the screen.
;
; OUTPUTS:
;  The dither values are updated in the STR structure.
;  =error   The error message if one occurred.
;
; USAGE:
;  IDL>apl_dither,str,allstr,prevstr
;
; By D.Nidever  2010
;-
pro apql_dither,str,allstr,prevstr,silent=silent,error=error


; Initialize all values to BAD
str.dither_prevexp_measured = 99.00
str.dither_prevexp_header = 99.00
str.dither_relative = 99.00
str.dither_status = 1


; Error Handling
;------------------
; Establish error handler. When errors occur, the index of the
; error is returned in the variable Error_status:  
; CATCH, Error_status 

;This statement begins the error handler:  
; if (Error_status ne 0) then begin 
;    error = !ERROR_STATE.MSG  
;    if not keyword_set(silent) then print,error
;    CATCH, /CANCEL 
;    return
; endif

; Not enough inputs
if n_elements(str) eq 0 then begin
  error = 'Not enough inputs'
  if not keyword_set(silent) then $
    print,'Syntax - apql_dither,str,allstr,prevstr,silent=silent,error=error'
  return
endif

; make sure the system variable exists before we proceed
DEFSYSV,'!apql',exists=apql_exists
if not apql_exists then begin
  error = 'APQL_DITHER: system variable !apql is not defined'
  if not keyword_set(silent) then print,error
  return
endif

; Check that plugmap exists
if not PTR_VALID(!apql.plugmap.datastr) then begin
  error = 'APQL_DITHER: no plugmap data available'
  if not keyword_set(silent) then print,error
  return
endif


npix = 2048L

frameid = str.frameid
readnum = str.readnum

nprevstr = n_elements(prevstr)

; Exposure types
CASE str.exptype of

  ; OBJECT Exposure
  ;----------------
  'OBJECT': begin

    ; Checking previous exposures, same PLATE and both OBJECT
    if nprevstr gt 0 then begin

      ; Same PLATE and both OBJECT
      if prevstr[nprevstr-1].exptype eq 'OBJECT' and prevstr[nprevstr-1].plate eq str.plate then begin

        ; Measuring dither shift with median sky spectra
        medflux = fltarr(3*npix)
        for i=0,2 do medflux[i*npix:i*npix+npix-1]=str.medsky[i].flux
        ; Use the LAST previous exposure
        nprevstr = n_elements(prevstr)
        medflux1 = fltarr(3*npix)
        for i=0,2 do medflux1[i*npix:i*npix+npix-1]=prevstr[nprevstr-1].medsky[i].flux
        ;medflux1 = prevstr[nprevstr-1].medskyflux

        ; Cross-correlate them
        XCORLB,medflux,medflux1,5,lastsh,chisq,sherr,error=error
        if n_elements(error) gt 0 then lastsh=99.00

        ;print,'Measured Dither Shift relative to ',file_basename(lastfile,'.fits')
        ;print,'Dither Shift = ',stringize(lastsh,ndec=3),' pixels'

        str.dither_prevexp_measured = lastsh

        ; Using header dither position values
        ;-------------------------------------
        head = headfits(str.filename)
        ;head1 = headfits(lastfile)
        head1 = prevstr[nprevstr-1].header

        ; Header dither positions
        ;ditherkey = 'DITHER'
        ditherkey = 'DITHPIX'
        lastframe = prevstr[nprevstr-1].frameid
        lastread = prevstr[nprevstr-1].readnum
        ditherpos1 = sxpar(head1,ditherkey)

        ;print,'Header dither position for ',lastframe,': ',stringize(ditherpos1,ndec=3),' pixels'
        ditherpos = sxpar(head,ditherkey)
        ;print,'Header dither position for ',frameid,': ',stringize(ditherpos,ndec=3),' pixels'

        ; The header values are in microns,  one pixel is 18 microns
        ;header_ditherunit_to_pixels = 1.0/18.0
        ;dithershift_header = (float(ditherpos)-float(ditherpos1)) * header_ditherunit_to_pixels
        ; Header values are in PIXELS
        dithershift_header = float(ditherpos)-float(ditherpos1)
        ; if the detector shifts to the right, the spectrum shifts the
        ; opposite way, to the left.

        str.dither_prevexp_header = dithershift_header

      ; Previous exposure of different PLATE or not OBJECT
      endif else begin
        str.dither_prevexp_measured = 99.00
        str.dither_prevexp_header = 99.00
      endelse


    ; First exposure, no prevstr
    endif else begin
      str.dither_prevexp_measured = 99.00
      str.dither_prevexp_header = 99.00
      ; Load the median skyline spectrum
      ;medlinestr = (*str.medskylinestr)
      medflux = fltarr(3*npix)
      for i=0,2 do medflux[i*npix:i*npix+npix-1]=str.medsky[i].flux
    endelse


    ; Are there any variations within the exposure
    ;----------------------------------------------
    nallstr = n_elements(allstr)
    if nallstr gt 0 then begin

      ; We have sky lines
      if PTR_VALID(str.medskylinestr) and PTR_VALID(allstr[0].medskylinestr) then begin

        ; Load sky lines
        medlinestr = (*str.medskylinestr)

        ; Compare to the 2nd read

        readnumber = allstr[0].readnum
        imedlinestr = (*allstr[0].medskylinestr)
        imedflux = fltarr(3*npix)
        for j=0,2 do imedflux[j*npix:j*npix+npix-1]=allstr[0].medsky[j].flux

        ; Match to current list
        SRCOR2,imedlinestr.fiber,imedlinestr.gaussx,medlinestr.fiber,medlinestr.gaussx,0.3,ind1,ind2,opt=1,/silent
        dum = where(ind1 ne -1,nmatch)
        if nmatch eq 0 then $
           SRCOR2,imedlinestr.fiber,imedlinestr.gaussx,medlinestr.fiber,medlinestr.gaussx,0.8,ind1,ind2,opt=1,/silent
        dum = where(ind1 ne -1,nmatch)
        if nmatch eq 0 then $
          SRCOR2,imedlinestr.fiber,imedlinestr.gaussx,medlinestr.fiber,medlinestr.gaussx,1.3,ind1,ind2,opt=1,/silent
        dum = where(ind1 ne -1,nmatch)
        if nmatch eq 0 then $
          SRCOR2,imedlinestr.fiber,imedlinestr.gaussx,medlinestr.fiber,medlinestr.gaussx,1.8,ind1,ind2,opt=1,/silent
        dum = where(ind1 ne -1,nmatch)

        ; Measure median shift of lines
        if nmatch gt 1 then begin
           medshift = median([medlinestr[ind2].gaussx-imedlinestr[ind1].gaussx])
           ;readstr[i].medshift = medshift
        endif else begin
          medshift = 99.0
        endelse

        XCORLB,imedflux,medflux,5,sh,error=error
        if n_elements(error) gt 0 then medshift=99.00
        ;readstr[i].xmedshift = sh
        str.dither_relative = medshift

      endif

      ;stop

    end ; nallstr gt 0


    ; DITHER STATUS
    ;--------------
    ; 0-good, 1-bad

    ; Check the dither values
    ;-------------------------
    str.dither_status = 0  ; okay for now, 
    ; Measured dither shift wrt previous exposure
    if abs(abs(str.dither_prevexp_measured)-0.5) gt 0.1 then str.dither_status = 1

    ; Header dither position wrt previous exposure
    if abs(abs(str.dither_prevexp_header)-0.5) gt 0.1 then str.dither_status = 1

    ; Header vs. Measured dither shift
    dither_diff = abs(str.dither_prevexp_measured - str.dither_prevexp_header)
    if dither_diff gt 0.1 then str.dither_status = 1

    ; Relative positions
    if nallstr gt 0 then begin
      maxreldither = max(abs(str.dither_relative))
      if maxreldither gt 0.03 then str.dither_status = 1
    endif

  end ; exptye=object


  ; ARCLAMP Exposures
  ;-------------------
  'ARCLAMP': begin

    ; Checking previous exposures
    if n_elements(prevstr) gt 0 then begin

      ; Both ARCLAMP
      if prevstr[nprevstr-1].exptype eq 'ARCLAMP' then begin

        ; Check the headers
        ;--------------------
        head = headfits(str.filename)
        head1 = prevstr[nprevstr-1].header

        ; Header dither positions
        ditherkey = 'DITHPIX'
        lastframe = prevstr[nprevstr-1].frameid
        lastread = prevstr[nprevstr-1].readnum
        ditherpos1 = sxpar(head1,ditherkey)

        ditherpos = sxpar(head,ditherkey)

        ; The header values are in microns,  one pixel is 18 microns
        ;header_ditherunit_to_pixels = 1.0/18.0
        ;dithershift_header = (float(ditherpos)-float(ditherpos1)) * header_ditherunit_to_pixels
        ; Header values are in PIXELS
        dithershift_header = float(ditherpos)-float(ditherpos1)
        ; if the detector shifts to the right, the spectrum shifts the
        ; opposite way, to the left.

        str.dither_prevexp_header = dithershift_header

      endif ; both arclamp

    ; First exposure
    endif else begin
      str.dither_prevexp_measured = 99.00
      str.dither_prevexp_header = 99.00
    endelse


    ; DITHER STATUS
    ;--------------
    ; 0-good, 1-bad

    ; Check the dither values
    ;-------------------------
    str.dither_status = 0  ; okay for now, 

    ; Header dither position wrt previous exposure
    if abs(abs(str.dither_prevexp_header)-0.5) gt 0.1 then str.dither_status = 1

  end ; exptype=arclamp

  ; ASDAF
  'ASDAF': begin

    ; Check that the previous exposure was ASDAF and had the same ra/dec
    ; This is already checked in apql_wrapper.pro

    ; Checking previous exposures
    if n_elements(prevstr) gt 0 then begin

      ; Check the headers
      ;--------------------
      head = headfits(str.filename)
      head1 = prevstr[nprevstr-1].header

      ; Header dither positions
      ditherkey = 'DITHPIX'
      lastframe = prevstr[nprevstr-1].frameid
      lastread = prevstr[nprevstr-1].readnum
      ditherpos1 = sxpar(head1,ditherkey)

      ditherpos = sxpar(head,ditherkey)

      ; The header values are in microns,  one pixel is 18 microns
      ;header_ditherunit_to_pixels = 1.0/18.0
      ;dithershift_header = (float(ditherpos)-float(ditherpos1)) * header_ditherunit_to_pixels
      ; Header values are in PIXELS
      dithershift_header = float(ditherpos)-float(ditherpos1)
      ; if the detector shifts to the right, the spectrum shifts the
      ; opposite way, to the left.

      str.dither_prevexp_header = dithershift_header

    endif


    ; DITHER STATUS
    ;--------------
    ; 0-good, 1-bad

    ; Check the dither values
    ;-------------------------
    str.dither_status = 0  ; okay for now, 

    ; Header dither position wrt previous exposure
    if abs(abs(str.dither_prevexp_header)-0.5) gt 0.1 then str.dither_status = 1

  end ; exptype=asdaf

  ; Dome and Flat won't be dithered
  else: begin
    ; Set dither status to ok
    str.dither_status = 0
  end

ENDCASE


;stop

end
