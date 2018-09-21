;+
;
; APQL_SNRMAG
;
; This plots the S/N vs. Hmag for all the fibers
; Quicklook software
;
; What information to return:
;  -SNR and Hmag per object fiber
;  -fit to SNR vs. magnitude
;  -coefficients for fitted SNR at H=11.5
;  -coefficients for SNR "goal"
;  -(S/N)^2 and read number for all reads
;
; INPUTS:
;  str      The quicklook structure.
;  allstr   The quicklook structure of all previous reads of this exposure.
;  /debug   Makes some diagnostic plots.
;  /silent  Don't print anything to the screen.
;
; OUTPUTS:
;  The S/N values are updated in the STR structure.
;  =error   The error message if one occurred.
;
; USAGE:
;  IDL>apl_snrmag,str,allstr
;
; By D.Nidever  2010
;-
pro apql_snrmag,str,allstr,debug=debug,silent=silent,error=error

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

   ; Initialize all values to BAD
   str.medsnr = -1
   str.snr_standard = -1
   str.logsnr_hmag_coef = 0
   str.logsnr_hmag_coef_goal = 0
   str.snr2_time_coef = 0
   str.expected_total_readnum = -1
   str.delta_snr2_standard = -999
   str.snr2_time_recent_coef = 0
   str.expected_total_readnum_recent = -1

   ; Not enough inputs
   if n_elements(str) eq 0 then begin
     error = 'Not enough inputs'
     if not keyword_set(silent) then print,'Syntax - apql_snrmag,str,allstr,debug=debug,silent=silent,error=error'
     return
   endif

   ; Check that the extracted spectra are there
   if not PTR_VALID(str.frame) then begin
     error = 'APQL_SNRMAG: NO extracted spectra'
     if not keyword_set(silent) then print,error
     return
   endif

   timestep = 10.6  ; seconds
   npix = 2048L
   nfibers = 300L
   ;nfibers = n_elements(plugmap.fiberdata)

   ; Calculate the median S/N in the middle chip
   flux = (*str.frame)[1].flux > 0.0
   err = (*str.frame)[1].err > 1.0
   snr = MEDIAN(flux/(err>1),dim=1)
   str.medsnr = snr

   CASE str.exptype of

   ; OBJECT Exposures
   'OBJECT': begin

     if PTR_VALID(!apql.plugmap.datastr) then begin

       ; Match to the magnitudes
       ;fiberid = plugmap.fiberdata.fiberid
       good_fibers = where( (*!apql.plugmap.datastr).fiberdata.spectrographid eq 2 and $
                            (*!apql.plugmap.datastr).fiberdata.holetype eq 'OBJECT' AND $
                            (*!apql.plugmap.datastr).fiberdata.fiberid ge 1 and $
                            (*!apql.plugmap.datastr).fiberdata.fiberid le 300,ngood_fibers)
       if ngood_fibers eq 0 then return
       fiberdata = (*!apql.plugmap.datastr).fiberdata[good_fibers]
       ;fiberind = fiberid-1
       ; fiberid=1 is at the top of the detector or index=299
       ; index = 300-fiberid
       ADD_TAG,fiberdata,'FIBERIND',-1,fiberdata
       fiberdata.fiberind = 300-fiberdata.fiberid
       ;fiberid = (*!apql.plugmap.datastr).fiberdata.fiberid  ; 1-indexed
       ;;fiberind = fiberid-1
       ;; fiberid=1 is at the top of the detector or index=299
       ;; index = 300-fiberid
       ;fiberind = 300-fiberid
       ;;hmag = plugmap.fiberdata.mag[1] ; J, H, Ks
       ;hmag = (*!apql.plugmap.datastr).fiberdata.mag[1] ; J, H, Ks
       ;snr_match = snr[fiberind]   ; S/N matched to the fiber data
       ;str.hmag[fiberind] = hmag
       str.hmag[fiberdata.fiberind] = fiberdata.mag[1]
       str.objtype[fiberdata.fiberind] = fiberdata.objtype

       ; Maybe only select OBJECT fibers?
       ;;objind = where(plugmap.fiberdata.objtype eq 'STAR',nobjind)
       ;;objind = where(objtype eq 'STAR',nobjind)
       ;;objind = where(objtype eq 'SCIENCE',nobjind)
       ;objtype = (*!apql.plugmap.datastr).fiberdata.objtype
       ;objind = where(objtype eq 'STAR' or objtype eq 'EXTOBJ',nobjind)
       objind = where( str.objtype eq 'STAR' or str.objtype eq 'EXTOBJ',nobjind)
       if nobjind eq 0 then begin
          print, 'No STAR or EXTOBJ found in plugmap for exptype=OBJECT'
          str.snr_standard = median([snr])
       endif else begin
          ;str.objtype[fiberind] = objtype

          obj_hmag = str.hmag[objind]
          obj_snr = str.medsnr[objind]
          ;obj_hmag = hmag[objind]
          ;obj_snr = snr_match[objind]

          ; Sort them by Hmag
          si = sort(obj_hmag)
          obj_hmag = obj_hmag[si]
          obj_snr = obj_snr[si]

          ; Only want ones with decent S/N
          gdobj = where(obj_snr gt 1.0,ngdobj)

          ; Fit log(S/N) vs. Hmag
          ;   should the slope be fixed?? clustering around -0.2
          ;   It's -0.5/2.512 = -0.199045
          ;p = where(obj_snr[si] gt 1.0,count)
          if ngdobj lt 2 then begin
              ; skip if we have nothing (Any-Star-Down-Any-Fiber)
              str.snr_standard = median([snr])
          endif else begin
              obj_hmag = obj_hmag[gdobj]
              obj_snr = obj_snr[gdobj]

              ;icoef = poly_fit(obj_hmag[si[p]],alog10(obj_snr[si[p]]),1,/double)
              icoef = poly_fit(obj_hmag,alog10(obj_snr>1),1,/double)
              str.logsnr_hmag_coef = reform(icoef)
              ; S/N at "standard" Hmag
              logsnr_standard = poly(str.hmag_standard,icoef)
              str.snr_standard = 10.^logsnr_standard

              ; plotting
              ;debug=1
              if keyword_set(debug) then begin
                y0 = 0.06
                dy = 0.32
                dyoff = 0.07
                setdisp,/silent
                yr = minmax(alog10(obj_snr[si]>1))
                yr = [yr[0]-range(yr)*0.05-0.5,yr[1]+range(yr)*0.05]
                plot,obj_hmag,alog10(obj_snr>1),ps=1,xr=[5,13],yr=yr,xs=1,ys=1,$
                     xtit='H',ytit='log(S/N)',tit=str.frameid+' Read='+str.readnum,$
                     position=[0.08,y0+2*dy,0.98,y0+3*dy-dyoff]
                ;     position=[0.08,0.56,0.98,0.96]
                oplot,[5,13],poly([5,13],icoef),co=250
                oplot,[str.hmag_standard],[logsnr_standard],ps=4,co=150,sym=2
                xyouts,5.5,yr[0]+0.3,strtrim(icoef[0],2)+' '+strtrim(icoef[1],2),align=0,charsize=1.3
                xyouts,5.5,yr[0]+0.1,'S/N standard = '+stringize(str.snr_standard,ndec=2),align=0,charsize=1.3
              endif

              ; S/N of standard H magnitude for CURRENT read
              ;curind = where(readstr.readnum eq long(str.readnum),ncurind)
              ;str.snr_standard = readstr[curind].snr_standard

              ; Overplot the S/N "goal"
              ;  use the same slope as other fits
              ;  Should use fixed value of -0.199045
              ;med_slope = median([readstr.logsnr_hmag_coef[1]])
              med_slope = str.logsnr_hmag_coef[1]
              ; y = mx+b
              ; b = y-mx
              yoffset = alog10(str.snr_standard_goal>1) - med_slope*str.hmag_standard
              goal_coef = [yoffset,med_slope]
              ;oplot,x,10.^(poly(x,goal_coef)),linestyle=0,co=cogreen,thick=thick*pbin
              str.logsnr_hmag_coef_goal = goal_coef
          endelse
       endelse

     endif ; plugmap

   end ; object

   ; ASDAF
   'ASDAF': begin

     ; Find fibers with above average S/N
     ;medsnr = median(str.medsnr)
     ;sigsnr = MAD(str.medsnr) > 0.5
     ;highind = where(str.medsnr gt 5*sigsnr+medsnr,nhighind)
     ;if nhighind gt 0 then str.snr_standard=median([str.medsnr[highind]]) else $
     ;  str.snr_standard = max(str.medsnr)

     ; Just use the maximum S/N, that should be of the star
     str.snr_standard = max(str.medsnr)

   end ; asdaf


   ; NON-OBJECT Exposures
   ;----------------------
   ; flat, dome
   else: begin

     ; For non-object exposures SNR_STANDARD means SNR_MEDIAN
     str.snr_standard = median([snr])

     ; Check for missing fibers
     if str.exptype eq 'QUARTZFLAT' or str.exptype eq 'DOMEFLAT' then begin
       fiberflux = (*str.frame)[1].flux
       medfiberflux = median(fiberflux[1600:2000,*],dim=1)        ; median flux in the green, all 300 fibers
       medall_medfiberflux = median(medfiberflux)            ; median flux of all 300 fibers
       missingfiber_thresh = 0.05*medall_medfiberflux > 50   ; threshold
       bdfibers = where(medfiberflux lt missingfiber_thresh,nbdfibers)

       ; Some missing fibers
       if nbdfibers gt 0 then begin
         missing_fibers = 300-bdfibers
         missing_fibers = missing_fibers[sort(missing_fibers)]
         print,strtrim(nbdfibers,2),' missing fibers: ',strjoin(strtrim(missing_fibers,2),', ')
         str.nmissingfibers = n_elements(missing_fibers)
         str.missingfibers = PTR_NEW(missing_fibers)
       endif

       ;stop

     endif ; check for missing fibers

   end  ; non-object exposures

   ENDCASE

print,'snr_standard = ',str.snr_standard
print,'logsnr_hmag_coef = ',str.logsnr_hmag_coef

   ; SNR2 vs. time
   ;-------------------

   ; S/N vs. time

   ; Plotting
   if keyword_set(debug) and n_elements(allstr) gt 0 then begin
     xr = [min(long(allstr.readnum))-1,max(long(allstr.readnum))+1]
     yr = [0.0,max(allstr.snr_standard^2)*1.05]
     plot,[allstr.readnum],[allstr.snr_standard^2],ps=8,xr=xr,yr=yr,xs=1,ys=1,$
          xtit='Read Number',ytit='(S/N standard)^2',tit='(S/N standard)^2 vs. Readnum',$
          position=[0.08,y0+dy,0.98,y0+2*dy-dyoff],/noerase
     ;     position=[0.08,0.08,0.98,0.48],/noerase
     x = scale_vector(findgen(100),xr[0],xr[1])
     if n_elements(allstr) ge 2 then begin
       coef = poly_fit(allstr.readnum,allstr.snr_standard^2,1)
       oplot,x,poly(x,coef),co=250 ;linestyle=2
     endif
   endif


   if n_elements(allstr) gt 0 then begin

     readnum = long([allstr.readnum,str.readnum])
     snr_standard = [allstr.snr_standard,str.snr_standard]
     coef = poly_fit(readnum*timestep,snr_standard^2,1)
     str.snr2_time_coef = coef
     ;oplot,x,poly(x,coef),linestyle=dashlstyle,thick=thick*pbin,color=color ; =2

     ; When do we expect to reach our goal
     ;  y=mx+b
     ;  x = (y-b)/m
     ; limit the number of reads to avoid Infinte (Inf) which causes lots of problems later
     expected_total_readnum = ((str.snr_standard_goal^2 - coef[0])/(coef[1]*timestep)) < 9999.0
     str.expected_total_readnum = expected_total_readnum
     ; Overplot expected limit
     ;oplot,[0,0]+expected_total_readnum,yr,linestyle=dotlstyle,thick=thick*pbin,color=color

     ; Measure Delta SNR2 with respect to previous read
     nallstr = n_elements(allstr)
     str.delta_snr2_standard = str.snr_standard^2 - allstr[nallstr-1].snr_standard^2

   ; Read=2
   endif else begin

     ;coef = poly_fit([1,readstr.readnum],[0.0,readstr.snr_standard^2],1)
     coef = poly_fit([1,str.readnum]*timestep,[0.0,str.snr_standard^2],1)
     str.snr2_time_coef = coef

     ; When do we expect to reach our goal
     ; limit the number of reads to avoid Infinte (Inf) which causes lots of problems later
     expected_total_readnum = ((str.snr_standard_goal^2 - coef[0])/(coef[1]*timestep)) < 9999.0
     str.expected_total_readnum = expected_total_readnum
     ; Overplot expected limit
     ;oplot,[0,0]+expected_total_readnum,yr,linestyle=dotlstyle,thick=thick*pbin,color=color

   endelse


   ; plotting
   if keyword_set(debug) and n_elements(allstr) gt 0 then begin
     xr = [min(long(allstr.readnum))-1,max(long(allstr.readnum))+1]
     ;yr = [0.0,max(readstr.avgobjsnr)*1.05]
     ;yr = [0.0,40]
     yr = [0.0,60]
     plot,[allstr.readnum],[allstr.delta_snr2_standard],ps=8,xr=xr,yr=yr,xs=1,ys=1,$
          xtit='Read Number',ytit='delta (S/N standard)^2',tit='delta (S/N standard)^2 vs. Readnum',$
          position=[0.08,y0,0.98,y0+dy-dyoff],/noerase
     if n_elements(allstr) ge 2 then begin
       med_delta = median([allstr.delta_snr2_standard])
       oplot,xr,[0,0]+med_delta,co=250
     endif
   endif


   ; Linear fits of (S/N standard)^2 vs. time
   ;------------------------------------------

   ; Linear fit based on the last few points
   if n_elements(allstr) gt 0 then begin
     ;recent_ind = where(readstr.readnum eq long(str.readnum) or readstr.readnum eq long(str.readnum)-1,nrecent_ind)
     ;lastind = first_el(maxloc(allstr.readnum))

     ; Use the last 6 reads
     nallstr = n_elements(allstr)
     indall = indgen(nallstr)
     lastind = indall[(nallstr-5)>0:nallstr-1]

     readnum = [allstr[lastind].readnum, str.readnum]
     snr_standard = [allstr[lastind].snr_standard, str.snr_standard]
     recent_coef = AP_ROBUST_POLY_FIT(readnum*timestep,snr_standard^2,1)
     ;print,'Recent coef=',recent_coef
     ;print,snr_standard
     str.snr2_time_recent_coef = recent_coef
     ; When do we expect to reach our goal
     ;  y=mx+b
     ;  x = (y-b)/m
     ; limit the number of reads to avoid Infinte (Inf) which causes lots of problems later
     expected_total_readnum_recent = ((str.snr_standard_goal^2 - recent_coef[0])/(recent_coef[1]*timestep)) < 9999.0
     str.expected_total_readnum_recent = expected_total_readnum_recent

   endif

   ; SNR STATUS
   ;------------------
   ;if str.snr_standard ge str.snr_standard_goal then str.continue_status=1
   ; THIS IS NOW SUPERCEDED BY APQL_STATUSCHECK.PRO

  ;wait,0.5
  ;stop

end
