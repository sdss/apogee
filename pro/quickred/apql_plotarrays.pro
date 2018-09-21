;+
;
; APQL_PLOTARRAYS
;
; This plots an image of the three arrays for the
; Quicklook software
;
; INPUTS:
;  str       The quicklook structure
;  /silent   Don't print anything to the screen.
;
; OUTPUTS:
;  The binned data arrays are put into the STR structure.
;  =error    The error message if one occurred.
;
; USAGE:
; IDL>apql_plotarrays,str
;
; By D.Nidever  2010
;-
pro apql_plotarrays,str,silent=silent,error=error

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
     if not keyword_set(silent) then print,'Syntax - apql_plotarrays,str,silent=silent,error=error'
     return
   endif

   ; Check that the CDS image exists
   if not PTR_VALID(str.cds_image) then begin
     error = 'APQL_PLOTARRAYS: NO CDS image'
     if not keyword_set(silent) then print,error
     return
   endif

   ; Default parameters
   npix = 2048L
   chipgap = 300
   nchips = 3
   ;npick = 3


   ; Bin the entire array
   nbin = str.arraydisplay_nbin
   im = (*(str.cds_image))[0:nchips*npix-1,*]
   npixbin = npix/nbin
   ;binim = REBIN(im,npixbin,npixbin*nchips)  ; averages values
   binim = REBIN(im,npixbin*nchips,npixbin)  ; averages values
   ;binim = transpose(binim)          ; flip it
   zscale,binim,z1,z2,contrast=0.10 ;0.25
   ; scale it
   minim = min(binim)
   maxim = max(binim)
   ; real_values = byte_values * bscale + bzero
   bzero = minim
   bscale = (maxim-minim)/255. > 1  ; don't scale up
   byte_binim = ( binim - bzero )/bscale
   byte_binim = byte( round( byte_binim ) ) ; round and make byte type
   ; stuff into the structure
   str.arraydisplay.data = byte_binim
   str.arraydisplay.bscale = bscale
   str.arraydisplay.bzero = bzero
   str.arraydisplay.zscale = [z1,z2]


   ; Loop through the arrays
   yloarr = [500,1000,1500]
   yhiarr = yloarr+99
   nsub = n_elements(yloarr)

   ; Loop through the subranges
   FOR i=0,nsub-1 DO BEGIN

     ylo = yloarr[i]
     yhi = yhiarr[i]
     nypix = yhi-ylo+1

     im = (*(str.cds_image))[0:nchips*npix-1,ylo:yhi]
     
     ; Only bin in Y-dimension
     npixbin = npix/nbin
     nypixbin = nypix/nbin
     ;binim = REBIN(im,nxpix,npixbin*nchips)  ; averages values
     binim = REBIN(im,npixbin*nchips,nypix)  ; averages values
     zscale,binim,z1,z2,contrast=0.10 ;0.25
     ; scale it
     minim = min(binim)
     maxim = max(binim)
     ; real_values = byte_values * bscale + bzero
     bzero = minim
     bscale = (maxim-minim)/255. > 1  ; don't scale up
     byte_binim = ( binim - bzero )/bscale
     byte_binim = byte( round( byte_binim ) ) ; round and make byte type
     ; stuff into the structure
     if PTR_VALID(str.arraydisplay_sub[i].data) then PTR_FREE,str.arraydisplay_sub[i].data
     str.arraydisplay_sub[i].data = PTR_NEW(byte_binim,/no_copy)
     str.arraydisplay_sub[i].bscale = bscale
     str.arraydisplay_sub[i].bzero = bzero
     str.arraydisplay_sub[i].yrange = [ylo,yhi]
     str.arraydisplay_sub[i].zscale = [z1,z2]

     ;stop

   ENDFOR

end
