;+
; NAME:
;   traceset2xy
;
; PURPOSE:
;   Convert from a trace set to an array of x,y positions
;
; CALLING SEQUENCE:
;   traceset2xy, tset, xpos, ypos
;
; INPUTS:
;   tset       - Structure containing trace set
;
; OPTIONAL KEYWORDS:
;   xpos       - Input positions to evaluate YPOS; if not specified (or 0),
;                then generate an [NX,NTRACE] array of each pixel position
;
; OUTPUTS:
;   xpos       - X positions corresponding to YPOS
;   ypos       - Y centers as an [nx,nTrace] array
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   djs_laxisgen()
;   flegendre()
;   fpoly()
;
; REVISION HISTORY:
;   19-May-1999  Written by David Schlegel, Princeton.
;   01-Dec-2000  Handle scalar xpos correctly - D. Finkbeiner
;   10-Jul-2001  Added fpoly- S.Burles
;-
;------------------------------------------------------------------------------
;pro traceset2xy, tset, xpos, ypos
function trace_func,xpos,coeff,func_name=func,xmin=xmin,xmax=xmax

   ; Need 3 parameters
   ;if (N_params() LT 3) then begin
   ;   print, 'Syntax - traceset2xy, tset, xpos, ypos'
   ;   return,-1
   ;endif


   if (func EQ 'legendre' OR $
       func EQ 'chebyshev' OR $
       func EQ 'poly') then begin

      ndim = size(coeff, /n_dim)
      dims = size(coeff, /dim)

      if (ndim EQ 1) then begin
         ncoeff = dims[0]
         nTrace = 1
      endif else if (ndim EQ 2) then begin
         ncoeff = dims[0]
         nTrace = dims[1]
      endif else begin
         message, 'TSET.COEFF contains invalid number of dimensions'
      endelse

      nx = long(xmax - xmin + 1)

      xmid = 0.5 * (xmin + xmax)
      xrange = xmax - xmin

      if (NOT keyword_set(xpos)) then $
        xpos = djs_laxisgen([nx, nTrace], iaxis=0) + xmin

      ypos = xpos*0.0
      for iTrace=0, nTrace-1 do begin
         xvec = 2.0 * (xpos[*,iTrace]-xmid)/xrange
         if (func EQ 'poly') then legarr = fpoly(xvec, ncoeff)
         if (func EQ 'legendre') then legarr = flegendre(xvec, ncoeff)
         if (func EQ 'chebyshev') then legarr = fchebyshev(xvec, ncoeff)
         ypos[*,iTrace] = legarr # [coeff[*,iTrace]]
      endfor

   endif else begin
      ;error, 'Unknown function' + func
      print,'Unknown function,' + func
   endelse

   ;stop

   return,ypos
end
;------------------------------------------------------------------------------
