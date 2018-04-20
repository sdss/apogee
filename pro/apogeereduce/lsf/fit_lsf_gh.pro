function fit_lsf_gh,x,par,dp,binsize=binsize,nlines=nlines,loarr=loarr,hiarr=hiarr,$
                    porder=porder,forcepositive=forcepositive,$
                    wproftype=wproftype,wporder=wporder

;+
;
; FIT_LSF_GH
;
; This gets the fitted LSF for a number of observed emission
; lines.  This is to be used with MPFITFUN.PRO in APLSF.PRO
;
; INPUTS:
;  x         The X-array for which we need to evaluate the function
;  par       The input paramters:.  They are:
;    2*Nlines parameters for the Height and Center of each line.
;               (Height1, Center1, Height2, Center2, ...)
;    X0       An additive x-offset.  This is only used to
;               evaluate the GH parameters that vary globally
;               with X.  This is not supposed to vary.
;    GHcoefs  The polynomial coefficients for sigma and the
;               Hermite parameters.  There are Porder[i]+1
;               coefficients for parameter i.  The Hermite parameters
;               start with H1 since we fix H0=1.
;    Wcoefs      The polynomial coefficients for the wings parameters. 
;                  There are WPorder[i]+1 coefficients for wing
;                  parameters i.
;
;  =binsize   The binsize to use in X-units.  If binsize is greater
;               than zero then a "binned" function will be used.  If
;               binsize=0 then a "normal, unbinned" function is used.
;  =nlines    The number of lines to fit.
;  =loarr     An array of LOWER indices for X to be used for the
;               lines.  For example,  line 0 should only be fit for
;               X[loarr[0]:hiarr[0]].
;  =hiarr     An array of UPPER indices for X to be used for the lines.
;  =porder    This array gives the polynomial order for the
;               global variation (in X) of each LSF parameter.
;               That includes sigma and the Hermite
;               coefficients (starting with H1 because we fix H0=1)
;               There will be Porder[i]+1 coefficients for parameter i.
;  =Wproftype  The Wing profile type
;  =WPorder    An array similar to Porder to give the polynomial
;                 order for the global varition (in X) of each wing parameter.
;
; OUTPUTS:
;  y          The sum of Nlines LSFs for X.
;  dp         The derivative function [Nx,Npar]
;
; USAGE:
;  IDL>y = fit_lsf_gh(x,par,binsize=binsize,nlines=nlines,loarr=loarr,hiarr=hiarr,porder=porder)
;
; By D.Nidever  March/April 2010
;-

; This runs gausshermitebin.pro for lots of lines
; and adds them up

nloarr = n_elements(loarr)
nhiarr = n_elements(hiarr)

; Horder
Horder = n_elements(Porder)-1
nGHcoefs = total(Porder+1)
nWpar = n_elements(WPorder)
if nWpar gt 0 then nWcoefs = total(WPorder+1)

; Breaking up the parameters
npar = n_elements(par)
height = par[0:3*nlines-3:3]
center = par[1:3*nlines-2:3]
yoffset = par[2:3*nlines-1:3]
X0 = par[3*nlines]           ; Xoffset
;GHcoefs = par[3*nlines+1:*]  ; GHcoefs
GHcoefs = par[3*nlines+1:3*nlines+nGHcoefs]  ; GHcoefs
if nWpar gt 0 then $
  Wcoefs = par[3*nlines+nGHcoefs+1:3*nlines+nGHcoefs+nWcoefs]  ; Wcoefs


; Each GH coefficient can itself be a function of X
;  need a 2D array of parameters
;ghpar = fltarr(6,5)   ; 6 GH coefficients x Ncoef
;ccount = 2*nlines
;for i=0,5 do begin
;  ncoef = porder[i]+1
;  ghpar[i,0:nocef-1] = par[ccount:ccount+ncoef-1]
;  ccount += ncoef
;end
;;ghpar = par[2*nlines:*]


nx = n_elements(x)
y = x*0.0

;dp = fltarr(nx,npar)

; Construct the input parameters
;  lsfpar = [binsize, X0, Horder, Porder, GHcoefs]
;  We are assuming here that the number of parameters
;  in GHcoefs is consistent with Porder

; Wings
if nWpar gt 0 then begin
  lsfpar = [binsize, X0, Horder, Porder, GHcoefs, Wproftype, nWpar, Wcoefs]
; No Wings
endif else begin
  lsfpar = [binsize, X0, Horder, Porder, GHcoefs]
endelse

; Loop through the separate lines
for i=0,nlines-1 do begin
  ; Only use certain indices
  if nloarr gt 0 and nhiarr gt 0 then begin
    ilo = loarr[i]
    ihi = hiarr[i]
  endif else begin
    ilo = 0
    ihi = nx-1
  endelse

  ;ipar = [height[i],center[i],ghpar]
  lsf = LSF_GH(x[ilo:ihi],center[i],lsfpar)
  if keyword_set(forcepositive) then begin
    bdind = where(lsf lt 0.0,nbdind)
    if nbdind gt 0 then lsf[bdind] = 999999.
    if nbdind gt 0 then print,'negative lsf'
  end

  y[ilo:ihi] += height[i] * lsf


  ; Yoffset
  y[ilo:ihi] += yoffset[i]

  ;stop

  ;ipar = [height[i],center[i],ghpar]
  ;y[ilo:ihi] +=  GAUSSHERMITEBIN(x[ilo:ihi],ipar,binsize=binsize)
  ;;y[ilo:ihi] +=  GAUSSHERMITEBIN(x[ilo:ihi],ipar,binsize=binsize,dp=idp)

  ;dp[ilo:ihi,i*2] += idp[*,0]           ; height
  ;dp[ilo:ihi,i*2+1] += idp[*,1]         ; center
  ;dp[ilo:ihi,2*nlines:*] += idp[*,2:*]  ; GH parameters
end

;stop

return,y

end
