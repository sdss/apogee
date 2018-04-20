function skyfit_lsf_gh_dev,par,dp,x=x,y=y,err=err,binsize=binsize,nlines=nlines,loarr=loarr,hiarr=hiarr,$
                    porder=porder,forcepositive=forcepositive,doublet=doublet,dbl_sep=dbl_sep,$
                    wproftype=wproftype,wporder=wporder

;+
;
; SKYFIT_LSF_GH_DEV
;
; This gets the fitted LSF for a number of observed emission
; lines.  This is to be used with MPFITFUN.PRO in APLSF.PRO
;
; INPUTS:
;  x          The X-array for which we need to evaluate the function
;  par        The input paramters:.  They are:
;    2*Nlines  parameters for the Height and Center of each line.
;                (Height1, Center1, Height2, Center2, ...)
;    X0        An additive x-offset.  This is only used to
;                evaluate the GH parameters that vary globally
;                with X.  This is not supposed to vary.
;    GHcoefs   The polynomial coefficients for sigma and the
;                Hermite parameters.  There are Porder[i]+1
;                coefficients for parameter i.  The Hermite parameters
;                start with H1 since we fix H0=1.
;    Wcoefs      The polynomial coefficients for the wings parameters. 
;                  There are WPorder[i]+1 coefficients for wing
;                  parameters i.
;
;  =binsize    The binsize to use in X-units.  If binsize is greater
;                than zero then a "binned" function will be used.  If
;                binsize=0 then a "normal, unbinned" function is used.
;  =nlines     The number of lines to fit.
;  =loarr      An array of LOWER indices for X to be used for the
;                lines.  For example,  line 0 should only be fit for
;                X[loarr[0]:hiarr[0]].
;  =hiarr      An array of UPPER indices for X to be used for the lines.
;  =porder     This array gives the polynomial order for the
;                global variation (in X) of each GH LSF parameter.
;                That includes sigma and the Hermite
;                coefficients (starting with H1 because we fix H0=1)
;                There will be Porder[i]+1 coefficients for parameter i.
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
;nGHcoefs = n_elements(GHcoefs)
GHcoefs = par[3*nlines+1:3*nlines+nGHcoefs]  ; GHcoefs
if nWpar gt 0 then $
  Wcoefs = par[3*nlines+nGHcoefs+1:3*nlines+nGHcoefs+nWcoefs]  ; Wcoefs


nx = n_elements(x)
model = y*0.0

;dp = fltarr(nx,npar)

; Construct the input parameters
;  lsfpar = [binsize, X0, Horder, Porder, GHcoefs]
;  We are assuming here that the number of parameters
;  in GHcoefs is consistent with Porder

; Wings
if nWpar gt 0 then begin
  lsfpar = [binsize, X0, Horder, Porder, GHcoefs, Wproftype, nWpar, WPorder, Wcoefs]
; No Wings
endif else begin
  lsfpar = [binsize, X0, Horder, Porder, GHcoefs]
endelse

if n_elements(doublet) eq 0 then doublet=height*0

; Deal with the doublets
;  put them at the end
doublet_ind = where(doublet eq 1,ndoublet_ind)
nefflines = total(doublet+1)
if ndoublet_ind gt 0 then begin
  ; put the extra doublet lines at the end
  effheight = [height,height[doublet_ind]]
  effcenter = center
  effcenter[doublet_ind]-=0.5*dbl_sep[doublet_ind]
  effcenter = [effcenter, center[doublet_ind]+0.5*dbl_sep[doublet_ind]]
  ;effyoffset = [yoffset,yoffset[doublet_ind]]
  effyoffset = yoffset   ; divide yoffset for doublet by 2
  effyoffset[doublet_ind]*=0.5
  effyoffset = [effyoffset, yoffset[doublet_ind]*0.5]
  effloarr = [loarr,loarr[doublet_ind]]
  effhiarr = [hiarr,hiarr[doublet_ind]]
  origind = [indgen(nlines),doublet_ind]   ; back to the original index

; no doublets
endif else begin
  effheight = height
  effcenter = center
  effyoffset = yoffset   ; divide yoffset for doublet by 2
  effloarr = loarr
  effhiarr = hiarr
  origind = indgen(nlines)
endelse

npixarr = effhiarr-effloarr+1
ntotpix = total(npixarr)

; Get the indices that go from X->X2
dum = lonarr(ntotpix)+1
dum[0] = effloarr[0]
indhi = total(npixarr,/cum)-1
indlo = indhi-npixarr+1
dum[indlo[1:*]] = effloarr[1:*]-effhiarr[0:nefflines-2]
ind = long(total(dum,/cum))
x2 = x[ind]

; Get the indices that go from center->center2
dum2 = lonarr(ntotpix)
dum2[indlo[1:*]] = 1
ind2 = long(total(dum2,/cum))
center2 = effcenter[ind2]

; Get the LSF
lsf1d = LSF_GH(x2,center2,lsfpar,dlsfpar2,/globalderiv) ;,double=0)

; Multiply the profiles by the height
height2 = effheight[ind2]
model1d = lsf1d*height2

; Add the Y-offset
model1d += effyoffset[ind2]

model = y*0
for i=0,nefflines-1 do model[effloarr[i]:effhiarr[i]]+=model1d[indlo[i]:indhi[i]]


; Get fixed parameters
;  fixed params have 0
;  free params have 1
if arg_present(dp) then begin
  fixed = ( dp eq 0 )
  dp = make_array(n_elements(x), n_elements(par), value=0.0d0)  ; initialize
endif else fixed=par*0

; Load up the derivative array
if arg_present(dp) then begin

  ; Loop through the lines
  for i=0,nefflines-1 do begin

    ilo = effloarr[i]
    ihi = effhiarr[i]

    ilo2 = indlo[i]
    ihi2 = indhi[i]

    ; Which parameters are fixed for this line
    ;fixed1 = lonarr(3+nGHcoefs)
    fixed1 = lonarr(3+nGHcoefs+nWcoefs)
    fixed1[0:2] = fixed[3*origind[i]:3*origind[i]+2]
    fixed1[3:*] = fixed[3*nlines+1:*]

    ; Do we need /globalderiv
    ;if total(fixed1[3:*]) eq 0 then globalderiv=1 else globalderiv=0
    if total(1-fixed1[3:*]) gt 0 then globalderiv=1 else globalderiv=0


    ; Height
    ;  first param in dlsf is the height (just the lsf)
    if fixed1[0] eq 0 then $
      dp[ilo:ihi,3*origind[i]] += lsf1d[ilo2:ihi2] / err[ilo:ihi]
    ; Center
    ;  second param in dlsf is the center
    if fixed1[1] eq 0 then $
      dp[ilo:ihi,3*origind[i]+1] += effheight[i]*dlsfpar2[ilo2:ihi2,1] / err[ilo:ihi]
    ; Yoffset
    ;   only once for a doublet
    if fixed1[2] eq 0 and i le (nlines-1) then $
      dp[ilo:ihi,3*origind[i]+2] += 1 / err[ilo:ihi]

    ; Add to the derivatives for the GHcoefs
    if keyword_set(globalderiv) then $
      dp[ilo:ihi,3*nlines+1:*] += effheight[i]*dlsfpar2[ilo2:ihi2,2:*] / (err[ilo:ihi]#replicate(1,n_elements(par)-2))

  end ; line loop

end  ; dp input

dev = (y-model)/err

return,dev



;  THE REGULAR OLD METHOD
;------------------------

OLDMETHOD:

model = y*0

;t0 = systime(1)
; Get fixed parameters
;  fixed params have 0
;  free params have 1
if arg_present(dp) then begin
  fixed = ( dp eq 0 )
  ; using double-precision seems to make it faster
  dp = make_array(n_elements(x), n_elements(par), value=0.0d0)  ; initialize
  ;dp = make_array(n_elements(x), n_elements(par), value=x[0]*0)  ; initialize
endif else fixed=par*0

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

  ; Which parameters are fixed for this line
  fixed1 = lonarr(3+nGHcoefs)
  fixed1[0:2] = fixed[3*i:3*i+2]
  fixed1[3:*] = fixed[3*nlines+1:*]

  ; Do we need /globalderiv
  if total(fixed1[3:*]) eq 0 then globalderiv=1 else globalderiv=0

  ; SINGLE Line
  ;--------------
  if doublet[i] eq 0 then begin

    ;ipar = [height[i],center[i],ghpar]
    ; Get LSF with derivative
    ;if arg_present(dp) and fixed[3*i+1] eq 0 then begin
    if arg_present(dp) then begin
      lsf = LSF_GH(x[ilo:ihi],center[i],lsfpar,dlsf,globalderiv=globalderiv)
    ; Just the LSF
    endif else begin
      lsf = LSF_GH(x[ilo:ihi],center[i],lsfpar)
    endelse

    if keyword_set(forcepositive) then begin
      bdind = where(lsf lt 0.0,nbdind)
      if nbdind gt 0 then lsf[bdind] = 999999.
      if nbdind gt 0 then print,'negative lsf'
    end

    model[ilo:ihi] += height[i] * lsf


    ; Add to the dp parameters
    if arg_present(dp) then begin

      ; Height
      ;  first param in dlsf is the height (just the lsf)
      if fixed[3*i] eq 0 then $
        dp[ilo:ihi,3*i] += lsf
      ; Center
      ;  second param in dlsf is the center
      if fixed[3*i+1] eq 0 then $
        dp[ilo:ihi,3*i+1] += height[i]*dlsf[*,1]
      ; Yoffset
      if fixed[3*i+2] eq 0 then $
        dp[ilo:ihi,3*i+2] += 1

      ; Add to the derivatives for the GHcoefs
      if keyword_set(globalderiv) then $
        dp[ilo:ihi,3*nlines+1:*] += height[i]*dlsf[*,2:*]

      ;stop

    endif


  ; DOUBLET
  ;--------------
  endif else begin


    ; First line
    center1 = center[i]-0.5*dbl_sep[i]

    ; Get LSF with derivative
    ;if arg_present(dp) and fixed[3*i+1] eq 0 then begin
    if arg_present(dp) then begin
      lsf1 = LSF_GH(x[ilo:ihi],center1,lsfpar,dlsf1,globalderiv=globalderiv)
    ; Just the LSF
    endif else begin
      lsf1 = LSF_GH(x[ilo:ihi],center1,lsfpar)
    endelse

    if keyword_set(forcepositive) then begin
      bdind = where(lsf1 lt 0.0,nbdind)
      if nbdind gt 0 then lsf1[bdind] = 999999.
      if nbdind gt 0 then print,'negative lsf'
    end

    model[ilo:ihi] += height[i] * lsf1

    ; Second line
    center2 = center[i]+0.5*dbl_sep[i]

    ; Get LSF with derivative
    ;if arg_present(dp) and fixed[3*i+1] eq 0 then begin
    if arg_present(dp) then begin
      lsf2 = LSF_GH(x[ilo:ihi],center2,lsfpar,dlsf2,globalderiv=globalderiv)
    ; Just the LSF
    endif else begin
      lsf2 = LSF_GH(x[ilo:ihi],center2,lsfpar)
    endelse

    if keyword_set(forcepositive) then begin
      bdind = where(lsf2 lt 0.0,nbdind)
      if nbdind gt 0 then lsf2[bdind] = 999999.
      if nbdind gt 0 then print,'negative lsf'
    end

    model[ilo:ihi] += height[i] * lsf2


    ; Add to the dp parameters
    if arg_present(dp) then begin

      ; Height
      ;  first param in dlsf is the height (just the lsf)
      if fixed[3*i] eq 0 then $
        dp[ilo:ihi,3*i] += lsf1+lsf2
      ; Center
      ;  second param in dlsf is the center
      if fixed[3*i+1] eq 0 then $
        dp[ilo:ihi,3*i+1] += height[i] * (dlsf1[*,1] + dlsf2[*,1])
      ; Yoffset
      ;  same for both
      if fixed[3*i+2] eq 0 then $
        dp[ilo:ihi,3*i+2] += 1

      ; Add to the derivatives for the GHcoefs
      if keyword_set(globalderiv) then $
        dp[ilo:ihi,3*nlines+1:*] += height[i] * (dlsf1[*,2:*] + dlsf2[*,2:*])

      ;stop

    endif

    ;stop

  endelse ; doublet

  ; Yoffset
  model[ilo:ihi] += yoffset[i]

  ;stop

  ;ipar = [height[i],center[i],ghpar]
  ;y[ilo:ihi] +=  GAUSSHERMITEBIN(x[ilo:ihi],ipar,binsize=binsize)
  ;;y[ilo:ihi] +=  GAUSSHERMITEBIN(x[ilo:ihi],ipar,binsize=binsize,dp=idp)

  ;dp[ilo:ihi,i*2] += idp[*,0]           ; height
  ;dp[ilo:ihi,i*2+1] += idp[*,1]         ; center
  ;dp[ilo:ihi,2*nlines:*] += idp[*,2:*]  ; GH parameters

  ;stop

end


; Add to the dp parameters
;  For the GH params we need to just call skyfit_lsf_gh.pro
if arg_present(dp) then begin

  ;; Loop through the GH coefs
  ;nGHcoefs = n_elements(GHcoefs)
  ;for i=0,nGHcoefs-1 do begin
  ;
  ;  ind = 3*nlines+1+i  ; parameter index
  ;
  ;  ; Free parameter
  ;  if fixed[ind] eq 0 then begin
  ;
  ;    ; Vary that one parameter slightly
  ;    par1 = par
  ;    parstep =  ( abs(par[ind]*0.1) > 0.001)
  ;    par1[ind] += parstep
  ;
  ;    model1 = SKYFIT_LSF_GH(x,par1,binsize=binsize,nlines=nlines,loarr=loarr,hiarr=hiarr,$
  ;                       porder=porder,forcepositive=forcepositive,doublet=doublet,dbl_sep=dbl_sep)
  ;
  ;    dp[*,ind] = (model1-model)/parstep
  ;
  ;    stop
  ;
  ;  endif
  ;
  ;end


  ; DIVIDE BY THE ERRORS
  ;  mpfit wants the derivative of the deviates NOT the model
  ;  the observed values cancel out in the calculation, but
  ;  the errors don't
  dp /= err # replicate(1.0,npar)

end

dev = (y-model)/err

return,dev

end
