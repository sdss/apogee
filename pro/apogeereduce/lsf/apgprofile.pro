function apgprofile,x,xcenter,par,dp=dp,stp=stp

;+
;
; APGPROFILE
;
; This program returns an APOGEE profile (LSF or PSF) with wings
; It calls both gausshermite.pro and profwings.pro to get the
; profile.  The profile is always normalized.
;
; INPUTS:
;  x       The array of X-values for which to compute the function values.
;            This needs to be a 1-dimensional array [Nx]
;  xcenter The position of the LSF center (in X-units). This needs
;            to be a 1-dimensional array [Nx]
;  par     The parameters. [binsize, Horder, GHparams, Wproftype, Wparams]
;            This needs to be 2-dimensional [Nx, Npar]
;  /stp    Stop at the end of the program
;
; OUTPUTS:
;  y       The output array of the function values for X
;  dp      The derivative array [Nx, Npar]
;     The derivs are: [height, center, sigma, H0, H1, H2, H3, H4, Wing pars]
; USAGE:
;  IDL>y = apgprofile(x,xcenter,par)
;
; By D.Nidever April 2011
;-

nx = n_elements(x)
nxcenter = n_elements(xcenter)
npar = n_elements(par)

; Not enough inputs
if nx eq 0 or nxcenter eq 0 or npar eq 0 then begin
  print,'lsf = apgprofile(x,xcenter,par,dp=dp,stp=stp)'
  return,-1
endif

sz1 = size(x)
sz2 = size(par)
nx1 = sz1[1]
npar = sz2[2]

; Parameters:
; binsize
; Horder     the GH order
; GHparams   the GH parameters: sigma, H1, H2, etc.
; Wproftype   the wing profile type
; Wparams    the wing parameters

binsize = par[0,0]  ; needs to be scalar
Horder = par[0,1]   ; needs to be scalar
GHparams = par[*,2:Horder+2]
; Wings
if (npar gt Horder+3) then begin
  Wproftype = par[0,Horder+3]  ; needs to be scalar
  Wparams = par[*,Horder+4:*]
  nWparams = n_elements(Wparams[0,*])  ; needs to be scalar
endif else nWparams=0

; Derivatives requested
;----------------------
if arg_present(dp) then begin

  ; Construct the GH input parameters
  ;  scale and H0 are fixed to 1 so it is normalized
  ;GHinpar = [1.0d0, Xcenter, GHparams[0], 1.0d0]
  ;if Horder gt 0 then GHinpar=[GHinpar, GHparams[1:*]]  ; add Hermite parameters
  GHinpar = dblarr(nx1,Horder+4)
  GHinpar[*,0] = 1
  GHinpar[*,1] = Xcenter
  GHinpar[*,2] = GHparams[*,0]
  GHinpar[*,3] = 1
  if Horder gt 0 then for i=0,Horder-1 do GHinpar[*,4+i] = GHparams[*,1+i]

  ; Correct the GH normalization for the wings
  if nWparams gt 0 then GHinpar[*,0] = 1.0-Wparams[*,0]

  ; Wings
  ;-------
  if nWparams gt 0 then begin
    Winpar = dblarr(nx1,nWparams+1)
    Winpar[*,0] = Wparams[*,0]
    Winpar[*,1] = Xcenter
    for i=1,nWparams-1 do Winpar[*,1+i] = Wparams[*,i]
    ;Winpar = [Wparams[0],Xcenter,wparams[1:*]]

    wlsf = PROFWINGSBIN(x,wproftype,Winpar,dp=dwings,binsize=binsize)
  endif

  ; Gauss-Hermite
  ;----------------
  ghlsf = GAUSSHERMITEBIN(x,GHinpar,binsize=binsize,dp=dlsf)

  ; Combine the profile
  lsf = ghlsf
  if nWparams gt 0 then lsf += wlsf

  ; Combine the derivatives
  ;  The gausshermitebin output derivatives are: height, center, sigma, GHpars
  ;    this includes the derivate for H0.
  ;  The PROFWINGS output derivatives: Wnorm, center, Wparams
  ;  The derivatives for "center" are a combination of GH and Wing, but
  ;  the others are independent of each other
  ; The derivs are: [height, center, sigma, H0, H1, H2, H3, H4, Wing pars]
  dp = dblarr(nx1,Horder+4+nWparams)
  dp[*,0:Horder+3] = dlsf
  if nWparams gt 0 then begin
    dp[*,1] += dwings[*,1]             ; add "center" deriv from Wings
    dp[*,Horder+4] = dwings[*,0]       ; Wing normalization
    ; this also changes the GH normalization, so we need to add this in
    dp[*,Horder+4] += -ghlsf/GHinpar[*,0]
    dp[*,Horder+5:*] = dwings[*,2:*]   ; Other wing params
  endif

; NO derivatives requested
;--------------------------
endif else begin

  ; Construct the GH input parameters
  ;  scale and H0 are fixed to 1 so it is normalized
  ;GHinpar = [1.0d0, Xcenter, GHparams[0], 1.0d0]
  ;if Horder gt 0 then GHinpar=[GHinpar, GHparams[1:*]]  ; add Hermite parameters
  GHinpar = dblarr(nx1,Horder+4)
  GHinpar[*,0] = 1
  GHinpar[*,1] = Xcenter
  GHinpar[*,2] = GHparams[*,0]
  GHinpar[*,3] = 1
  if Horder gt 0 then for i=0,Horder-1 do GHinpar[*,4+i] = GHparams[*,1+i]

  ; Correct the GH normalization for the wings
  if nWparams gt 0 then GHinpar[*,0] = 1.0-Wparams[*,0]

  ; Wings
  ;-------
  if n_elements(Wparams) gt 0 then begin
    Winpar = dblarr(nx1,nWparams+1)
    Winpar[*,0] = Wparams[*,0]
    Winpar[*,1] = Xcenter
    for i=1,nWparams-1 do Winpar[*,1+i] = Wparams[*,i]
    ;Winpar = [Wparams[0],xcenter,wparams[1:*]]

    wlsf = PROFWINGSBIN(x,wproftype,Winpar,binsize=binsize)
  endif

  ; Gauss-Hermite
  ;----------------
  ghlsf = GAUSSHERMITEBIN(x,GHinpar,binsize=binsize)

  ; Combine the profile
  lsf = ghlsf
  if nWparams gt 0 then lsf += wlsf

endelse


if keyword_set(stp) then stop

return,lsf

end
