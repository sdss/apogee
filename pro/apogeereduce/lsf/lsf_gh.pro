function lsf_gh,x,xcenter,par,dlsfgh,globalderiv=globalderiv,double=dbl,stp=stp,nogauss=nogauss,nowings=nowings

;+
;
; LSF_GH
;
; This returns the LSF using a superposition of 
; Gauss-Hermite functions.  The LSF is *always* normalized.
; Therefore, the coefficient for H0 is fixed to 1.
;
; NOTE: There are three ways to input the X/Xcenter values.
;         See below.
;
; INPUTS:
;  x       The array of X-values for which to compute the LSF.
;            This must be a 1 or 2-dimensional array (see below).
;  xcenter The position of the LSF center (in X-units).  This
;            must be a scalar (and the X-values 1-dim) or a
;            1-dimensional array (see below).
;
;   NOTE: There are three ways to input X/Xcenter:
;     1) X as a 1-dimensional array, and Xcenter as a scalar
;          so the same center for all X.
;          Used by LSF_GH.PRO previousy.
;     2) X as a 2-dimensional array [Ncenter,Nlsf] and X-center
;          as a 1-dimensional array [Ncenter].  Used by
;          LSF_GH2D.PRO previously
;     3) X as a 1-dimensional array and X-center as a 1-dimensional
;          array with the same number of elements as X (i.e. there
;          is an Xcenter for every X.  Used by LSF_GH1D.PRO
;          previously.
;
;  par     The input parameters:
;
;    binsize  The width of a pixel in X-units.  If this is non-zero
;               then a "binned" Gauss-Hermite function is used.  If
;               binsize=0 then a "normal, unbinned" Gauss-Hermite
;               function is used.
;    X0       An additive x-offset.  This is only used to
;               evaluate the GH parameters that vary globally
;               with X.
;    Horder   The highest Hermite order, Horder=0 means
;               only a constant term (i.e. only Gaussian).
;               There are Horder Hermite coefficients (since we fix H0=1).
;    Porder   This array gives the polynomial order for the
;               global variation (in X) of each LSF parameter.
;               That includes sigma and the Horder Hermite
;               coefficients (starting with H1 because we fix H0=1)
;               There will be Porder[i]+1 coefficients for
;               parameter i.
;    GHcoefs  The polynomial coefficients for sigma and the
;               Horder Hermite parameters.  There are Porder[i]+1
;               coefficients for parameter i.  The Hermite parameters
;               start with H1 since we fix H0=1.
;    Wproftype  The Wing profile type
;    nWpar      The number of wing parameters.
;    WPorder    An array similar to Porder to give the polynomial
;                 order for the global varition (in X) of each wing parameter.
;    Wcoefs     The polynomial coefficients for the wings parameters. 
;                 There are WPorder[i]+1 coefficients for wing
;                 parameters i.
;
;  /nogauss  This sets H0=0 which removes the constant Gaussian
;               term and makes the sum=0.
;  /stp      Stop at the end of the program.
;
; OUTPUTS:
;  lsf       The LSF function
;  dlsf      The partial derivative of lsf [Nx,Npar].
;
; USAGE:
;  IDL>lsf=lsf_gh(x,xcenter,par)
;
; BY D. Nidever  March/April 2010
;    June 2012, merged lsf_gh, lsf_gh1d, and lsf_gh2d into one
;-

nx = n_elements(x)
nxcenter = n_elements(xcenter)
npar = n_elements(par)

; Not enough input parameters
if nx eq 0 or nxcenter eq 0 or npar eq 0 then begin
  print,'Syntax - lsf = lsf_gh(x,xcenter,par,dp)'
  return,-1
endif

if n_elements(dbl) eq 0 then dbl=1  ; double-precision by default

; Determine the type of X/Xcenter input
inptype = 0
x_sz = size(x)
xcenter_sz = size(xcenter)
if x_sz[0] eq 1 and xcenter_sz[0] eq 0 then inptype=1
if x_sz[0] eq 2 and xcenter_sz[0] eq 1 then inptype=2
if x_sz[0] eq 1 and xcenter_sz[0] eq 1 then inptype=3
if inptype eq 0 then begin
  print,'X/XCENTER dimensions not compatible'
  print,'1) X [Nlsf] and Xcenter (scalar)'
  print,'2) X [Ncenter,Nlsf] and Xcenter [Ncenter]'
  print,'3) X [N] and Xcenter [N]'
  return,-1
endif

; Make internal X/Xcenter arrays depending on the
; input type
case inptype of
; normal
1: begin
  XX = x
  Xcen = dblarr(x_sz[1])+Xcenter
end
; 2D
;  X [Ncenter,Nlsf]
;  Xcenter [Ncenter]
2: begin
  XX = (x)(*)
  Xcen = ( Xcenter#replicate(1.0d0,x_sz[2]) )(*)
end
; 1D
3: begin
  XX = x
  Xcen = Xcenter
end
else:
endcase

; Breaking up the parameters
binsize = par[0]
Xoffset = par[1]   ; Additive Xoffset
Horder = par[2]
Porder = par[3:Horder+3]   ; Horder+1 array
nGHcoefs = total(Porder+1)

; Getting the GH parameters that vary globally
cpar = par[Horder+4:Horder+4+nGHcoefs-1]
if keyword_set(dbl) then coefarr = dblarr(Horder+1,max(Porder)+1) else $
  coefarr = fltarr(Horder+1,max(Porder)+1)
cstart = [0,TOTAL(Porder+1,/cum)]  ; extra one at the end
; Coefarr might have extra zeros at the end, but it shouldn't
;  make a difference.
for i=0,Horder do $
  coefarr[i,0:Porder[i]] = cpar[cstart[i]:cstart[i]+Porder[i]]

; XX and Xcen are now a 1D array
; There is a separate center for each X value.
; Now we need parameters for EACH X-element
sz1 = size(XX)
if keyword_set(dbl) then GHpar = dblarr(sz1[1],Horder+1) else $
  GHpar = fltarr(sz1[1],Horder+1)
nGHpar = Horder+1

; To evaluate the GH parmeters that vary with X we evaluate
; them at Center+Xoffset
Xcenter1 = Xcen+Xoffset
for i=0,Horder do GHpar[*,i] = POLY(Xcenter1,reform(coefarr[i,*]))


; Wing parameters
if npar gt (3+Horder+1+nGHcoefs) then begin
  wcpar = par[3+Horder+1+nGHcoefs:*] ; the wing part of "par"

  ; Nwpar     number of W parameters
  ; WPorder   the polynomial order for each
  ; Wing coefficients
  wproftype = wcpar[0]
  nWpar = wcpar[1]
  wPorder = wcpar[2:2+nWpar-1]
  nWcoefs = total(wPorder+1)

  ; Getting the Wing parameters that vary globally
  wcoef = wcpar[nWpar+2:*]
  wcoefarr = dblarr(nWpar,max(wPorder)+1)
  wcstart = [0,TOTAL(wPorder+1,/cum)]  ; extra one at the end
  ; wcoefarr might have extra zeros at the end, but it shouldn't
  ;  make a difference.
  for i=0,nWpar-1 do $
    wcoefarr[i,0:wPorder[i]] = wcoef[wcstart[i]:wcstart[i]+wPorder[i]]

  ; To evaluate the Wing parmeters that vary with X we evalute
  ; them at the Center+Xoffset
  Wparam = dblarr(sz1[1],nWpar)
  for i=0,nWpar-1 do Wparam[*,i] = POLY(Xcenter1,reform(wcoefarr[i,*]))

end ; wing parameters


if keyword_set(nowings) then undefine,wparam

; "Binned" Gauss-Hermite function
if (binsize gt 0) then begin

  ; Derivatives requested
  if arg_present(dlsfgh) then begin

    ; Wings
    if n_elements(wparam) gt 0 then begin

      inpar = dblarr(sz1[1],nGHpar+nWpar+3)
      inpar[*,0] = binsize
      inpar[*,1] = Horder
      inpar[*,2:nGHpar+1] = GHpar
      inpar[*,nGHpar+2] = wproftype
      inpar[*,nGHpar+3:nGHpar+nWpar+2] = Wparam

      lsf = APGPROFILE(XX,Xcen,inpar,dp=dlsf)

    ; NO wings
    endif else begin

      inpar = dblarr(sz1[1],nGHpar+2)
      inpar[*,0] = binsize
      inpar[*,1] = Horder
      inpar[*,2:nGHpar+1] = GHpar

      lsf = APGPROFILE(XX,Xcen,inpar,dp=dlsf)
    endelse


  ; No derivatives
  endif else begin

    ; Wings
    if n_elements(wparam) gt 0 then begin

      inpar = dblarr(sz1[1],nGHpar+nWpar+3)
      inpar[*,0] = binsize
      inpar[*,1] = Horder
      inpar[*,2:nGHpar+1] = GHpar
      inpar[*,nGHpar+2] = wproftype
      inpar[*,nGHpar+3:nGHpar+nWpar+2] = Wparam

      lsf = APGPROFILE(XX,Xcen,inpar)

    ; NO wings
    endif else begin

      inpar = dblarr(sz1[1],nGHpar+2)
      inpar[*,0] = binsize
      inpar[*,1] = Horder
      inpar[*,2:nGHpar+1] = GHpar

      lsf = APGPROFILE(XX,Xcen,inpar)
    endelse

  endelse ; no derivative


; No binning
endif else begin

  print,'No-binning not supported yet'
  return,-1

  ;; Derivatives requested
  ;if arg_present(dlsfgh) then begin
  ;  lsf = GAUSSHERMITE(x,inpar,dp=dlsf)
  ;endif else begin
  ;  lsf = GAUSSHERMITE(x,inpar)
  ;endelse

endelse

;if keyword_set(nogauss) then inpar[*,3]=0   ; No constant Gaussian term, H0

;stop

; Derivatives requested
if arg_present(dlsfgh) then begin

  ; Convert to global input parameter derivatives
  if keyword_set(globalderiv) then begin

    ; Breaking up the parameters
    ; binsize - 0  FIXED
    ; Xoffset - 1  FIXED
    ; Horder  - 2  FIXED
    ; Porder  - 3:Horder+3
    ; GHcoefs - Horder+4:npar-1


    ; We will return the derivatives of the height, center
    ;  and the GHcoefs.
    nGHcoefs = total(Porder+1)
    nWcoefs = total(WPorder+1)

    ; Initialize the output array
    dlsfgh = make_array(sz1[1],nGHcoefs+2+nWcoefs,value=lsf[0]*0)

    ; The first one is the height
    dlsfgh[*,0] = dlsf[*,0]
    ; The first one is the center
    dlsfgh[*,1] = dlsf[*,1]

    derivparcntr = 2   ; counter, where the next parameters will start in the array

    ; The parameters in dlsf are: [height, center, sigma, H0, H1, H2, H3, H4]
    ; H0 is always fixed at 1. we want H1-H4
    ; so this is now the derivative in [sigma, H1, H2, H3, H4]
    szlsf = size(dlsf)
    ;dlsf_ghpar = make_array(szlsf[1],szlsf[2]-3,value=lsf[0]*0)
    dlsf_ghpar = make_array(szlsf[1],Horder+1,value=lsf[0]*0)
    dlsf_ghpar[*,0] = dlsf[*,2]
    if Horder gt 0 then dlsf_ghpar[*,1:Horder] = dlsf[*,4:Horder+3]

    ; Loop through the GH functions
    For i=0,Horder do begin

      GHparind = indgen(Porder[i]+1) + derivparcntr
      numGHpar = long(Porder[i]+1)
      GHparind0 = derivparcntr
      GHparind1 = GHparind0 + numGHpar-1

      ; We have the derivative wrt the GH coefficient
      ; but we want it wrt the polynomial coefficient
      ; GHpar = c_0 + c_1*x + c_2*x^2
      ; So we need to multiply by the derivative of GHpar wrt
      ; that parameter, e.g. d GHpar/ d c_i =  x^i
      ; All that matters it the center, Xcenter1
     
      ; dlsf_ghpar is [Nlines,Nlsfpix,5]
      ; Xcenter1 is [Nlines,Nlsfpix]
      ;num = [szlsf[1], szlsf[2], numGHpar]
      ;dlsfgh[*,*,GHparind0:GHparind1] = REBIN(dlsf_ghpar[*,*,i],num) * $
      ;                                  REBIN(Xcenter1,num)^REBIN(indgen(1,1,numGHpar),num)

      ; This is 2x faster
      for j=0,numGHpar-1 do $
         dlsfgh[*,GHparind0+j] = dlsf_ghpar[*,i] * Xcenter1^j


      ; increment the counter
      derivparcntr += numGHpar

    End

    ; The parameters in dlsf are: [height, center, sigma, H0, H1, H2, H3, H4,
    ;                                W1, W2, ...]
    dlsf_wpar = dlsf[*,Horder+4:*]

    ; Loop through the Wing parameters
    For i=0,nWpar-1 do begin

      Wparind = indgen(WPorder[i]+1) + derivparcntr
      numWpar = long(WPorder[i]+1)
      Wparind0 = derivparcntr
      Wparind1 = Wparind0 + numWpar-1

      ; We have the derivative wrt the GH coefficient
      ; but we want it wrt the polynomial coefficient
      ; GHpar = c_0 + c_1*x + c_2*x^2
      ; So we need to multiply by the derivative of GHpar wrt
      ; that parameter, e.g. d GHpar/ d c_i =  x^i
      ; All that matters it the center, Xcenter1

      ;dlsfgh[*,Wparind0:Wparind1] = dlsf_wpar[*,i] # (Xcenter1^indgen(WPorder[i]+1))

      ; This is 2x faster
      for j=0,numWpar-1 do $
         dlsfgh[*,Wparind0+j] = dlsf_wpar[*,i] * Xcenter1^j

      ; increment the counter
      derivparcntr += long( WPorder[i]+1 )

    End

  ; Return normal LSF derivatives
  ;  This is the derivative of the LSF+wing parameters evaluated at
  ;  each Xcenter.  This is *not* of the global polynomial coefficients
  endif else begin
    dlsfgh = dlsf
  endelse ; normal derivatives


  ; Format differently for INPTYPE=2
  if inptype eq 2 then begin
    dlsfgh_temp = dlsfgh
    dlsfgh_sz = size(dlsfgh)
    dlsfgh = make_array(x_sz[1],x_sz[2],dlsfgh_sz[2],value=dlsf[0]*0)
    for i=0,dlsfgh_sz[2]-1 do dlsfgh[*,*,i]=dlsfgh_temp[*,i]
  endif

endif  ; derivates requested


; Change format to 2D for INPTYPE=2
if inptype eq 2 then begin
  lsf_temp = lsf
  lsf = dblarr(x_sz[1],x_sz[2])
  lsf[*,*] = lsf_temp
endif


if keyword_set(stp) then stop

return,lsf

end
