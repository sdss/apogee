function profwingsbin,x,proftype,par,dp=dp,binsize=binsize,stp=stp

;+
;
; PROFWINGSBIN
;
; This function returns profile wings with the APOGEE LSF or PSF.
;
; The normalization is such that the integral/sum over all pixels is
; equal to par[0]*par[3].
;
; To just scale a profile that has an integral of unity, keep par[3]=1.
;
; INPUTS:
;  X         The array of X-values for which to compute the function values.
;              This needs to be 1-dimensional.
;  proftype  The profile type:
;                 1-Gaussian
;                 2-Lorentzian
;                 3-Exponential
;                 4-1/r^2
;                 5-1/r^3
;                 6-Moffat  not supported yet
;                 7-Voigt   not supported yet
;  par       The parameters.  The 1st gives the "area under the curve", the
;               normalization.  The 2nd gives the center of the
;               profile.  The meaning of the other paramters depend on
;               the profile type.
;               This needs to be 2-dimensional [Nx,Npar].
;  =binsize  The size of each "pixel" in X units.  The default is 1.0.
;  /stp      Stop at the end of the program.
;
; OUTPUTS:
;  y         The output array of the function values for X
;  dp        The derivative array [Nx, Nprofile, Npar]
;
; USAGE:
;  IDL>y = profwingsbin(x,1,par,binsize=1)
;
; By D. Nidever  April 2011
;-

apgundef,y,dp

nx = n_elements(x)
nproftype = n_elements(proftype)
npar = n_elements(par)

; Not enough inputs
if nx eq 0 or nproftype eq 0 or npar eq 0 then begin
  print,'y = profwingsbin(x,proftype,par,dp=dp,binsize=binsize,stp=stp)'
  return,-1
endif

; Binsize
if n_elements(binsize) eq 0 then binsize=1

sz1 = size(x)
sz2 = size(par)
nx1 = sz1[1]
npar = sz2[2]

; The profile types:
;   1-Gaussian
;   2-Lorentzian
;   3-Exponential
;   4-1/r^2
;   5-1/r^3
;   6-Moffat  not supported yet
;   7-Voigt   not supported yet

; Parameters
; The 1st parameter is always the "area under the curve"
; and gives the total normalization.
; The 2nd paramter is always the CENTER.
; The meaning of the other parameters depend on
; profile type.

; PROFILE TYPE
;-------------
CASE proftype of

; 1. Gaussian
;--------------
1: begin

  if npar lt 3 then begin
    print,'Need 3 parameters for Gaussian'
    return,-1
  endif

  ; y = A*exp(-0.5*(x-center)^2/sigma^2)

  ; Parameters are: [normalization, center, sigma]
  ; Gaussian area = ht*wid*sqrt(2*pi)
  center = par[*,1]
  sigma = par[*,2]
  ;ht = par[0]/(sigma*sqrt(2d*!dpi))
  ;y = ht * exp(-0.5d*(x-center)^2/sigma^2)

  ; Rescale the X-values using Gaussian center and sigma
  ;  have two X-values: (1) left side of pixel, (2) right side of pixel
  w = ( X - par[*,1])/par[*,2]
  w1 = ( X-0.5d*binsize - par[*,1])/par[*,2]
  w2 = ( X+0.5d*binsize - par[*,1])/par[*,2]
  sqr2 = sqrt(2d)  ; some sqrt constants

  y = par[0]*0.5d*(erf(w2/sqr2) - erf(w1/sqr2))
  ; erf returns a normalized gaussian (2x) so we didn't
  ; need to scale the amplitude

  ; Derivative, only if requested
  if arg_present(dp) then begin
    ; Gaussian values at the edges
    eexp1 = exp(-0.5*w1^2.)
    eexp2 = exp(-0.5*w2^2.)
    sqr2pi = sqrt(2d*!dpi)
    dydw = (eexp2-eexp1)*par[*,0]/sqr2pi

    dp = dblarr(nx1,npar)  ; initialize
    dp[*,0] = y/par[*,0]
    dp[*,1] = dydw * (-1d/par[*,2])
    dp[*,2] = dydw * (-w/par[*,2]) - y/par[*,2]
    ;dp[*,2] *= binsize
  end ; derivative

end  ; gaussian

; 2. Lorentzian
;--------------
2: begin

  if npar lt 3 then begin
    print,'Need 3 parameters for Lorentzian'
    return,-1
  endif

  ; y = A / ( ((x-center)/sigma)^2 + 1)

  ; Parameters are: [normalization, center, sigma]
  ; Area = ht*sigma*pi
  center = par[*,1]
  sigma = par[*,2]
  ;ht = par[0]/(sigma*!dpi)
  ;y = ht / ( ((x-center)/sigma)^2 + 1d )

  ; Rescale the X-values using Gaussian center and sigma
  ;  have two X-values: (1) left side of pixel, (2) right side of pixel
  w = ( X - par[*,1])/par[*,2]
  w1 = ( X-0.5d*binsize - par[*,1])/par[*,2]
  w2 = ( X+0.5d*binsize - par[*,1])/par[*,2]
  sqr2 = sqrt(2d)  ; some sqrt constants

  y = par[*,0]*(atan(w2) - atan(w1))/!dpi

  ; Derivative, only if requested
  if arg_present(dp) then begin
    ; Lorentzian values at the edges
    sqr2pi = sqrt(2d*!dpi)
    dydw = ( 1.0/(w2^2.+1) - 1.0/(w1^2.+1) )*par[*,0]/!dpi

    dp = dblarr(nx1,npar)  ; initialize
    dp[*,0] = y/par[*,0]
    dp[*,1] = dydw * (-1d/par[*,2])
    dp[*,2] = dydw * (-w/par[*,2]) - y/par[*,2]
  end ; derivative

end

; 3. Exponential
;---------------
3: begin

  if npar lt 3 then begin
    print,'Need 3 parameters for Exponential'
    return,-1
  endif

  ; y = A*exp(-abs(x-center)/scale)
  ;
  ; Parameters are: [normalization, center, scale]
  ; Area = 2*A*scale
  center = par[*,1]
  scale = par[*,2]
  ;ht = 2*par[0]*scale
  ;y = ht * exp(-abs(x-center)/scale)

  ; Rescale the X-values using Gaussian center and sigma
  ;  have two X-values: (1) left side of pixel, (2) right side of pixel
  w = ( X - par[*,1])/par[*,2]
  w1 = ( X-0.5d*binsize - par[*,1])/par[*,2]
  w2 = ( X+0.5d*binsize - par[*,1])/par[*,2]
  sqr2 = sqrt(2d)  ; some sqrt constants

  y = 2*par[*,0]*scale^2*(-exp(-abs(w2))*signs(w2) + exp(-abs(w1))*signs(w1))

  ; this has a problem at x=xcenter
  ; Need to treat any points across x=xcenter separately
  cross = where( (w1 lt 0 and w2 gt 0) or (w1 eq 0) or (w2 eq 0),ncross)
  if ncross gt 0 then begin
      y[cross] = 2*par[*,0]*scale^2*(1 - exp(-abs(w2[cross]))) + $
                 2*par[*,0]*scale^2*(1 - exp(-abs(w1[cross])))
  endif

  ; Derivative, only if requested
  if arg_present(dp) then begin
    ; Exponential values at the edges
    eexp1 = exp(-abs(w1))
    eexp2 = exp(-abs(w2))
    sqr2pi = sqrt(2d*!dpi)
    dydw = (eexp2-eexp1)*2*par[*,0]*scale^2

    dp = dblarr(nx1,npar)  ; initialize
    dp[*,0] = y/par[*,0]
    dp[*,1] = dydw * (-1d/par[*,2])
    dp[*,2] = dydw * (-w/par[*,2]) + y/par[*,2]
    ; This is slightly low for x=center
  end ; derivative

end

; 4. 1/r^2
;----------
4: begin

  if npar lt 2 then begin
    print,'Need 2 parameters for 1/r^2'
    return,-1
  endif

  ; y = A / ( (x-center)^2 + 1)
  ; very similar to lorentzian, with sigma=1

  ; Parameters are: [nxormalization, center]
  ; Area = ht*pi
  center = par[*,1]
  ;ht = par[0]/!dpi
  ;y = ht / ( (x-center)^2 + 1d )

  ; Rescale the X-values using Gaussian center and sigma
  ;  have two X-values: (1) left side of pixel, (2) right side of pixel
  w = ( X - par[*,1])
  w1 = ( X-0.5d*binsize - par[*,1])
  w2 = ( X+0.5d*binsize - par[*,1])
  sqr2 = sqrt(2d)  ; some sqrt constants

  y = par[*,0]*(atan(w2) - atan(w1))/!dpi


  ; Derivative, only if requested
  if arg_present(dp) then begin
    ; 1/r^2 values at the edges
    sqr2pi = sqrt(2d*!dpi)
    dydw = ( 1.0/(w2^2.+1) - 1.0/(w1^2.+1) )*par[*,0]/!dpi

    dp = dblarr(nx1,npar)  ; initialize
    dp[*,0] = y/par[*,0]
    dp[*,1] = -dydw
  end ; derivative

end

; 5. 1/r^3
;----------
5: begin

  if npar lt 2 then begin
    print,'Need 2 parameter for 1/r^3'
    return,-1
  endif

  ; y = A / ( abs(x-center)^3 + 1)
  ; very similar to lorentzian, with sigma=1 and now power=3

  ; Parameters are: [normalization, center]
  ; Area = 4*ht*pi/(3*sqrt(3))
  center = par[*,1]
  ht = par[*,0]*3d*sqrt(3d)/(4d*!dpi)
  ;y = ht / ( abs(x-center)^3 + 1d )

  ; Rescale the X-values using Gaussian center and sigma
  ;  have two X-values: (1) left side of pixel, (2) right side of pixel
  w = ( X - par[*,1])
  w1 = ( X-0.5d*binsize - par[*,1])
  w2 = ( X+0.5d*binsize - par[*,1])
  sqr2 = sqrt(2d)  ; some sqrt constants

  ; Integral( dx/(x^3+1) ) = (1/6)*ln( (x+1)^2/(x^2-x+1) )
  ;                            + (1/sqrt(3))*atan( (2x-1)/sqrt(3) )
  y = ht*( ( (1/6d)*alog( (abs(w2)+1)^2/(w2^2-abs(w2)+1) ) + (1/sqrt(3d))*atan( (2*abs(w2)-1)/sqrt(3d)) )*signs(w2) - $
           ( (1/6d)*alog( (abs(w1)+1)^2/(w1^2-abs(w1)+1) ) + (1/sqrt(3d))*atan( (2*abs(w1)-1)/sqrt(3d)) )*signs(w1) )

  ; Need to treat any points across x=xcenter separately
  cross = where( (w1 lt 0 and w2 gt 0) or (w1 eq 0) or (w2 eq 0),ncross)
  if ncross gt 0 then begin
      y[cross] = ht*( ( (1/6d)*alog( (abs(w2[cross])+1)^2/(w2[cross]^2-abs(w2[cross])+1) ) + (1/sqrt(3d))*atan( (2*abs(w2[cross])-1)/sqrt(3d)) ) - $
                      ( (1/6d)*alog( (0+1)^2/(0-0+1) ) + (1/sqrt(3d))*atan( (0-1)/sqrt(3d)) ) ) + $
                 ht*( -( (1/6d)*alog( (0+1)^2/(0-0+1) ) + (1/sqrt(3d))*atan( (2*0-1)/sqrt(3d)) ) - $
                      ( (1/6d)*alog( (abs(w1[cross])+1)^2/(w1[cross]^2-abs(w1[cross])+1) ) + (1/sqrt(3d))*atan( (2*abs(w1[cross])-1)/sqrt(3d)) )*signs(w1) )
  endif

  ; Derivative, only if requested
  if arg_present(dp) then begin
    ; 1/r^3 values at the edges
    sqr2pi = sqrt(2d*!dpi)
    dydw = ( 1.0/(abs(w2)^3.+1) - 1.0/(abs(w1)^3.+1) )*ht

    dp = dblarr(nx1,npar)  ; initialize
    dp[*,0] = y/par[*,0]
    dp[*,1] = -dydw
  end ; derivative

end

; 6. Moffat
;----------
6: begin

  ;                 GAUSSIAN#          Lorentzian#         Moffat#
  ;
  ;   Model     A[0]*exp(-0.5*u^2)   A[0]/(u^2 + 1)   A[0]/(u^2 + 1)^A[3]
  ;
  ;   A[0]         Peak Value          Peak Value        Peak Value
  ;   A[1]        Peak Centroid       Peak Centroid     Peak Centroid
  ;   A[2]       Gaussian Sigma           HWHM%             HWHM%
  ;   A[3]         + A[3]    *          + A[3]   *      Moffat Index
  ;   A[4]         + A[4]*x  *          + A[4]*x *         + A[4]   *
  ;   A[5]                                                 + A[5]*x *
  ;
  ;   Notes: # u = (x - A[1])/A[2]
  ;          % Half-width at half maximum
  ;          * Optional depending on NTERMS

  print,'Moffat profile not supported yet'
  return,-1

  ; y = A / ( ((x-center)/sigma)^2 + 1)^pow

  ; Parameters are: [normalization, center, sigma]
  ; Gaussian area = ht*sigma*pi
  center = par[*,1]
  sigma = par[*,2]
  pow = par[*,3]
  ;ht = par[0]/(sigma*!dpi)
  ht = par[*,0]  ; for now
  y = ht / ( ((x-center)/sigma)^2 + 1d )^pow

  ; Don't know that area for a Moffat function!!!!

  ; Derivative, only if requested
  if arg_present(dp) then begin
    dp = dblarr(nx1,npar)  ; initialize
    dp[*,0] = y/par[*,0]
    ;dp[*,*,1] = 
    ;dp[*,*,2] = 
    ;dp[*,*,3] = 
  end ; derivative

end

; 7. Voigt
;----------
7: begin

  ; We need a convolution for this, so maybe this isn't good
  ; Also the Voigt is a Gaussian convolved with a Lorentzian so
  ; not sure this is needed if we already have a Gaussian for the main LSF

  print,'Voigt profile not supported yet'
  return,-1

end

else: begin
  print,'Profile ',strtrim(ptype,2),' is NOT suported'
  return,-1
end

ENDCASE


if keyword_set(stp) then stop

return,y

end

