function gausshermitebin,x,par,dp=dGH,binsize=binsize,average=average,double=dbl,stp=stp

;+
;
; GAUSSHERMITEBIN
;
; This function returns the value for a superposition of "binned"
; Gauss-Hermite polynomials, or more correctly "Hermite functions".
; These are orthogonal functions.  The "unbinned" version of this
; program is called gausshermite.pro.
;
; The normalization is such that the integral/sum over all pixels is
; equal to par[0]*par[3]*binsize (par[0]*par[3] if /average is set).
;
; To just scale a profile that has an integral of unity, keep par[3]=1
; and use /average. 
;
; Currently this only goes up to H4.
;
; INPUTS:
;  X         The array of X-values for which to compute the function values.
;              This needs to be 1-dimensional.
;  par       The Gauss-Hermite coefficients: [height, center, sigma,
;              H0, H1, H2, H3, H4, ...] up to H9.  This needs to be
;              2-dimensional [Nx,Npar].
;  =binsize  The size of each "pixel" in X units.  The default is 1.0.
;  /average  Output the average value in each bin, not the sum.  This
;              divides the normal output by "binsize".  The is set by
;              default!!!  So use average=0 is you don't want the average.
;  /stp      Stop at the end of the program.
;
; OUTPUTS:
;  GHfunc    The output array of the function values for X
;  dGH       The derivative array [Nx, Npar]
;
; USAGE:
;  IDL>y = gausshermitebin(x,par)
;
; By D. Nidever  March/April 2010
;   much of it copied from GHM_int.pro by Ana Perez
;-

apgundef,GHfunc,dGH

nx = n_elements(x)
npar = n_elements(par)

; Not enough parameters
if nx eq 0 or npar eq 0 then begin
  print,'GHfunc = gausshermitebin(x,par,dp=dGH,binsize=binsize,average=average,stp=stp)'
  return,-1
endif

; AVERAGE BY DEFAULT!!!
;if n_elements(average) eq 0 then average=1  ; NO !!
if n_elements(average) eq 0 then average=0

if n_elements(dbl) eq 0 then dbl=1   ; double-precision by default

; This function gives the "binned" flux in a pixel.

; Parameters are:
;  3 Gaussian parameters: height, center, sigma
;  5 Hermite polynomial coefficients

sz1 = size(x)
sz2 = size(par)
nx1 = sz1[1]
npar = sz2[2]

; HH hermite polynomials
; GHfunc  the output value
; dGH  derivative of GHfunc
;nherm = 5       ; number of hermite coefficients to use
nherm = npar-3   ; 5 is the "standard" way to run it
;hpar = dblarr(nherm)
;if npar gt 3 then $
;  hpar[0:npar-4] = par[3:*]  ; hermite coefficients
;hpar = dblarr(5,nx1,nx2)
;for i=0,nherm-1 do hpar[i,*,*]=par[*,*,i+3]
hpar = par[*,3:*]
hpar = dblarr(nx1,10)
if npar gt 3 then $
  hpar[*,0:nherm-1] = par[*,3:*]  ; hermite coefficients

; Initializing arrays
;;dw = dblarr(npar,nx)
;integ = dblarr(nherm,nx)
;hh = dblarr(nherm)
;integ = dblarr(nx1,npar)
integ = dblarr(nx1,nherm)
;hh = dblarr(nx1,nherm)
if n_elements(binsize) eq 0 then binsize=1.0
       ; how wide is one pixel in X (length to integrate over)

sqrpi = sqrt(!dpi)
sqr2pi = sqrt(2d*!dpi)


; Rescale the X-values using Gaussian center and sigma
;  have two X-values: (1) left side of pixe, (2) right side of pixel
w = ( X - par[*,1])/par[*,2]
w1 = ( X-0.5*binsize - par[*,1])/par[*,2]  ; previously called "w"
w2 = ( X+0.5*binsize - par[*,1])/par[*,2]  ; previously called "wup"

; Hermite polynomials are (probabilist's):
; H0 = 1
; H1 = x
; H2 = x^2-1
; H3 = x^3-3x
; H4 = x^4-6x^2+3
; H5 = x^5-10x^3+15x
; H6 = x^6-15x^4+45x^2-15
; H7 = x^7-21x^5+105x^3-105x
; H8 = x^8-28x^6+210x^4-420x^2+105
; H9 = x^9-36x^7+378x^5-1260x^3+945x
;
; So in terms of the powers of X with coefficients
; Ci for each Hermite polynomial:
; 0th: C0-C2+3*C4
; 1st: C1-3*C3
; 2nd: C2-6*C4
; 3rd: C3
; 4th: C4
; 5th: C5-21C7+378C9
; 6th: C6-28C8
; 7th: C7-36C9
; 8th: C8

; The Hermite function is:
;  Psi_n(x) = 1/sqrt(n!*2*pi) * exp(-x^2/2) * H_n(x)

; So we are doing:
;  = Integral( P[0]*exp(-0.5*w^2) * SUM_i( P[i+3]*H[i,w] / sqrt(i!*2*pi) ) )
;  where H are the hermite polynomials, and the Sum over i is from 0->4
;  the integral is from x1->x2, w=(x-P[1])/P[2]

; The normalization is sqrt(n!)
; 0th: 1
; 1st: 1
; 2nd: sqrt(2)
; 3rd: sqrt(6)
; 4th: sqrt(24)
; 5th: sqrt(120)
; 6th: sqrt(720)
; 7th: sqrt(5040)
; 8th: sqrt(40320)
; 9th: sqrt(362880)
sqr2 = sqrt(2d)  ; some sqrt constants
sqr3 = sqrt(3d)
sqr6 = sqrt(6d)
sqr24 = sqrt(24d)
sqr120 = sqrt(120d)
sqr720 = sqrt(720d)
sqr5040 = sqrt(5040d)
sqr40320 = sqrt(40320d)
sqr362880 = sqrt(362880d)

; HH are the coefficients for the powers of X
;   hpar will be zero for orders higher than are actually desired
hh = dblarr(nx1,10)      ; the coefficients for the powers of X
hh[*,0] = hpar[*,0] - hpar[*,2]/sqr2 + hpar[*,4]*3d/sqr24 - hpar[*,6]*15d/sqr720 + hpar[*,8]*105d/sqr40320
hh[*,1] = hpar[*,1] - hpar[*,3]*3d/sqrt(6d) + hpar[*,5]*15d/sqr120 - hpar[*,7]*105d/sqr5040 + hpar[*,9]*945d/sqr362880
hh[*,2] = hpar[*,2]/sqr2 - hpar[*,4]*(6d/sqr24) + hpar[*,6]*45d/sqr720 - hpar[*,8]*420d/sqr40320
hh[*,3] = hpar[*,3]/sqrt(6d) - hpar[*,5]*10d/sqr120 + hpar[*,7]*105d/sqr5040 - hpar[*,9]*1260d/sqr362880
hh[*,4] = hpar[*,4]/sqr24 - hpar[*,6]*15d/sqr720 + hpar[*,8]*210d/sqr40320
hh[*,5] = hpar[*,5]/sqr120 - hpar[*,7]*21d/sqr5040 + hpar[*,9]*378d/sqr362880
hh[*,6] = hpar[*,6]/sqr720 - hpar[*,8]*28d/sqr40320
hh[*,7] = hpar[*,7]/sqr5040 - hpar[*,9]*36d/sqr362880
hh[*,8] = hpar[*,8]/sqr40320
hh[*,9] = hpar[*,9]/sqr362880

; Gaussian values at the edges
eexp1 = exp(-0.5*w1^2.)
eexp2 = exp(-0.5*w2^2.)

; Integrals of exp(-w^2/2)*w^i:
;  i=0   sqrt(pi/2) * erf(w^2/sqrt(2)) 
;  i=1   -exp(-w^2/w)
;  i=2   -w*exp(-w^2/2) * Integral( exp(-w^2/2) )
;  i=3   -w^2*exp(-w^2/2) * 2*Integral( w*exp(-w^2/2) )
;  i=4   -w^3*exp(-w^2/2) * 3*Integral( w^2*exp(-w^2/2) )
;  i=5   and so on
;  i=6   
;  i=7   
;  i=8   
;  after i=1 we can use a recursion relation
integ[*,0] = sqrt(!dpi/2)*(erf(w2/sqr2) - erf(w1/sqr2))
if nherm gt 1 then integ[*,1] = -eexp2 + eexp1
if nherm gt 2 then begin
  for i=2,nherm-1 do integ[*,i] = ( -w2^(i-1)*eexp2 + w1^(i-1)*eexp1 ) + (i-1)*integ[*,i-2]
end

;  Now multiply times the polynomial coefficients and sum
;if keyword_set(dbl) then GHfunc = dblarr(nx1) else GHfunc = fltarr(nx1)
GHfunc = dblarr(nx1)
for i=0,nherm-1 do GHfunc += hh[*,i]*integ[*,i]
GHfunc *= par[*,0]/sqr2pi


; Using the "average" in each pixel
if keyword_set(average) then GHfunc /= binsize


; Derivative, only if requested
if arg_present(dGH) then begin

  ; initialize
  ;if keyword_set(dbl) then dGH = dblarr(nx1,npar) else $
  ;  dGH = fltarr(nx1,npar)  ; initialize
  dGH = dblarr(nx1,npar)   ; initialize

  dGH[*,0] = GHfunc/par[*,0]

  ; The derivative of GHfunc wrt to w, dGH/dw
  ;  Since GH is the integral from w1->w2 taking the derivative
  ;  gets rid of the integral and we just evaluate the function at w1/w2
  dGHdw = eexp2*(hh[*,0] + hh[*,1]*w2 + hh[*,2]*w2^2 + hh[*,3]*w2^3 + hh[*,4]*w2^4 + hh[*,5]*w2^5 + hh[*,6]*w2^6 + hh[*,7]*w2^7 + hh[*,8]*w2^8)-$
          eexp1*(hh[*,0] + hh[*,1]*w1 + hh[*,2]*w1^2 + hh[*,3]*w1^3 + hh[*,4]*w1^4 + hh[*,5]*w1^5 + hh[*,6]*w1^6 + hh[*,7]*w1^7 + hh[*,8]*w1^8)
  dGHdw *= par[*,0]/sqr2pi
  if keyword_set(average) then dGHdw /= binsize

  dGH[*,1] = dGHdw * (-1d/par[*,2])
  dGH[*,2] = dGHdw * (-w/par[*,2]) - GHfunc/par[*,2]


  dGH[*,3] = par[*,0] * integ[*,0] / sqr2pi
  SWITCH npar of
    13: dGH[*,12] = (par[*,0]/sqr2pi) * ( integ[*,9] - 36*integ[*,7] + 378*integ[*,5] - 1260*integ[*,3] + 945*integ[*,1] )/sqr362880 ; H9
    12: dGH[*,11] = (par[*,0]/sqr2pi) * ( integ[*,8] - 28*integ[*,6] + 210*integ[*,4] - 420*integ[*,2] + 105*integ[*,0] )/sqr40320 ; H8
    11: dGH[*,10] = (par[*,0]/sqr2pi) * ( integ[*,7] - 21*integ[*,5] + 105*integ[*,3] - 105*integ[*,1] )/sqr5040   ; H7
    10: dGH[*,9] =  (par[*,0]/sqr2pi) * ( integ[*,6] - 15*integ[*,4] + 45*integ[*,2] - 15*integ[*,0] )/sqr720   ; H6
    9:  dGH[*,8] =  (par[*,0]/sqr2pi) * ( integ[*,5] - 10*integ[*,3] + 15*integ[*,1] )/sqr120    ; H5
    8:  dGH[*,7] =  (par[*,0]/sqr2pi) * ( integ[*,4] - 6*integ[*,2] + 3*integ[*,0] )/sqr24  ; H4
    7:  dGH[*,6] =  (par[*,0]/sqr2pi) * ( integ[*,3] - 3*integ[*,1] )/sqr6  ; H3
    6:  dGH[*,5] =  (par[*,0]/sqr2pi) * ( integ[*,2] - integ[*,0] )/sqr2  ; H2
    5:  dGH[*,4] =  (par[*,0]/sqr2pi) * integ[*,1]     ; H1
  ENDSWITCH
end

if keyword_set(stp) then stop

return,GHfunc

end

