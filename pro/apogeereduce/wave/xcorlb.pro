pro  xcorlb,star,temp,range,shift,chisq,shifterr,errspec=errspec0,mask=mask0,plot=pl,out=ot,error=error,stp=stp

;+
;
; XCORLB
;
; Cross-correlate two spectra.
;
; INPUTS:
;  star      Stellar object spectrum
;  temp      Stellar template spectrum
;  range     The range of pixels to search for.
;  =errspec  The error spectrum for "star"
;  =mask     The mask to use with "star".  1-for points to use,
;              0-for points to ignore.
;              This can also be an array of weights.
;  /plot     Plot
;  /out      Verbose output
;  /stp      Stop at the end of the program.
;
; OUTPUTS:
;  shift     The best shift of "star" relative to "temp", i.e. a shift
;             of +10 means that "star" is shifted 10 pixels to the
;             right of "temp", OR in other words you need to shift
;             temp 10 pixels to the right in order for it to match "star".
;  chisq     The reduced chi squared of the best shift
;  shifterr  The 1-sigma uncertainty (confidence interval) for
;              "shift".  This is only meaningful if a "real" errspec is input.
;
; USAGE:
;  IDL>xcorlb,spec1,spec2,10,shift,chisq,shifterr
;
; By G.Marcy 1988
;  Modified by D.Nidever  Oct 2008
;-

apgundef,error,shift,chisq,shifterr

nstar = n_elements(star)
ntemp = n_elements(temp)
nrange = n_elements(range)

; Not enough inputs
if nstar eq 0 or ntemp eq 0 or nrange eq 0 then begin
   print,'Syntax - xcorlb, star,temp,range,shift,chisq,shifterr,errspec=errspec,mask=mask,plot=pl,out=ot,error=error'
   print, 'Uses spline to find extremum'   
   error = 'Not enough inputs'
   return
endif

; Error spectrum not input
if n_elements(errspec0) gt 0 then errspec=errspec0 else errspec=fltarr(nstar)+1.0

; Error spectrum not of the right length
if n_elements(errspec) ne nstar then begin
  error = 'ERRSPEC must have same number of elements as STAR'
  print,error
  return
endif

; Mask input
nmask = n_elements(mask0)
if nmask gt 0 then begin

  mask = mask0

  ; Not the right size
  if nmask ne nstar then begin
    error = 'MASK and STAR not of the same length'
    print,error
    return
  endif

  ; No good points
  gdmask = where(mask gt 0.0,ngdmask)
  if ngdmask lt 1 then begin
  ;if total(mask) lt 1 then begin
    error = 'NO good points in mask'
    print,error
    return
  endif

  ;bd = where(mask ne 1 and mask ne 0,nbd)
  ;if nbd gt 0 then begin
  ;  error = 'MASK must be 0s and 1s ONLY'
  ;  print,error
  ;  return
  ;endif

endif else begin
  mask = star*0.+1.0
endelse

; The mask might be an array of weights, not just 1s and 0s.
; Normalize to N, where N is the number of non-zero points
mask = mask > 0.0                ; no negative weights allowed!
dum = where(mask gt 0.0,nused)
mask = (mask/total(mask))*nused

;Measures the shift of temp. relative to star (a shift to rt. is +)
;Accuracy is typically 0.05pxl.   
;G. Marcy 12/88

ln = n_elements(temp)
ls = n_elements(star)
len = min([ln,ls])
newln = len - 2*range    ; Leave "RANGE" on ends for overhang.
; range too large
if 2*range ge len then begin
  error = 'RANGE too large for arrays'
  print,error
  return
endif
;te = temp/(total(temp)/ln)    ; NOTE!  This normalization can cause problems
;st = star/(total(star)/ls)    ;         if there are some bad pixels
;te = temp/(total(temp>0.0)/ln)
;st = star/(total(star>0.0)/ls)
te = temp       ; internal arrays
st = star
newend = range + newln - 1

; THIS IS CHOPPING OFF TOO MANY PIXELS AT THE ENDS
; MINIMIZE THIS!!!!!

; Degrees of freedom
dof = nused-1

; Shift relative to each other
x = findgen(2 * range+1) - range
chi = fltarr(2 * range+1)
for j = -range,range do begin     ; Goose step, baby.

  ; Pearson's chi square test
  ; chisq = total( (O-E)^2/E^2 )

  ; Chi square statistic
  ; chisq = total( (O-E)^2/ sig^2 )

  ; Move the TEMPLATE not the star
  ; Do NOT move the mask, otherwise the edges will be "trimmed"
  ;  as it's being shifted and the number of unmasked opints
  ;  might change with the shift and artificially affect chisq
  dif = te[range+j:newend+j] - st[range:newend]
  ;chi[j+range] = total( (dif^2) / errspec[range:newend]^2 )
  chi[j+range] = total( ((dif*mask[range:newend])^2) / errspec[range:newend]^2 )

  ;dif = te[range:newend] - st[range+j:newend+j]
  ;chi[j+range] = total( (dif^2) / errspec[range+j:newend+j]^2 )
  ;;chi[j+range] = total( ((dif^2)*mask[range+j:newend+j]) / errspec[range+j:newend+j]^2 )
  ;;dif = st[range:newend] - te[range+j:newend+j]
  ;;chi[j+range] = total( ((dif^2)*mask[range:newend]) / errspec[range:newend]^2 )
  ;chi[j+range] = total( ((dif*mask[range+j:newend+j])^2) / errspec[range+j:newend+j]^2 )
  ;;chi[j+range] = total((dif^2)/te[range:newend])  ;Too bad sdev. doesn't work.

  ;stop

endfor


; WE CHANGED THE ORDER
; flip to keep it consistent
chi = reverse(chi)

; Interpolate to get the best shift
len = n_elements(x) * 100
xl = findgen(len)
xl = xl/100. - range
xp = xl[0:len-100]
;cp = fspline(x,chi,xp)
;cp = spline(x,chi,xp)
gd = where(finite(chi) eq 1,ngd)  ; only want finite points
if ngd eq 0 then begin
  error = 'NO good chisq points'
  print,error
  return
endif
cp = cspline(x[gd],chi[gd],xp)   ; cspline is ~25x faster

; Check for problems
mm0 = first_el(where(chi eq min(chi)))
mm = first_el(where(cp eq min(cp)))
if min(cp) lt 0.0 or abs(xp[mm]-x[mm0]) ge 1.0 then begin
  ;dum = interpolate(chi[gd],xp)
  cp = interpol(chi[gd],x[gd],xp)  ; linear interpolation
endif
mm = where(cp eq min(cp))
shift = xp[mm[0]]

; Best chisq
chisq = min(cp)

; Confidence interval, min(chisq)+1
lo = (mm[0]-100) > 0
hi = (mm[0]+100) < (n_elements(xp)-1)
coef = poly_fit(xp[lo:hi],cp[lo:hi],2)
; coef = [c,b,a]
; (y-k)=a*(x-h)^2
;  k=4ac-b^2/4a  h=-b/2a (axis of symmetry)
;  The shape is determined by "a" alone.
;  Want x (relative to h) where y=k+1
;   y=a*x^2=1 -> x=1/sqrt(a)
shifterr = 1.0/sqrt(coef[2])


; Plot the spectra and chisq distribution
if keyword_set(pl) then begin
  !p.multi=[0,1,2]
  plot,findgen(ls),star,/xsty
  oplot,findgen(ln)+shift,temp,co=250,linestyle=2
  plot, x, chi,ps=1,/ynozero,/xsty,xtit='X Shift',ytit='Chi',tit='Shift='+strtrim(shift,2)+' R.Chisq='+strtrim(chisq/dof,2)
  oplot,xp,cp
  oplot,[shift],[chisq],ps=1,co=250,sym=2
  !p.multi=0
endif


; Verbose output     
if keyword_set(ot) then begin
  print,'shift = ',shift
  print,'R.Chisq = ',chisq/dof
endif

if keyword_set(stp) then stop

end

