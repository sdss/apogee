function aspcap_elemfit,str,calstr,el,flag,noext=noext

; returns a _correction_ to the uncalibrated elemental abundances

je=where(strtrim(calstr.elem,2) eq el,nje)
je=je[0]
; get independent variables, setting to edge of range if outside of range
;if el eq 'M' then flag=str.paramflag[3] else if el eq 'alpha' then flag=str.paramflag[6] $
;else flag=str.elemflag[je]
teff=aspcap_xlim(str.fparam[0],calstr[je].temin,calstr[je].temax,flag,paramflagval('CALRANGE_WARN'))
mh=aspcap_xlim(str.fparam[3],calstr[je].femin,calstr[je].femax,flag,paramflagval('CALRANGE_WARN'))
;if el eq 'M' then str.paramflag[3]=flag else if el eq 'alpha' then str.paramflag[6]=flag $
;else str.elemflag[je]=flag

; internal calibration
x=teff-calstr[je].te0
par=calstr[je].par
elemfit=calstr[je].elemfit
if elemfit eq 4 then begin
  fit=par[-1]*x+par[-4]*mh
endif
if elemfit eq 2 then begin
  fit=par[-1]*x+par[-2]*x^2
endif
if elemfit eq 5 then begin
  fit=par[-1]*x+par[-4]*mh*x+par[-2]*x^2+par[-5]*x^2*mh
endif
if elemfit eq 1 then begin
  fit=par[-1]*x
endif
if elemfit eq 3 then begin
  fit=par[-1]*x+par[-2]*x^2+par[-3]*x^3
endif
if elemfit eq 0 then begin
  fit=0
endif
; external calibration
extpar=calstr[je].extpar
if ~keyword_set(noext) then begin
  if calstr[je].extfit eq 11 then begin
    mhclip=mh
    j=where(mhclip lt -1)
    mhclip[j]=-1
    j=where(mhclip gt -0.5)
    mhclip[j]=-0.5
    fit += extpar[0] + (mhclip-(-1.))*(extpar[1]-extpar[0])/0.5
  endif else fit+=(extpar[0]+extpar[1]*str.fparam[3]+extpar[2]*str.fparam[3]^2)
endif

return,fit

end
