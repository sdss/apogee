;+
;
; APG_BINARYCHECK
;
; This checks if the APOGEE RV variations are consistent with the
; RV measurement uncertainties (binary=0) or if there is RV
; variability above this level due to companions or pulsations.
;
; INPUTS:
;  rvstr   Structure with individual visit RVs.
;  =prob   The probability to use for the binary threshold.  Default prob=0.001.
;
; OUTPUTS:
;  binary       The binarity flag:
;                -1  Not enough good points.  Need at least 3 with S/N>5
;                 0  RV variations consistent with RV measurement uncertainties
;                 1  Higher RV variations due to companions, pulsations, etc.
;  staberv_chisq      Chi squared value of RV variations being "stable".
;  staberv_rchisq     Reduced Chi squared value.
;  chisq_threshold    Chi squared threshold based on the number of
;                       "good" RV measurements and the threshold probability.
;  stablerv_chisq_prob The probability of getting this chisq value.
;
; USAGE:
;  IDL>rvstr,binary,stablerv_rchisq,stablerv_chisq,chisq_threshold,stablerv_chisq_prob
;
; By D.Nidever  May 2012
;-

pro apg_binarycheck,rvstr,binary,stablerv_chisq,stablerv_rchisq,chisq_threshold,stablerv_chisq_prob,$
                    prob=prob,synth=synth

apgundef,stablerv_chisq,stablerv_rchisq,chisq_threshold,stablerv_chisq_prob
binary = -1

; Not enough inputs
if n_elements(rvstr) eq 0 then begin
  print,'Syntax - rvstr,binary,stablerv_chisq,stablerv_rchisq,chisq_threshold,stablerv_chisq_prob'
  return
endif

; Are there enough good points
if keyword_set(synth) then $
gdvisits = where(rvstr.synthvrelerr lt 1e5 and rvstr.snr gt 5,ngdvisits) else $
gdvisits = where(rvstr.vrelerr lt 1e5 and rvstr.snr gt 5,ngdvisits)
if ngdvisits lt 3 then return

; Threshold probability
if n_elements(prob) eq 0 then prob = 0.001

; Compute mean velocities
;  maybe weight by 1/verr^2 instead of SNR
if keyword_set(synth) then begin
  vrelerr = rvstr[gdvisits].synthvrelerr > 0.010
  vhelio = rvstr[gdvisits].synthvhelio
endif else begin
  vrelerr = rvstr[gdvisits].vrelerr > 0.010
  vhelio = rvstr[gdvisits].vhelio
endelse
snr = rvstr[gdvisits].snr
mnvhelio = TOTAL(vhelio*snr)/TOTAL(snr)

; Check for binarity/RV variability
;-----------------------------------
  
; Compute chisq for the velocity being constant
stablerv_rchisq = sqrt( TOTAL( ((vhelio-mnvhelio)/vrelerr)^2 ) / (ngdvisits-1) )  ; reduced chisq
stablerv_chisq = TOTAL( ((vhelio-mnvhelio)/vrelerr)^2 )     ; normal chisq
; Calculate the chisq threshold (with DOF) that P=0.05
chisq_threshold = chisqr_cvf(prob,ngdvisits-1)
; The probability of having this chisq
stablerv_chisq_prob = 1-chisqr_pdf(stablerv_chisq,ngdvisits-1)

; Binary threshold
if stablerv_chisq gt chisq_threshold then binary=1 else binary=0
end
