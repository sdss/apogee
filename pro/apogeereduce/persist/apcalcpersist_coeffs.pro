FUNCTION APCALCPERSIST_COEFFS, cc_data, flux_data, domask=domask

; This calculates the five double-exponential parameters given
; the "master table" of coefficients and a 2D stimulus image.
  
flux_data_neg_i = WHERE(flux_data LT 0D)
IF (flux_data_neg_i[0] NE -1) THEN BEGIN
    flux_data[flux_data_neg_i] = 0D
ENDIF

out_coeffs = DBLARR(2048L, 2048L,5L)

; constant parameter
;FOR i = 0L, 2047L DO BEGIN
;    FOR j = 0L, 2047L DO BEGIN
;        out_coeffs[i, j, 0] = FUNC_PERSISTENCE_EXPPOW(flux_data[i, j], reform(cc_data[i, j, 1:4]))
;    ENDFOR
;ENDFOR
;y = p[0] + p[1] * EXP(p[2] * x) * x^p[3]
out_coeffs[*,*,0] = cc_data[*,*,1] + cc_data[*,*,2] * EXP(cc_data[*,*,3] * flux_data) * flux_data^cc_data[*,*,4]

; long-term amplitude
;FOR i = 0L, 2047L DO BEGIN
;    FOR j = 0L, 2047L DO BEGIN
;        out_coeffs[i, j, 1] = FUNC_PERSISTENCE_EXPPOW(flux_data[i, j], reform(cc_data[i, j, 5:8]))
;    ENDFOR
;ENDFOR
out_coeffs[*,*,1] = cc_data[*,*,5] + cc_data[*,*,6] * EXP(cc_data[*,*,7] * flux_data) * flux_data^cc_data[*,*,8]

; long-term timescale
;FOR i = 0L, 2047L DO BEGIN
;    FOR j = 0L, 2047L DO BEGIN
;        out_coeffs[i, j, 2] = FUNC_PERSISTENCE_EXPPOW(flux_data[i, j], reform(cc_data[i, j, 9:12]))
;    ENDFOR
;ENDFOR
out_coeffs[*,*,2] = cc_data[*,*,9] + cc_data[*,*,10] * EXP(cc_data[*,*,11] * flux_data) * flux_data^cc_data[*,*,12]

; short-term amplitude
;FOR i = 0L, 2047L DO BEGIN
;    FOR j = 0L, 2047L DO BEGIN
;        out_coeffs[i, j, 3] = FUNC_PERSISTENCE_EXPPOW(flux_data[i, j], reform(cc_data[i, j, 13:16]))
;    ENDFOR
;ENDFOR
out_coeffs[*,*,3] = cc_data[*,*,13] + cc_data[*,*,14] * EXP(cc_data[*,*,15] * flux_data) * flux_data^cc_data[*,*,16]

; short-term timescale
;FOR i = 0L, 2047L DO BEGIN
;    FOR j = 0L, 2047L DO BEGIN
;        out_coeffs[i, j, 4] = FUNC_PERSISTENCE_EXPPOW(flux_data[i, j], reform(cc_data[i, j, 17:20]))
;    ENDFOR
;ENDFOR
out_coeffs[*,*,4] = cc_data[*,*,17] + cc_data[*,*,18] * EXP(cc_data[*,*,19] * flux_data) * flux_data^cc_data[*,*,20]

; Multiply by the mask
if keyword_set(domask) then for i=0,4 do out_coeffs[*,*,i]*=cc_data[*,*,0]

RETURN, out_coeffs

END
