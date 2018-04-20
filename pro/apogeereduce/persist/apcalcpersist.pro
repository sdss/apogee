FUNCTION APCALCPERSIST $
, coeffs $
, in_time

in_reads = in_time / 10.6D

;out_flux = DBLARR(2048L, 2048L)
;FOR i = 0L, 2047L DO BEGIN
;    FOR j = 0L, 2047L DO BEGIN
;        out_flux[i, j] = FUNC_PERSISTENCE_EXP_DOUBLE(in_reads, reform(coeffs[i, j, *]))
;    ENDFOR
;ENDFOR
;y = p[0] + p[1] * EXP(p[2] * x) + p[3] * EXP(p[4] * x)
out_flux = coeffs[*,*,0] + coeffs[*,*,1] * EXP(coeffs[*,*,2] * in_reads) + coeffs[*,*,3] * EXP(coeffs[*,*,4] * in_reads)

RETURN, out_flux

END
