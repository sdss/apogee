FUNCTION APCALCPERSIST_LOADTABLE, cc_file

; This loads the persistence model "master table"
  
cc_data = DBLARR(2048L, 2048L, 21L)

; There are five parameters for the double-exponential fit.
; The two timescales are constant (for a given pixel), but the
; other three are a function of stimulus (exp*power+constant, 4 coefficients)

; mask
cc_data[*, *, 0] = MRDFITS(cc_file, 0, /SILENT)

; constant term
cc_data[*, *, 1] = MRDFITS(cc_file, 1, /SILENT)
cc_data[*, *, 2] = MRDFITS(cc_file, 2, /SILENT)
cc_data[*, *, 3] = MRDFITS(cc_file, 3, /SILENT)
cc_data[*, *, 4] = MRDFITS(cc_file, 4, /SILENT)

; long-term amplitude
cc_data[*, *, 5] = MRDFITS(cc_file, 5, /SILENT)
cc_data[*, *, 6] = MRDFITS(cc_file, 6, /SILENT)
cc_data[*, *, 7] = MRDFITS(cc_file, 7, /SILENT)
cc_data[*, *, 8] = MRDFITS(cc_file, 8, /SILENT)

; long-term timescale, constant
cc_data[*, *, 9] = MRDFITS(cc_file, 9, /SILENT)
;cc_data[*, *, 10] = DBLARR(2048L, 2048L)
;cc_data[*, *, 11] = DBLARR(2048L, 2048L)
;cc_data[*, *, 12] = DBLARR(2048L, 2048L)

; short-term amplitude
cc_data[*, *, 13] = MRDFITS(cc_file, 10, /SILENT)
cc_data[*, *, 14] = MRDFITS(cc_file, 11, /SILENT)
cc_data[*, *, 15] = MRDFITS(cc_file, 12, /SILENT)
cc_data[*, *, 16] = MRDFITS(cc_file, 13, /SILENT)

; short-term timescale constant
cc_data[*, *, 17] = MRDFITS(cc_file, 14, /SILENT)
;cc_data[*, *, 18] = DBLARR(2048L, 2048L)
;cc_data[*, *, 19] = DBLARR(2048L, 2048L)
;cc_data[*, *, 20] = DBLARR(2048L, 2048L)

RETURN, cc_data

END
