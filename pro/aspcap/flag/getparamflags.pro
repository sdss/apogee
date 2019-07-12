pro getparamflags,flag,badflag,descrip

flag=['GRIDEDGE_BAD','CALRANGE_BAD','OTHER_BAD','FERRE_BAD','PARAM_MISMATCH_BAD','','','',$
      'GRIDEDGE_WARN','CALRANGE_WARN','OTHER_WARN','FERRE_WARN','PARAM_MISMATCH_WARN','OPTICAL_WARN','ERR_WARN','FAINT_WARN',$
      'PARAM_FIXED','RV_WARN','','','','','','',$
      'LOGG_CAL_RC','LOGG_CAL_RGB','LOGG_CAL_MS','LOGG_CAL_RGB_MS','','','','']
badflag=[1,1,1,1,1,0,0,0,$
         0,0,0,0,0,0,0,0,$
         0,0,0,0,0,0,0,0,$
         0,0,0,0,0,0,0,0]

descrip=[$
 'Parameter within 1/8 grid spacing of grid edge ',$
 'Parameter outside valid range of calibration determination ',$
 'Other error condition ',$
 'Failed solution in FERRE ',$
 'Elemental abundance from window differs significantly from parameter abundance ',$
 ' ',$
 ' ',$
 ' ',$
 'Parameter within 1/2 grid spacing of grid edge ',$
 'Parameter in possibly unreliable range of calibration determination ',$
 'Other warning condition ',$
 'FERRE warning (not implemented?) ',$
 'Elemental abundance from window differs from parameter abundance ',$
 'Comparison with optical abundances suggests problem ',$
 'Large expected uncertainty or upper limit based on location in parameter space (Teff, [M/H], S/N) ',$
 'Warning based on faint star/RV combination ',$
 'Parameter set at fixed value, not fit',$
 'RV puts important line off of chip ',$
 ' ',$
 ' ',$
 ' ',$
 ' ',$
 ' ',$
 ' ',$
 'Use RC gravity calibration ',$
 'Use RGB gravity calibration ',$
 'Use MS gravity calibration ',$
 'Use RBG/MS transition gravity calibration ',$
 ' ',$
 ' ',$
 ' ',$
 ' ']
end
