pro getparamflags,flag,badflag,descrip

flag=['GRIDEDGE_BAD','CALRANGE_BAD','OTHER_BAD','FERRE_BAD','PARAM_MISMATCH_BAD','','','',$
      'GRIDEDGE_WARN','CALRANGE_WARN','OTHER_WARN','FERRE_WARN','PARAM_MISMATCH_WARN','','','',$
      'PARAM_FIXED','','','','','','','',$
      '','','','','','','','']
badflag=[1,1,1,1,1,0,0,0,$
         0,0,0,0,0,0,0,0,$
         0,0,0,0,0,0,0,0,$
         0,0,0,0,0,0,0,0]

descrip=[$
 'Parameter within 1/8 grid spacing of grid edge ',$
 'Parameter outside valid range of calibration determination ',$
 'Other error condition ',$
 'Failed solution in FERRE ',$
 'Elemental abundance from window differs from parameter abundance ',$
 ' ',$
 ' ',$
 ' ',$
 'Parameter within 1/2 grid spacing of grid edge ',$
 'Parameter in possibly unreliable range of calibration determination ',$
 'Other warning condition ',$
 ' ',$
 'Elemental abundance from window differs significantly from parameter abundance ',$
 ' ',$
 ' ',$
 ' ',$
 'Parameter set at fixed value, not fit',$
 ' ',$
 ' ',$
 ' ',$
 ' ',$
 ' ',$
 ' ',$
 ' ',$
 ' ',$
 ' ',$
 ' ',$
 ' ',$
 ' ',$
 ' ',$
 ' ',$
 ' ']
end
