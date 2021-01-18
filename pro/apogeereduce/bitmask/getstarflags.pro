pro getstarflags,flag,badflag,descrip

flag=['BAD_PIXELS','COMMISSIONING','BRIGHT_NEIGHBOR','VERY_BRIGHT_NEIGHBOR','LOW_SNR','','','',$
      'LOW_MTPFLUX','PERSIST_HIGH','PERSIST_MED','PERSIST_LOW','PERSIST_JUMP_POS','PERSIST_JUMP_NEG','','',$
      'SUSPECT_RV_COMBINATION','SUSPECT_BROAD_LINES','BAD_RV_COMBINATION','','','','','',$
      '','','','','','','','']
badflag=[1,0,0,1,0,0,0,0,$   ;0-7
         0,0,0,0,0,0,0,0,$   ;8-15
         0,0,1,0,0,0,0,0,$   ;16-23
         0,0,0,0,0,0,0,0]    ;24-31

descrip=[$
 'Spectrum has many bad pixels (>20%):  BAD',$                                                                  ;0
 'Commissioning data (MJD<55761), non-standard configuration, poor LSF: WARN',$                                 ;1
 'Star has neighbor more than 10 times brighter: WARN',$
 'Star has neighbor more than 100 times brighter: BAD',$
 'Spectrum has low S/N (S/N<5)',$                                                                               ;4
 '',$
 '',$
 '',$
 'Spectrum falls on fiber in MTP block with low (<0.5) relative flux',$
 'Spectrum has significant number (>20%) of pixels in high persistence region: WARN',$                          ;9
 'Spectrum has significant number (>20%) of pixels in medium persistence region: WARN',$
 'Spectrum has significant number (>20%) of pixels in low persistence region: WARN',$
 'Spectrum show obvious positive jump in blue chip: WARN',$
 'Spectrum show obvious negative jump in blue chip: WARN',$                                                     ;13
 '',$
 '',$
 'RVs from synthetic template differ significantly (~2 km/s) from those from combined template: WARN',$         ;16
 'Cross-correlation peak with template significantly broader than autocorrelation of template: WARN',$
 'RVs from synthetic template differ very significatly (~10 km/s) from those from combined template: BAD',$     ;18
 '',$
 '',$
 '',$
 '',$
 '',$
 '',$
 '',$
 '',$
 '',$
 '',$
 '',$
 '',$
 ''$
]
end
