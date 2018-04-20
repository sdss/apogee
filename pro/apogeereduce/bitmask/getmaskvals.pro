pro getmaskvals,flag,badflag,maskcontrib,descrip

flag=['BADPIX','CRPIX','SATPIX','UNFIXABLE','BADDARK','BADFLAT','BADERR','NOSKY',$
      'LITTROW_GHOST','PERSIST_HIGH','PERSIST_MED','PERSIST_LOW','SIG_SKYLINE','SIG_TELLURIC','NOT_ENOUGH_PSF','']

badflag=[1,1,1,1,1,1,1,1,$
         0,0,0,0,0,0,1,0]

maskcontrib=[0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,$
             0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.0]

descrip=[$
 'Pixel marked as BAD in bad pixel mask or from strong persistence jump',$
 'Pixel marked as cosmic ray in ap3d',$
 'Pixel marked as saturated in ap3d',$
 'Pixel marked as unfixable in ap3d',$
 'Pixel marked as bad as determined from dark frame',$
 'Pixel marked as bad as determined from flat frame',$
 'Pixel set to have very high error (not used)',$
 'No sky available for this pixel from sky fibers',$
 'Pixel falls in Littrow ghost, may be affected',$
 'Pixel falls in high persistence region, may be affected',$
 'Pixel falls in medium persistence region, may be affected',$
 'Pixel falls in low persistence region, may be affected',$
 'Pixel falls near sky line that has significant flux compared with object',$
 'Pixel falls near telluric line that has significant absorption',$
 'Less than 50 percent PSF in good pixels',$
 ''$
]
end
