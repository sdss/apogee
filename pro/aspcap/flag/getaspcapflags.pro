pro getaspcapflags,flag,badflag,descrip

flag=['TEFF_WARN','LOGG_WARN','VMICRO_WARN','M_H_WARN','ALPHA_M_WARN','C_M_WARN','N_M_WARN','STAR_WARN',$
      'CHI2_WARN','COLORTE_WARN','ROTATION_WARN','SN_WARN','SPEC_HOLE_WARN','ATMOS_HOLE_WARN','VSINI_WARN','',$
      'TEFF_BAD','LOGG_BAD','VMICRO_BAD','M_H_BAD','ALPHA_M_BAD','C_M_BAD','N_M_BAD','STAR_BAD',$
      'CHI2_BAD','COLORTE_BAD','ROTATION_BAD','SN_BAD','SPEC_HOLE_BAD','ATMOS_HOLE_BAD','VSINI_BAD','NO_ASPCAP_RESULT']
badflag=[2,2,0,0,0,0,0,2,$
         2,2,2,2,2,2,0,0,$
         1,1,1,0,0,0,0,1,$
         1,1,1,1,1,2,1,1]

descrip=[$
 'WARNING on effective temperature (see PARAMFLAG[0] for details) ',$
 'WARNING on log g (see PARAMFLAG[1] for details) ',$
 'WARNING on vmicro (see PARAMFLAG[2] for details) ',$
 'WARNING on metals (see PARAMFLAG[3] for details) ',$
 'WARNING on [alpha/M] (see PARAMFLAG[4] for details) ',$
 'WARNING on [C/M] (see PARAMFLAG[5] for details) ',$
 'WARNING on [N/M] (see PARAMFLAG[6] for details) ',$
 'WARNING overall for star: set if any of TEFF, LOGG, CHI2, COLORTE, ROTATION, SN warn are set ',$
 'high chi^2 (> 2*median at ASPCAP temperature (WARN)',$
 'effective temperature more than 500K from photometric temperature for dereddened color (WARN)',$
 'Spectrum has broad lines, with possible bad effects: FWHM of cross-correlation of spectrum with best RV template relative to auto-correltion of template > 1.5 (WARN)',$
 'S/N<70 (WARN)',$
 'Grid point within 2 grid steps of hole-filled synthesis ',$
 'Grid point within 2 grid steps of hole-filled atmosphere ',$
 ' ',$
 ' ',$
 'BAD effective temperature (see PARAMFLAG[0] for details) ',$
 'BAD log g (see PARAMFLAG[1] for details) ',$
 'BAD vmicro (see PARAMFLAG[2] for details) ',$
 'BAD metals (see PARAMFLAG[3] for details) ',$
 'BAD [alpha/M] (see PARAMFLAG[4] for details) ',$
 'BAD [C/M] (see PARAMFLAG[5] for details) ',$
 'BAD [N/M] (see PARAMFLAG[6] for details) ',$
 'BAD overall for star: set if any of TEFF, LOGG, CHI2, COLORTE, ROTATION, SN error are set, or any GRIDEDGE_BAD ',$
 'high chi^2 (> 5*median at ASPCAP temperature (BAD)',$
 'effective temperature more than 1000K from photometric temperature for dereddened color (BAD)',$
 'Spectrum has broad lines, with possible bad effects: FWHM of cross-correlation of spectrum with best RV template relative to auto-correltion of template > 2 (BAD)',$
 'S/N<50 (BAD)',$
 'Grid point within 1 grid steps of hole-filled synthesis ',$
 'Grid point within 1 grid steps of hole-filled atmosphere ',$
 ' ',$
 ' '$
]
end
