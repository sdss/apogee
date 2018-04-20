;+
;
; aptrimgrid(grid,visitstr)
;
; This program trims the grid to a temperature range based on the dereddend J-K color of the star
;
; INPUTS:
;  grid        The original RV synth grid
;  visitstr    The visitstr structure with information on the
;                 visit spectra for this star.
;
; OUTPUTS:
;  newgrid     The RV synth grid with the entries for incompatable temperatures remove
;
; USAGE:
;  IDL> newgrid = aptrimgrid(grid,visitstr)
;
; Written by Nicholas Troup 9/1/2015
;-


FUNCTION aptrimgrid, grid, visitstr
 
 ; Calculate dereddened J-K color
 jk = visitstr.J-visitstr.K
 trimTEFF = 1
 IF (visitstr.J GT 30) OR (visitstr.k GT 30) THEN trimTEFF = 0
 
 ; set mildly negative ak_targ to zero (but not -99!)
 IF visitstr.ak_targ lt 0 and visitstr.ak_targ gt -0.5 then visitstr.ak_targ=0.
 IF visitstr.ak_targ ge 0 THEN BEGIN
    ak = visitstr.ak_targ
 ENDIF ELSE BEGIN
    IF visitstr.ak_wise ge 0 THEN BEGIN
        ak = visitstr.ak_wise
    ENDIF ELSE BEGIN
        ak = !VALUES.F_NAN
        trimTeff = 0
    ENDELSE
 ENDELSE

 ;To do: change to using SFD_EBV for high-lattitude stars.  

 ejk = 1.5*ak
 jk0 = jk-ejk
 IF jk0 LT 0 THEN trimTeff = 0 
 ;Stars with jk0 <0 are more likely to have bad J/K values than they are to be stars hotter than Vega (10,000 K), so don't trim the teff range in case it is an error.
 
 ;Calculate Teff from color
 ;coef = [ 0.538671, 1.31645, -2.05162, 3.69036, -3.35314, 1.15256 ] ; Coefficients from David N's commissioning work

 
 ; Coefficients from Gonzalez Hernandez and Bonifacio 2009 
 ; Teff = b0 + b1*X + b2*X^2 + b3*X*Z + b4*Z + b5*Z^2
 ; X = (J-K)_0 ; Z = [Fe/H]
 dcoef = [0.6524, 0.5813, 0.1225, 0.0646, 0.0370, 0.0016]
 gcoef = [0.6517, 0.6312, 0.0168, 0.0381, 0.0256, 0.0013]
 coef = [] ;Coefficeints from assuming [Fe/H] = 0
 feh = [-2.5,-2,-1.5,-1,-0.5,0,0.5]
 dInvTemp = (dcoef[0] + dcoef[1]*jk0 + dcoef[2]*jk0^2 + dcoef[3]*jk0*feh + dcoef[4]*feh^2 + dcoef[5]*feh^2 > 0.168) < 1.87 
 gInvTemp = (gcoef[0] + gcoef[1]*jk0 + gcoef[2]*jk0^2 + gcoef[3]*jk0*feh + gcoef[4]*feh^2 + gcoef[5]*feh^2 > 0.168) < 1.87 
 invcoltemp = [dInvTemp, gInvTemp];( poly(jk0,coef) > 0.168) < 1.87
 colorTeff = 5040./invcoltemp 
 
 IF trimTeff THEN BEGIN
   ; Establish acceptable effective temperature range to search
   Ttol = 750
   Tmin = min(colorTeff) - Ttol;colorTeff-Trange
   Tmax = max(colorTeff) + Ttol;colorTeff+Trange
 ENDIF ELSE BEGIN
   Tmin = min(grid.Teff)
   Tmax = max(grid.Teff)
 ENDELSE
 print, 'Trimming grid to Teff = ('+string(Tmin)+','+string(Tmax)+')'
 
 ; Establish acceptable logg range to search based on washington photometry flags
 ; NOTE THESE MIGHT NEED TO CHANGE FOR APOGEE2 TARGFLAGS!
 loggmin = -1.
 loggmax = 6.
 ;APOGEE_TARGET1 Bit7 = Washington Giant
 IF (visitstr.apogee_target1 AND 2L^7) THEN loggmax = 4.
 ;APOGEE_TARGET1 Bit8 = Washington Dwarf
 IF (visitstr.apogee_target1 AND 2L^8) THEN loggmin = 3.
 print, 'Trimming grid to logg = ('+string(loggmin)+','+string(loggmax)+')'
 ind = where((grid.TEFF GE Tmin) AND (grid.TEFF LE Tmax) AND (grid.LOGG GE loggmin) AND (grid.LOGG LE loggmax),nind)
 IF nind LT 1 THEN ind = where((grid.LOGG GE loggmin) AND (grid.LOGG LE loggmax)) ;IF the restrictions remove all grid options, give up the temperature restriction.
 
 
 origdata = grid.origdata[ind,*]
 data = grid.data[ind,*]
 ndata = grid.ndata[ind,*]
 metals = grid.metals[ind]
 teff = grid.teff[ind]
 logg = grid.logg[ind]  
 
 numdata = n_elements(data[0,*])
 
 remove_tags, grid, ['origdata','data','ndata','metals','teff','logg'], newgrid
 
 add_tag, newgrid, 'origdata', make_array(nind,numdata),newgrid
 add_tag, newgrid, 'data', make_array(nind,numdata),newgrid
 add_tag, newgrid, 'ndata', make_array(nind,numdata),newgrid  
 add_tag, newgrid, 'metals', make_array(nind),newgrid  
 add_tag, newgrid, 'teff', make_array(nind),newgrid  
 add_tag, newgrid, 'logg', make_array(nind),newgrid 
 
 newgrid.origdata = origdata
 newgrid.data = data
 newgrid.ndata = ndata
 newgrid.metals = metals
 newgrid.teff = teff
 newgrid.logg = logg
 return, newgrid
 
END
