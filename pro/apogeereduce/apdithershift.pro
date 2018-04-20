pro apdithershift,frame1,frame2,shift,shifterr,xcorr=xcorr,lines=lines,object=object,plot=plot,pfile=pfile,stp=stp,shiftarr=shiftarr,plugmap=plugmap,nofit=nofit

;+
;
; APDITHERSHIFT
;
; This program measured the SHIFT in two dithered images
;
; INPUTS:
;  frame1    A structure giving 1D extracted data and headers for
;              all three chips for the FIRST dithered frame.
;  frame2    The same as "frame1" but for the SECOND dithered frame.
;  /xcorr    Use cross-correlation to measure the shift.
;  /lines    Use emission lines to measure the shift.  This is the default.
;  /object   This is an object spectrum (stellar spectra).  The
;              spectra will be normalized by the stellar continuum.
;  /pl       Show plots of the fits.
;  /stp      Stop at the end of the program.
;
; OUTPUTS:
;  shift     The shift in pixels between frame2 and frame1.  A positive
;              shift means that frame2 is to the RIGHT of frame1.
;  shifterr  The uncertainty of the shift measurement.
;
; USAGE:
;  IDL>apdithershift,frame1,frame2,shift
;
; By D. Nidever  March 2010
;-

apgundef,shift,shifterr

nframe1 = n_elements(frame1)
nframe2 = n_elements(frame2)

; Not enough inputs
if nframe1 eq 0 or nframe2 eq 0 then begin
  print,'Syntax - apdithershift,frame1,frame2,shift,shifterr,xcorr=xcorr,lines=lines,object=object,pl=pl,stp=stp'
  return
endif


; Checking the tags of the input structure in FRAME1
tags = tag_names(frame1)
needtags1 = ['CHIPA','CHIPB','CHIPC']
for i=0,n_elements(needtags1)-1 do begin
  if (where(tags eq needtags1[i]))[0] eq -1 then begin
    print,'TAG ',needtags1[i],' NOT FOUND in input structure'
    return
  end
end
needtags2 = ['HEADER','FLUX','ERR','MASK']
for i=0,2 do begin
  tags2 = tag_names(frame1.(i))
  for j=0,n_elements(needtags2)-1 do begin
    if (where(tags2 eq needtags2[j]))[0] eq -1 then begin
      print,'TAG ',needtags2[j],' NOT FOUND in input structure'
      return
    end
  end
end

; Checking the tags of the input structure in FRAME2
tags = tag_names(frame2)
needtags1 = ['CHIPA','CHIPB','CHIPC']
for i=0,n_elements(needtags1)-1 do begin
  if (where(tags eq needtags1[i]))[0] eq -1 then begin
    print,'TAG ',needtags1[i],' NOT FOUND in input structure'
    return
  end
end
needtags2 = ['HEADER','FLUX','ERR','MASK']
for i=0,2 do begin
  tags2 = tag_names(frame2.(i))
  for j=0,n_elements(needtags2)-1 do begin
    if (where(tags2 eq needtags2[j]))[0] eq -1 then begin
      print,'TAG ',needtags2[j],' NOT FOUND in input structure'
      return
    end
  end
end


; Temporary versions of the data
f1 = frame1
f2 = frame2

sz = size(f1.chipa.flux[*,*])
npix = sz[1]
nfibers = sz[2]

;-------------------------
; Using CROSS-CORRELATION
;-------------------------
if keyword_set(xcorr) then begin

  print,'Using Cross-correlation to measure the dither shift'

  if keyword_set(plugmap) then begin
    iplugind = where(plugmap.fiberdata.spectrographid eq 2,nplugged) 
    iplugind = 300-plugmap.fiberdata[iplugind].fiberid
  endif else iplugind=indgen(nfibers)
  f1chipaflux=f1.chipa.flux[*,iplugind]
  f1chipbflux=f1.chipb.flux[*,iplugind]
  f1chipcflux=f1.chipc.flux[*,iplugind]
  f2chipaflux=f2.chipa.flux[*,iplugind]
  f2chipbflux=f2.chipb.flux[*,iplugind]
  f2chipcflux=f2.chipc.flux[*,iplugind]
  sz = size(f1chipaflux[*,*])
  nfibers = sz[2]

  ; Should use the PlugMap to only pick object spectra??

  ; The input arrays are [2048,300,8]
  ;  the planes are: [spec, wave, error, flag, sky, errsky,
  ;      telluric, error_telluric]

  ; Object spectra, normalize
  ;  DON'T want to do this for SKY FIBERS!!!!
  if keyword_set(object) then begin
    print,'Median filtering the spectra'
    f1chipaflux = f1chipaflux / ( MEDFILT2D(f1chipaflux,100,dim=1,/edge_copy) > 1)
    f1chipbflux = f1chipbflux / ( MEDFILT2D(f1chipbflux,100,dim=1,/edge_copy) > 1)
    f1chipcflux = f1chipcflux / ( MEDFILT2D(f1chipcflux,100,dim=1,/edge_copy) > 1)
    f2chipaflux = f2chipaflux / ( MEDFILT2D(f2chipaflux,100,dim=1,/edge_copy) > 1)
    f2chipbflux = f2chipbflux / ( MEDFILT2D(f2chipbflux,100,dim=1,/edge_copy) > 1)
    f2chipcflux = f2chipcflux / ( MEDFILT2D(f2chipcflux,100,dim=1,/edge_copy) > 1)
  endif

  ; Do the cross-correlation
  nlags = 21
  lags = lindgen(nlags)-nlags/2
  lo = nlags
  hi = npix-nlags-1

  fiber=fltarr(nfibers,3)
  chip=fltarr(nfibers,3)
  xshiftarr=fltarr(nfibers,3)
  for ichip=0,2 do begin
    xcorr = fltarr(nlags,nfibers)
    for i=0,nlags-1 do begin
      if ichip eq 0 then xcorr[i,*] = TOTAL( f1chipaflux[lo:hi,*] * f2chipaflux[lo+lags[i]:hi+lags[i],*], 1 )
      if ichip eq 1 then xcorr[i,*] = TOTAL( f1chipbflux[lo:hi,*] * f2chipbflux[lo+lags[i]:hi+lags[i],*], 1 )
      if ichip eq 2 then xcorr[i,*] = TOTAL( f1chipcflux[lo:hi,*] * f2chipcflux[lo+lags[i]:hi+lags[i],*], 1 )
    endfor
    for i=0,nfibers-1 do begin
      ;fiber[i,ichip] = i
      fiber[i,ichip] = iplugind[i]
      chip[i,ichip] = ichip
      xshiftarr[i,ichip] = -100.

      xcorr0 = xcorr[*,i]
      if total(xcorr0) eq 0 then goto, BOMB1
      ; Now fit a Gaussian to it
      coef1 = ap_robust_poly_fit(lags,xcorr0,1)
      estimates = [max(xcorr0), 0.0, 2.0, coef1[0], coef1[1]]
      yfit1 = MPFITPEAK(lags,xcorr0,par1,nterms=5,estimates=estimates,/gaussian,/positive,$
                      perror=perror1,chisq=chisq1,dof=dof1,yerror=yerror1,status=status1)

      estimates = [par1[0], 0.0, par1[2], coef1[0], coef1[1], 0.0, 0.0]
      par2 = MPFITFUN('gausspoly',lags,xcorr0,xcorr0*0+1,estimates,perror=perror2,chisq=chisq2,$
                      dof=dof2,yerror=yerror2,status=status2,yfit=yfit2,/quiet)
  
      gd = where(abs(lags-par2[1]) lt 3.0*par2[2],ngd)
      if ngd eq 0 then goto,BOMB1
      par3 = MPFITFUN('gausspoly',lags[gd],xcorr0[gd],xcorr0[gd]*0+1,par2,perror=perror3,chisq=chisq3,$
                    dof=dof3,yerror=yerror3,status=status3,yfit=yfit3,/quiet)

      xshiftarr[i,ichip] = par3[1]
      BOMB1:
    endfor
  endfor

;  xcorr1 = fltarr(nlags,nfibers)
;  xcorr2 = fltarr(nlags,nfibers)
;  xcorr3 = fltarr(nlags,nfibers)
;  lo = nlags
;  hi = npix-nlags-1
;  for i=0,nlags-1 do begin
;    xcorr1[i,*] = TOTAL( f1.chipa.flux[lo:hi,*] * f2.chipa.flux[lo+lags[i]:hi+lags[i],*], 1 )
;    xcorr2[i,*] = TOTAL( f1.chipb.flux[lo:hi,*] * f2.chipb.flux[lo+lags[i]:hi+lags[i],*], 1 )
;    xcorr3[i,*] = TOTAL( f1.chipc.flux[lo:hi,*] * f2.chipc.flux[lo+lags[i]:hi+lags[i],*], 1 )
;  end
;  xcorr = [[xcorr1],[xcorr2],[xcorr3]]  ; combine the three chip xcorr arrays
;  xshiftarr = fltarr(nfibers*3)
;  print,'Fitting Cross-correlation peaks'
;
;  for i=0,nfibers*3-1 do begin
;
;    xcorr0 = xcorr[*,i]
;
;    ; Now fit a Gaussian to it
;    coef1 = ap_robust_poly_fit(lags,xcorr0,1)
;    estimates = [max(xcorr0), 0.0, 2.0, coef1[0], coef1[1]]
;    yfit1 = MPFITPEAK(lags,xcorr0,par1,nterms=5,estimates=estimates,/gaussian,/positive,$
;                      perror=perror1,chisq=chisq1,dof=dof1,yerror=yerror1,status=status1)
;
;    ;gd = where(abs(lags-a[1]) lt 3.0*a[2],ngd)
;    ;yfit = MPFITPEAK(lags[gd],xcorr0[gd],par,nterms=5,estimates=estimates,/gaussian,/positive,$
;    ;                  perror=perror,chisq=chisq,dof=dof,yerror=yerror,status=status)
;
;    estimates = [par1[0], 0.0, par1[2], coef1[0], coef1[1], 0.0, 0.0]
;    par2 = MPFITFUN('gausspoly',lags,xcorr0,xcorr0*0+1,estimates,perror=perror2,chisq=chisq2,$
;                    dof=dof2,yerror=yerror2,status=status2,yfit=yfit2,/quiet)
;
;    gd = where(abs(lags-par2[1]) lt 3.0*par2[2],ngd)
;    if ngd eq 0 then goto,BOMB1
;    par3 = MPFITFUN('gausspoly',lags[gd],xcorr0[gd],xcorr0[gd]*0+1,par2,perror=perror3,chisq=chisq3,$
;                    dof=dof3,yerror=yerror3,status=status3,yfit=yfit3,/quiet)
;
;
;    xshift = par3[1]
;    xshiftarr[i] = xshift
;
;    ; Plot the cross-correlation peak and fit
;    if keyword_set(pl) then begin
;      plot,lags,xcorr0,ps=-1,/ysty
;      ;oplot,lags[gd],yfit,ps=-1,co=250
;      oplot,lags[gd],yfit3,ps=-1,co=250
;      print,strtrim(i+1,2),' Shift = ',strtrim(xshift,2),' pixels'
;      wait,0.5
;      ;stop
;    end
;
;    ;stop
;    BOMB1:
;
;  end

  if keyword_set(shiftarr) then shiftarr=xshiftarr
  ; Measure final shift
  gd = where(xshiftarr gt -99)
  ROBUST_MEAN,xshiftarr[gd],shift,shiftsig
  shifterr = shiftsig/sqrt(nfibers*3)
; Printing the results
  print,'Shift = ',strtrim(shift,2),'+/-',strtrim(shifterr,2),' pixels'

  ; do a linear fit to the shifts
  ; don't use blue fibers in superpersistence region
  if keyword_set(nofit) then begin
   shift=[shift,0.]
   goto, fitend
  endif

  bad=where(xshiftarr lt -99 or (chip eq 2 and fiber gt 200),complement=gd)
  shift = AP_ROBUST_POLY_FIT(fiber[gd],xshiftarr[gd],1)
  print,'Fit coefficients: ', shift
  ; Plot all the xcorr shifts
  ;plot = 1 ;1
  if keyword_set(plot) then begin
    if keyword_set(pfile) then begin
      set_plot,'PS'
      file_mkdir,file_dirname(pfile)
      device,file=pfile+'.eps',/encap,/color,xsize=16,ysize=16
      smcolor,/ps
    endif else smcolor
    xr = [0,nfibers*3]
    ;yr = [min(xshiftarr),max(xshiftarr)]
    yr = [-3,3]*mad(xshiftarr)+median(xshiftarr)
  xr=[0,nfibers]
  yr=[-0.6,0.6]
    plot,indgen(nfibers),POLY(findgen(nfibers),shift),xtit='Spectrum #',ytit='Pixel Shift',xr=xr,yr=yr,xs=1,ys=1,thick=3
    oplot,xshiftarr[0:nfibers-1],color=2,ps=1
    oplot,xshiftarr[nfibers:2*nfibers-1],color=3,ps=1
    oplot,xshiftarr[2*nfibers:3*nfibers-1],color=4,ps=1
    ;oplot,[300,300],[-10,10],linestyle=2
    ;xyouts,150,-0.12,'Chip a',align=0.5,charsize=1.3
    ;oplot,[600,600],[-10,10],linestyle=2
    ;xyouts,450,-0.12,'Chip b',align=0.5,charsize=1.3
    ;xyouts,750,-0.12,'Chip c',align=0.5,charsize=1.3
    ;xyouts,mean(xr),yr[1]-0.1*(yr[1]-yr[0]),'Shift = '+strtrim(shift,2)+'+/-'+strtrim(shifterr,2)+' pixels',$
    ;       align=0.5,charsize=1.5
  legend,['Zero '+string(format='(f8.3)',shift[0]),'Slope '+string(format='(e10.2)',shift[1])],textcolor=[1,1],/top,/left

    if keyword_set(pfile) then begin
      device,/close
      ps2gif,pfile+'.eps',/delete,/eps,chmod='664'o
    endif

  endif

  fitend:

;----------------------
; Using EMISSION LINES
;----------------------
Endif else begin

  print,'Using Emission Lines to measure the dither shift'

  chiptag = ['a','b','c']
  apgundef,allmatchstr

  ; Loop through the chips
  For i=0,2 do begin

    print,'Fitting lines for Chip ',chiptag[i]

    str1 = f1.(i)
    str2 = f2.(i)
    ;cube1 = f1.(i).data
    ;cube2 = f2.(i).data

    ; Fit the peaks
    ;   this can find peaks in ThAr or object frames
    print,' Frame 1'
    APPEAKFIT,str1,linestr1
    print,' Frame 2'
    APPEAKFIT,str2,linestr2

    ; Add chip numbers
    ADD_TAG,linestr1,'CHIP',i+1,linestr1
    ADD_TAG,linestr2,'CHIP',i+1,linestr2

    ; Find where the fibers switch
    ;nlines1 = n_elements(linestr1)
    ;diff1 = linestr1.fiber-[0,linestr1[0:nlines1-2].fiber]
    ;break1 = [0,where(diff1 eq 1)]
    ;nlines2 = n_elements(linestr2)
    ;diff2 = linestr2.fiber-[0,linestr2[0:nlines2-2].fiber]
    ;break2 = [0,where(diff2 eq 1)]
    ; if one fiber is missing ALL lines then this won't work

    oneline = linestr1[0]
    STRUCT_ASSIGN,{dummy:0},oneline   ; zero out the contents
    matchstr = replicate({f1:oneline,f2:oneline},50000L)

    ; Now match the lines
    cntlines = 0
    for j=0,nfibers-1 do begin
      gd1 = where(linestr1.fiber eq j,ngd1)
      gd2 = where(linestr2.fiber eq j,ngd2)

      if ngd1 gt 0 and ngd2 gt 0 then begin
        ifiber1 = linestr1[gd1]
        ifiber2 = linestr2[gd2]
        thresh = 1.0
        SRCOR2,ifiber1.fiber,ifiber1.gaussx,ifiber2.fiber,ifiber2.gaussx,thresh,ind1,ind2,opt=1,/silent
        dum = where(ind1 ne -1,nmatch)

        if nmatch gt 0 then begin
          matchstr[cntlines:cntlines+nmatch-1].f1 = ifiber1[ind1]
          matchstr[cntlines:cntlines+nmatch-1].f2 = ifiber2[ind2]
          cntlines += nmatch
        endif
      endif

    end

    ; Trim trailing elements
    matchstr = matchstr[0:cntlines-1]

    ; Add to the total structure
    PUSH,allmatchstr,matchstr

    ;stop

  End  ; chip loop

  ; Calculate the shift
  diffy = allmatchstr.f2.gaussx-allmatchstr.f1.gaussx
  nlines = n_elements(allmatchstr)
  ; don't want lines in superpersistence region, give different shifts
  gdlines = where( (allmatchstr.f2.chip lt 3 or (allmatchstr.f2.chip eq 3 and allmatchstr.f2.fiber lt 200)) and $
                    allmatchstr.f1.gerror[1] lt 2 and allmatchstr.f2.gerror[1] lt 2,ngdlines)
  ROBUST_MEAN,diffy[gdlines],shift,shiftsig
  shifterr = shiftsig/sqrt(ngdlines)

  ; use height and error to weight the values
  ; don't use persistence region for blue chip
  err = sqrt( allmatchstr.f1.gerror[2]^2 + allmatchstr.f2.gerror[2]^2 ) > 0.001
  ; use height > 300 and err lt 0.15 (or some percentile)
  ; for blue chip use fiber < 200

  ; Plot all of the shifts
  if keyword_set(pl) then begin
    xr = [0,n_elements(allmatchstr)]
    yr = [min(diffy),max(diffy)]
    plot,diffy[gdlines],ps=8,sym=0.3,xtit='Emission Line #',ytit='Pixel Shift',xr=xr,yr=yr,xs=1,ys=1
    oplot,xr,[0,0]+shift ,co=250
    xyouts,mean(xr),yr[1]-0.1*(yr[1]-yr[0]),'Shift = '+strtrim(shift,2)+'+/-'+strtrim(shifterr,2)+' pixels',$
           align=0.5,charsize=1.5
  endif

  ; Printing the results
  print,'Shift = ',strtrim(shift,2),'+/-',strtrim(shifterr,2),' pixels'
  ; This should be accurate to ~0.001 if there are ~18 lines per fiber

  ;stop

Endelse  ; using emission lines


if keyword_set(stp) then stop

end
