pro apfindtrace,str,tracestr,npoly=npoly,pl=pl,stp=stp,sigkern=sigkern,$
                nthreshsig=nthreshsig,thresh=thresh0,peakcol=peakcol0

;+
;
; APFINDTRACE
;
; Find aperture traces in a flat field frame
;
; INPUTS: 
;  str         A structure with the 2D image information
;                in tags FLUX, ERR and MASK.  This is
;                the structure loaded by APLOADFRAME but
;                for a single chip.
;  =npoly      The polynomial order to use for the fitting.
;  =sigkern    The Gaussian sigma to use (in pixels) for the convolution
;  =nthreshsig The threshold to use for peak detection in sigma of the
;                 median spectrum
;  =thresh     Threshold to use.
;  =peakcol    The column to use for detection of the fiber peaks.
;                Normally this is the central column.
;  /pl         Do some diagnostic plotting.
;  stp         Stop at the end of the program.
;
; OUTPUTS:
;  tracestr    A structure of trace information for each
;                fiber.
; 
; USAGE:
;  IDL>apfindtrace,str,tracestr
;
; By D.Nidever  Feb. 2010
;-

;file = '/net/stream/apogee/sp2d/ap2D-a-00000003.fits'

apgundef,tracestr,cube,head,message

; Not enough inputs
if n_elements(str) eq 0 then begin
  print,'Syntax - apfindtrace,str,tracestr,npoly=npoly,pl=pl,stp=stp,sigkern=sigkern,'
  print,'                     nthreshsig=nthreshsig,thresh=thresh,peakcol=peakcol'
  return
endif

; Checking the tags of the input structure
needtags = ['FILENAME','HEADER','FLUX','ERR','MASK']
tags = tag_names(str)
for i=0,n_elements(needtags)-1 do begin
  if (where(tags eq needtags[i]))[0] eq -1 then begin
    print,'TAG ',needtags[i],' NOT FOUND in input structure'
    return
  end
end

im = float(str.flux)

;im = reform(cube[*,*,0])
sz = size(im)
nx = sz[1]
ny = sz[2]

print,'Fitting traces for ',str.filename

; Take a sum half-way up
if n_elements(peakcol0) gt 0 then peakcol=(peakcol0>0)<(nx-1) else peakcol=nx/2
halfup = nx/2  ;ny/2
bin = 51  ;101
;mid = TOTAL(im[*,halfup-bin/2:halfup+bin/2],2)/bin
;mid = TOTAL(im[halfup-bin/2:halfup+bin/2,*],1)/bin
;mid = MEDIAN(im[peakcol-bin/2:halfup+bin/2,*],dim=1)
xlo = (peakcol-bin/2) > 0
xhi = (peakcol+bin/2) < (nx-1)
mid = MEDIAN(im[xlo:xhi,*],dim=1)

; Get Sigma of MID
;  sigma by looking at differences of neighboring pixels
;sig0 = MAD( im[halfup-bin/2-1:halfup+bin/2-1,*] - im[halfup-bin/2:halfup+bin/2,*] )
sig0 = MAD( shift(im[xlo:xhi,*],-1) - im[xlo:xhi,*] )
;  get sigma from ERR array
gdpix = where((str.mask and badmask()) eq 0,ngdpix)
if ngdpix gt 1000 then mederr = median(str.err[gdpix]) else mederr=median(str.err)
sig = mederr/sqrt(bin)  ; sigma for MID

; Convolve with a symmetric kernel
;kernel = [0.02,0.1,0.3,1.0,0.3,0.1,0.02]
xx = findgen(21)-10.
if n_elements(sigkern) eq 0 then sigkern=1.5
kernel = exp(-0.5*(xx^2/sigkern^2))      ; Gaussian kernel with sigma=1.5 pixels
;kernel = gfunc(findgen(7)-3,[1.0,0.0,0.5])
mid2 = convol(mid,kernel,/center,/norm)

; Find peaks, must be greater than neighbors and above threshold
;thresh = median(mid2) > mean(mid2)
;thresh = ( median(mid2) > mean(mid2) ) * 0.2 > 7*mad(mid2)
;thresh = ( median(mid2) > mean(mid2) ) * 0.1 > 4*sig   ; this includes the broken fiber
if n_elements(nthreshsig) gt 0 then thresh=nthreshsig*sig else $
  thresh = ( median(mid2) > mean(mid2) ) * 0.5 > 7*sig   ; this excludes the broken fiber
if n_elements(thresh0) gt 0 then thresh=thresh0  ; using input threshold
lft = [0.0,mid2[0:ny-2]]
rgt = [mid2[1:ny-1],0.0]
peak = where(mid2 gt lft AND mid2 gt rgt and mid2 gt thresh,npeak)

print,strtrim(npeak,2),' Peaks found'
if npeak gt 300 then stop

if npeak eq 0 then begin
  print,'No peaks found'
  return
endif

;if npeak lt 300 then stop

;stop

; Get centroids along Y

; Starting in middle and working towards the END
;------------------------------------------------
apgundef,allx,ally,height
npts = 20
;xshift = fltarr(npts)
cenarr = fltarr(npeak,npts)-1e9
matcharr = lonarr(npeak,npts)
sumarr = fltarr(ny,npts)
;step = 2048/(npts+1)
step = nx/(npts+1)
offarr = [ (findgen(npts/2)+1)*step, (findgen(npts/2)+1)*(-step)]+halfup
for i=0,npts-1 do begin
    
  ;off = halfup+(i+1)*step
  off = offarr[i]
  bin = 51
  ;sum = TOTAL(im[*,off-bin/2:off+bin/2],2)/bin
  ;sum = MEDIAN(im[*,off-bin/2:off+bin/2],dim=2)
  sum = MEDIAN(im[off-bin/2:off+bin/2,*],dim=1)
  sum2 = convol(sum,kernel,/center,/norm)
  sumarr[*,i] = sum2

  ; Find peaks, must be greater than neighbors and above threshold
  ;thresh = median(sum) > mean(sum)
  ;thresh = ( median(sum2) > mean(sum2) ) * 0.2 > 7*sig
  if n_elements(nthreshsig) gt 0 then thresh=nthreshsig*sig else $
    thresh = ( median(sum2) > mean(sum2) ) * 0.2 > 7*sig   ; this excludes the broken fiber
  if n_elements(thresh0) gt 0 then thresh=thresh0     ; using input threshold
  sumlft = [0.0,sum2[0:ny-2]]
  sumrgt = [sum2[1:ny-1],0.0]
  ;peak2 = where(sum gt sumlft AND sum gt sumrgt and sum gt thresh,npeak)
  peak2 = where(sum2 gt sumlft AND sum2 ge sumrgt and sum2 gt thresh,npeak2)

  ; Compare to last peaks
  ;-----------------------

  ; First ones, compare to MID
  if abs(off-halfup) eq step then begin

    srcor2,peak,peak*0,peak2,peak2*0,2,ind1,ind2,opt=1,/silent
    dum = where(ind1 ne -1,nmatch)
      if nmatch gt 0 then ind1_match = ind1 else ind1_match=-1

    lastpeak = peak   ; initialize lastpeak with good peaks

  ; Compare to last one
  endif else begin

    gd = where(lastpeak ge 0,ngd)

    ; We had some matches last time
    if ngd gt 0 then begin
      srcor2,lastpeak[gd],lastpeak[gd]*0,peak2,peak2*0,2,ind1,ind2,opt=1,/silent
      dum = where(ind1 ne -1,nmatch)
      if nmatch gt 0 then ind1_match = gd[ind1] else ind1_match=-1
    ; No matches last time, nothing to match it to
    end else begin
      ind1=-1 & ind1_match=-1 & ind2=-1 & nmatch=0
    endelse

  endelse

  ; We have matches
  if nmatch gt 0 then begin
    cenarr[ind1_match,i] = peak2[ind2]
    matcharr[ind1_match,i] = 1
  endif else begin
    ;print,'No matches'
    ;stop
  endelse

  ;if max(abs(peak2[ind2])) gt 3000 then stop
   ;if i eq 14 then stop

  ; Update last peaka for the matches
  if nmatch gt 0 then $
    lastpeak[ind1_match] = peak2[ind2]

  ;print,nmatch

  ;stop

end

;stop

; Now put them all in order
;  the points with off<halfup are reversed, need to flip them
indhi = where(offarr gt halfup,nindhi)
indlo = where(offarr lt halfup,nindlo)
cenarr = [[reverse(cenarr[*,indlo],2)],[cenarr[*,indhi]]]
sumarr = [[reverse(sumarr[*,indlo],2)],[sumarr[*,indhi]]]
matcharr = [[reverse(matcharr[*,indlo],2)],[matcharr[*,indhi]]]

;stop

; Need npoly=4, otherwise the wings are NOT fit well.
if n_elements(npoly) eq 0 then npoly=4
tracestr = REPLICATE({peaky:0.0,wty:0.0,gaussy:0.0,fwhm:0.0,height:0.0,gcoef:fltarr(3),fitcoef:dblarr(npoly+1),$
                      sigcoef:dblarr(npoly+1),rms:-1.0,chisq:-1d0,modelcoef:fltarr(npoly+1),fixed:0,$
                      coef:fltarr(npoly+1)},npeak)

;y = lindgen(ny)
x = lindgen(nx)
halfwid = 3L 

; Loop through each fiber/peak
FOR i=0,npeak-1 do begin

  if (i+1) mod 100 eq 0 then print,i+1,'/',npeak,format='(I3,A1,I3)'

  ; Loop through sub-regions
  ;wtx = fltarr(ny)-1e9
  sum_wty = fltarr(nx)
  num_wty = lonarr(nx)
  For j=0,npts-1 do begin

    if cenarr[i,j] ge 0 then begin

      ; X indices
      ;;lo = peak[i]-3
      ;;hi = peak[i]+3
      ;xlo = cenarr[i,j]-3
      ;xhi = cenarr[i,j]+3
      ;ylo = (j+1)*step-step/2
      ;if j eq 0 then ylo=0
      ;yhi = (j+1)*step+step/2
      ;if j eq npts-1 then yhi=nx-1
      ylo = cenarr[i,j]-halfwid
      yhi = cenarr[i,j]+halfwid
      ;xlo = (j+1)*step-step/2
      xlo = round( (j+1)*step-step*1.5 ) > 0
      if j eq 0 then xlo=0
      ;xhi = (j+1)*step+step/2
      xhi = round( (j+1)*step+step*1.5 )
      if j eq npts-1 then xhi=nx-1

      ; Each fiber is about 4 pixels wide

      ; Fit Gaussian to MID peak
      ;  a = [Height, X, Sigma]
      ;estimates = [max(sumarr[xlo:xhi,j]),cenarr[i,j],1.0]
      ;yfit = MPFITPEAK(lindgen(7)+xlo,sumarr[xlo:xhi,j],a,nterms=3,estimates=estimates,/gaussian,/positive)
      ; this isn't actually used???

      ; Find a flux-weighted X position
      flux = im[xlo:xhi,ylo:yhi] > 1
      ;pos = (lindgen(7)+xlo)#replicate(1.0,yhi-ylo+1)
      ;totwt = TOTAL(flux,1)
      ;totwtx = TOTAL(flux*pos,1)
      ;wtx1 = totwtx/totwt
      ;wtx[ylo:yhi] = wtx1
      pos = replicate(1.0,xhi-xlo+1)#(lindgen(2*halfwid+1)+ylo)
      totwt = TOTAL(flux,2)
      totwty = TOTAL(flux*pos,2)
      wty1 = totwty/totwt
      sum_wty[xlo:xhi] += wty1
      num_wty[xlo:xhi]++

      ;stop

    end

  End ; loop through sub-regions

  ; Take average of multiple measurements
  wty = sum_wty/(num_wty > 1)
  gd = where(num_wty ge 1,ngd)


  ; Fit Gaussian to MID peak
  ;  a = [Height, X, Sigma]
  ylo = cenarr[i,npts/2]-halfwid
  yhi = cenarr[i,npts/2]+halfwid
  if abs(cenarr[i,npts/2]) gt 1e6 then begin
    gdcen = where(abs(cenarr[i,*]) lt 1e5,ngdcen)
    if ngdcen eq 0 then begin
      if not keyword_set(silent) then print,'Cannot trace fiber '+strtrim(i+1,2)
      goto,BOMB
    endif
    ylo = median(cenarr[i,gdcen])-halfwid
    yhi = median(cenarr[i,gdcen])+halfwid
  endif
  estimates = [mid[peak[i]],peak[i],1.0]
  parinfo = replicate({limited:[1,1],limits:[0.0,0.0]},3)
  parinfo[0].limits = [0.5,1.5]*mid[peak[i]]
  parinfo[1].limits = [-0.5,0.5]+peak[i]
  parinfo[2].limits = [0.2,5.0]
  measure_error = sqrt(mid[ylo:yhi]>1)
  yfit = MPFITPEAK(lindgen(2*halfwid+1)+ylo,mid[ylo:yhi],par,nterms=3,estimates=estimates,$
                   measure_error=measure_error,/gaussian,/positive,parinfo=parinfo,$
                   status=status,chisq=chisq,dof=dof)
  ;rchisq = chisq/dof

  ;plot,lindgen(2*halfwid+1)+ylo,mid[ylo:yhi]
  ;oplot,lindgen(2*halfwid+1)+ylo,yfit,co=250

  ; Fit it with a polynomial
  coef = ROBUST_POLY_FIT(x[gd],wty[gd],npoly)
  
  ; Fit CHEBYSHEV function
  ;cheby_coef = CHEBY_FIT(y,wtx,3,nsig=0,smooth=0)
  ; THIS GIVES *TERRIBLE* SOLUTIONS!!!!

  ; Get "good" points
  diff = wty-poly(x,coef)
  rms = mad(diff,/zero)
  ;gd = where(abs(diff)/rms lt 3.0,ngd)
  gd = where(abs(diff)/rms lt 3.0,ngd)


  ; Normal POLY_FIT to get uncertainties
  coef2 = POLY_FIT(x[gd],wty[gd],npoly,sigma=sigma,chisq=chisq)

  ; Stuff it in the trace structure
  tracestr[i].peaky = peak[i]
  tracestr[i].wty = wty[i]
  tracestr[i].gaussy = par[1]
  tracestr[i].fwhm = par[2]*2.35482
  tracestr[i].height = par[0]
  tracestr[i].gcoef = par
  tracestr[i].fitcoef = coef
  tracestr[i].sigcoef = sigma
  tracestr[i].rms = rms
  tracestr[i].chisq = chisq


  ; Plotting
  ;pl=1
  if keyword_set(pl) then begin

    ; Plot the points/fit over the image
    ;lo2 = (lo-5)>0
    ;hi2 = (hi+5)<(nx-1)
    ;x2 = lindgen(hi2-lo2+1)+lo2
    ;displayc,im[lo2:hi2,*],x2,tit='Fiber '+strtrim(i+1,2),xtit='X',ytit='Y'
    ;oplot,wtx,y,ps=1
    ;oplot,poly(y,coef),y,co=250
    ;
    ;wait,1

    !p.multi=[0,1,2]
    !p.charsize=1.5
    yr = [-5,5]+median(wty)
    plot,x,wty,ps=1,xtit='X',ytit='Y',xs=1,ys=1,yr=yr,tit='Fiber '+strtrim(i+1,2)
    oplot,x[gd],wty[gd],ps=1,co=200
    oplot,x,poly(x,coef),co=250,thick=2
    ;oplot,x,cheby(x,cheby_coef,[min(x),max(x)]),co=150
    plot,x,wty-poly(x,coef),ps=1,xtit='X',ytit='Y',yr=[-5.0*rms,5.0*rms],xs=1,ys=1
    oplot,[0,nx],[0.0,0.0],co=250,thick=2
    xyouts,100,3.5*rms,'RMS = '+strtrim(rms,2),charsize=1.5,align=0
    !p.multi=0
    !p.charsize=1.5

    ;print,range(poly(x,coef))

    wait,0.2 ;0.4
    ;stop

  end

  BOMB:

  ;stop

END


; We could fit the polynomial coefficients as a function X/Peak now
; to get it more accurately

; This is where you might want to do a Cholesky decomposition or
; something like that

; The POLY FITS are already QUITE GOOD.  Probably don't need to fix anything

; Plot all of the fits
;if keyword_set(pl) then begin
;
;  plot,[0],[0],/nodata,xr=[0,2047],yr=[0,2047],xs=1,ys=1,xtit='X',ytit='Y'
;  for i=0,npeak-1 do oplot,x,poly(x,tracestr[i].coef)
;  wait,1
;
;end

; Fit the coefficients as a function of Y
;-----------------------------------------
tracestr.coef = tracestr.fitcoef   ; use fitted coeffients for all to start
if n_elements(tracestr) gt 4 then begin
  npolycoef = 2
  coefstr = REPLICATE({coef:fltarr(npolycoef+1),rms:0.0,nbd:0},npoly+1)
  print,'Checking coefficients'
  for i=0,npoly do begin
    if keyword_set(verbose) then $
      print,'Checking COEF=',strtrim(i+1,2)

    cof = ROBUST_POLY_FIT(tracestr.peaky,tracestr.fitcoef[i],npolycoef)
    model = poly(tracestr.peaky,cof)
    rms = MAD(tracestr.fitcoef[i]-model,/zero)

    coefstr[i].coef = cof
    coefstr[i].rms = rms
    tracestr.modelcoef[i] = model  ; save the model

    bd = where(abs(tracestr.fitcoef[i]-model) gt 4*rms,nbd)

    ; Plotting
    ;pl = 1
    if keyword_set(pl) then begin
      plot,tracestr.peaky,tracestr.fitcoef[i],ps=1
      oplot,tracestr.peaky,model,co=250
      if nbd gt 0 then oplot,[tracestr[bd].peaky],[tracestr[bd].fitcoef[i]],ps=1,co=150,sym=2
    endif

    ; Fix bad values
    if nbd gt 0 then begin
      if keyword_set(verbose) then $
        print,'Fixing ',strtrim(nbd,2),' bad value(s)'
      tracestr[bd].coef[i] = model[bd]
      tracestr[bd].fixed = 1
     coefstr[i].nbd = nbd
    end

    ; If one coefficient was fixed for a fiber, then the rest should
    ; probably be as well

  end

  ; For fibers with "fixed" coefs we should use the model values for ALL
  ;  coefficients
  bd = where(tracestr.fixed eq 1,nbd)
  if nbd gt 0 then begin
    for i=0,nbd-1 do tracestr[bd[i]].coef = tracestr[bd[i]].modelcoef
  endif

endif ; fit coefficients

; Remove bad traces
bdtrace = where(tracestr.chisq lt 0.0,nbdtrace)
if nbdtrace gt 0 then begin
  if nbdtrace lt n_elements(tracestr) then REMOVE,bdtrace,tracestr else $
    apgundef,tracestr
endif

; Plot the set of coefficients
if keyword_set(pl) then begin
  !p.multi=[0,2,3]
  !p.charsize=1.5
  plot,tracestr.peaky,tracestr.coef[0],ps=1,xtit='Y',ytit='Coef 0',tit='Coef 0'
  plot,tracestr.peaky,tracestr.coef[1],ps=1,xtit='Y',ytit='Coef 1',tit='Coef 1'
  plot,tracestr.peaky,tracestr.coef[2],ps=1,xtit='Y',ytit='Coef 2',tit='Coef 2'
  plot,tracestr.peaky,tracestr.coef[3],ps=1,xtit='Y',ytit='Coef 3',tit='Coef 3'
  plot,tracestr.peaky,tracestr.coef[4],ps=1,xtit='Y',ytit='Coef 4',tit='Coef 4'
  !p.multi=0
  !p.charsize=1.0
endif

;stop

if keyword_set(stp) then stop

end
