function appeakfit_p2wcen,xcen,wave,dw=dw

; convert center in pixel units to wavelength
; given WAVE array

apgundef,dw,wcen

npix = n_elements(wave)
xcen_int = floor(xcen)
if xcen_int lt 0 then xcen_int=0
if xcen_int gt (npix-1) then xcen_int=npix-1
xcen_dec = xcen-xcen_int

if xcen_int lt (npix-1) then dw=wave[xcen_int+1]-wave[xcen_int] else $
  dw=wave[npix-1]-wave[npix-2]

; Wavelength at pixel value + decimal * DW
wcen = wave[xcen_int] + dw*xcen_dec

return,wcen

end

;------------------------------------------------------------------

pro appeakfit,str,linestr,wave=wave,nogauss=nogauss,nocont=nocont,error=error,fibers=fibers,smooth=smooth,$
              nsigthresh=nsigthresh,silent=silent,verbose=verbose,pl=pl,count=count

;+
;
; APPEAKFIT
;
; This program finds peaks in a 2Kx300 extracted 2D spectrum
; and fits Gaussians to each one.  The spectrum can be ThAr
; or an object spectrum with airglow lines.
;
; INPUTS:
;  str          A structure that contains the EXTRACTED flux,
;                 variance, and mask 2D arrays with tags of FLUX,
;                 ERR and MASK.
;  =wave        The wavelength array to use with "cube"
;  /nogauss     Don't do any Gaussian fitting.
;  /nocont      Don't remove the continuum.  It was previously removed.
;  =nsigthresh  The sigma threshold for line detection.
;                 The default is 4.
;  /pl          Show plots of the fits.
;  /silent      Don't print anything to the screen.
;  /verbose     Print lots of information to the screen.
;  /stp         Stop at the end of the program.
;
; OUTPUTS;
;  linestr      The line structure.
;  =count       The number of lines found
;  =error       The error message if one occurred.
;
; USAGE:
;  IDL>appeakfit,str,linestr
;
; By D.Nidever  March 2010
;-

; Error Handling
;------------------
; Establish error handler. When errors occur, the index of the  
; error is returned in the variable Error_status:  
CATCH, Error_status 

;This statement begins the error handler:  
if (Error_status ne 0) then begin 
   error = !ERROR_STATE.MSG  
   if not keyword_set(silent) then print,error
   CATCH, /CANCEL 
   return
endif

apgundef,linestr
count = 0

; Not enough inputs
if n_elements(str) eq 0 then begin
  error = 'Not enough inputs'
  if not keyword_set(silent) then begin
    print,'Syntax - appeakfit,str,linestr,wave=wave,nogauss=nogauss,nocont=nocont,silent=silent,'
    print,'                   count=count,nsigthresh=nsigthres,verbose=verbose,error=error,stp=stp'
  endif
  return
endif

; Check that STR is a structure
if size(str,/type) ne 8 then begin
  error = 'STR must be a structure'
  if not keyword_set(silent) then print,error
  return
endif

; Check that we have all the tags we need
needtags2 = ['FLUX','ERR','MASK']
for i=0,2 do begin
  tags2 = tag_names(str)
  for j=0,n_elements(needtags2)-1 do begin
    if (where(tags2 eq needtags2[j]))[0] eq -1 then begin
      if not keyword_set(silent) then $
        print,'TAG ',needtags2[j],' NOT FOUND in input structure'
      return
    end
  end
end

if not keyword_set(smooth) then smooth=0

sz = size(str.flux)
npix = sz[1]
nfibers = sz[2]
if sz[0] eq 1 then nfibers=1

; Starting LINE structure
dumstr = {fiber:-1L,peakx:0.0,gaussx:0.0d0,height:0.0,par0:dblarr(4),gpar:dblarr(4),gpar0:dblarr(4),$
           gerror:dblarr(4),chisq:0.0,rchisq:0.0,dof:0L,rms:0.0,sumflux:0.0,status:0}
linestr = REPLICATE(dumstr,100*nfibers)

cntlines = 0

; if not requesting specific fibers, do them all!
if not keyword_set(fibers) then fibers=indgen(nfibers)

; Find the peaks in the fibers
For ii=0L,n_elements(fibers)-1 do begin
  i=fibers[ii]

  if nfibers gt 1 then begin
    fiber = str.flux[*,i]
    err = str.err[*,i]
  endif else begin
    fiber = str.flux
    err = str.err
  endelse
  ;fiber = reform(cube[*,i,0])
  ;var = reform(cube[*,i,1]) > 0.1  ; variance be > 0
  resid = fiber            ; lines will be removed

  x = findgen(npix)

  ; Subtract out a median-filtered spectrum
  if not keyword_set(nocont) then begin

    ;;BINDATA,y,fiber,xbin1,ybin1,binsize=512,/mode
    ;;coef = ROBUST_POLY_FIT(xbin1,ybin1,4)
    ;coef = ROBUST_POLY_FIT(x,fiber,4)
    ;medfiber = POLY(x,coef)
    ;;medfiber = MEDFILT1D(fiber,300,/edge_copy)
    medfiber = MEDFILT1D(fiber,101,/edge_copy)
    fibersub = fiber-medfiber
    ;fibersub = fibersub-MEDIAN(fibersub)
    ;pos = where(fibersub gt 0.0,npos)
    ;fake = [fibersub[pos],-fibersub[pos]]
    ;sig = MAD(fake)

    ; Fit a skewed Gaussian to the histogram
    ;  to get a better estimate of the continuum
    ;----------------------------------------
    sigfiber = MAD(fibersub) > 1   ; make sure it isn't too small
    bin = round(sigfiber)/10. > 1  ; round to nearest 10th
    ; Get histogram
    hist = HISTOGRAM(fibersub,bin=bin,locations=xbin,min=-5.0*sigfiber,max=5.0*sigfiber)
    xbin = xbin+0.5*bin
    nbin = n_elements(xbin)
    ; Fit regular Gaussian
    ;estimates = [max(hist),xbin[first_el(maxloc(hist))],mad(fibersub)]
    estimates = [max(hist),xbin[first_el(maxloc(hist))],sigfiber]
    gfit1 = MPFITPEAK(xbin,hist,par1,estimates=estimates,nterms=3,/gaussian,/positive)
    ; Calculate the skewness of the distribution
    gd = where(abs(fibersub) lt 5.0*sigfiber,ngd)  ; get the main distribution, better for skew
    skew = SKEWNESS(fibersub[gd])
    if finite(skew) eq 0 then skew=0.0
    ; Fit a skewed Gaussian
    par2 = MPFITFUN('skewgauss',xbin,hist,hist*0+1,[par1,skew],yfit=yfit2,/quiet)
    yfit2 = skewgauss(xbin,par2)
    nb = 5000
    xb = scale_vector(findgen(nb),min(xbin),max(xbin))
    yfit3 = skewgauss(xb,par2)
    maxind = first_el(maxloc(yfit3))
    skewmax = xb[maxind]  ; position of maximum
    halfind = first_el(where(lindgen(nb) gt maxind and yfit3 lt 0.5*max(yfit2)))
    if halfind eq -1 then halfind=nb/2
    skewhalf = xb[halfind]          ; position half way down from maximum on positive side
    ; Now subtract this from "fibersub"
    ;fibersub = fibersub-skewhalf
    ;medfiber = medfiber + skewhalf
    fibersub = fibersub-skewmax
    medfiber = medfiber + skewmax

  ; No continuum
  endif else begin
    medfiber = fiber*0.0
    fibersub = fiber
  endelse

  ; Find peaks, must be greater than neighbors and above threshold
  ;sig = mad(fibersub)
  if smooth gt 0 then smfibersub = GSMOOTH(fibersub,smooth) else smfibersub=fibersub
  ;sig = mad(smfibersub-shift(smfibersub,1))
  sig = mad(fibersub-shift(fibersub,1))
  ;;sig = median(sqrt(var))  ; use the error spectrum
  if n_elements(nsigthresh) eq 0 then nsigthresh=4.0
  thresh = (nsigthresh>0.5)*sig
  lft = [0.0,smfibersub[0:npix-2]]
  rgt = [smfibersub[1:npix-1],0.0]
  xpix = lindgen(npix)
  peak = where(smfibersub gt lft AND smfibersub gt rgt and fibersub gt thresh and fibersub gt 5.0*err and $
               (xpix gt 3 and xpix lt npix-4),npeak)
  ;lft = [0.0,fibersub[0:npix-2]]
  ;rgt = [fibersub[1:npix-1],0.0]
  ;xpix = lindgen(npix)
  ;peak = where(fibersub gt lft AND fibersub gt rgt and fibersub gt thresh and fiber gt 5.0*err and $
  ;             (xpix gt 3 and xpix lt npix-4),npeak)

  if keyword_set(verbose) then $
    print,'Fiber ',strtrim(i+1,2),' - ',strtrim(npeak,2),' Lines found'
  ; Add new elements to LINESTR if necessary
  if cntlines+npeak gt n_elements(linestr) then begin
    linestr = [linestr, REPLICATE(dumstr,(cntlines+npeak)-n_elements(linestr))]
  endif

  ; X-array
  xall = lindgen(npix)

  ; Loop through the peaks
  For j=0L,npeak-1 do begin

    ; Fit with Gaussians
    
    ; X indices
    npixwide0 =  3  ;10 ; 5
    lo0 = (peak[j]-npixwide0)>0
    hi0 = (peak[j]+npixwide0)<(npix-1)
    num0 = hi0-lo0+1


    xx0 = xall[lo0:hi0]
    xx0cen = peak[j]
    dx = 1

    ; Get quantitative estimates of height, center, sigma
    flux = resid[lo0:hi0]-medfiber[lo0:hi0]
    flux -= median(flux)   ; put the median at zero
    flux =  flux > 0       ; don't want negative pixels
    ht0 = max(flux)
    totflux = total(flux)
    linestr[cntlines].sumflux = totflux
    ;  Gaussian area is A = ht*wid*sqrt(2*pi)
    sigma0 = (totflux*dx)/(ht0*sqrt(2*!dpi)) > 0.01
    cen0 = TOTAL(flux*xx0)/totflux
    cen0 = ( (xx0cen-dx*0.5) > cen0 ) < (xx0cen+dx*0.5)   ; constrain the center
    ; Use linear-least squares to calculate height and sigma
    psf1 = exp(-0.5d*(xx0-cen0)^2/sigma0^2) ; normalized Gaussian
    wtht1 = total(flux*psf1)/total(psf1*psf1) ; linear least squares
    ; Second iteration
    sigma1 = (totflux*dx)/(wtht1*sqrt(2*!dpi))
    psf2 = exp(-0.5d*(xx0-cen0)^2/sigma1^2) ; normalized Gaussian
    wtht2 = total(flux*psf2)/total(psf2*psf2)
    ;par0 = [ht0, cen0, sigma0, median(medfiber[lo:hi])]
    par0 = [wtht2, cen0, sigma1, median(medfiber[lo0:hi0])]
    ;par0 = [wtht2, cen0, sigma1, median(resid[lo0:hi0])]

    ; Now get more pixels to fit if necessary
    npixwide =  ceil((2*sigma1)/dx) > 3
    lo = (peak[j]-npixwide)>0
    hi = (peak[j]+npixwide)<(npix-1)
    num = hi-lo+1
    xx = xall[lo:hi]

    ; Gaussian Fitting
    ;-----------------
    if not keyword_set(nogauss) then begin

      estimates = double(par0)
      estimates[0] = fibersub[peak[j]]
      estimates[1] = xx0cen

      ; Fit Gaussian to MID peak
      ;  a = [Height, X, Sigma]
      parinfo = replicate({fixed:0,limits:[0.0,0.0],limited:[1,1]},4)
      parinfo[0].limits = [ 0.5*fibersub[peak[j]], 2.0*fibersub[peak[j]] ]
      ;parinfo[0].limits = [ 0.8, 1.2 ]*estimates[0]
      parinfo[1].limits = [-1.0,1.0]*dx+estimates[1]
      parinfo[2].limits = [0.2*dx, 5.0]
      parinfo[3].limits = [-1, 1]*(3*sig > 0.3*abs(estimates[3])) + estimates[3]

      ;if keyword_set(wave) then begin
      ;  ;estimates = [fibersub[peak[j]], wave[peak[j],i], 0.25, medfiber[peak[j]]]
      ;  ;parinfo[1].limits = [wave[peak[j],i]-1.0,wave[peak[j],i]+1.0]
      ;  parinfo[1].limits = [-0.5,0.5]*dx+estimates[1]
      ;  ;parinfo[2].limits = [0.0, 2]
      ;  parinfo[2].limits = [0.2*dx, 4]
      ;endif else begin 
      ;  ;estimates = [fibersub[peak[j]], peak[j], 1.0, medfiber[peak[j]]]
      ;  ;parinfo[1].limits = [peak[j]-1.0,peak[j]+1.0]
      ;  parinfo[1].limits = [-1.0,1.0]*dx+estimates[1]
      ;  ;parinfo[2].limits = [0.2, 5.0]
      ;  parinfo[2].limits = [0.2*dx, 5.0]
      ;endelse    
      ;parinfo[3].limits = [0.0, 0.2*abs(medfiber[peak[j]])]+medfiber[peak[j]]
      ;parinfo[3].limits = [-0.01, 0.2*abs(medfiber[peak[j]])]+medfiber[peak[j]]
      ;parinfo[3].limits = [-0.4, 0.4]*abs(estimates[3])+estimates[3]
      ;parinfo[3].limits = [-1, 1]*abs(estimates[3])+estimates[3]


      ; Regular Gaussian Fitting
      yfit1 = MPFITPEAK(xx,resid[lo:hi],a,nterms=4,estimates=estimates,error=err[lo:hi],parinfo=parinfo,$
                       /gaussian,/positive,perror=perror1,chisq=chisq1,dof=dof1,yerror=yerror1,status=status1)
      if status1 lt 1 then goto,BOMB

      ; This can sometimes put the peak BETWEEN the pixels, totally wrong!!!!

      ; Now fit with a Binned Gaussian
      func = 'gaussbin'
      pars = MPFITFUN(func,xx,resid[lo:hi],err[lo:hi],a,parinfo=parinfo,yfit=yfit,/quiet,$
                      perror=perror,bestnorm=chisq,dof=dof,status=status)
      if status lt 1 then goto,BOMB

      rms = sqrt(mean((resid[lo:hi]-yfit)^2))
      pcerror = perror * sqrt(chisq/dof)

      ; Remove line from spectrum
      resid_withline = resid  ; resid with this line
      resid[lo:hi] -= yfit

      ; Put values in LINESTR
      ;lineind = cntlines + j
      linestr[cntlines].fiber = i
      linestr[cntlines].chisq = chisq
      linestr[cntlines].rchisq = chisq/dof
      linestr[cntlines].dof = dof
      linestr[cntlines].rms = rms
      linestr[cntlines].status = status1
      linestr[cntlines].height = pars[0]

      ; in pixels
      if not keyword_set(wave) then begin
        linestr[cntlines].peakx = peak[j]
        linestr[cntlines].gaussx = pars[1]
        linestr[cntlines].par0 = par0
        linestr[cntlines].gpar = pars
        linestr[cntlines].gpar0 = a
        linestr[cntlines].gerror = pcerror

      ; in wavelengths
      endif else begin
        wpar0 = par0
        wpar0[1] = appeakfit_p2wcen(par0[1],wave[*,i],dw=dw0)
        wpar0[2] *= abs(dw0)

        wa = a
        wa[1] = appeakfit_p2wcen(a[1],wave[*,i],dw=dw2)
        wa[2] *= abs(dw2)

        wpars = pars
        wpars[1] = appeakfit_p2wcen(pars[1],wave[*,i],dw=dw1)
        wpars[2] *= abs(dw1)

        wpcerror = pcerror
        wpcerror[1] *= abs(dw1)
        wpcerror[2] *= abs(dw1)

        linestr[cntlines].peakx = wave[peak[j]]
        linestr[cntlines].gaussx = wpars[1]
        linestr[cntlines].par0 = wpar0
        linestr[cntlines].gpar = wpars
        linestr[cntlines].gpar0 = wa
        linestr[cntlines].gerror = wpcerror
      endelse

    ; No Gaussian Fitting
    ;---------------------
    endif else begin

      yfit = par0[0]*exp(-0.5*(xx-par0[1])^2/par0[2]^2) + par0[3]

      rms = sqrt(mean((resid[lo:hi]-yfit)^2))

      ; Remove line from spectrum
      resid_withline = resid  ; resid with this line
      resid[lo:hi] -= yfit

      ; Put in LINESTR
      ;lineind = cntlines + j
      linestr[cntlines].fiber = i
      linestr[cntlines].height = par0[0]
      linestr[cntlines].rms = rms

      ; in pixels
      if not keyword_set(wave) then begin
        linestr[cntlines].peakx = peak[j]
        linestr[cntlines].gaussx = par0[1]
        linestr[cntlines].par0 = par0
      ; in wavelengths
      endif else begin
        wpar0 = par0
        wpar0[1] = appeakfit_p2wcen(par0[1],wave[*,i],dw=dw0)
        wpar0[2] *= abs(dw0)

        linestr[cntlines].peakx = wave[peak[j]]
        linestr[cntlines].gaussx = wpar0[1]
        linestr[cntlines].par0 = wpar0
      endelse

    endelse


    ; Plotting
    ;pl=1 ;0 ;1 ;1
    if keyword_set(pl) then begin
      loadct,39
      plot,xx,fiber[lo:hi],yr=[0,max(fiber[lo:hi])],xtit='X',ytit='Counts',$
           tit='Fiber '+strtrim(i+1,2)+' Peak '+strtrim(j+1,2),xs=1,ys=1
      oplot,xx,medfiber[lo:hi],co=80
      oplot,xx,err[lo:hi],co=150
      oplot,xx,resid_withline[lo:hi],co=200
      oplot,xx,resid[lo:hi],co=100
      oplot,xx,yfit,co=250
      legend,['Original spectrum','Continuum','Error','Resid With Line','After Line Removed','Fit'],$
             textcolor=[255,80,150,200,100,250],/top,/left
      wait,0.5
    endif

    ; Increment the counter
    cntlines++

    ;stop

    BOMB:

  End

  ;stop

End

; Trim extra lines
linestr = linestr[0:cntlines-1]

; How many lines per fiber
;nlines = histogram(linestr.fiber,bin=1,min=0,max=nfibers-1)
nlines = lonarr(n_elements(fibers))
for i=0,n_elements(fibers)-1 do begin
  dum = where(linestr.fiber eq fibers[i],nline)
  nlines[i] = nline
endfor

if not keyword_set(silent) then begin
  print,'Median NLINES per fiber = ',strtrim(long(median(nlines)),2)
  ;print,''
endif

count = n_elements(linestr)

;stop

if keyword_set(stp) then stop

end
