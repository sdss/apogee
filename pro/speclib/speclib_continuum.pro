pro speclib_continuum,order,niter,lowrej,highrej,y,ycont,mask=mask,interpolate=interpolate,x=x,ivar=ivar

;+
;	Rough non-interactive spectral normalization fitting a order-th order
;	polynomial to the highest points in vector y in an iterative manner.
;
;	IN: order - integer - order of the polynomial to fit the continuum
;							-1<order<0 returns ycont equal to y processed
;								through a -order*100 percentile filter 
;								(median filter when order=-0.5)
;							order=-1 indicates returns ycont with all elements 
;								equal to the median 
;							order=-2 indicates that we are using pem.pro
;							order=-3 indicates we are using a running average\
;						        -200<order<-100. indicates we're using pem with 
;								percentile=-(100.+order) (order=-150. for median)
;	    niter - integer - number of iterations (cont >=0)
;					    - width of filter for median/percentile (-1<order<0)
;						- no meaning (order=-1)
;						- number of pieces (order=-2, i.e. pem) 
;						- width of the pieces (order = -3)
;	   lowrej - float   - a point is rejected when goes lowrej times sigma 
;				under the fit
;	  highrej - float   - a point is rejected when goes highrej times sigma
;				above the fit
;		y - fltarr  - spectrum to fit
;
;	OUT:ycont - fltarr  - continuum
;
;	KEYWORDS: mask -- integer array of the same size as y, flagging with 0s
;					the positions that should be masked (not used for fitting)
;
;			 interpolate -- when a mask is present this keyword will cause
;					the masked data to be replaced by linear interpolation rather
;					than simply ignored
;			
;			 x -- this can be used for entering an array of abscissae,
;					when other that the default ([0,1,2...]) is desired. Its
;					dimension must match that of y. This keyword can also be
;					used to output the actual abscissae used in the fitting.
;
;			ivar -  this can be used for entering an array with inverse variance values
;				for y. When present, these will be used in the fittings when order>=0
;
;	NOTE: It is assumed that the frequency/energy/wavelength axis is
;		evenly spaced.
;	
;	C. Allende Prieto, UT, Aug 1999
;			 , UT, Oct 2006, modified to keep using the previous
;					poly_fit, as that with a newer IDL 
;					version fails with with continuum.pro
;			 , IAC, Feb 2011, adapted it so that order<0 indicates the median
;			 , IAC, April 2011, added 'mask' and 'interpolate' keywords
;			 , IAC, May 2011, added median/percentile filter
;			 , IAC, July 2011, added pem-anomaly and smooth normalizations and x keyword
;			 , IAC, February 2012, force /interpolate when using a mask in pem mode
;			 , IAC, March 2016, added ivar keyword
;-
if N_params() LT 5 then begin
      print,'% continuum: continuum,order,niter,lowrej,highrej,y,ycont[,mask=mask,interpolate=interpolate,x=x,ivar=ivar]'
      return
endif
n=n_elements(y)

;errors
if n_elements(ivar) ne 0 then begin
if n_elements(ivar) ne n_elements(y) then begin
    print,'% CONTINUUM: warning -- ivar does not have the same size as y -- ivar will be ignored'
    undefine,ivar
endif
endif

if n_elements(ivar) gt 0 then begin
  macho=machar()
  if isa(y,'Double') then macho=machar(/double)
  tiny=macho.eps
  wbad=where(ivar le tiny)
  if max(wbad) gt -1 then ivar[wbad]=2.*tiny
  err=sqrt(1.d0/ivar)
endif

;deal with mask 
if (n_elements(mask) eq 0) then mask=replicate(1,n) else begin
	if n_elements(mask) ne n then begin
		print,'% continuum: the arrays mask and y have not the same number of elements!'
		return
	endif
	if max(where(mask ne 0)) lt 0 then begin
		print,'% continuum: the array mask has no elements different from zero!'
		return	
	endif
	; must interpolate when -100>order>=-10 (pem mode) 
	if (order gt -100. and order le -10.) then interpolate=1 
endelse
if (n_elements(x) gt 0) then begin
	if n_elements(x) ne n then begin
		print,'% continuum: the arrays x and y have not the same number of elements!'
		return
	endif
endif else begin
	x=findgen(n)
endelse 

xx=x[where(mask ne 0)]
yy=y[where(mask ne 0)] 

if n_elements(ivar) gt 0 then err=err[where(mask ne 0)]
if keyword_set(interpolate) then begin
	yy=interpol(yy,xx,x)
	if n_elements(ivar) gt 0 then err=interpol(err,xx,x)
	xx=x
endif
m=n_elements(yy)
if (m le n/10.) then begin
	print,'% CONTINUUM: Warning-- the data points not masked are less than 10% of all available!'
endif
if (m lt order+1) then begin
	print,'% CONTINUUM: Error-- the data points not masked are less than order+1, I quit!'
	return
endif

for i=0,niter-1 do begin
	if order ge 0 then begin
		if n_elements(ivar) gt 0 then $
		  coef=poly_fit(xx,yy,order,yfit,yband,sigma,measure_errors=err) else $
		  coef=poly_fit(xx,yy,order,yfit,yband,sigma) 
		if n_elements(yfit) eq 1 and order eq 0 then yfit=replicate(yfit,n_elements(yy))
	endif else begin
	        if order gt -1. then yfit=cmedian(yy,niter,percentile=-order) else $
		case order of 
		-1: begin
			yfit=yy	
			yfit[*]=median(yy)
		    end
		-2: pem,yy,niter,yfit 
		-3: yfit=smooth(yy,niter,/edge_truncate) 
		else: begin
		         if (order gt -200.) and (order lt -100.) then begin
			   pem,yy,niter,yfit,percentile=-(100.+order)/100.
		         endif else begin
			   print,'% CONTINUUM: Error invalid value for order=',order
		           return
		         endelse
		      end
		endcase
		break
	endelse
	for j=0l,m-1 do begin		
		if (yy(j) lt yfit(j)-lowrej*sigma or $
		yy(j) gt yfit(j)+highrej*sigma) then $
		yy(j)=yfit(j)
	endfor
	;oplot,findgen(n),y
endfor

if order ge 0 then begin
	;coef=poly_fit(xx,yy,order,yfit,yband,sigma)
	ycont=poly(x,coef)
endif else begin
	ycont=yfit
endelse

end
