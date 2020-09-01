function histograma,y,nbins=nbins,binsize=binsize,x=x,cum=cum,ind=ind,forcex=forcex

;+
;
;	Calculates an histogram of the vector y.
;
;	IN: y	- fltarr	vector
;
;	KEYWORDS:	nbins - number of bins to be used (default is 10)
;			    This value is ignored if binsize is provided.
;
;		binsize	- size of the bin to use (otherwise chosen to span the
;					range of y with 10 elements)
;			x	- vector with the centers of the intervals
;				  used to calculate the histogram. An element of
;				  y is included in bin i when 
;
;				x(i)-binsize/2. <= y(i) < x(i)+binsize/2.
;
;			This array is internally created, unless forcex is set,
;				in which case the input array is adopted.
;
;		cum - when switched on, a cumulative distrib. is returned
;
;		ind - an optional array with the indices of the elements
;			in the input array in each bin (padded with -1s).
;
;			forcex -- set this keyword to impose an input x array,
;		otherwise x is internally created and replaces the input one.
;
;	Carlos Allende Prieto, UT, Sep 1999
;			     , UT, April 2005, changed to use long integers
;				 , A&M, October 2010, added cum keyword
;			 , Radazul, June 2011, added nbins and ind keywords, 
;			changed nbins default and allowed x to be an input
;				even if it doesn't contain all elements of y,
;			as long as x is equidistant, when setting xforce.
;-
;

if keyword_set(forcex) then begin
	nbins=n_elements(x)
	if nbins lt 2 then begin
		print,'% HISTOGRAMA: The intput x array must have at least 2 elements. I quit!'
		return,-1	
	endif
	binsize=x[1]-x[0]
	delta=x-shift(x,1)
	delta=delta[1:n_elements(delta)-1]
	if max(delta)-min(delta) gt 1d-3*mean(delta) then begin
		print,'% HISTOGRAMA: The intput x array is not equidistantly spaced. I quit!'
		stop
		return,-1
	endif
endif else begin
	if not keyword_set(nbins) then nbins=10
	if not keyword_set(binsize) then binsize=(max(y)-min(y))/float(nbins) else $
		nbins=ceil((max(y)-min(y))/binsize)
	left=nbins*binsize-(max(y)-min(y))
	if (nbins eq 0) then nbins=1
	x=findgen(nbins)*binsize+min(y)-left/2.+binsize/2.		
endelse

hist=lonarr(nbins)
if arg_present(ind) then begin
	ind=lonarr(nbins,n_elements(y),/nozero)
	ind[*]=-1
endif

for i=0l,nbins-1 do begin
	pos=where(y ge x(i)-binsize/2. and y lt x(i)+binsize/2.)
	if (i eq nbins-1) then $
	pos=where(y ge x(i)-binsize/2. and y le x(i)+binsize/2.)
	if (max(pos) ne -1) then hist(i)=n_elements(pos)
	if keyword_set(cum) and i gt 0 then hist(i)=hist(i-1)+hist(i)
	if n_elements(ind) gt 0 then ind[i,0:n_elements(pos)-1]=pos
endfor
return,hist
end
