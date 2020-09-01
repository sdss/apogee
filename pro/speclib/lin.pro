pro 	lin,grid,x,y,normalized=normalized,indi=indi

;+
;	Multilinear interpolation in N dimensions
;
;	INPUT:	grid	-fltarr		grid to interpolate on 
;		x	-fltarr		vector of desired coordinates [indices!]
;
;	OUTPUT: y	-float		result
;
;	KEYWORD:	normalized	signals that the x vector refers to
;					normalized coords. (0-1)
;
;			indi		vector with the order for the 
;					interpolations (default is
;					indi=[ndim,ndim-1 ...])
;
;	C. Allende Prieto, March 2005
;-

npar = n_params()
if (npar eq 0) then begin
	print,'% LIN: lin,grid,x,y[,normalized=normalized,indi=indi]'
	return
endif 


s=size(grid)
ndim=s[0]
np=s[1:ndim]

if not keyword_set(indi) then begin
	indi=intarr(ndim)
	for i=1,ndim do begin
		indi[i-1]=ndim-i+1
	endfor
endif

base=2
ee=intarr(ndim,base^ndim)
for i=1,ndim do begin
	nr=base^(i-1)
	c1=1
	for j=1,nr do begin
		c2=1
		for k=1,base^ndim/nr do begin
			;ee(i-1,c1-1)=(c2-1)/base^(ndim-i)
			ee(ndim-indi[i-1],c1-1)=(c2-1)/base^(ndim-i)
			c1=c1+1
			c2=c2+1	
		endfor
	endfor
endfor

;print,'indi=',indi
;print,'ee=',ee

if keyword_set(normalized) then begin
	t=x*(np-1.)	    ; 0-np(i)	indices)
endif else begin
	t=x
endelse


ntimes=lonarr(ndim)
ntot=long(np(ndim-1))
ntimes(ndim-1)=1
for i=2,ndim do begin
	ntimes(ndim-i)=ntot
	ntot=ntot*long(np(ndim-i))
endfor


wrk=dblarr(2^ndim)

factor=lonarr(ndim)
factor[*]=1
for i=1,ndim-1 do begin
	for j=0,i-1 do begin
		factor[i]=factor[i]*long(np[j])
	endfor
endfor

	
for i=1,2^ndim	do begin
	indices=[floor(t)+ee(*,i-1)]
	w=where(indices eq np)	; taking care of points on end edges
	if (max(w) gt -1) then indices[w]=indices[w]-1 
	index=long(total(factor*indices))
	wrk(i-1)=grid[index]		
endfor

;interpolate
for i=1,ndim	do begin
	delta=(t(indi(i-1)-1)-floor(t(indi(i-1)-1)))
	for j=1,2^(ndim-i) 	do begin
		wrk(j-1)=(wrk(2*j-1)-wrk(2*j-1-1))* $
			;(t(ndim-i)-floor(t(ndim-i)))+ wrk(2*j-1-1)
			;(t(indi(i-1)-1)-floor(t(indi(i-1)-1)))+ wrk(2*j-1-1)
			delta + wrk(2*j-1-1)
	endfor
endfor
y=wrk[0]

;stop

end		