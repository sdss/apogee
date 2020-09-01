pro	getaa,ndim,np,aa,ntot=ntot

;+
; Generates a matrix with ndim columns and ntot [=sum_{i=1}^{ndim} np(i)]
; rows which can be used to transform ndim nested loops into a single loop
;
;	IN:	ndim	- integer	Number of dimensions
;		np	- intarr	Elements in each dimension
;
;	OUT:    aa	- intarr	2xndim array with indices
;
;	KEYWORDS: ntot  - integer	output computed ntot on demand
;
; Example
;
; for i=0,np[0]-1 do begin
;	for j=0,np[1]-1 do begin
;		for k=0,np[2]-1 do begin
;			print,i,j,k
;		endfor
;	endfor
; endfor
;
; is equivalent to
;
; getaa,3,np,aa
; for i=0,np[0]*np[1]*np[2]-1 do begin
;	print,aa[0,i],aa[1,i],aa[2,i]
; endfor
; 
; C. Allende Prieto, September 2008, adapted from the fortran90 
;			version that I wrote some time ago.
; 				"  , May 2010, upgraded integers to long
;-

if N_Params() lt 2 then begin
	print,'% GETAA: use -- getaa,ndim,np,aa'
	return
endif

np=[0,np]
ntot=long(np(ndim))
for j=2,ndim do begin
	ntot=ntot*long(np(ndim-j+1))
endfor

aa=intarr(ndim,ntot)

;print,'ntot=',ntot

;working arrays
v=lonarr(ndim+1)
c=lonarr(ndim+1)
nr=lonarr(ndim+1)

;initialize counters and values for the indices
;set number of repetitions 
for i=1,ndim do begin
   v(i)=1
   c(i)=0
   nr(i)=1
   for j=i+1,ndim do begin
   	nr(i)=nr(i)*np(j)
   endfor
endfor

for i=1l,ntot do begin
	for j=1,ndim do begin
		c(j)=c(j)+1
		if (c(j) gt nr(j)) then begin
			v(j)=v(j)+1
			c(j)=1
			if (v(j) gt np(j)) then begin
				v(j)=1
			endif
		endif
		aa(j-1,i-1)=v(j)-1
	endfor
	;print,aa(*,i-1)		
endfor

np=np[1:ndim]

end


