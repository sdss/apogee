function cmedian,array,width,percentile=percentile
;+
;	Extension of the IDL intrinsic median function (median filter) 
;	to use other percentiles. By default percentile=0.5 and
;	calculates the median. 
;	
;	IN: array - 			Input data array
;		width -  integer 	Width over which the median filter is applied
;
;	OUT: array after applying the median filter
;	
;	KEYWORD: percentile		Desired percentile (0.5 for median)
;
;
;	NOTE: Unlike the IDL intrinsic median
;	even if an array with an even number of elements is passed
;	the result is still an element of the array - after sorting
;	the returned element is that with the index 
;	floor(number_of_data_points*percentile). The behaviour in the
;	edges is also different from the intrinsic 'median'.
;
;	C. Allende Prieto, IAC,  May 2011
;	
;-

if N_params() lt 1 then begin
	 print,'% CMEDIAN: use -- result=cmedian(array,[width,percentile=percentile])'
	 return,0
endif

nel=n_elements(array)

if not keyword_set(width) then width=nel
if not keyword_set(percentile) then percentile=0.5
if width gt nel then begin
	print,'% CMEDIAN: Error -- width must be smaller than the number of elements in array!'
	return,0
endif

if width eq nel then begin
	np=n_elements(array)
	wdata=sort(array)
	result=array[wdata[floor(np*percentile)]]
endif else begin
	result=array
	result[*]=0.
	for i=0,nel-1 do begin
		;result[i]=median(array[max([i-width/2,0]):min([i+width/2,nel-1])])
		data=array[max([i-width/2,0]):min([i+width/2,nel-1])]
		np=n_elements(data)
		wdata=sort(data)
		result[i]=data[wdata[floor(np*percentile)]]
	endfor
endelse

return,result
end



