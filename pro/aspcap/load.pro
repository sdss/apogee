pro load,file,d,skip=skip,tell=tell,reject=reject,error=error

;
;	Reads data organized in columns in a file
;
;	IN: file		- string	name of input file
;	OUT: d 			- str array	ncols x nrows elements
;
;
;	KEYWORD: skip		- integer	if >0 this is the number
;						of lines to skip before starting
;						to read data
;
;		tell		- switch	when on, it tells what's read
;
;		reject		- string	set to a string which, when
;						present in a line flags it to
;						be ignored
;
;		error		- integer	0 is returned when the file
;						was read ok, a -1 when there
;						was a problem
;	
;
;	C. Allende Prieto, UT, March 2003
;		"	 , HET, June 2003 - 'reject' added
;			 , UT, Nov 2004 - strbrk substituted by intrinsic 
;					  strsplit
;			 , UT, July 2005 - 'error' added
;			 , IAC, May 2010 -- promoted j loop variable to long 
;
npar = n_params()
if (npar eq 0) then begin
	print,'load,file,d1,d2,d3, ...'
	return
endif 

error=0
;n=wc(file)
n=file_lines(file)
get_lun,lun
openr,lun,file
if keyword_set(skip) then begin
	if (floor(skip) eq skip and skip ge 1) then begin
		header=strarr(skip)
		readf,lun,header
	endif else begin
		print,'% LOAD: Skip has an illegal value'
		error=-1
		return
	endelse
endif else skip=0
s=strarr(n-skip)
readf,lun,s
close,lun
free_lun,lun


;strbrk,s(0),out
out=strsplit(s(0),/extract)
ncols=n_elements(out)
d=strarr(ncols,n)
for j=0l,ncols-1 do begin
	d(j,0)=out(j)
endfor
i=1l
while (i le n-1-skip) do begin
;for i=1l,n-1-skip do begin
	pass=0
	while not pass do begin
		pass=1
		if keyword_set(reject) then begin
			w=strpos(s(i),reject)
			if (max(w) gt -1) then begin
				i=i+1	
				pass=0
			endif
		endif	
	endwhile
	;strbrk,s(i),out
	out=strsplit(s(i),/extract)
	ncols2=n_elements(out)
	if keyword_set(tell) then begin
		print,out
		help,out
	endif
	if (ncols2 ne ncols and strlen(strcompress(out[0],/remove_all)) eq 0) then begin
		print,'% LOAD: the number of columns is NOT constant'
		error=-1
		;stop
		return
	endif 
;	for j=0,ncols-1 do begin
;		d(j,i)=out(j)
;	endfor
        if ncols2 ne ncols then d[*,i]=0. else d[*,i]=out
;endfor
i=i+1
endwhile

;cleanup last blank line(s)
w=-1
while (max(w) eq -1) do begin 
	i=i-1
	d=d[*,0:i]
	w=where(d[*,i] ne '')
endwhile

end
