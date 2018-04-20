pro aplincorr,slice_in,lindata,slice_out,stp=stp

sz = size(slice_in)
nreads = sz[2]

szlin = size(lindata)

; a separate coefficient for each output (512 columns)
corr = fltarr(2048,nreads)
npar=szlin[2]

; loop over quadrants
slice_out = slice_in
for i=0,3 do begin
    corr[512*i:512*i+511,*]=lindata[i,0]
    ; use difference with second read as count level
    ;x=alog10(slice_in[512*i:512*i+511,*])
    x=slice_in[512*i:512*i+511,*]
    ;for j=2,nreads-1 do x[*,j]=alog10((x[*,j]-x[*,1])*(j+1.)/(j-1.))
    for j=2,nreads-1 do x[*,j]=(x[*,j]-x[*,1])*(j+1.)/(j-1.)
    bd=where(finite(x) eq 0,nbd)
    if nbd gt 0 then x[bd]=0
    term=x
    for j=1,npar-1 do begin
      corr[512*i:512*i+511,*]+=lindata[i,j]*term
      term*=x
    endfor
    ; set first read correction equal to second
    for j=0,1 do corr[*,j]=corr[*,2]     ;+(corr[*,2]-corr[*,3]) needs to work if only 3 reads!
    ;for j=0,1 do corr[*,j]=1.
endfor
slice_out[0:2047] = slice_in[0:2047]/corr

if keyword_set(stp) then stop
end
