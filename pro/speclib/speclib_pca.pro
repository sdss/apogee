pro speclib_pca,file,nvar,addnoise=addnoise

ntot=file_lines(file)
openr,lun,file,/get_lun
line=''
readf,lun,line
npix=n_elements(strsplit(line))
data=fltarr(npix,ntot)
point_lun,lun,0
readf,lun,data
free_lun,lun

;remove the mean from each variable
means=total(data,2)/ntot
data=data-rebin(means,npix,ntot)
;f2 = f2 - rebin (means2,npix2,ntot)

if keyword_set(addnoise) then for i=0,npix-1 do data[i,*]+=randomn(12345,ntot)*addnoise
;compute derived variables based upon the principal components
t0=systime(1)
result=PCOMP(data, COEFFICIENTS = coefficients, nvariables= nvar, $
		EIGENVALUES=v, VARIANCES=variances, /COVARIANCE, /double)
print,'time -> ',systime(1)-t0,' seconds to perform PCA'
	
w = coefficients[*,0:nvar-1]/REBIN(v, npix, npix) 
c = data ## transpose(w*rebin(v,npix,npix))
p = c ## w
		
print, 'Reconstruction error: ', TOTAL((p - data)^2) 
print, 'Energy conservation: ', TOTAL(data^2),TOTAL(v)*(ntot-1) 
print, '     Mode   Eigenvalue  PercentVariance' 
for j=0,nvar-1 do print,j+1,v[j],variances[j]*100.

outfile=file
strput,outfile,'p',strpos(outfile,'n_')
print,'outfile: ',outfile
openw,lun,outfile,/get_lun
printf,lun,means,format='(100000(1x,e20.10))'
printf,lun,v,format='(100000(1x,e20.10))'
for i=0,nvar-1 do printf,lun,w[*,i],format='(100000(1x,e20.10))'
for i=0l,ntot-1 do printf,lun,c[*,i],format='(100000(1x,e20.10))'
free_lun,lun

end
