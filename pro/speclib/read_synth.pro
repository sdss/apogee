pro read_synth,synthfile,x,wav,flux,grid=grid,ax=ax,hd=hd,npix=npix,ntot=ntot,ndim=ndim,$
	lambda=lambda,nnpix=nnpix,continuum=continuum,wave=wave,logw=logw,resolution=resolution,$
	npca=npca,means=means,vv=vv,ww=ww,constant=constant,vacuum=vacuum,$
	precontinuum=precontinuum,param=param

;+
;	Reads and interpolates in a synth grid
;
;	IN: synthfile	- string	Name of the grid
;	    x	 	- fltarr	Vector of the parameters for which
;					 we wish to obtain interpolated fluxes
;
;	OUT: wav	- fltarr	Vector of wavelengths
;	     flux	- fltarr	Vector of interpolated fluxes
;	
;	KEYWORDS: grid	- fltarr	extracts the full grid through this 
;					variable; if provided from a pparamrevious
;					call, the reading will be avoided
;		  ax	- struct	axes info. Tags are: 
;					label: axis label (string)
;					n : axis number of elements (integer)
;					llimit:lower limit (first element)(string)
;					step: step (string) 
;
;					-- for multi-grids the information of the first header 
;						with information is adopted
;							
;
;		  hd	- strarr	full header
;
;		  npix  - long  	number of pixels (spectral direction)
;
;		  ntot  - long  	total of elements in other dimensions 
;							[= product(np)]
;
;		  ndim  - long		number of dimensions 
;
;		lambda  - dblarr	Array with wavelengths
;
;		vacuum - integer    0/1 for vacuum/std. air scales
;
;		param - fltarr		2D array with parameters: as many columns as
;					parameters (ndim) and many rows as models (ntot)
;
;		 MULTI keywords  ------------------------------------------------------
;
;			nnpix  intarr	number of pixels (for each header in MULTI headers)
;		continuum  fltarr	continuum normalization info (MULTI x 4)
;		precontinuum  fltarr	precontinuum normalization info (MULTI x 4) 
;					(applied to the data before the continuum correction)
;		    wave   fltarr	wavelength information (MULTI x 2)
;		    logw   intarr      wavelength(0), log10(wavelength), 
;					or ln(wavelength) scale (MULTI)
;		resolution fltarr	resolving power (MULTI)
;
;
;		 NPCA keywords  ------------------------------------------------------
;
;		constant -double	constant that has been added to data
;							(it is not subtracted here, only read)
;		npca	 intarr		array with numbers of pca components per section
;		means				cean flux array (npix)	-- only for npca files
;			vv				array of eigenvalues (of the covariance matrix) 
;							(nvar) -- only npca files
;			ww				array of first nvar eigenvectors (npix x nvar)
;										-- only npca files
;
;
;
;	NOTE: An axis can be composed by 
;	axis1=findgen(ax[0].n)*ax[0].step+ax[0].llimit
;
;
;	C. Allende Prieto, UT, March 2005
;					 , IAC, August 2010 -- using a transposed grid for efficiency
;					 , IAC, November 2010 -- adapted to handle multi grids
;					 , IAC, November 2010 -- adapted to handle npca grids
;					 , IAC, August 2011 -- altered to deal with grids without CONTINUUM
;					 , IAC, Nov 2011 -- changed to handle TRANSPOSED grids
;					 , IAC, Dec 2011 -- added multi keyword continuum
;					 , IAC, Feb 2012 -- changed ax.llimit/ax.step to string variables
;					 , Radazul, July 2012, introduce logw=2
;					 , Radazul, Sept 2012, add precontinuum keyword	
;					 , IAC, April 2016, added param keyword
;-

npar = n_params()
if (npar eq 0) then begin
	print,'% READ_SYNTH: read_synth,synthfile,x,wav,flux[,grid=grid,ax=ax,hd=hd,$'
	print,'% READ_SYNTH:   npix=npix,ntot=ntot,nnpix=nnpix,continuum=continuum,$'
	print,'% READ_SYNTH:   npca=npca,means=means,vv=vv,ww=ww,constant=constant,$'
	print,'% READ_SYNTH:   precontinuum=precontinuum,param=param]'
	return
endif

openr,lun,synthfile,/get_lun

;reading the header
hola='hola'
header='HEADER'
label='label'

multi=0
npix=0
vacuum=-1
while not eof(lun) and strmid(hola,1,1) ne '/' do begin

	readf,lun,hola
	header=[header,hola]
	;print,hola
	
	if (strpos(hola,'NPCA') gt -1) then $
		npca=long(strsplit(strmid(hola,7,strlen(hola)-7),/extract))
		
	if (strpos(hola,'CONSTANT') gt -1) then $
		constant=double(strmid(hola,11,strlen(hola)-11))
		
	if (strpos(hola,'MULTI') gt -1) then multi=fix(strmid(hola,8,strlen(hola)-8))
	
	if (strpos(hola,'ID') gt -1) then id=strmid(hola,5,strlen(hola)-5)
	
	if (strpos(hola,'N_OF_DIM') gt -1) then $
		ndim=fix(strmid(hola,11,strlen(hola)-11))
	
	if (strpos(hola,'N_P') gt -1) then $
		np=fix(strsplit(strmid(hola,6,strlen(hola)-6),/extract))
	
	if (strpos(hola,'NPIX') gt -1) then begin
		npix=long(strmid(hola,7,strlen(hola)-7))
		nnpix=npix
	endif
	
	if (strpos(hola,'LABEL') gt -1) then begin
		label=[label,strmid(hola,12,strlen(hola)-12)]
	endif
	
	if (strpos(hola,'LLIMITS') gt -1) then $
		llimits=double(strsplit(strmid(hola,10,strlen(hola)-10),/extract))
		
	if (strpos(hola,'STEPS') gt -1) then $
		steps=double(strsplit(strmid(hola,8,strlen(hola)-8),/extract))
		
	if (strpos(hola,'WAVE') gt -1) then $
		wave=double(strsplit(strmid(hola,7,strlen(hola)-7),/extract))

	if (strpos(hola,'LOGW') gt -1) then $
		logw=fix(strmid(hola,7,strlen(hola)-7))

	if (strpos(hola,'VACUUM') gt -1) then $
		vacuum=fix(strmid(hola,9,strlen(hola)-9))

	if (strpos(hola,'RESOLUTION') gt -1) then $
		resolution=double(strmid(hola,13,strlen(hola)-13))

	if (strpos(hola,'PRECONTINUUM') gt -1) then begin
		precontinuum=double(strsplit(strmid(hola,15,strlen(hola)-15),/extract))
	endif else begin
		if (strpos(hola,'CONTINUUM') gt -1) then $
			continuum=double(strsplit(strmid(hola,12,strlen(hola)-12),/extract))
	endelse
		
endwhile

;check that we dont' mix npca with multi
if n_elements(npca) gt 0 and multi gt 1 then begin
	print,'% READ_SYNTH: Warning      --  This is an NPCA grid -- '
	print,'% READ_SYNTH: --> npix will be adopted from the first header'
endif

;reset npix if this is a non-npca multi header
;(should not be there, but to be safe)
if n_elements(npca) eq 0 and multi gt 1 then npix=0


;reading multiple headers if present
for i=1,multi do begin
  hola='hola'
  label2='label'
  npix2=0
  vacuum2=-1
  while not eof(lun) and strmid(hola,1,1) ne '/' do begin
	readf,lun,hola
	
	header=[header,hola]
	
	if (strpos(hola,'ID') gt -1) then id2=strmid(hola,5,strlen(hola)-5)
	
	if (strpos(hola,'N_OF_DIM') gt -1) then $
		ndim2=fix(strmid(hola,11,strlen(hola)-11))
	
	if (strpos(hola,'N_P') gt -1) then $
		np2=fix(strsplit(strmid(hola,6,strlen(hola)-6),/extract))
	
	if (strpos(hola,'NPIX') gt -1) then $
		npix2=long(strmid(hola,7,strlen(hola)-7))
		
	if (strpos(hola,'LABEL') gt -1) then begin
		label2=[label2,strmid(hola,12,strlen(hola)-12)]
	endif
	
	if (strpos(hola,'LLIMITS') gt -1) then $
		llimits2=double(strsplit(strmid(hola,10,strlen(hola)-10),/extract))
		
	if (strpos(hola,'STEPS') gt -1) then $
		steps2=double(strsplit(strmid(hola,8,strlen(hola)-8),/extract))
		
	if (strpos(hola,'WAVE') gt -1) then $
		wave2=double(strsplit(strmid(hola,7,strlen(hola)-7),/extract))

	if (strpos(hola,'LOGW') gt -1) then $
		logw2=fix(strmid(hola,7,strlen(hola)-7))
		
	if (strpos(hola,'VACUUM') gt -1) then $
		vacuum2=fix(strmid(hola,9,strlen(hola)-9))		
	
	if (strpos(hola,'RESOLUTION') gt -1) then $
		resolution2=double(strmid(hola,13,strlen(hola)-13))

	if (strpos(hola,'PRECONTINUUM') gt -1) then begin
		precontinuum2=double(strsplit(strmid(hola,15,strlen(hola)-15),/extract))
	endif else begin
		if (strpos(hola,'CONTINUUM') gt -1) then $
			continuum2=double(strsplit(strmid(hola,12,strlen(hola)-12),/extract)) 
	endelse
	
  endwhile
  if i eq 1 then begin
  	;track info
	label=label2
	ndim=ndim2
	np=np2
	llimits=llimits2
	steps=steps2
	nnpix=npix2
	if n_elements(continuum2) gt 0 then continuum=continuum2
	if n_elements(precontinuum2) gt 0 then precontinuum=precontinuum2
	if n_elements(npca) eq 0 then npix=npix2
	wave=wave2
	logw=logw2
	vacuum=vacuum2
	resolution=resolution2
  endif else begin
  	;check and track
  	if max(ndim2 - ndim) gt 0 or $
  		max(np2-np) gt 0 or $
  		max(llimits2-llimits) gt 0 or $
  		max(steps2-steps) gt 0 then begin
  			print,'% READ_SYNTH: ERROR - inconsistent MULTI header!'
  			return
  	endif
  	nnpix=[nnpix,npix2]
  	if n_elements(continuum2) gt 0 then begin
  		if n_elements(continuum) gt 0 then continuum=[[continuum],[continuum2]] else begin
	  		continuum=continuum2
	  		print,'% READ_SYNTH: WARNING - The  CONTINUUM keyword is not consistently included in all MULTI headers!'
	  	endelse
	endif
  	if n_elements(precontinuum2) gt 0 then begin
  		if n_elements(precontinuum) gt 0 then precontinuum=[[precontinuum],[precontinuum2]] else begin
	  		precontinuum=precontinuum2
	  		print,'% READ_SYNTH: WARNING - The  PRECONTINUUM keyword is not consistently included in all MULTI headers!'
	  	endelse
	endif
	if n_elements(npca) eq 0 then npix=npix+npix2
  	wave=[[wave],[wave2]]
  	logw=[logw,logw2]
  	vacuum=[vacuum,vacuum2]
  	resolution=[resolution,resolution2]
  endelse
endfor

;build wavelength array for output
if n_elements(nnpix) eq 1 then begin
	lambda=findgen(nnpix[0])*wave[1]+wave[0]
	if logw eq 1 then lambda=10.d0^lambda
	if logw eq 2 then lambda=exp(lambda)
endif else begin
	lambda=findgen(nnpix[0])*wave[1,0]+wave[0,0]
	if logw[0] eq 1 then lambda=10.d0^lambda
	if logw[0] eq 2 then lambda=exp(lambda)
	for i=1,n_elements(nnpix)-1 do begin
		lambda2=findgen(nnpix[i])*wave[1,i]+wave[0,i]
		if logw[0] eq 1 then lambda2=10.d0^lambda2
		if logw[0] eq 2 then lambda2=exp(lambda2)
		lambda=[lambda,lambda2]
	endfor
endelse

;remove first element of the label vector and the header
hd=header[1:n_elements(header)-1]
label=label[1:n_elements(label)-1]

;report to user
print,'% READ_SYNTH: Grid id is ',id
print,'% READ_SYNTH: --> NP= ',np
print,'% READ_SYNTH: --> NPIX = ',npix

ntot=long(product(np))

;define axes
ax=replicate({label:label[0],n:np[0],llimit:llimits[0],step:steps[0]},ndim)
for i=1,ndim-1 do begin
	ax[i].label=label[i]
	ax[i].n=np[i]
	ax[i].llimit=llimits[i]
	ax[i].step=steps[i]
endfor

;define param array
getaa,ndim,np,aa
param=fltarr(ndim,ntot)
for i=0,ndim-1 do param[i,*]=aa[i,*]*ax[i].step+ax[i].llimit

;if npca grid, then read supplementary data
if n_elements(npca) gt 0 then begin
	means=dblarr(total(npca))
	readf,lun,means
	vv=dblarr(total(npca))
	readf,lun,vv
	vv=transpose(vv)
	ww=dblarr(total(npca),npix/n_elements(npca))
	readf,lun,ww
endif

;if grid is provided, no need to read the data again
if (not keyword_set(grid) or n_elements(grid) eq 0) then begin

	;transposed?
	transposed=0
	wt=strpos(header,'TRANSPOSED')
	tpos=where(wt gt -1)
	if n_elements(tpos) gt 1 then begin
		print,'% READ_SYNTH: Error, more than one TRANSPOSED keyword in the header!'
		return
	endif 
	if tpos[0] gt -1 then begin
		tval=header[tpos]
		wt=strpos(header[tpos],'=')
		tval=strmid(header[tpos],wt+1,100)
		if fix(tval) eq 1 then transposed=1
	endif	

	;define grid
	if transposed then begin
		print,'% READ_SYNTH: Found TRANSPOSED=1 in the header of the grid'
		grid=fltarr(ntot,npix)
	endif else begin
		if ndim lt 8 then grid=fltarr([npix,reverse(np)]) else begin
			if n_elements(x) gt 0 then begin
				print,'% READ_SYNTH: IDL cannot handle arrays with more than 8 dimensions'
				print,'%             You cannot interpolate in the grid, but can still '
				print,'% 			 read the data into a 2D array, omitting the array x'
				return
			endif
			grid=fltarr(npix,ntot)
			print,'% READ_SYNTH: IDL cannot handle arrays with more than 8 dimensions'
			print,'% 				A 2D array is returned (npix,ntot)=(',npix,',',ntot,')'			
		endelse
	endelse

	print,'% READ_SYNTH: Reading data ...'
	;now get the numbers
	readf,lun,grid
	
	;transpose if transposed
	if transposed then begin
		grid=transpose(grid)
		rnp=reverse(np)
		case ndim of
			1: grid=reform(grid,npix,ntot,/overwrite) ; do nothing!
			2: grid=reform(grid,npix,rnp[0],rnp[1],/overwrite)
			3: grid=reform(grid,npix,rnp[0],rnp[1],rnp[2],/overwrite)
			4: grid=reform(grid,npix,rnp[0],rnp[1],rnp[2],rnp[3],/overwrite)
			5: grid=reform(grid,npix,rnp[0],rnp[1],rnp[2],rnp[3],rnp[4],/overwrite)
			6: grid=reform(grid,npix,rnp[0],rnp[1],rnp[2],rnp[3],rnp[4],rnp[5],/overwrite)
			7: grid=reform(grid,npix,rnp[0],rnp[1],rnp[2],rnp[3],rnp[4],rnp[5],rnp[6],/overwrite)
			else: begin
				if n_elements(x) gt 0 then begin
					print,'% READ_SYNTH: IDL cannot handle arrays with more than 8 dimensions'
					print,'%             You cannot interpolate in the grid, but can still '
					print,'% 			 read the data into a 2D array, omitting the array x'
					return
				endif
				print,'% READ_SYNTH: IDL cannot handle arrays with more than 8 dimensions'
				print,'% 				A 2D array is returned (npix,ntot)=(',npix,',',ntot,')'
			end
		endcase	
	endif

endif  ;reading or not
close,lun
free_lun,lun,/force

wav=-1
flux=-1
;interpolate linearly, if x is present and has a valid size
if n_elements(x) ne ndim then begin

	if n_elements(x) eq 0 then begin
		print,'% READ_SYNTH: x not present -- grid read, fluxes are not interpolated'
	endif else begin
			print,'% READ_SYNTH: x has not ',ndim,' elements -- fluxes are not interpolated'
	endelse

endif else begin

	ulimits=llimits+steps*(np-1)			    ; from x's to p's and t's
	p=0.5d0+(2.d0*x-llimits-ulimits)/(2.d0*(ulimits-llimits)) ; 0-1	normalized
	t=p*(np-1.d0)	            	    ; 0-np(i)	indices)

	if multi eq 0 then begin
		wav=wave[0]+dindgen(npix)*wave[1]
	endif else begin
		wav=wave[0,0]+dindgen(nnpix[0])*wave[1,0]	
		for i=1,multi-1 do begin
			wav=[wav,wave[0,i]+dindgen(nnpix[i])*wave[1,i]]
		endfor
	endelse

	flux=findgen(npix)
	
	;ti=systime(1)
	for l=1l,npix	do begin	;loop over wavelengths
		;lin,grid[*,*,*,*,*,l-1],t,y ;fast, but not general -- need as many *s as dimensions
		;lin,grid,[t,l-1],y      ;slower though general, but breaks down for large arrays
		;lin,transpose(reform((transpose(grid))[l-1,*,*,*,*,*,*,*])),t,y ;also general, but superslow
		lin,reform(grid[l-1,*,*,*,*,*,*,*]),reverse(t),y ;this works, just need to keep grid transposed
		flux(l-1)=y	
	endfor	
	;print,systime(1)-ti

endelse

end
