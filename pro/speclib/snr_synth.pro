pro snr_synth,synthfile,snr,nruns=nruns,edges=edges,noise=noise,renorm=renorm,snr_sys=snr_sys,skip=skip
;+
;	Creates input files for FERRE using the spectra of given
;	synth grid, adding noise.
;
;	IN: synthfile    - string        Name of the grid
;	    snr          - float	 S/N for the output spectra
;					 (can be an array)
;
;	OUT: FERRE input files (ipf, frd, err)
;
;	KEYWORDS: nruns  - float	 number of realizations per
;					 entry in the grid (default=1)
;
;		  edges - 		 when on, the grid edges are 
;					 included (default is off)
;
;		  noise -		 when set to >0 an alternative
;					noise model is used. Possible values
;
;			0 (default)	 Gaussian: IDL's randomn intrinsic
;			1 		 Gaussian: IDL's randomu on vars. x and y
;					 and b=sqrt(-2.d0*alog(x))*cos(2.d0*!pi*y)
;			2		 white noise (-1/snr <b< 1/snr)
;
;
;		  renorm		renormalize each spectrum by calling 'continuum'
;						with this input array (if renorm=[0,1,anything,anything]
;						then this is simply a division by the mean.
;
;		  snr_sys		when set, this involves an additional systematic
;					(Gaussian) error introduced to the spectrum as a whole
;
;			skip		when > 1, then only 1 out of 'skip' spectra are printed
;						to the output files
;
;	C. Allende Prieto, MSSL-UCL, September 2008
;				"	 , IAC, April 2010 -- enhanced the renorm option
;				"    , IAC, February 2011 -- added skip keyword
;-

if N_params() lt 2 then begin
	print,'% SNR_SYNTH: use -- snr_synth,synthfile,snr[,nruns=nruns,edges=edges,noise=noise,renorm=renorm,snr_sys=snr_sys]'
	return
endif


;get grid info
;synthkeys,synthfile,ndim,np,wave,npix,loggw,resolution,continuum,llimit,step,label
read_synth,synthfile,ndim=ndim,npix=npix,logw=loggw,resolution=resolution,$
	continuum=continuum,ax=ax,/grid
np=ax.n
llimit=ax.llimit
step=ax.step
label=ax.label


if not keyword_set(skip) then skip=1

if keyword_set(renorm) then begin
if n_elements(renorm) ne 4 then begin
	print,'% SNR_SYNTH: Error -- renorm must have 4 elements! (see continuum.pro)'
	return
endif
endif

;get loop indices
getaa,ndim,np,aa,ntot=ntot

;set nruns, edges, and noise
if not keyword_set(nruns) then nruns=1
if not keyword_set(noise) then noise=0
if not keyword_set(edges) then edges=0

flux=dblarr(npix)
root=strmid(synthfile,2,strpos(synthfile,'.')-2)
openw,1,strcompress(root+'.ipf',/rem)
openw,2,strcompress(root+'.frd',/rem)
openw,3,strcompress(root+'.err',/rem)

for j=0,n_elements(snr)-1 do begin

	print,'-> Working on S/N=',snr[j]

	openr,10,synthfile

	;skip header
	line='hola'
	multi=0
	while not eof(10) and strmid(line,1) ne '/' do begin
		readf,10,line
		print,line
		if (strpos(line,'MULTI') gt -1) then multi=fix(strmid(line,8,strlen(line)-8))
	endwhile
	;multigrid
	for i=1,multi do begin
		  readf,10,line
		  while not eof(10) and strmid(line,1) ne '/' do begin
		  	readf,10,line
		  	print,line
		  endwhile
	endfor

	single=0.0
	for i=1l,ntot do begin
	
	  ; retain only 1/skip 
	  if abs(float(i/skip)-(float(i)/skip)) gt 1e-3 then begin 
	  
	    readf,10,single
	    
	  endif else begin
	
		obj=strcompress('snr_synth_'+string(snr[j])+'_'+string(i),/rem)
		par=dblarr(ndim)
		borderline=0
		for k=0,ndim-1 do begin
			par[k]=aa[k,i-1]*step[k]+llimit[k]
			if aa[k,i-1] eq 0 or aa[k,i-1] eq np[k]-1 then borderline=1
		endfor

		print,'lote ',j,' --> ',i,' of ',ntot
			
		readf,10,flux
		
		;print,'snr=',snr[j],'  spec#=',i,'  coor=',aa[*,i]
		;plot,flux

		;skip edges
		if not edges and borderline then begin
		endif else begin
		for k=0,nruns-1 do begin	
			
			case noise of
			0: begin
				b=randomn(seed,npix,/double)
			end
			1: begin
				x=randomu(seed,npix,/double)
				y=randomu(seed,npix,/double)
				b=sqrt(-2.d0*alog(x))*cos(2.d0*!pi*y)
			end
			2: begin
				b=randomn(seed,npix,/double)*2.d0-1.0d0
			end
			else: begin
				print,'%snr_synth: noise must be 0,1 or 2'
				return
			end
			endcase

			flux2=flux*(1.d0+ b /double(snr[j]))
			error=flux2/double(snr[j])
			;oplot,flux2,col=140
			;stop

			if keyword_set(snr_sys) then begin
				b2=randomn(seed,1,/double)
				flux2=flux2*(1.d0+b2[0]/snr_sys)
				error=sqrt(error^2+1.d0/snr_sys^2)
			endif

			
			if keyword_set(renorm) then begin ;flux2=flux2/mean(flux2)
				continuum,renorm[0],renorm[1],renorm[2],renorm[3],flux2,conti,$
				  ivar=1.d0/error^2
				flux2=flux2/conti
				error=error/conti
			endif
			
			if nruns eq 1 then obj2=obj else $
			obj2=strcompress(obj+'_'+string(k+1,format='(i)'),/rem)
	
			print,'obj2 is',obj2

			printf,1,obj2,par,format='(1x,a35,'+string(ndim)+'(f11.2,1x))'
			printf,2,flux2,$
			  format='(1x,'+string(npix)+'(e15.8,1x))'
			printf,3,error,$
			  format='(1x,'+string(npix)+'(e15.8,1x))'
			  
		endfor
		endelse
	  endelse	; skip 
	endfor
	close,10
endfor
close,1
close,2
close,3

end
