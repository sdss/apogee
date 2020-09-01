pro	inter_synth,synthfile
;+
;	Generates files for testing FERRE interpolations
;
;	IN: synthfile	-string		Name of the input synth file
;
;	The routine outputs two ipf files and a header for a ferre grid.
;
;	The first ipf file (run1.ipf) gives the coordinates for a grid
;	that is one step smaller in each dimension and shifted by half
;	a step (e.g. | . | . |, where | are the original positions and .
;	the new ones). The second ipf file (run2.ipf) recovers subset of 
;	the original grid contained within the 2nd one (e.g. the  central |
;	in the previous 1D example).
;	
;	The output header can be attached to the ferre spectra produced
;	by running ferre with the original grid and run1.ipf in nov=0
;	(interpolation) mode, to create a new ferre synth module. The 
;	output fluxes interpolating in that new grid with run2.ipf can
;	be compared with the results for run2.ipf from the original grid
;	to estimate interpolation errors.
;
;	C. Allende Prieto, El Portezuelo, April 2016
;
;-

if N_params() lt 1 then begin
	print,'% INTER_SYNTH: -- usage -- inter_synth,synthfile'
	return
endif

;reading synth file info
read_synth,synthfile,/g,param=param,npix=npix,ntot=ntot,ndim=ndim,ax=ax,hd=hd

;creating ipf files
np=ax.n
ntot2=product(np-1)
ntot3=product(np-2)
getaa,ndim,np-1,aa2
getaa,ndim,np-2,aa3
param2=fltarr(ndim,ntot2)
param3=fltarr(ndim,ntot3)
for i=0,ndim-1 do param2[i,*]=aa2[i,*]*ax[i].step+ax[i].llimit+ax[i].step/2.0
for i=0,ndim-1 do param3[i,*]=aa3[i,*]*ax[i].step+ax[i].llimit+ax[i].step

openw,lun,'run1.ipf',/get_lun
for i=0l,ntot2-1 do printf,lun,strcompress('run1_entry'+string(i,format='(i09)'),/rem),param2[*,i],format='(a,99(x,f10.4))'
close,lun
openw,lun,'run2.ipf'
for i=0l,ntot3-1 do printf,lun,strcompress('run2_entry'+string(i,format='(i09)'),/rem),param3[*,i],format='(a,99(x,f10.4))'
close,lun

;creating intermediate header for intermediate synthfile
root=strmid(synthfile,2,strpos(synthfile,'.dat')-2)
print,'creating run1.ipf, run2.ipf and '+root+'.hd'

openw,lun,root+'.hd'
for i=0l,n_elements(hd)-1 do begin
	if strpos(hd[i],'N_P') gt -1 then begin
		printf,lun,' N_P =',np-1,format='(a,99(x,i8))' 
		continue
	endif
	if strpos(hd[i],'LLIMITS') gt -1 then begin
		printf,lun,' LLIMITS =',ax.llimit+ax.step/2.0,format='(a,99(1x,f12.7))'
		continue
	endif

	printf,lun,hd[i]	
endfor
close,lun
free_lun,lun

end

