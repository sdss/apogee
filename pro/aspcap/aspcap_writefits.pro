pro aspcap_writefits,str,file,aspcapstar=aspcapstar,aspcap_vers=aspcap_vers,apred_vers=apred_vers,mjddir=mjddir

; routine for writing output ASPCAP FITS structures given IDL ASPCAP structure
; with aspcapstar keyword, also writes individual aspcapStar files for each object

; Jon Holtzman, 2012-13

if keyword_set(mjddir) then $ 
  outfile=file_dirname(file)+'/'+mjddir+'/'+file_basename(file) else outfile=file

; make intro header with software versions
mkhdr,hdr,0
sxaddpar,hdr,'idlver',!version.release+' '+!version.os+' '+!version.arch,'IDL Version'
;sxaddpar,hdr,'idlutver',idlutils_version(),'Software version of IDLUTILS'
;sxaddpar,hdr,'idlwrap',idlwrap_version(),'Software version of IDLWRAP'
;sxaddpar,hdr,'ferre',ferre_version(),'Software version of FERRE'
;sxaddpar,hdr,'speclib',speclib_version(),'Software version of SPECLIB'
sxaddpar,hdr,'apogee',getenv('APOGEE_VER'),'Software version of APOGEE'
leadstr='ASPCAP'
sxaddhist,leadstr+'The data are in separate extensions:',hdr
sxaddhist,leadstr+' HDU0 = Header only',hdr
sxaddhist,leadstr+' HDU1 - Table with stellar data and ASPCAP parameters',hdr
sxaddhist,leadstr+' HDU2 - Table with ASPCAP spectra, etc.',hdr
sxaddhist,leadstr+' HDU3 - Table with parameter names, wavelength array',hdr
mwrfits,0,outfile,hdr,/create

; HDU1 has parameter table with entries for each object
file_mkdir,file_dirname(outfile)
mwrfits,str.param,outfile

; HDU2 has table with spectra, errors, and best fits
mwrfits,str.spec,outfile

; HDU3 has library information common to all objects
mwrfits,str.lib,outfile

if keyword_set(mjddir) then begin
  file_delete,file,/allow
  file_link,outfile,file
endif

; write individual aspcapStar files if requested
if keyword_set(aspcapstar) then begin
  w0=4.179d0
  npix=8575
  dw=6.d-6
  pixmin=(str.lib.wavemin-w0)/dw
  pixmax=(str.lib.wavemax-w0)/dw
  ipixmin=nint(pixmin)
  ipixmax=nint(pixmax)
  if max([abs(pixmin-ipixmin),abs(pixmax-ipixmax)]) ge 0.001 then begin
    print,pixmin,pixmax
    print,str.lib.wavemin,str.lib.wavemax
    print,'Library wavelengths do not seem to match?'
    stop
  endif
  ;idlutils_vers=idlutils_version()
  ;idlwrap_vers=idlwrap_version() 
  ;ferre_vers=ferre_version()
  ;speclib_vers=speclib_version()
  apogee_vers=getenv('APOGEE_VER')
  for i=0,n_elements(str.param)-1 do begin
    flux=fltarr(npix)
    err=fltarr(npix)
    spec_bestfit=fltarr(npix)
    ipix=0
    for j=0,2 do begin
      n=ipixmax[j]-ipixmin[j]+1
      flux[ipixmin[j]:ipixmax[j]]=str.spec[i].spec[ipix:ipix+n-1]
      err[ipixmin[j]:ipixmax[j]]=str.spec[i].err[ipix:ipix+n-1]
      spec_bestfit[ipixmin[j]:ipixmax[j]]=str.spec[i].spec_bestfit[ipix:ipix+n-1]
      ipix+=n
    endfor
    outdir=file_dirname(file)
    outfile=outdir+'/aspcapStar-'
    if keyword_set(apred_vers) then outfile=outfile+apred_vers+'-'
    if keyword_set(aspcap_vers) then outfile=outfile+aspcap_vers+'-'
    outfile=outfile+strtrim(str.param[i].apogee_id,2)+'.fits'
    if keyword_set(mjddir) then $ 
      tmpfile=file_dirname(outfile)+'/'+mjddir+'/'+file_basename(outfile) else tmpfile=outfile
    mkhdr,hdr,0
    sxaddpar,hdr,'idlver',!version.release+' '+!version.os+' '+!version.arch,'IDL Version'
    ;sxaddpar,hdr,'idlutver',idlutils_vers,'Software version of IDLUTILS'
    ;sxaddpar,hdr,'idlwrap',idlwrap_vers,'Software version of IDLWRAP'
    ;sxaddpar,hdr,'ferre',ferre_vers,'Software version of FERRE'
    ;sxaddpar,hdr,'speclib',speclib_vers,'Software version of SPECLIB'
    sxaddpar,hdr,'speclib',apogee_vers,'Software version of SPECLIB'
    sxaddpar,hdr,'OBJ',str.param[i].apogee_id,' Object ID'
    sxaddpar,hdr,'RA',str.param[i].ra,' RA (J2000)'
    sxaddpar,hdr,'DEC',str.param[i].dec,' DEC (J2000)'
    sxaddpar,hdr,'GLON',str.param[i].glon,' Galactic longitude'
    sxaddpar,hdr,'GLAT',str.param[i].glat,' Galactic latitude'
    sxaddpar,hdr,'J',str.param[i].j,' 2MASS J magnitude'
    sxaddpar,hdr,'H',str.param[i].h,' 2MASS H magnitude'
    sxaddpar,hdr,'K',str.param[i].k,' 2MASS Ks magnitude'
    sxaddpar,hdr,'LOCID',str.param[i].location_id,' LOCATION ID'
    sxaddpar,hdr,'FIELD',str.param[i].field,' Field'
    sxaddpar,hdr,'TARG1',str.param[i].apogee_target1,' APOGEE_TARGET1 target flag'
    sxaddpar,hdr,'TARG2',str.param[i].apogee_target2,' APOGEE_TARGET2 target flag'
    sxaddpar,hdr,'AKTARG',str.param[i].ak_targ,' Targeting extinction'
    sxaddpar,hdr,'AKMETHOD',str.param[i].ak_targ_method,' extinction method used for targeting'
    sxaddpar,hdr,'AKWISE',str.param[i].ak_wise,' WISE extinction'
    sxaddpar,hdr,'SFD_EBV',str.param[i].sfd_ebv,' SFD E(B-V)'
    sxaddpar,hdr,'SNR',str.param[i].snr,' signal-to-noise'
    sxaddpar,hdr,'VRAD',str.param[i].vrad,' radial velocity of this spectrum'
    sxaddpar,hdr,'VHELIO',str.param[i].vhelio,' heliocentric velocity (km/s)'
    sxaddpar,hdr,'VSCATTER',str.param[i].vscatter,' STDEV of VHELIO'
    sxaddpar,hdr,'VERR',str.param[i].verr,' weighted error of VHELIO'
    ;sxaddpar,hdr,'VLSR',str.param[i].vlsr,' mean LSR velocity'
    ;sxaddpar,hdr,'VGSR',str.param[i].vgsr,' mean GSR velocity'
    sxaddpar,hdr,'NVISITS',str.param[i].nvisits,' number of visits'

    ;mwrfits,0,outfile,hdr
    file_mkdir,file_dirname(tmpfile)
    FITS_WRITE,tmpfile,0,hdr
    mkhdr,hdr1,flux,/image
    sxaddpar,hdr1,'CRVAL1',w0
    sxaddpar,hdr1,'CDELT1',dw
    sxaddpar,hdr1,'CRPIX1',1
    sxaddpar,hdr1,'CTYPE1','LOG-LINEAR'
    sxaddpar,hdr1,'DC-FLAG',1
    mwrfits,flux,tmpfile,hdr1
    mwrfits,err,tmpfile,hdr1
    mwrfits,spec_bestfit,tmpfile,hdr1
    mwrfits,str.param[i],tmpfile

    if keyword_set(mjddir) then begin
      file_delete,outfile,/allow
      file_link,tmpfile,outfile
    endif
  endfor
endif

end

