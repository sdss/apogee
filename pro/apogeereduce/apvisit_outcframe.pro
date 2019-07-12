pro apvisit_outcframe,frame,plugmap,outfiles,silent=silent,stp=stp

;+
;
; AP1DVISIT_OUTCFRAME
;
; This outputs a sky-corrected frame, apCframe
; This is called from ap1dvisit.pro
;
; INPUTS:
;  frame     The structure that contains the sky-corrected
;              frame with all three chips.
;  plugmap   The Plug Map structure for this plate
;  /silent   Don't print anything to the screen
;  /stp      Stop at the end of the program.
;
; OUTPUTS:
;  The frame is written to the "cframe" directory.
;
; USAGE:
;  IDL>ap1dvisit_outcframe,frame
;
; By D.Nidever  May 2010
; Modifications J. Holtzman 2011+
;-

dirs=getdir()

nframe = n_elements(frame)
nplugmap = n_elements(plugmap)

; Not enough inputs
if nframe eq 0 or nplugmap eq 0 then begin
  print,'Syntax - ap1dvisit_outcframe,frame,plugmap,silent=silent,stp=stp'
  return
endif

; Checking the tags of the input structure
tags = tag_names(frame)
needtags1 = ['CHIPA','CHIPB','CHIPC','SHIFT','TELLSTAR']
for i=0,n_elements(needtags1)-1 do begin
  if (where(tags eq needtags1[i]))[0] eq -1 then begin
    print,'TAG ',needtags1[i],' NOT FOUND in input structure'
    return
  end
end
needtags2 = ['HEADER','FLUX','ERR','MASK','WAVELENGTH','SKY','SKYERR','TELLURIC','TELLURICERR','LSFCOEF','WCOEF']
for i=0,2 do begin
  tags2 = tag_names(frame.(i))
  for j=0,n_elements(needtags2)-1 do begin
    if (where(tags2 eq needtags2[j]))[0] eq -1 then begin
      print,'TAG ',needtags2[j],' NOT FOUND in input structure'
      return
    end
  end
end


; Get the frame information
rawname = frame.chipa.filename
info = apfileinfo(rawname,/silent)
id8 = info.fid8
if id8 eq '' then id8=info.suffix

; apCframe files contain the following:
;    * HDU #0 = Header only
;    * HDU #1 = Flux in units of ADU [FLOAT]
;    * HDU #2 = Error in units of ADU [FLOAT]
;    * HDU #3 = Pixel mask [32-bit INT]
;    * HDU #4 = Wavelength in units of A [DOUBLE]
;    * HDU #5 = Sky flux in units of ADU [FLOAT]
;    * HDU #6 = Sky error in units of ADU [FLOAT]
;    * HDU #7 = Telluric absorption flux in units of [FLOAT]
;    * HDU #8 = Telluric error [FLOAT]
;    * HDU #9 = Wavelength solution coefficients [BINARY FITS TABLE or DOUBLE]
;    * HDU #10 = LSF coefficients [BINARY FITS TABLE or DOUBLE]
;    * HDU #11 = Plug-map structure from plPlugMapM file [BINARY FITS TABLE]
;    * HDU #12 = Plugmap header
;    * HDU #13 = Telluric scaling table
;    * HDU #14 = Shift information table
;
; There is a separate file for each chip - [abc]


chiptag = ['a','b','c']

; Loop through the three chips
For i=0,2 do begin

  ; Update the header:
  ;-------------------
  header = frame.(i).header

  ; Remove the trailing blank lines
  indend = where(stregex(header,'^END',/boolean) eq 1,nindend)
  if indend[0] eq -1 then indend=n_elements(header)-1
  header = header[0:indend[0]]

  ; Add extension explanations
  ;----------------------------
  leadstr = 'AP1DVISIT: '
  sxaddhist,leadstr+systime(0),header
  info = GET_LOGIN_INFO()
  sxaddhist,leadstr+info.user_name+' on '+info.machine_name,header
  sxaddhist,leadstr+'IDL '+!version.release+' '+!version.os+' '+!version.arch,header
  sxaddhist,leadstr+' APOGEE Reduction Pipeline Version: '+getvers(),header
  sxaddhist,leadstr+'Output File:',header
  sxaddhist,leadstr+' HDU0 - Header only',header
  sxaddhist,leadstr+' HDU1 - Flux (ADU)',header
  sxaddhist,leadstr+' HDU2 - Error (ADU)',header
  sxaddhist,leadstr+' HDU3 - flag mask (bitwise OR combined)',header
  sxaddhist,leadstr+'        1 - bad pixels',header
  sxaddhist,leadstr+'        2 - cosmic ray',header
  sxaddhist,leadstr+'        4 - saturated',header
  sxaddhist,leadstr+'        8 - unfixable',header
  sxaddhist,leadstr+' HDU4 - Wavelength (Ang)',header
  sxaddhist,leadstr+' HDU5 - Sky (ADU)',header
  sxaddhist,leadstr+' HDU6 - Sky Error (ADU)',header
  sxaddhist,leadstr+' HDU7 - Telluric',header
  sxaddhist,leadstr+' HDU8 - Telluric Error',header
  sxaddhist,leadstr+' HDU9 - Wavelength coefficients',header
  sxaddhist,leadstr+' HDU10 - LSF coefficients',header
  sxaddhist,leadstr+' HDU11 - Plugmap structure',header
  sxaddhist,leadstr+' HDU12 - Plugmap header',header
  sxaddhist,leadstr+' HDU13 - Telluric structure',header
  sxaddhist,leadstr+' HDU14 - Shift structure',header


  ; Create filename
  ;   apCframe-[abc]-ID8.fits 
  outfile= outfiles[i]
  if not keyword_set(silent) then print,'Writing Cframe to ',outfile

  ; HDU #0 = Header only
  ;----------------------
  FITS_WRITE,outfile,0,header

  ; HDU #1 = Flux in units of ADU [FLOAT]
  ;---------------------------------------
  flux = float(frame.(i).flux)
  MKHDR,header1,flux,/image
  sxaddpar,header1,'CTYPE1','Pixel'
  sxaddpar,header1,'CTYPE2','Fiber'
  sxaddpar,header1,'BUNIT','Flux (ADU)'
  MWRFITS,flux,outfile,header1,/silent

  ; HDU #2 = Flux Error in ADU [FLOAT]
  ;------------------------------------
  bderr=where(frame.(i).err eq baderr(),nbd)
  err = float(frame.(i).err)
  if nbd gt 0 then err[bderr] = baderr()
  MKHDR,header2,err,/image
  sxaddpar,header2,'CTYPE1','Pixel'
  sxaddpar,header2,'CTYPE2','Fiber'
  sxaddpar,header2,'BUNIT','Flux Error (ADU)'
  MWRFITS,errout(err),outfile,header2,/silent

  ; HDU #3 = Pixel mask [32-bit INT]
  ;---------------------------------
  mask = fix(frame.(i).mask)
  MKHDR,header3,mask,/image
  sxaddpar,header3,'CTYPE1','Pixel'
  sxaddpar,header3,'CTYPE2','Fiber'
  sxaddpar,header3,'BUNIT','Flag Mask (bitwise)'
  sxaddhist,'Explanation of BITWISE flag mask (OR combined)',header3
  sxaddhist,' 1 - bad pixels',header3
  sxaddhist,' 2 - cosmic ray',header3
  sxaddhist,' 4 - saturated',header3
  sxaddhist,' 8 - unfixable',header3
  MWRFITS,mask,outfile,header3,/silent

  ; HDU #4 = Wavelength in units of Ang [DOUBLE]
  ;----------------------------------------------
  wave = double(frame.(i).wavelength)
  MKHDR,header4,wave,/image
  sxaddpar,header4,'CTYPE1','Pixel'
  sxaddpar,header4,'CTYPE2','Fiber'
  sxaddpar,header4,'BUNIT','Wavelength (Ang)'
  MWRFITS,wave,outfile,header4,/silent

  ; HDU #5 = Sky flux in units of ADU [FLOAT]
  ;-------------------------------------------
  sky = float(frame.(i).sky)
  MKHDR,header5,sky,/image
  sxaddpar,header5,'CTYPE1','Pixel'
  sxaddpar,header5,'CTYPE2','Fiber'
  sxaddpar,header5,'BUNIT','Sky (ADU)'
  MWRFITS,sky,outfile,header5,/silent

  ; HDU #6 = Sky error in units of ADU [FLOAT]
  ;--------------------------------------------
  skyerr = float(frame.(i).skyerr)
  MKHDR,header6,skyerr,/image
  sxaddpar,header6,'CTYPE1','Pixel'
  sxaddpar,header6,'CTYPE2','Fiber'
  sxaddpar,header6,'BUNIT','Sky Error (ADU)'
  MWRFITS,skyerr,outfile,header6,/silent

  ; HDU #7 = Telluric absorption flux in units of [FLOAT]
  ;--------------------------------------------------------
  telluric = float(frame.(i).telluric)
  MKHDR,header7,telluric,/image
  sxaddpar,header7,'CTYPE1','Pixel'
  sxaddpar,header7,'CTYPE2','Fiber'
  sxaddpar,header7,'BUNIT','Telluric'
  MWRFITS,telluric,outfile,header7,/silent

  ; HDU #8 = Telluric error [FLOAT]
  ;-------------------------------------
  telerr = float(frame.(i).telluricerr)
  MKHDR,header8,telerr,/image
  sxaddpar,header8,'CTYPE1','Pixel'
  sxaddpar,header8,'CTYPE2','Fiber'
  sxaddpar,header8,'BUNIT','Telluric Error'
  MWRFITS,telerr,outfile,header8,/silent

  ; HDU #9 = Wavelength solution coefficients [DOUBLE]
  ;-----------------------------------------------------
  wcoef = double(frame.(i).wcoef)
  MKHDR,header9,wcoef,/image
  sxaddpar,header9,'CTYPE1','Fiber'
  sxaddpar,header9,'CTYPE2','Parameters'
  sxaddpar,header9,'BUNIT','Wavelength Coefficients'
  sxaddhist,'Wavelength Coefficients to be used with PIX2WAVE.PRO:',header9
  sxaddhist,' 1 Global additive pixel offset',header9
  sxaddhist,' 4 Sine Parameters',header9
  sxaddhist,' 7 Polynomial parameters (first is a zero-point offset',header9
  sxaddhist,'                     in addition to the pixel offset)',header9
  MWRFITS,wcoef,outfile,header9,/silent

  ; HDU #10 = LSF coefficients [DOUBLE]
  ;-------------------------------------
  lsfcoef = double(frame.(i).lsfcoef)
  MKHDR,header10,lsfcoef,/image
  sxaddpar,header10,'CTYPE1','Fiber'
  sxaddpar,header10,'CTYPE2','Parameters'
  sxaddpar,header10,'BUNIT','LSF Coefficients'
  sxaddhist,'LSF Coefficients to be used with LSF_GH.PRO:',header10
  sxaddhist,'  binsize  The width of a pixel in X-units.  If this is non-zero',header10
  sxaddhist,'             then a "binned" Gauss-Hermite function is used.  If',header10
  sxaddhist,'             binsize=0 then a "normal, unbinned" Gauss-Hermite',header10
  sxaddhist,'             function is used.',header10
  sxaddhist,'  X0       An additive x-offset.  This is only used to',header10
  sxaddhist,'             evaluate the GH parameters that vary globally',header10
  sxaddhist,'             with X.',header10
  sxaddhist,'  Horder   The highest Hermite order, Horder=0 means',header10
  sxaddhist,'             only a constant term (i.e. only Gaussian).',header10
  sxaddhist,'             There are Horder Hermite coefficients (since we fix H0=1).',header10
  sxaddhist,'  Porder   This array gives the polynomial order for the',header10
  sxaddhist,'             global variation (in X) of each LSF parameter.',header10
  sxaddhist,'             That includes sigma and the Horder Hermite',header10
  sxaddhist,'             coefficients (starting with H1 because we fix H0=1)',header10
  sxaddhist,'             There will be Porder[i]+1 coefficients for',header10
  sxaddhist,'             parameter i.',header10
  sxaddhist,'  GHcoefs  The polynomial coefficients for sigma and the',header10
  sxaddhist,'             Horder Hermite parameters.  There are Porder[i]+1',header10
  sxaddhist,'             coefficients for parameter i.  The Hermite parameters',header10
  sxaddhist,'             start with H1 since we fix H0=1.',header10
  MWRFITS,lsfcoef,outfile,header10,/silent

  ; HDU #11 = Plug-map structure from plPlugMapM file [BINARY FITS TABLE]
  ;----------------------------------------------------------------------
  plugdata = plugmap.fiberdata
  MWRFITS,plugdata,outfile,/silent ; first write the data with no header
                                   ; MWRFITS will add the necessary info

  ; HDU # 12 = Plug-map header values
  ;-------------------------------------
  ; remove FIBERDATA and GUIDEDATA
  pltags = tag_names(plugmap)
  plind = indgen(n_elements(pltags))
  bd = where(pltags eq 'FIBERDATA',nbd)
  tmp = where(pltags eq 'GUIDEDATA',nbd)
  if nbd gt 0 then bd=[bd,tmp]
  REMOVE,bd,plind
  newplug = CREATE_STRUCT(pltags[plind[0]],plugmap.(plind[0]))
  for k=1,n_elements(plind)-1 do newplug = CREATE_STRUCT(newplug,pltags[plind[k]],plugmap.(plind[k]))
  MWRFITS,newplug,outfile,/silent

  ; HDU # 13 = Telluric table
  MWRFITS,frame.tellstar,outfile,/silent

  ; HDU # 14 = Telluric table
  MWRFITS,frame.shift,outfile,/silent

  ;; Now modify the header
  ;header11 = HEADFITS(outfile,exten=11)
  ;pltags = tag_names(plugmap)
  ;; Add plate information to output header
  ;for j=0,n_elements(pltags)-1 do begin
  ;  sz = size(plugmap.(j))
  ;  type = size(plugmap.(j),/type)
  ;  if pltags[j] ne 'HDR' and pltags[j] ne 'FIBERDATA' and sz[0] lt 2 then begin
  ;    ;sxaddpar,header11,pltags[j],plugmap.(j),' PLUGMAPTAG'
  ;    line = 'PLMAPTG: '+strtrim(pltags[j],2)+'='
  ;    if type eq 7 then line=line+"'"+strtrim(plugmap.(j),2)+"'" else $
  ;      line=line+strtrim(plugmap.(j),2)
  ;    sxaddhist,line,header11
  ;  endif
  ;end
  ;; Add original plPlugMap header as history lines
  ;for j=0,n_elements(plugmap.hdr)-1 do begin
  ;  if strtrim(plugmap.hdr[j],2) ne '' then $
  ;    sxaddhist,'PLMAPHD: '+plugmap.hdr[j],header11
  ;  end
  ;; Now modify the header in the file
  ;MODFITS,outfile,0,header11,exten_no=11

  ;stop

end

if keyword_set(stp) then stop

end
