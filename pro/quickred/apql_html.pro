pro apql_html,str

; This makes the HTML file for the quicklook software

cowhite = '#FFFFFF'
cogreen = '#44CC44'
cored = '#FF2222'
colblue = '#6666FF'
copink = 'pink'
colred = '#FF6666'
angstrom = '&#8491;'             ; HTML Angstrom symbol

filename = str.lastread

;dir = file_dirname(filename)
base = file_basename(filename,'.fits')
dum = strsplit(base,'-',/extract)  ; apRaw-00000118-005
;frameid = dum[1]
;readnum = dum[2]

qldir = str.outdir+'/ql/'
;name = frameid+'-'+readnum
name = str.frameid+'-'+str.readnum



; Make HTML page showing representative spectra and display arrays
;-------------------------------------------------------------------
name = str.frameid+'-'+str.readnum
apgundef,lines
PUSH,lines,'<HTML>'
PUSH,lines,'<HEAD>'
PUSH,lines,'<TITLE>'
PUSH,lines,'APOGEE Exposure and Spectra Display for '+name+' Read='+strtrim(str.readnum,2)
PUSH,lines,'</TITLE>'
PUSH,lines,'</HEAD>'
PUSH,lines,''
PUSH,lines,'<BODY>'
PUSH,lines,''
PUSH,lines,'<center>'
PUSH,lines,''
PUSH,lines,'<h2>Array Display</h2>'
PUSH,lines,'<table border=0>'
PUSH,lines,'<tr><td><center><h4>All Fibers</h4></center></td></tr>'
PUSH,lines,'<tr><td><img src="'+file_basename(str.plotarrays_all_figfile)+'.'+str.figext+'"></td></tr>'
PUSH,lines,'<tr><td><center><h4>Fiber XX-XX (Rows 1500-1600)</h4></center></td></tr>'
PUSH,lines,'<tr><td><img src="'+file_basename(str.plotarrays_sub_figfiles[0])+'.'+str.figext+'"></td></tr>'
PUSH,lines,'<tr><td><center><h4>Fiber XX-XX (Rows 1000-1100)</h4></center></td></tr>'
PUSH,lines,'<tr><td><img src="'+file_basename(str.plotarrays_sub_figfiles[1])+'.'+str.figext+'"></td></tr>'
PUSH,lines,'<tr><td><center><h4>Fiber XX-XX (Rows 500-600)</h4></center></td></tr>'
PUSH,lines,'<tr><td><img src="'+file_basename(str.plotarrays_sub_figfiles[2])+'.'+str.figext+'"></td></tr>'
PUSH,lines,'</table>'
PUSH,lines,''
PUSH,lines,'<p>'
PUSH,lines,'<h2>Representative Spectra</h2>'
PUSH,lines,'<img src="'+file_basename(str.plotspec_figfile)+'.'+str.figext+'">'
PUSH,lines,''
PUSH,lines,'</center>'
PUSH,lines,''
PUSH,lines,'</BODY>'
PUSH,lines,'</HTML>'
PUSH,lines,''
; Write the HTML file
dataview_htmlfile = str.outdir+'apql_dataview_'+name+'.html'
WRITELINE,dataview_htmlfile,lines
str.dataview_htmlfile = dataview_htmlfile



; Make SKY HTML Page
;---------------------
name = str.frameid+'-'+str.readnum
apgundef,lines
PUSH,lines,'<HTML>'
PUSH,lines,'<HEAD>'
PUSH,lines,'<TITLE>'
PUSH,lines,'APOGEE Sky for '+name+' Read='+strtrim(str.readnum,2)
PUSH,lines,'</TITLE>'
PUSH,lines,'</HEAD>'
PUSH,lines,''
PUSH,lines,'<BODY>'
PUSH,lines,''
PUSH,lines,'<center>'
PUSH,lines,'<h2>Sky Flux Levels</h2>'
PUSH,lines,'<table>'

str.sky_status = 1  ; okay for now

; Average line flux
;  don't want it to saturate
if str.skyvar_avglineflux gt 60000L then begin
  str.sky_status = 0
  color = cored
endif else begin
  color = cogreen
endelse
PUSH,lines,'<tr><td>Average Line Flux</td><td bgcolor="'+color+'">'+stringize(str.skyvar_avglineflux,ndec=0)+' counts</td>'+$
           '<td><60,000 counts</td></tr>'
; Average line flux RATE
;  don't want the lines to saturate in 600sec
if str.skyvar_avglineflux_rate gt 100L then begin
  str.sky_status = 0
  color = cored
endif else begin
  color = cogreen
endelse
PUSH,lines,'<tr><td>Average Line Flux Rate</td><td bgcolor="'+color+'">'+stringize(str.skyvar_avglineflux_rate,ndec=1)+$
           ' counts/sec</td><td><100 counts/sec</td></tr>'
; Continuum Flux
if str.skyvar_contflux gt 10000 then begin
  str.sky_status = 0
  color = cored
endif else begin
  color = cogreen
endelse
PUSH,lines,'<tr><td>Continuum Flux</td><td bgcolor="'+color+'">'+stringize(str.skyvar_contflux,ndec=0)+' counts</td>'+$
           '<td><10,000 counts</td></tr>'
; Continuum Flux RATE
if str.skyvar_contflux_rate gt 20 then begin
  str.sky_status = 0
  color = cored
endif else begin
  color = cogreen
endelse
PUSH,lines,'<tr><td>Continuum Flux Rate</td><td bgcolor="'+color+'">'+stringize(str.skyvar_contflux_rate,ndec=1)+$
           ' counts/sec</td><td><20 counts/sec</td></tr>'
PUSH,lines,'</table>'
PUSH,lines,'<p>'
PUSH,lines,''

PUSH,lines,'<h2>Sky Variation Across All Sky Fibers</h2>'
PUSH,lines,'<table>'

; Median Flux Deviation
if str.skyvar_meddev_perc gt 5 then begin
  str.sky_status = 0
  color = cored
endif else begin
  color = cogreen
endelse
PUSH,lines,'<tr><td>Median flux deviation</td><td bgcolor="'+color+'">'+stringize(str.skyvar_meddev_perc,ndec=3)+' %</td>'+$
           '<td><5%</td></tr>'

; Stddev Flux Deviation
if str.skyvar_stddev_perc gt 5 then begin
  str.sky_status = 0
  color = cored
endif else begin
  color = cogreen
endelse
PUSH,lines,'<tr><td>Stddev flux deviation</td><td bgcolor="'+color+'">'+stringize(str.skyvar_stddev_perc,ndec=3)+' %</td>'+$
           '<td><5%</td></tr>'
PUSH,lines,'</table>'
PUSH,lines,'<p>'

PUSH,lines,'<img src="'+file_basename(str.skyvar_figfile)+'.'+str.figext+'">'
PUSH,lines,'<p>'
PUSH,lines,'<h2>Airglow Line Position</h2>'
PUSH,lines,'<img src="'+file_basename(str.skycheck_figfile)+'.'+str.figext+'" width=800>'
PUSH,lines,'<p>'
PUSH,lines,'</center>'
PUSH,lines,''
PUSH,lines,'</BODY>'
PUSH,lines,'</HTML>'
PUSH,lines,''
; Write the HTML file
sky_htmlfile = str.outdir+'apql_sky_'+name+'.html'
WRITELINE,sky_htmlfile,lines
str.sky_htmlfile = sky_htmlfile




; Make MAIN HTML PAGE
;---------------------

apgundef,lines
PUSH,lines,'<HTML>'
PUSH,lines,'<HEAD>'
PUSH,lines,'<TITLE>'
PUSH,lines,'APOGEE QUICKLOOK Diagnostics for apRaw-'+name
PUSH,lines,'</TITLE>'
PUSH,lines,'</HEAD>'
PUSH,lines,''
PUSH,lines,'<BODY>'
PUSH,lines,''
;PUSH,lines,'<h2>APOGEE QUICKLOOK</h2>'
;PUSH,lines,'<hr>'
PUSH,lines,''
; This is a two-element "placement" table to keep the two tables next to each ohter
PUSH,lines,'<table>'
PUSH,lines,'<tr><td colspan=2 align=center><h2>APOGEE QUICKLOOK</h2></td></tr>'
PUSH,lines,'<tr><td align=center><font size=4><b>Current Exposure</b></font></td>'
PUSH,lines,'<td align=center><font size=4><b>Exposure List for Plate '+strtrim(str.plate,2)+'</b></font></td></tr>'
PUSH,lines,'<tr><td align=center valign=top>'

; Filename
; plate name
; current exposure time
; current read
;     link to see figure of the array
; is fits header okay?       link to see fits header
; dither position okay?      link to see dither plots
; airglow lines okay?        link to see airglow line positions
; airglow variability okay?  link to see airglow line variabilty plots
; wavelengths okay?          link to see wavelength solution plot
; overall status:
; show SNR plot
; estimated total exposure time needed to obtain required S/N level
; show values and color the background in red/green

; for S/N plot show S/N at H=12 (or someething) as a function of time

; maybe make the "keyord" a link

; Might need separate HTML page for the popup links

PUSH,lines,'<table border=1>'
; Exposure name
;----------------
;  with link to display arrays
;PUSH,lines,'<tr><td><a href="'+file_basename(str.dataview_htmlfile)+'" target="_blank" height=300 width=400>Exposure</a></td>'
PUSH,lines,'<tr><td><A HREF="javascript:void(0)"onclick="window.open('+"'"+file_basename(str.dataview_htmlfile)+"','dataview','"+$
           "height=800, width=900,scrollbars=yes'"+')"><b>Exposure</b></a></td>'
PUSH,lines,'<td align="center">apRaw-'+str.frameid+'</td></tr>'
; Read number
;-------------
PUSH,lines,'<tr><td><b>Read Number</b></td>'
PUSH,lines,'<td align="center">'+strtrim(long(str.readnum),2)+'</td></tr>'
; Exposure time
;---------------
PUSH,lines,'<tr><td><b>Exp.Time</b></td>'
PUSH,lines,'<td align="center">'+strtrim(str.exptime,2)+' sec</td></tr>'
; Plate
;---------------
PUSH,lines,'<tr><td><b>Plate</b></td>'
PUSH,lines,'<td align="center">'+strtrim(str.plate,2)+'</td></tr>'

; Fits header
;-------------
;PUSH,lines,'<tr><td><a href="'+file_basename(str.fitsheader_htmlfile)+'" target="_blank" height=300 width=400>FITS header</a></td>'
PUSH,lines,'<tr><td><A HREF="javascript:void(0)"onclick="window.open('+"'"+file_basename(str.fitsheader_htmlfile)+"','fitsheader','"+$
           "height=380, width=600,scrollbars=yes'"+')"><b>FITS header</b></a></td>'

if str.fitsheader_status eq 1 then begin
  status = 'OKAY'
  bgcolor = cogreen
endif else begin
  status = 'ERROR'
  bgcolor = cored
endelse
PUSH,lines,'<td bgcolor="'+bgcolor+'" align="center">'+status+'</td></tr>'

; Dither position okay
;----------------------
;PUSH,lines,'<tr><td><a href="'+file_basename(str.dither_htmlfile)+'" target="_blank" height=300 width=400>Dither Position</a></td>'
PUSH,lines,'<tr><td><A HREF="javascript:void(0)"onclick="window.open('+"'"+file_basename(str.dither_htmlfile)+"','dither','"+$
           "height=600, width=700,scrollbars=yes'"+')"><b>Dither Position</b></a></td>'
if str.dither_status eq 1 then begin
  ;status = 'OKAY'
  bgcolor = cogreen
endif else begin
  ;status = 'ERROR'
  bgcolor = cored
endelse
sdither = stringize(str.dither_prevexp_measured,ndec=2)+'/'+$
           stringize(str.dither_prevexp_header,ndec=2)+' pix'
; First exposure
if str.dither_prevexp_measured gt 90. and str.dither_prevexp_header gt 90 then begin
  bgcolor = cogreen
  sdither = 'N.A.'
endif
;PUSH,lines,'<td bgcolor="'+bgcolor+'" align="center">'+status+'</td></tr>'
PUSH,lines,'<td bgcolor="'+bgcolor+'" align="center">'+sdither+'</td></tr>'
; Link to see airglow line positions

; Sky 
;-----
;   -airglow line positions
;   -airglow flux
;   -airglow variability
;   -continuum flux
;PUSH,lines,'<tr><td><a href="'+file_basename(str.sky_htmlfile)+'" target="_blank" height=300 width=400>Sky</a></td>'
PUSH,lines,'<tr><td><A HREF="javascript:void(0)"onclick="window.open('+"'"+file_basename(str.sky_htmlfile)+"','sky','"+$
           "height=800, width=900,scrollbars=yes'"+')"><b>Sky</b></a></td>'
if str.sky_status eq 1 then begin
  status = 'OKAY'
  bgcolor = cogreen
endif else begin
  status = 'ERROR'
  bgcolor = cored
endelse
PUSH,lines,'<td bgcolor="'+bgcolor+'" align="center">'+status+'</td></tr>'

; Wavelengths
;-------------
;   -wavelength range
;   -solution coefficients?
;PUSH,lines,'<tr><td><a href="'+file_basename(str.wavelength_htmlfile)+'" target="_blank" height=300 width=400>Wavelengths</a></td>'
PUSH,lines,'<tr><td><A HREF="javascript:void(0)"onclick="window.open('+"'"+file_basename(str.wavelength_htmlfile)+"','wavelengths','"+$
           "height=800, width=800,scrollbars=yes'"+')"><b>Wavelengths</b></a></td>'
if str.wavelength_status eq 1 then begin
  ;status = 'OKAY'
  bgcolor = cogreen
endif else begin
  ;status = 'ERROR'
  bgcolor = cored
endelse
;PUSH,lines,'<td bgcolor="'+bgcolor+'" align="center">'+status+'</td></tr>'
PUSH,lines,'<td bgcolor="'+bgcolor+'" align="center">'+stringize(median(str.wavelength_diff),ndec=2)+' '+angstrom+'</td></tr>'

; SNR
;------------------
if str.snr_standard ge str.snr_standard_goal then str.continue_status=1
PUSH,lines,'<tr><td><b>S/N at H='+stringize(str.hmag_standard,ndec=1)+'</b></td>'
if str.continue_status eq 1 then begin
  ;status = 'CONTINUE'
  bgcolor = cogreen
endif else begin
  ;status = 'STOP'
  bgcolor = cored
endelse
PUSH,lines,'<td bgcolor="'+bgcolor+'" align="center">'+stringize(str.snr_standard,ndec=1)+'/'+stringize(str.snr_standard^2,ndec=1)+'</td></tr>'

; Estimated exptime
;-------------------
time_per_read = str.exptime/long(str.readnum)
expected_exptime = str.expected_total_readnum * time_per_read
if str.continue_status eq 1 then begin
  ;status = 'CONTINUE'
  bgcolor = cogreen
endif else begin
  ;status = 'STOP'
  bgcolor = cored
endelse
PUSH,lines,'<tr><td><b>Estimated Exp.Time</b></td>'
PUSH,lines,'<td bgcolor="'+bgcolor+'" align="center">'+stringize(expected_exptime,ndec=1)+' sec</td></tr>'
PUSH,lines,'</table>'
PUSH,lines,''
PUSH,lines,'</td><td align=center valign=top>'  ; new cell in placement table
PUSH,lines,''


; Plate Exposures Table
;------------------------------------------------------------
; Check all of the files on disk for the same night, but BEFORE
;  this exposure.  Then get the exposures for the same plate.
; Maybe save the STR structure and look at those.

; Start the HTML table
PUSH,lines,'<table border=1>'
; Nexposure for this plate
;-------------------------
PUSH,lines,'<tr><td align="center"><b>Num</b></td><td align="center"><b>Exposure</b></td>'+$
           '<td align="center"><b>Exp.Time</b></td><td align="center"><b>Nreads</b></td>'+$
           '<td align="center"><b>Dither</b></td><td><b>S/N at H='+stringize(str.hmag_standard,ndec=1)+'</b></td></tr>'

total_exptime = 0.0
total_snr2 = 0.0
numcount = 1

; Check the files
strfiles = file_search(str.outdir+'apRaw-????????-???_str.fits',count=nstrfiles)
if nstrfiles gt 0 then begin

  ; Get unique exposure/frame numbers
  strbase = file_basename(strfiles,'_str.fits')  ; apRaw-00000113-002 format
  dum = strsplitter(strbase,'-',/extract)
  frameids = reform(dum[1,*])
  rnum = long(reform(dum[2,*]))
  ui_frameids = uniq(frameids,sort(frameids))
  uframeids = frameids[ui_frameids]
  nframeids = n_elements(uframeids)
  platearr = strarr(nframeids)   ; plates for the exposures
  laststrfiles = strarr(nframeids)

  ; Loop through the frames
  for i=0,nframeids-1 do begin
    ind = where(frameids eq uframeids[i],nind)  ; reads for this frame
    irnum = rnum[ind]
    maxirnum = max(irnum)
    maxind = (where(irnum eq maxirnum))[0]    ; want the last read
    lastind = ind[maxind]
    laststrfiles[i] = strfiles[lastind]
    ; Get associated FITS file
    fitsfile = str.datadir+'apRaw-'+uframeids[i]+'-'+string(maxirnum,format='(I03)')+'.fits'
    if file_test(fitsfile) eq 1 then begin
      head = headfits(fitsfile)
      plate = sxpar(head,'PLATE',count=nplate)
      if nplate gt 0 then platearr[i] = strtrim(plate,2)
    endif
  end ; frames loop



  ; Other exposures with the same plate
  sameplate_ind = where(platearr eq str.plate and uframeids ne str.frameid,nsameplate_ind)
  if nsameplate_ind gt 0 then begin

    ; Loop through the exposures
    for i=0,nsameplate_ind-1 do begin

      ; Load the STR file
      istr = MRDFITS(laststrfiles[sameplate_ind[i]],1,/silent)

      ; Only use exposures BEFORE this one
      ;if istr.jd lt str.jd and istr.jd gt 0 then begin
      if long(istr.frameid) lt long(str.frameid) then begin

        ; Color-code by dither-position
        bgcolor=colblue
        if istr.dither_prevexp_measured ge 0.0 then bgcolor=colred
        if i eq 0 then bgcolor=cowhite  ; first is always white
        sdither = stringize(istr.dither_prevexp_measured,ndec=1)
        if istr.dither_prevexp_measured gt 0.0 then sdither='+'+sdither
        if istr.dither_prevexp_measured gt 50. then sdither='0.0'   ; first exposure

        ; exposure name, exptime, nreads, dither position, snr
        PUSH,lines,'<tr bgcolor="'+bgcolor+'"><td align="center">'+strtrim(numcount,2)+'</td><td align="center">'+istr.frameid+'</td>'+$
                   '<td align="center">'+stringize(float(istr.exptime),ndec=0)+' sec</td>'+$
                   '<td align="center">'+strtrim(long(istr.readnum),2)+'</td><td align="center">'+sdither+'</td>'+$
                   '<td align="center">'+stringize(istr.snr_standard,ndec=1)+'</td></tr>'
        total_exptime += float(istr.exptime)
        total_snr2 += istr.snr_standard^2
        numcount++

      end ; earlier exposure

    end ; exposure loop

  ; First exposure for this plate
  end else begin
    ;PUSH,lines,'<td align="center">XX</td></tr>'
  endelse

endif

; Now add the current exposure
;------------------------------

; Color-code by dither-position
bgcolor=colblue
if str.dither_prevexp_measured ge 0.0 then bgcolor=colred
if i eq 0 then bgcolor=cowhite  ; first is always white
sdither = stringize(str.dither_prevexp_measured,ndec=1)
if str.dither_prevexp_measured gt 0.0 then sdither='+'+sdither
if str.dither_prevexp_measured gt 50. then sdither='0.0'   ; first exposure

PUSH,lines,'<tr bgcolor="'+bgcolor+'"><td align="center">'+strtrim(numcount,2)+'</td><td align="center">'+str.frameid+'</td>'+$
           '<td align="center">'+stringize(float(str.exptime),ndec=0)+' sec</td>'+$
           '<td align="center">'+strtrim(long(str.readnum),2)+'</td><td align="center">'+sdither+'</td>'+$
           '<td align="center">'+stringize(str.snr_standard,ndec=1)+'</td></tr>'
total_exptime += float(str.exptime)
total_snr2 += str.snr_standard^2

; Total EXPTIME for this plate
;------------------------------
PUSH,lines,'<tr><td colspan=2><b>Total Exptime</b></td>'
PUSH,lines,'<td align="center">'+stringize(total_exptime,ndec=0)+' sec</td><td></td><td></td>'+$
           '<td align="center">'+stringize(sqrt(total_snr2),ndec=1)+'</td></tr>'
PUSH,lines,'</table>'


; End of placement table
;PUSH,lines,'</tr></table>'



; SNR plot
;----------
;PUSH,lines,'<p>'
;PUSH,lines,'<h2>Signal to Noise</h2>'
PUSH,lines,'<tr><td colspan=2 align=center><p>'
PUSH,lines,'<img src="'+file_basename(str.snrmag_figfile)+'.'+str.figext+'"></td></tr>'
PUSH,lines,'</table>'  ; end of placement table
PUSH,lines,''
PUSH,lines,''
PUSH,lines,'</BODY>'
PUSH,lines,'</HTML>'
PUSH,lines,''



;PUSH,lines,'<h2>Dither Position</h2>'
;PUSH,lines,'<img src="'+file_basename(str.dither_figfile)+'.'+str.figext+'">'
;PUSH,lines,'<p>'
;PUSH,lines,'<h2>Representative Spectra</h2>'
;PUSH,lines,'<img src="'+file_basename(str.plotspec_figfile)+'.'+str.figext+'">'
;PUSH,lines,'<p>'
;PUSH,lines,'<h2>Array Display</h2>'
;PUSH,lines,'<table border=0>'
;PUSH,lines,'<tr><td><center><h4>All Rows</h4></center></td></tr>'
;PUSH,lines,'<tr><td><img src="'+file_basename(str.plotarrays_all_figfile)+'.'+str.figext+'"></td></tr>'
;PUSH,lines,'<tr><td><center><h4>Rows 1500-1600</h4></center></td></tr>'
;PUSH,lines,'<tr><td><img src="'+file_basename(str.plotarrays_sub_figfiles[0])+'.'+str.figext+'"></td></tr>'
;PUSH,lines,'<tr><td><center><h4>Rows 1000-1100</h4></center></td></tr>'
;PUSH,lines,'<tr><td><img src="'+file_basename(str.plotarrays_sub_figfiles[1])+'.'+str.figext+'"></td></tr>'
;PUSH,lines,'<tr><td><center><h4>Rows 500-600</h4></center></td></tr>'
;PUSH,lines,'<tr><td><img src="'+file_basename(str.plotarrays_sub_figfiles[2])+'.'+str.figext+'"></td></tr>'
;PUSH,lines,'</table>'
;PUSH,lines,'<p>'
;PUSH,lines,'<h2>SNR versus Magnitude</h2>'
;PUSH,lines,'<img src="'+file_basename(str.snrmag_figfile)+'.'+str.figext+'">'
;PUSH,lines,'<p>'
;PUSH,lines,'<h2>Airglow Line Position</h2>'
;PUSH,lines,'<img src="'+file_basename(str.skycheck_figfile)+'.'+str.figext+'">'
;PUSH,lines,'<p>'
;PUSH,lines,'<h2>Sky Variability</h2>'
;PUSH,lines,'<img src="'+file_basename(str.skyvar_figfile)+'.'+str.figext+'">'
;PUSH,lines,'<p>'
;PUSH,lines,'<h2>Wavelengths</h2>'
;PUSH,lines,'<img src="'+file_basename(str.wavesol_figfile)+'.'+str.figext+'">'
;PUSH,lines,'<p>'
;PUSH,lines,''
;PUSH,lines,''
;PUSH,lines,''
;PUSH,lines,'</BODY>'
;PUSH,lines,'</HTML>'
;PUSH,lines,''

; Write the HTML file
outfile = str.outdir+'apql_'+name+'.html'
WRITELINE,outfile,lines


;stop

end
