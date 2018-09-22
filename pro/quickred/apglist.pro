pro apglist,input,fields0,outname=outname,stp=stp,mjd=mjd,comfile=comfile

;+
;
; APGLIST
;
; This program lists important header parameters for APOGEE data
;
; INPUTS:
;  input     An input list of files to print information for.
;              (e.g. apR-a-*.fits).
;  fields    An array of fields to print.  The default is
;              ['DATE-OBS','NREAD','OBJECT','FILENAME','COMMENT']
;  =outname  Write the output to this file
;  =mjd      Specify the MJD date for a header
;  /stp      Stop at the end of the program.
;
; OUTPUTS:
;  The outputs are printed to the screen
;
; USAGE:
;  IDL>apglist,'apR-a*.apz'
;
; By D.Nidever February 2011
;-

; Not enough inputs
if n_elements(input) eq 0 then begin
  print,'Syntax - apglist,input,fields,outname=outname,stp=stp,mjd=mjd'
  return
endif

; Default parameters
if n_elements(input) eq 0 then input='apR-a-*.fits'

; Load the files
LOADINPUT,input,files,count=nfiles,/exist
if nfiles eq 0 then begin
  print,'No files'
  return
endif

; Where are we writing
if n_elements(outname) gt 0 then begin
  openw,unit,outname[0],/get_lun
  openw,htmlunit,outname[0]+'.html',/get_lun
  print,'Writing output to file ',outname[0]
  printf,htmlunit,'<HTML><BODY> <H2><CENTER> APOGEE LOG'
  daycnv,mjd+2400000.5,yy,mm,dd
  date=string(mm)+'/'+strtrim(string(dd),2)+'/'+strtrim(string(yy),2)
  if keyword_set(mjd) then printf,htmlunit,date,'(MJD=',mjd,')'
  printf,htmlunit,'</CENTER></H2><TABLE BORDER=2><TR>'
endif else begin
  unit = -1    ; write to the screen
endelse

; get comments from a comment file if specified
ncom=0
if n_elements(comfile) gt 0 then begin
  readcol,comfile,format='(a)',com,delimiter=';'
  ncom=n_elements(com)
  icom=intarr(ncom)
  comment=strarr(ncom)
  scomment=strarr(ncom)
  for i=0,ncom-1 do begin
    catch,error_status
    if error_status ne 0 then i+=1
    reads,com[i],num
    if error_status eq 0 then begin
      icom[i]=abs(num)
      c=strtrim(com[i])
      scomment[i]=strmid(c,0,1)
      comment[i]=strmid(c,strpos(c,' ')+1)
    endif else i+=1
    catch,/cancel
  endfor
endif

if n_elements(fields0) eq 0 then begin
  fields = ['DATE-OBS','NREAD','FILENAME','COMMENT']
endif else begin
  fields = fields0
endelse
fields = strupcase(strtrim(fields,2))
nfields = n_elements(fields)

; Base filename lengths
bases = file_basename(files)
lenbases = strlen(bases)+2

; Get the array of format keys and header line
headline = string('NUM',format='(A-4)')+' '+string('FILENAME',format='(A-'+strtrim(max(lenbases)+2,2)+')')
fmtarr = strarr(nfields)
printf,htmlunit,'<TD>FILENAME'
for j=0,nfields-1 do begin

  ; Format statement
  case fields[j] of
    'DATE-OBS': fmt='(A-25)'
    'OBJECT': fmt='(A-10)'
    'NREAD': fmt='(I-6)'
    'NFRAMES': fmt='(I-8)'
    'EXPTIME': fmt='F-10.2)'
    'ORIGNAME': fmt='(A-25)'
    'FILENAME': fmt='(A-20)'   ;'(A-25)'
    'FILTER': fmt='(A-10)'
    'FILTER1': fmt='(A-10)'
    'EXPTYPE': fmt='(A-14)'
    'PLATEID': fmt='(A-10)'
    'LAMPQRTZ': fmt='(I-1)'
    'LAMPTHAR': fmt='(I-1)'
    'LAMPUNE': fmt='(I-1)'
    'COMMENT': fmt='(A-60)'
    else: fmt='(A-10)'
  endcase
  fmtarr[j] = fmt

  ; Get the number of characters
  fmtlen = strlen(fmt)
  ;len = strmid(fmt,2,fmtlen-3)
  len = strmid(fmt,3,fmtlen-4)
  if strpos(fields[j],'LAMP') ge 0 then f=strmid(fields[j],4) else f=fields[j]
  if strpos(f,'DATE-OBS') ge 0 then f='UT' 
  headline += string(f,format='(A-'+strtrim(len,2)+')')
  printf,htmlunit,'<TD>'+f
end
headlen = strlen(headline)
printf,unit,strjoin(replicate('-',headlen))
;printf,unit,'  '+headline
printf,unit,headline
printf,unit,strjoin(replicate('-',headlen))

oldcomment = ''
; Loop through the files
for i=0,nfiles-1 do begin

  file = files[i]
  dum = strsplit(file_basename(file),'.',/extract)
  exten = first_el(dum,/last)
  if exten eq 'apz' then num_exten=1 else num_exten=0

  head = headfits(files[i],exten=num_exten)
  head1 = headfits(files[i],exten=num_exten+1)
  base = file_basename(files[i])

  line = string(i+1,format='(I-4)')+' '+string(base,format='(A-'+strtrim(max(lenbases)+2,2)+')')

  basecomponents=strsplit(base,'-',/extract)
  apid=file_basename(basecomponents[2],'.apz')
  apnum=0L
  reads,apid,apnum
  apnum=apnum mod 1000

  ; add pre-exposure comment if provided
  kcom=-1
  if ncom gt 0 then begin
    kcom=where(icom eq apnum)
    for kk=0,n_elements(kcom)-1 do $
    if kcom[kk] ge 0 then if scomment[kcom[kk]] eq '-' then printf,htmlunit,'<TR><TD COLSPAN=',nfields+1,' BGCOLOR=lightyellow>'+comment[kcom[kk]]
  endif

  printf,htmlunit,'<TR><TD>'+apid

  ; Loop over the fields/keywords to print out
  for j=0,nfields-1 do begin

    catch,error_status
    if error_status ne 0 then j=j+1
    value = sxpar(head,fields[j],count=count)

    ; Get DATE-OBS from first read
    if fields[j] eq 'DATE-OBS' then begin
      value=sxpar(head1,fields[j],count=count)
      comp=strsplit(value,'T',/extract)
      value=comp[1]
    endif

    ; get airmass from SECZ, or ARMASS if SECZ doesn't exist
    if fields[j] eq 'SECZ' then begin
      value=sxpar(head1,'ALT',count=count)
      sz=size(value)
      if count gt 0 and (sz[1] eq 4 or sz[1] eq 5) then $
        value=1./cos((90-value)*3.14159/180.) $
      else begin
        value=sxpar(head1,'ARMASS',count=count)
        if count eq 0 then value=''
      endelse
    endif

    if fields[j] eq 'OBSCMNT' then begin
      value=sxpar(head1,fields[j],count=count)
      if value eq oldcomment then value='' else oldcomment=value
    endif

    if count eq 0 then value=''

    ; Trim .fits from FILENAME
    if fields[j] eq 'FILENAME' then value=file_basename(value,'.fits')

    ; sxpar has problems with the COMMENT keyword
    ; because it interprets it like a HISTORY line
    if fields[j] eq 'COMMENT' then begin
      ind_comment = where(strmid(head,0,9) eq 'COMMENT =',count)
      if count gt 0 then begin
        txt = strmid(head[ind_comment],9)
        tmp_head = [head[0:1],'COMMENT1='+txt,'END'+replicate(' ',77)]  ; make temporary header
        value = sxpar(tmp_head,'COMMENT1')
      endif else begin
        value = ''
      endelse
      for kk=0,n_elements(kcom)-1 do $
        if kcom[kk] ge 0 then begin
          if scomment[kcom[kk]] ne '+' and scomment[kcom[kk]] ne '-' then $
            value+=comment[kcom[kk]]
      endif
    endif
    if strpos(fields[j],'LAMP') ge 0 then value=fix(value)
    svalue=strtrim(string(value),2)
    line += string(svalue,format=fmtarr[j])
    printf,htmlunit,'<TD>'+svalue
  end

  printf,unit,line
  ; add post-exposure comment if provided
  for kk=0,n_elements(kcom)-1 do $
  if kcom[kk] ge 0 then if scomment[kcom[kk]] eq '+' then printf,htmlunit,'<TR><TD COLSPAN=',nfields+1,' BGCOLOR=lightyellow>'+comment[kcom[kk]]

end

printf,unit,strjoin(replicate('-',headlen))

; Close the file
if n_elements(outname) gt 0 then begin
  close,unit
  free_lun,unit
  printf,htmlunit,'</table></body></html>'
  free_lun,htmlunit
endif

if keyword_set(stp) then stop

end
