; pro rdlibhead,libfile,libstr0,libstr

; Read the header of a synthetic library file into  structures
;  


; rdsinglehead reads a single header block
pro rdsinglehead,lun,libstr

; read the header until we find a '/' line
hola='hola'
multi=0
nlabel=0
while not eof(lun) and strmid(hola,1,1) ne '/' do begin
  readf,lun,hola
  ; see if we have a keyword = value pair by looking for multiple words
  words=strsplit(hola,/extract)
  ; card name is first word
  card=words[0]
  nwords=n_elements(words)
  if strpos(card,'=') ge 0 then begin
    card=strmid(card,0,strpos(card,'='))
    nwords+=1
  endif
;  if strpos(card,'=') ge 0 then strput,card,' ',strpos(card,'=')
  if nwords gt 2 then begin
    ; value are words after = sign
    wordpos=strsplit(hola,'=')
    value=strtrim(strmid(hola,wordpos[1]),2)
    if strpos(value,"'") ge 0 then begin
      ; if there's a quote in the value, then it's a string
      quote=strpos(value,"'")
      while quote ge 0 do begin
        strput,value,' ',quote
        quote=strpos(value,"'")
      endwhile 
      ; Special handling for LABEL(i) cards to put them into an array, for
      ;   loading later, else load into structure
      if strpos(card,'LABEL') ge 0 then begin
        paren=strpos(card,'(')
        reads,strmid(card,paren+1),ilabel
        label[ilabel-1]=strtrim(value,2)
        nlabel+=1
      endif else $
      add_tag,libstr,card,value,libstr
    endif else begin
      ; we have numeric value(s). Find how many
      sdata=strsplit(value,/extract)
      ndata=n_elements(sdata)
      ; if there's a decimal point call it a float, else integer
      if strpos(value,'.') ge 0 then data=dblarr(ndata) else data=intarr(ndata)
      tmp=0.d0
      for i=0,ndata-1 do begin
        reads,sdata[i],tmp
        data[i]=tmp
      endfor
      ; load into structure
     ; if multiple cards exist, take last one (needed for CONTINUUM card in apStar libraries!)
     if tag_exist(libstr,card) then begin
        remove_tags,libstr,card,newstr
        libstr=newstr
     endif
      add_tag,libstr,card,data,libstr
      ; Keep track of two special cards
      if card eq 'N_OF_DIM' then label=strarr(data)
      if card eq 'MULTI' then multi=data
    endelse
  endif

endwhile
if nlabel gt 0 then add_tag,libstr,'LABEL',label,libstr

end

; main routine
pro rdlibhead,libfile,libstr0,libstr

print,'reading ',libfile
; open the file for reading
if ~file_test(libfile+'.hdr') then stop,'File '+libfile+'.hdr does not exist!'
openr,lun,libfile+'.hdr',/get_lun

; start the output structure with the filename
libstr0={file: libfile}

; read the first header block
rdsinglehead,lun,libstr0
if tag_exist(libstr0,'multi') then multi=libstr0.multi else multi=0

; if multi>0, needo read multi more headers
if multi gt 0 then begin
 for imulti=0,multi[0]-1 do begin
   tmpstr={file: libfile}
   rdsinglehead,lun,tmpstr
   if imulti eq 0 then libstr=tmpstr else libstr=[libstr,tmpstr]
 endfor
endif

if not tag_exist(libstr0,'N_OF_DIM') then add_tag,libstr0,'N_OF_DIM',libstr[0].n_of_dim,libstr0
if not tag_exist(libstr0,'LLIMITS') then add_tag,libstr0,'LLIMITS',libstr[0].llimits,libstr0
if not tag_exist(libstr0,'STEPS') then add_tag,libstr0,'STEPS',libstr[0].steps,libstr0
if not tag_exist(libstr0,'N_P') then add_tag,libstr0,'N_P',libstr[0].n_p,libstr0
if not tag_exist(libstr0,'LABEL') then add_tag,libstr0,'LABEL',libstr[0].label,libstr0

free_lun, lun

end



