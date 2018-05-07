;==========================================================================
; getnums
;       parses a string with file names into array of file numbers
;         aaa-bbb      parses to list of numbers bewteen aaa and bbb (e.g., for consecutive frames)
;         aaa,bbb,ccc  parses to list of aaa bbb ccc (e.g. for non-consecutive frames)
;         aaa          parses to aaa (e.g., for single frames)
;
function getnums,list

  undefine,out
  if strpos(list,',') gt 0 then files=strsplit(list,',',/extract) else files=list
  for i=0,n_elements(files)-1 do begin
    if strpos(files[i],'-') gt 0 then begin
      elements=strsplit(files[i],'-',/extract)
      startval=0L & endval=0L
      reads,elements[0],startval
      reads,elements[1],endval
      n=endval-startval+1
      ims=startval+indgen(n)
    endif else begin
      ims=0L
      reads,files[i],ims
    endelse
    if n_elements(nums) eq 0 then nums=ims else nums=[nums,ims]
  endfor
  return,nums

  if strpos(files,'-') gt 0 then begin
    elements=strsplit(files,'-',/extract)
    startval=0L & endval=0L
    reads,elements[0],startval
    reads,elements[1],endval
    n=endval-startval+1
    return,startval+indgen(n)
  endif else if strpos(files,',') gt 0 then begin
    elements=strsplit(files,',',/extract)
    nums=lonarr(n_elements(elements))
    nnn=0L
    for i=0,n_elements(elements)-1 do begin
      reads,elements[i],nnn
      nums[i]=nnn
    endfor
    return,nums
  endif else begin
    nnn=0L
    reads,files,nnn
    return,nnn
  endelse
  
end
