pro aspcap_namedtags,allstar,labs

tagnames=tag_names(allstar)

; Assign values to the _named_ parameter tags
params=aspcap_params(paramtags)
nparam=n_elements(params)
index=intarr(nparam)
for ipar=0L, nparam-1L do index[ipar]= where(strtrim(labs.param_symbol,2) eq params[ipar])

for iparam=0,n_elements(params)-1 do begin
  itag=where(strtrim(tagnames,2) eq strtrim(paramtags[iparam],2),ntag)
  if ntag gt 0 then allstar.(itag)= -9999.99
  ierrtag=where(strtrim(tagnames,2) eq strtrim(paramtags[iparam],2)+'_ERR',ntag)
  if ntag gt 0 then allstar.(itag)= -999.99
  if ntag gt 0 and index[iparam] ge 0 then begin
    allstar.(itag)=allstar.param[index[iparam]]
    gd = where(allstar.param_cov[index[iparam],index[iparam]] gt 0,ngd)
    if ngd gt 0 then  allstar[gd].(ierrtag) = sqrt(allstar[gd].param_cov[index[iparam],index[iparam]])
  endif
endfor

; vmicro, vmacro/vsini different for dwarfs and giants
allstar.vmicro=-9999.99
gd=where(allstar.param[2] gt -999,ngd)
if ngd gt 0 then allstar[gd].vmicro=10.^allstar[gd].param[2]
dw=where(strpos(allstar.aspcap_class,'BA') ge 0 or strpos(allstar.aspcap_class,'GKd') ge 0  or $
         strpos(allstar.aspcap_class,'Fd') ge 0 or strpos(allstar.aspcap_class,'Md') ge 0,ndw, complement=giant)
if ndw gt 0 then begin
  allstar[dw].vmacro=10.^0.
  allstar[dw].vsini=10.^allstar[dw].param[7]
endif
allstar[giant].vmacro=10.^allstar[giant].param[7]

; Assign values to the _named_ element tags
elems=aspcap_elems(elemtags,elemtoh,nelem=nelem)
nelem=n_elements(elems)
eindex=intarr(nelem)
for ipar=0L, nelem-1L do eindex[ipar]= where(strtrim(labs.elem_symbol,2) eq elems[ipar])

ife = where(strtrim(labs.elem_symbol,2) eq 'Fe')
for ielem=0,n_elements(elems)-1 do begin
  itag=where(strtrim(tagnames,2) eq strtrim(elemtags[ielem],2),ntag)
  if ntag gt 0 then allstar.(itag)= -9999.99
  ierrtag=where(strtrim(tagnames,2) eq strtrim(elemtags[ielem],2)+'_ERR',ntag)
  if ntag gt 0 then allstar.(itag)= -999.99
  iflagtag=where(strtrim(tagnames,2) eq strtrim(elemtags[ielem],2)+'_FLAG',ntag)
  if ntag gt 0 then allstar.(itag)= 0L
  if eindex[ielem] ge 0 then begin
    if itag ge 0 then begin
      if ielem eq ife then begin
        ; Fe is special, since we don't want [Fe/Fe]!
        allstar.(itag)=allstar.x_h[eindex[ielem]]
        allstar.(ierrtag) = allstar.x_h_err[eindex[ielem]]
        allstar.(iflagtag) = allstar.elemflag[eindex[ielem]]
      endif else begin
        gd=where(allstar.x_h[eindex[ielem]] gt -9998. and allstar.x_h[ife] gt -9998.,ngd)
        if ngd gt 0 then begin
          allstar[gd].(itag)=allstar[gd].x_h[eindex[ielem]]-allstar[gd].x_h[ife]
          allstar[gd].(ierrtag) = allstar[gd].x_m_err[eindex[ielem]]
        endif
        allstar.(iflagtag) = allstar.elemflag[eindex[ielem]]
      endelse

    endif
  endif
endfor

end

