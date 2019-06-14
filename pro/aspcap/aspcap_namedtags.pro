pro aspcap_namedtags,allstar,labs

tagnames=tag_names(allstar)

; Assign values to the _named_ parameter tags
params=aspcap_params(paramtags)
nparam=n_elements(params)
index=intarr(nparam)
for ipar=0L, nparam-1L do index[ipar]= where(strtrim(labs.param_symbol,2) eq params[ipar])

; flag faint stars and don't populate named tags
faint=where(allstar.h gt 14.6,nfaint,complement=gd,ncomplement=ngd)
for iparam=0,n_elements(params)-1 do begin
  itag=where(strtrim(tagnames,2) eq strtrim(paramtags[iparam],2),ntag)
  if ntag gt 0 then allstar.(itag)= -9999.99
  ierrtag=where(strtrim(tagnames,2) eq strtrim(paramtags[iparam],2)+'_ERR',ntag)
  if ntag gt 0 then allstar.(ierrtag)= -999.99
  if ntag gt 0 and index[iparam] ge 0 then begin
    if nfaint gt 0 then $
      for ifaint=0,nfaint-1 do allstar[faint[ifaint]].paramflag[index[iparam]] = allstar[faint[ifaint]].paramflag[index[iparam]] or paramflagval('FAINT_WARN')
    if ngd gt 0 then begin
      allstar[gd].(itag)=allstar[gd].param[index[iparam]]
      gd1 = where(allstar[gd].param_cov[index[iparam],index[iparam]] gt 0,ngd1)
      if ngd1 gt 0 then  allstar[gd[gd1]].(ierrtag) = sqrt(allstar[gd[gd1]].param_cov[index[iparam],index[iparam]])
    endif
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
  if ntag gt 0 then allstar.(ierrtag)= -999.99
  iflagtag=where(strtrim(tagnames,2) eq strtrim(elemtags[ielem],2)+'_FLAG',ntag)
  if ntag gt 0 then allstar.(iflagtag)= 0L
  if eindex[ielem] ge 0 then begin
    if itag ge 0 then begin
     ; special handling for DR16 Ce given RV issue for line near edge of chip
     if elems[ielem] eq 'Ce' then begin
        bd=where(allstar.vhelio_avg gt 120,nbd)
        if nbd gt 0 then begin
          allstar[bd].x_h[ielem] = -9999.99
          allstar[bd].x_h_err[ielem] = -999.99
          allstar[bd].x_m[ielem] = -9999.99
          allstar[bd].x_m_err[ielem] = -999.99
          allstar[bd].elemflag[eindex[ielem]] = allstar[bd].elemflag[eindex[ielem]] or paramflagval('RV_WARN')
        endif
     endif
     ; set FAINT_WARN
     if nfaint gt 0 then $
       for ifaint=0,nfaint-1 do allstar[faint[ifaint]].elemflag[eindex[ielem]] = allstar[faint[ifaint]].elemflag[eindex[ielem]] or paramflagval('FAINT_WARN')
     allstar.(iflagtag) = allstar.elemflag[eindex[ielem]]
     ; decision for DR16: only populate named tags with zero warnings for this element AND for Fe
     named=where(allstar.elemflag[eindex[ielem]] eq 0 and allstar.elemflag[eindex[ife]] eq 0,n)
     if n gt 0 then begin
      if ielem eq ife then begin
        ; Fe is special, since we don't want [Fe/Fe]!
        allstar[named].(itag)=allstar[named].x_h[eindex[ielem]]
        allstar[named].(ierrtag) = allstar[named].x_h_err[eindex[ielem]]
      endif else begin
        gd=where(allstar[named].x_h[eindex[ielem]] gt -9998. and allstar[named].x_h[ife] gt -9998.,ngd)
        if ngd gt 0 then begin
          allstar[named[gd]].(itag)=allstar[named[gd]].x_h[eindex[ielem]]-allstar[named[gd]].x_h[ife]
          allstar[named[gd]].(ierrtag) = allstar[named[gd]].x_m_err[eindex[ielem]]
        endif
      endelse
     endif
    endif
  endif
endfor

end

