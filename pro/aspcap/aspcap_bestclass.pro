function aspcap_bestclass,allparam,allspec,alllib,nocoarse=nocoarse,goodclass=goodclass,classes=classes

; routine to choose the best fit from FERRE solutions that
;   are done with different libraries (classes)
; with keyword goodclass, only use specified "allowed" classes

if keyword_set(goodclass) then begin
  ntot=0
  for i=0,n_elements(goodclass)-1 do begin
    gd=where(strtrim(allparam.class,2) eq strtrim(goodclass[i],2),ngd)
    if ntot eq 0 then allgd=gd else allgd=[allgd,gd]
    ntot+=ngd
  endfor
  if ntot eq 0 then stop,'No stars from allowed classes: ',goodclass
endif else allgd=indgen(n_elements(allparam))

if ~keyword_set(classes) then begin
  iclass=uniq(allparam.class,sort(allparam.class))
  classes=allparam[iclass].class
  nclass=n_elements(iclass)
endif else nclass=n_elements(classes)
star=allparam[allgd].apogee_id
istar=uniq(star,sort(star))

nparam=n_elements(allparam[0].fparam)
newparam=struct_rename_tags(allparam,'','',addtag='FPARAM_CLASS',addval=fltarr(nparam,nclass)-9999.)
newparam=struct_rename_tags(newparam,'','',addtag='CHI2_CLASS',addval=fltarr(nclass)-9999.)
; assumes all .lib structures are the same!!!!
lib=struct_rename_tags(alllib[0],'','',addtag='CLASSES',addval=classes)

for i=0,n_elements(istar)-1 do begin
  chi2=fltarr(nclass)-1
  index=intarr(nclass)
  print,star[istar[i]]
  if keyword_set(nocoarse) then $
  j=where(strpos(allparam[allgd].class,'coarse') lt 0 and strtrim(allparam[allgd].apogee_id,2) eq strtrim(star[istar[i]],2),nj) else $
  j=where(strpos(allparam[allgd].class,'coarse') ge 0 and strtrim(allparam[allgd].apogee_id,2) eq strtrim(star[istar[i]],2),nj) 
  if nj gt 0 then begin
    for jj=0,n_elements(j)-1 do begin
      index[jj]=j[jj]
      chi2[jj]=allparam[allgd[index[jj]]].param_chi2
      ; Need to penalize at edge of grid a little bit, to allow for better solutions in adjacent grid with
      ;   slightly larger chi^2 (although not clear why!). But can't penalize too much, because it's
      ;   possible that it's at the edge correctly, but another grid in which the match is terrible won't
      ;   be at the edge of the grid (e.g., a  hot star in the Mcoarse grid!)
      ;if (allparam[allgd[index[jj]]].paramflag[0] and paramflagval('GRIDEDGE_WARN')) gt 0 then chi2[jj]*=1.25
      if (allparam[allgd[index[jj]]].paramflag[0] and paramflagval('GRIDEDGE_BAD')) gt 0 then chi2[jj]*=1.25
      ; penalize M grid above 3550 (l30e)
  ;    if strpos(allparam[allgd[index[jj]]].class,'M') ge 0 and allparam[allgd[index[jj]]].fparam[0] gt 3500.+250./8. then chi2[jj]*=10.
      ; penalize GK grid below 3750/3990 (l30f and l30g)
      if strpos(allparam[allgd[index[jj]]].class,'GK') ge 0 and allparam[allgd[index[jj]]].fparam[0] lt 3990. then chi2[jj]*=10.
    endfor
    gdchi2=where(chi2 gt 0)
    bestchi2=min(chi2[gdchi2],ibest)
    ibest=gdchi2[ibest]

    ; save the parameters from ALL classes for the final output, i.e. in the best entry, under fparam_class tag
    j=where(allparam.apogee_id eq star[istar[i]],nj)
    for jj=0,n_elements(j)-1 do begin
      iclass=where(strtrim(classes,2) eq strtrim(allparam[j[jj]].class),nc)
      if nc eq 0 then stop,'no matching class found in aspcap_bestclass!'
      ; save the parameters from this class for the final output
      newparam[allgd[index[ibest]]].fparam_class[*,iclass] = newparam[j[jj]].fparam
      newparam[allgd[index[ibest]]].chi2_class[iclass] = newparam[j[jj]].param_chi2
    endfor
    if i eq 0 then begin
      finalparam=newparam[allgd[index[ibest]]]
      finalspec=allspec[allgd[index[ibest]]]
    endif else begin
      finalparam=[finalparam,newparam[allgd[index[ibest]]] ]
      finalspec=[finalspec,allspec[allgd[index[ibest]]] ]
    endelse
  endif else stop,'Error in aspcap_bestclass, not finding match for star: ',star[istar[i]]
endfor

finalstr={param: finalparam, spec:finalspec, lib: lib}

return,finalstr
end
