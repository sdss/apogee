pro fixobject,objects

; fix any negative RAs
bd=where(objects.ra lt 0,nbd)
if nbd gt 0 then objects[bd].ra = 0.
if nbd gt 0 then objects[bd].dec = 0.

tags=['j',$
      'j_err', $
      'h', $
      'h_err', $
      'k', $
      'k_err', $
      'wash_m', $
      'wash_m_err', $
      'wash_t2', $
      'wash_t2_err', $
      'ddo51', $
      'ddo51_err', $
      'irac_3_6', $
      'irac_3_6_err', $
      'irac_4_5', $
      'irac_4_5_err', $
      'irac_5_8', $
      'irac_5_8_err', $
      'irac_8_0', $
      'irac_8_0_err', $
      'wise_4_5', $
      'wise_4_5_err', $
      'targ_4_5', $
      'targ_4_5_err', $
      'pmra', $
      'pmdec', $
      'ak_targ', $
      'ak_wise', $
      'sfd_ebv' $
      ]
tagnames=tag_names(objects)
for i=0,n_elements(tags)-1 do begin
  j=where(strtrim(tagnames,2) eq strtrim(strupcase(tags[i]),2),nj)
  if nj gt 0 then begin
    bd=where(finite(objects.(j)) eq 0,nbd)
    if nbd gt 0 then objects[bd].(j) = -9999.99
  endif
endfor
end
   
