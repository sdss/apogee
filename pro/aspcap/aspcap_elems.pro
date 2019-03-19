function aspcap_elems,tagnames,elemtoh,elemfitnames,nelem=nelem

if n_elements(nelem) eq 0 then nelem=0

;tagnames=strupcase(['Al_H','Ca_H','C_H','Fe_H','K_H','Mg_H','Mn_H','Na_H','Ni_H','N_H','O_H','Si_H','S_H','Ti_H','V_H'])
;elemtoh=[1,0,0,1,1,0,1,1,1,0,0,0,0,0,1]
;return,['Al','Ca','C','Fe','K','Mg','Mn','Na','Ni','N','O','Si','S','Ti','V']

; add Co, Cr, Cu, Ce, TiII?, CI?

tagnames=strupcase(['C_Fe','CI_Fe','N_Fe','O_Fe','Na_Fe','Mg_Fe','Al_Fe','Si_Fe','P_Fe','S_Fe','K_Fe','Ca_Fe','Ti_Fe','TiII_Fe','V_Fe','Cr_Fe', 'Mn_Fe','Fe_H','Co_Fe', 'Ni_Fe','Cu_Fe', 'Ge_Fe','Rb_Fe', 'Ce_Fe', 'Nd_Fe', 'Yb_Fe'])
;tagnames=strupcase(['C_Fe','CI_Fe','N_Fe','O_Fe','Na_Fe','Mg_Fe','Al_Fe','Si_Fe','P_Fe','S_Fe','K_Fe','Ca_Fe','Ti_Fe','TiII_Fe','V_Fe','Cr_Fe', 'Mn_Fe','Fe_H','Co_Fe', 'Ni_Fe','Cu_Fe', 'Ge_Fe','Ce_Fe', 'Rb_Fe', 'Y_Fe', 'Nd_Fe'])
;'Al1_Fe','Al2_Fe','Al3_Fe','Al4_Fe','Al5_Fe',$
;'Ti1_Fe','Ti2_Fe','Ti3_Fe','Ti4_Fe','Ti5_Fe','Ti6_Fe','Ti7_Fe','Ti8_Fe','Ti9_Fe'])
elemtoh=[0,0,0,0,1,0,1,0,1,0,1,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,$
1,1,1,1,1]
;0,0,0,0,0,0,0,0,0]
elemfitnames=['[C/M]','[CI/M]','[N/M]','[O/M]','[Na/H]','[Mg/M]','[Al/H]','[Si/M]','[P/H]','[S/M]','[K/H]','[Ca/M]','[Ti/M]','[TiII/M]','[V/H]','[Cr/H]','[Mn/H]','[Fe/H]','[Co/H]','[Ni/H]','[Cu/H]','[Ge/H]','[Rb/H]','[Ce/H]', '[Nd/H]', '[Yb/H]']
;elemfitnames=['[C/M]','[CI/M]','[N/M]','[O/M]','[Na/H]','[Mg/M]','[Al/H]','[Si/M]','[P/H]','[S/M]','[K/H]','[Ca/M]','[Ti/M]','[TiII/M]','[V/H]','[Cr/H]','[Mn/H]','[Fe/H]','[Co/H]','[Ni/H]','[Cu/H]','[Ge/H]','[Ce/H]','[Rb/H]', '[Y/H]', '[Nd/H]']
return,['C','CI','N','O','Na','Mg','Al','Si','P','S','K','Ca','Ti','TiII','V','Cr','Mn','Fe','Co','Ni','Cu','Ge','Rb','Ce','Nd','Yb']
;return,['C','CI','N','O','Na','Mg','Al','Si','P','S','K','Ca','Ti','TiII','V','Cr','Mn','Fe','Co','Ni','Cu','Ge','Ce','Rb','Y','Nd']
;'Al_15960', 'Al_15972', 'Al_16723', 'Al_16754', 'Al_16767',$
;'Ti_15190','Ti_15319','Ti_15339','Ti_15431','Ti_15547','Ti_15607','Ti_15720','Ti_16334','Ti_16639']



if nelem eq 16 then begin
  tagnames=strupcase(['C_H','N_H','O_H','Na_H','Mg_H','Al_H','Si_H','S_H','K_H','Ca_H','Ti_H','V_H','Mn_H','Fe_H','Ni_H','Nd_H'])
  elemtoh=[0,0,0,1,0,1,0,0,1,0,0,1,1,1,1,1]
  elemfitnames=['[C/M]','[N/M]','[O/M]','[Na/H]','[Mg/M]','[Al/H]','[Si/M]','[S/M]','[K/H]','[Ca/M]','[Ti/M]','[V/H]','[Mn/H]','[Fe/H]','[Ni/H]','[Nd/H]']
  return,['C','N','O','Na','Mg','Al','Si','S','K','Ca','Ti','V','Mn','Fe','Ni','Nd']
endif else begin
  tagnames=strupcase(['C_H','N_H','O_H','Na_H','Mg_H','Al_H','Si_H','S_H','K_H','Ca_H','Ti_H','V_H','Mn_H','Fe_H','Ni_H'])
  elemtoh=[0,0,0,1,0,1,0,0,1,0,0,1,1,1,1]
  elemfitnames=['[C/M]','[N/M]','[O/M]','[Na/H]','[Mg/M]','[Al/H]','[Si/M]','[S/M]','[K/H]','[Ca/M]','[Ti/M]','[V/H]','[Mn/H]','[Fe/H]','[Ni/H]']
  return,['C','N','O','Na','Mg','Al','Si','S','K','Ca','Ti','V','Mn','Fe','Ni']
endelse
end

