chubak=mrdfits('chubak.fits.gz',1)
restore,'table1.dat'
str={id: 'NONE', ra: 0.d0, dec: 0.d0, rv: 0.d0, rvsig: 0.d0, jd: 0., dt: 0., src: 'NONE'}
out=replicate(str,n_elements(chubak)+n_elements(str_table1))
for i=0,n_elements(chubak)-1 do begin
  out[i].id = chubak[i].name
  out[i].ra = chubak[i].ra
  out[i].dec =  chubak[i].dec
  out[i].rv = chubak[i].rv
  out[i].rvsig = chubak[i].sigmarv
  out[i].jd = chubak[i].jd
  out[i].dt = chubak[i].deltat
  out[i].src =  'Chubak & Marcy'
endfor

n=n_elements(chubak)
for i=0,n_Elements(str_table1)-1 do begin
  out[n+i].id = str_table1[i].simbad_name
  out[n+i].ra = 15*(str_table1[i].rah+str_table1[i].ram/60.+str_table1[i].ras/3600.)
  out[n+i].dec = abs(str_table1[i].decd)+abs(str_table1[i].decm/60.)+abs(str_table1[i].decs/3600.)
  if str_table1[i].decd lt 0 or str_table1[i].decm lt 0 or str_table1[i].decs lt 0 then out[n+i].dec=-1*out[n+i].dec
  out[n+i].rv = str_table1[i].crv/1000.
  out[n+i].rvsig = str_table1[i].rms/1000.
  out[n+i].jd =  str_table1[i].mean_t
  out[n+i].dt =  str_table1[i].del_t
  out[n+i].src =  'Nidever et al'
endfor

end
