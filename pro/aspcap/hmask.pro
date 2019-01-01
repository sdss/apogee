pro hmask,obs,wid=wid

if not keyword_set(wid) then wid=10

for i=11,21 do begin
 wmask=rydberg(4,i)
 bad=where(obs.wave gt wmask-wid and obs.wave lt wmask+wid,nbad)
print,'masking ', wmask-wid,wmask+wid,nbad
 if nbad gt 0 then begin
   obs.flux[bad,*]=0.
   obs.invar[bad,*]=0.
 endif
endfor

end
