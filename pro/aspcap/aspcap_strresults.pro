pro aspcap_strresults,file,spec

; This procedure reads the ferre output files and put 
; the information into an structure.

load,file+'.mdl',mdl
load,file+'.wav',wav
load,file+'.frd',frd
load,file+'.err',err
load,file+'.ipf',ipf
load,file+'.spm',spm

spec={wav:double(wav),frd:float(frd),err:double(err),syn:float(mdl),ipf:ipf,spm:spm}

end
