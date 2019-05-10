; rerun aspcap_correct on an allStsar file
function recorrect,file,caldir=caldir,nowrite=norwrite,suffix=suffix,rgblinear=rgblinear

if ~keyword_set(caldir) then caldir='cal/'
a=mrdfits(file+'.fits',1)
b=mrdfits(file+'.fits',2)
c=mrdfits(file+'.fits',3)
aspcap_correct,a,c.elem_symbol,caldir,rgblinear=rgblinear
for i=0,n_elements(a) do aspcapflags[i]=aspcapflag(aspcap[i].aspcapflag,0)

aspcap_namedtags,a,c

if ~keyword_set(norwrite) then begin
  if ~keyword_set(suffix) then suffix='cal'
  mwrfits,a,file+suffix+'.fits',/create
  mwrfits,b,file+suffix+'.fits'
  mwrfits,c,file+suffix+'.fits'
endif

return,a
end
