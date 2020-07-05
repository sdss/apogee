from apogee.aspcap import aspcap
from apogee.aspcap import mask

els=aspcap.elems()
for el in els[0]: mask.mkmask(el,globalmask='mask_v02_aspcap.txt')
