import numpy as np
import glob
from holtz.apogee import aspcap
import pdb

w=np.loadtxt('wave.dat')

for el in aspcap.elems()[0] :
  m=np.loadtxt(el+'.filt')
  f=open(el+'.wind','w')
  for i in range(1,len(m)) :
    if m[i]>0 and m[i-1]==0 :
      wstart=w[i]
    elif m[i]==0 and m[i-1]>0 :
      f.write(str(wstart)+' '+str(w[i])+'\n')

