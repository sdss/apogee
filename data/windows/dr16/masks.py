from __future__ import print_function
from apogee.aspcap import aspcap
import numpy as np
import pdb

out=open('masks.dat','w')
out.write('{:<12s}'.format('Wave'))
els=aspcap.elems()
masks=[]
for el in els[0]:
  masks.append(np.loadtxt(el+'.mask'))
  out.write('{:<8s}'.format('   '+el))
out.write('\n')

w=np.loadtxt('wave.dat')
for iline, line in enumerate(w) :
  out.write('{:12.3f}'.format(line[0]))
  for mask in masks :
    out.write('{:8.3f}'.format(mask[iline]))
  out.write('\n')
