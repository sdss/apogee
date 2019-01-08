from apogee.utils import apload
from pyvista import tv
import matplotlib.pyplot as plt
from tools import plots
import numpy as np
import pdb

def lsfsum(y) :
  for col in range(500,2500,500) :
    print( y[1].data[:,150,col].sum(),y[1].data[:,150,col].sum(),y[1].data[:,150,col].sum())

load=apload.ApLoad(apred='r11')
fig,ax=plots.multi(4,3,hspace=0.001,wspace=0.001)
pfig,pax=plots.multi(9,3)
for ichip,chip in enumerate(['a','b','c']) :
    #t=tv.TV()
    y1=load.apLSF(3430016)[chip]
    print('y1: ')
    lsfsum(y1)
    y2=load.apLSF(7510018)[chip]
    print('y2: ')
    lsfsum(y2)
    y3=load.apLSF(11130063)[chip]
    print('y3: ')
    lsfsum(y3)
    y4=load.apLSF(14600018)[chip]
    print('y4: ')
    lsfsum(y4)
    y5=load.apLSF(18430026)[chip]
    print('y5: ')
    lsfsum(y5)
    y6=load.apLSF(22330043)[chip]
    print('y6: ')
    lsfsum(y6)
    y7=load.apLSF(25560065)[chip]
    print('y7: ')
    lsfsum(y7)
    ii=0
    for ipar in range(26) :
        iy=ipar%3
        ix=ipar//3
        pax[iy,ix].plot(y1[0].data[ii,:])
        ii+=1
    for icol in range(4) :
      for fiber in [150] :
        col=500+icol*500
        n=y1[1].data.shape[0]
        plots.plotl(ax[ichip,icol],np.arange(21),y1[1].data[n/2-10:n/2+11,fiber,col])
        n=y2[1].data.shape[0]
        plots.plotl(ax[ichip,icol],np.arange(21),y2[1].data[n/2-10:n/2+11,fiber,col])
        n=y3[1].data.shape[0]
        plots.plotl(ax[ichip,icol],np.arange(21),y3[1].data[n/2-10:n/2+11,fiber,col])
        n=y4[1].data.shape[0]
        plots.plotl(ax[ichip,icol],np.arange(21),y4[1].data[n/2-10:n/2+11,fiber,col])
        n=y5[1].data.shape[0]
        plots.plotl(ax[ichip,icol],np.arange(21),y5[1].data[n/2-10:n/2+11,fiber,col])
        n=y6[1].data.shape[0]
        plots.plotl(ax[ichip,icol],np.arange(21),y6[1].data[n/2-10:n/2+11,fiber,col])
        n=y7[1].data.shape[0]
        plots.plotl(ax[ichip,icol],np.arange(21),y7[1].data[n/2-10:n/2+11,fiber,col])
    #t.tv(y1[1].data[:,:,1024])
    #t.tv(y1[1].data[:,:,1024])
    ##t.tv(y3[:,:,1024])
    #t.tv(y4[1].data[:,:,1024])
    #t.tv(y7[1].data[:,:,1024])
    plt.show()
    pdb.set_trace()
