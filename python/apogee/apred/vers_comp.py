import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from tools import plots
from sdss import yanny
from apogee.utils import apload
import os
import pdb

chips=['a','b','c']
colors=['r','g','b']

# APOGEE-N
fig,ax=plots.multi(1,3,hspace=0.001,sharex=True,sharey=True)
t11=apload.ApLoad(apred='t11')
b=t11.ap1D(3190056)
plug=yanny.yanny(os.environ['MAPPER_DATA']+'/55880/plPlugMapM-5585-55880-01.par')
objType=np.array(plug['PLUGMAPOBJ']['objType'])
fibers=np.array(plug['PLUGMAPOBJ']['fiberId'])
tel=np.where(objType == 'SPECTROPHOTO_STD')[0]
rows=300-fibers[tel]
amed=np.median(b['a'][1].data[rows,:],axis=1)
bmed=np.median(b['b'][1].data[rows,:],axis=1)
cmed=np.median(b['c'][1].data[rows,:],axis=1)
anorm=np.median(bmed/amed)
cnorm=np.median(bmed/cmed)
anorm=1.
cnorm=1.
npix=190
design=np.zeros([3*npix*len(rows),5+len(rows)])
y=np.zeros([3*npix*len(rows)])
for ichip,chip in enumerate(chips) :
  for irow,row in enumerate(rows):
    x=b[chip][4].data - 16000.
    design[irow*3*npix+ichip*npix:irow*3*npix+ichip*npix+npix,0] = x[row,100:2000:10]**3
    design[irow*3*npix+ichip*npix:irow*3*npix+ichip*npix+npix,1] = x[row,100:2000:10]**2
    design[irow*3*npix+ichip*npix:irow*3*npix+ichip*npix+npix,2] = x[row,100:2000:10]
    if ichip == 0 : design[irow*3*npix+ichip*npix:irow*3*npix+ichip*npix+npix,3] = 1
    elif ichip == 2 : design[irow*3*npix+ichip*npix:irow*3*npix+ichip*npix+npix,4] = 1
    design[irow*3*npix+ichip*npix:irow*3*npix+ichip*npix+npix,5+irow] = 1
    y[irow*3*npix+ichip*npix:irow*3*npix+ichip*npix+npix] = np.log10(b[chip][1].data[row,100:2000:10])
gd=np.where(np.isfinite(y))[0]
design=design[gd,:]
y=y[gd]
w = np.linalg.solve(np.dot(design.T,design), np.dot(design.T, y))
def norm(w,coef,ichip) :
    x = w-16000.
    logflux = coef[0]*x**3 + coef[1]*x**2 + coef[2]*x 
    if ichip == 0 : logflux += coef[3]
    elif ichip == 2 : logflux += coef[4]
    return 10.**logflux

for row in rows[0:30:5] :
    for ichip,chip in enumerate(chips) :
        x=b[chip][4].data[row,:]
        plots.plotl(ax[0],x,b[chip][1].data[row,:],color='k',semilogy=True)
        plots.plotl(ax[0],x,b[chip][1].data[row,:]/norm(x,w,ichip),color=colors[ichip],semilogy=True)

r10=apload.ApLoad(apred='r10')
b=r10.ap1D(3190056)
for row in rows[0:30:5] :
    plots.plotl(ax[1],b['c'][4].data[row,:],b['c'][1].data[row,:],color='b',semilogy=True)
    plots.plotl(ax[1],b['b'][4].data[row,:],b['b'][1].data[row,:],color='g',semilogy=True)
    plots.plotl(ax[1],b['a'][4].data[row,:],b['a'][1].data[row,:],color='r',yr=[0,10000],semilogy=True)

r8=apload.ApLoad(apred='r8')
c=fits.open('/uufs/chpc.utah.edu/common/home/sdss/apogeework/apogee/spectro/redux/r8/red/55881/ap1D-c-03190056.fits')
b=fits.open('/uufs/chpc.utah.edu/common/home/sdss/apogeework/apogee/spectro/redux/r8/red/55881/ap1D-b-03190056.fits')
a=fits.open('/uufs/chpc.utah.edu/common/home/sdss/apogeework/apogee/spectro/redux/r8/red/55881/ap1D-a-03190056.fits')
for row in rows[0:30:5] :
    plots.plotl(ax[2],c[4].data[row,:],c[1].data[row,:],color='b',semilogy=True)
    plots.plotl(ax[2],b[4].data[row,:],b[1].data[row,:],color='g',semilogy=True)
    plots.plotl(ax[2],a[4].data[row,:],a[1].data[row,:],color='r',yr=[0,10000],semilogy=True)

for w in range(15200,17000,100) :
    for i in range(3) :
        plots.plotl(ax[i],[w,w],[0,10000],color='k')
plt.show()

# APOGEE-S

fig,ax=plots.multi(1,3,hspace=0.001,sharex=True,sharey=True)
t11=apload.ApLoad(apred='t11',telescope='lco25m',verbose=True)
b=t11.ap1D(26180012)
plug=yanny.yanny(os.environ['MAPPER_DATA_2S']+'/58180/plPlugMapM-10426-58180-01.par')
objType=np.array(plug['PLUGMAPOBJ']['objType'])
fibers=np.array(plug['PLUGMAPOBJ']['fiberId'])
tel=np.where(objType == 'SPECTROPHOTO_STD')[0]
rows=300-fibers[tel]
for row in rows[0:5]:
    plots.plotl(ax[0],b['c'][4].data[row,:],b['c'][1].data[row,:],color='b',semilogy=True)
    plots.plotl(ax[0],b['b'][4].data[row,:],b['b'][1].data[row,:],color='g',semilogy=True)
    plots.plotl(ax[0],b['a'][4].data[row,:],b['a'][1].data[row,:],color='r',yr=[0,10000],semilogy=True)

r10=apload.ApLoad(apred='r10',telescope='lco25m')
b=r10.ap1D(26180012)
for row in rows[0:5]:
    plots.plotl(ax[1],b['c'][4].data[row,:],b['c'][1].data[row,:],color='b',semilogy=True)
    plots.plotl(ax[1],b['b'][4].data[row,:],b['b'][1].data[row,:],color='g',semilogy=True)
    plots.plotl(ax[1],b['a'][4].data[row,:],b['a'][1].data[row,:],color='r',yr=[0,10000],semilogy=True)

t9=apload.ApLoad(apred='t9',telescope='lco25m')
b=t9.ap1D(26180012)
for row in rows[0:5]:
    plots.plotl(ax[2],b['c'][4].data[row,:],b['c'][1].data[row,:],color='b',semilogy=True)
    plots.plotl(ax[2],b['b'][4].data[row,:],b['b'][1].data[row,:],color='g',semilogy=True)
    plots.plotl(ax[2],b['a'][4].data[row,:],b['a'][1].data[row,:],color='r',yr=[0,10000],semilogy=True)

#r8=apload.ApLoad(apred='r8')
#c=fits.open('/uufs/chpc.utah.edu/common/home/sdss/apogeework/apogee/spectro/redux/r8/red/55881/ap1D-c-03190056.fits')
#plots.plotl(ax[2],c[4].data[100,:],c[1].data[100,:],color='b')
#c=fits.open('/uufs/chpc.utah.edu/common/home/sdss/apogeework/apogee/spectro/redux/r8/red/55881/ap1D-b-03190056.fits')
#plots.plotl(ax[2],c[4].data[100,:],c[1].data[100,:],color='g')
#c=fits.open('/uufs/chpc.utah.edu/common/home/sdss/apogeework/apogee/spectro/redux/r8/red/55881/ap1D-a-03190056.fits')
#plots.plotl(ax[2],c[4].data[100,:],c[1].data[100,:],color='r',yr=[0,10000])

for w in range(15200,17000,100) :
    for i in range(3) :
        plots.plotl(ax[i],[w,w],[0,10000],color='k')
plt.show()
