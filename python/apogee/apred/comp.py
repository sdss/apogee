import os
import numpy as np
from tools import plots
from apogee.utils import apload
from apogee.utils import bitmask
from astropy.io import fits

chips=['a','b','c']

def comp1d(frame,apred='test',rows=range(300)) :
    load=apload.ApLoad(apred=apred)
    new=load.ap1D(frame)
    old={}
    mjd=55562+int(frame//10000)
    fig,ax = plots.multi(1,3,hspace=0.001)
    x=np.arange(2048)
    for ichip,chip in enumerate(chips) :
        old[chip]=fits.open(os.environ['APOGEE_REDUX']+'/r8/red/{:d}/ap1D-{:s}-{:d}.fits'.format(mjd,chip,frame))
        for row in rows :
            plots.plotl(ax[ichip],x,new[chip][1].data[row,:]/old[chip][1].data[row,:],yr=[0,1.5])

def compCframe(plate,frame,apred='test',ratio=True,rows=range(300),yr=None,hdu=1) :
    load=apload.ApLoad(apred=apred)
    mjd=55562+int(frame//10000)
    new=load.apCframe('M67',plate,mjd,frame)
    old={}
    fig,ax = plots.multi(1,3,hspace=0.001)
    x=np.arange(2048)
    for ichip,chip in enumerate(chips) :
        old[chip]=fits.open(os.environ['APOGEE_REDUX']+'/r8/apo25m/{:d}/{:d}/apCframe-{:s}-{:d}.fits'.format(plate,mjd,chip,frame))
        for row in rows :
            if ratio :
               plots.plotl(ax[ichip],x,new[chip][hdu].data[row,:]/old[chip][hdu].data[row,:],yr=[0,1.5])
            else :
               plots.plotl(ax[ichip],x,new[chip][hdu].data[row,:],yr=yr)
               plots.plotl(ax[ichip],x,old[chip][hdu].data[row,:],yr=yr)
               plots.plotl(ax[ichip],x,new[chip][hdu].data[row,:]-old[chip][hdu].data[row,:],yr=yr)


#def comp2d(frame) :
#
def comp(plate=7267,mjd=56654,fiber=150,frame=10920059,field='M67') :
    r11=apload.ApLoad(apred='r11')

    v=r11.apVisit(plate,mjd,fiber)
    a=r11.ap1D(frame)
    c=r11.apCframe(field,plate,mjd,frame)

    v14=fits.open(os.environ['APOGEE_REDUX']+'/r8/apo25m/{:d}/{:d}/apVisit-r8-{:d}-{:d}-{:03d}.fits'.format(plate,mjd,plate,mjd,fiber))
    a14={}
    c14={}
    for chip in chips:
        a14[chip]=fits.open(os.environ['APOGEE_REDUX']+'/r8/red/{:d}/ap1D-{:s}-{:d}.fits'.format(mjd,chip,frame))
        c14[chip]=fits.open(os.environ['APOGEE_REDUX']+'/r8/apo25m/{:d}/{:d}/apCframe-{:s}-{:08d}.fits'.format(plate,mjd,chip,frame))

    fig,ax=plots.multi(1,3,hspace=0.01)
    x=np.arange(4096)
    pixmask=bitmask.PixelBitMask()
    for ichip,chip in enumerate(chips) :
        y=v[1].data[ichip,:]
        plots.plotl(ax[ichip],x,v[1].data[ichip,:]/v14[1].data[ichip,:])
        bd = np.where( ((v[3].data[ichip,:] & pixmask.badval()) > 0) |
                       ((v[3].data[ichip,:] & pixmask.getval('SIG_SKYLINE')) > 0) ) [0]
        y[bd]=np.nan
        plots.plotl(ax[ichip],x,y/v14[1].data[ichip,:])


    fig,ax=plots.multi(3,3,hspace=0.01)
    x=np.arange(2048)
    for ichip,chip in enumerate(chips) :
        plots.plotl(ax[ichip,0],x,c[chip][1].data[300-fiber,:])
        plots.plotl(ax[ichip,0],x,c14[chip][1].data[300-fiber,:])
        plots.plotl(ax[ichip,1],x,c[chip][1].data[300-fiber,:]/c14[chip][1].data[300-fiber])
        plots.plotl(ax[ichip,2],x,a[chip][1].data[300-fiber,:]/a14[chip][1].data[300-fiber])
