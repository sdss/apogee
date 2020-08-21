import pdb
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from apogee.aspcap import aspcap
from tools import plots

def apinspect(load,aptype,xval,yval,zval,xr=None,yr=None,zr=None,
              fparam=False,param=False,**kwargs) :

    matplotlib.use('TkAgg')
    if aptype == 'apField' :
        data = load.apField(kwargs['field'])[1].data
        spec = load.apField(kwargs['field'])[2].data
        wave = np.hstack(aspcap.apStarWave())
        specplot = [['SPEC'],['ERR'],['MASK']]
    elif aptype == 'aspcapField' :
        data = load.aspcapField(kwargs['field'])[1].data
        spec = load.aspcapField(kwargs['field'])[2].data
        wave = np.hstack(aspcap.gridWave())
        specplot = [['SPEC','SPEC_BESTFIT'],['ERR'],['MASK']]

    if xr is None : xr = limits(xval,'x')
    if yr is None : yr = limits(yval,'y')
    if zr is None : zr = limits(zval,'z')

    xdata = getdata(data,xval,fparam,param)
    ydata = getdata(data,yval,fparam,param)
    zdata = getdata(data,zval,fparam,param)

    fig,ax=plots.multi(1,1)
    plots.plotc(ax,xdata,ydata,zdata,
                xr=xr,yr=yr,zr=zr,colorbar=True,
                xt=label(xval),yt=label(yval),zt=label(zval))
    plots._data = data
    plots._id_cols=['APOGEE_ID']
    plots.event(fig)

    if spec is not None :
        fig_spec,ax_spec=plots.multi(1,len(specplot),hspace=0.001,sharex=True)
        key=' '
        while key != 'e' :
            if key == 'o' :
                obj=input('Enter APOGEE_ID: ')
                index=np.where(data['APOGEE_ID'] == obj)[0][0]
                pdb.set_trace()
            else :
                x,y,key,index = plots.mark(fig,index=True)
            for iax,splot in enumerate(specplot) :
                ax_spec[iax].cla()
                for s in splot :
                    plots.plotl(ax_spec[iax],wave,spec[s][index])
            plt.figure(fig_spec.number)
            plt.draw()


def getdata(data,val,fparam=False,param=False) :

    if val in aspcap.params()[2] :
        if fparam :
            vals = np.squeeze(data['FPARAM'][:,index(val)])
        elif param :
            vals = np.squeeze(data['PARAM'][: index(val)])
        else :
            vals = data[val]
    else :
        vals = data[val]

    if vals is None : print("Can't find quantity: ",val)

    return vals


def label(val) :
    if val == 'TEFF' : return 'Teff'
    if val == 'LOGG' : return 'log g'
    if val == 'M_H' : return '[M/H]'
    if val == 'ALPHA_M' : return '[alpha/M]'

def limits(val,axis) :
    """ get default limits for input quantity
    """
    if val == 'TEFF' : 
        lim=[3000,8000]
        if axis=='x' : lim.reverse()
    if val == 'LOGG' : 
        lim=[-1,6]
        if axis=='y' : lim.reverse()
    if val == 'M_H' : lim= [-2,.5]
    if val == 'ALPHA_M' : lim=[-0.5,0.75]
    try : return lim
    except : return
    
def index(val) :
    """ Return index in param array for desired quantity
    """
    pars = aspcap.params()[2]
    return np.where(pars == val)[0]
