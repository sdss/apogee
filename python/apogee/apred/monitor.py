import matplotlib.pyplot as plt
import numpy as np
from tools import plots
import pdb

def fwhm(a,mjd=[0,99999]):

    thar = np.where((a['THAR'] == 1) )[0]
    gd = np.where((a['THAR'] == 1) & (a['JD']>mjd[0]+2400000) & (a['JD']<=mjd[1]+2400000) )[0]
    xr= [a['JD'][gd].min()-2400000-10,a['JD'][gd].max()-2400000+10]
    print(len(gd))
    colors = ['r','g','b','m','c']
    fibers=[290,220,150,80,10]
    fig,ax = plots.multi(2,3,hspace=0.001,wspace=0.001,sharex=True,figsize=(14,8))
    for iline in range(2) :
        if iline == 0 : yt='FWHM'
        else : yt=''
        for ichip in range(3) :
            # use same yr for both lines
            med = np.median(a['GAUSS'][gd,0,ichip,:,2]*2.354)
            yr= [med-1,med+1]
            for ifiber in range(5) :
                plots.plotp(ax[ichip,iline],a['JD'][gd]-2400000,a['GAUSS'][gd,iline,ichip,ifiber,2]*2.354,color=colors[ifiber],yr=yr,xr=xr,
                            label='Fiber {:d}'.format(fibers[ifiber]),xt='MJD',yt=yt)
                allgd = np.where((a['GAUSS'][thar,iline,ichip,ifiber,2]*2.354 > 0.5) & (a['GAUSS'][thar,iline,ichip,ifiber,2]*2.354 < 5) )[0]
                allmed = np.median(a['GAUSS'][thar[allgd],iline,ichip,ifiber,2]*2.354)
                print(iline,ichip,ifiber,allmed)
                plots.plotl(ax[ichip,iline],ax[ichip,iline].get_xlim(),[allmed,allmed],color=colors[ifiber],xr=xr,yr=yr)
            ax[ichip,iline].legend(fontsize='xx-small',loc='upper left')
        plt.draw()
    fig.suptitle('FWHM of two (left/right) ThArNe lines ')
    pdb.set_trace()

