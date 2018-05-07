import numpy as np
from astropy.io import fits
from tools import plots
import pdb
import matplotlib.pyplot as plt

fig,ax=plots.multi(3,1,figsize=(12,5),wspace=0.001)
plt.show()

darks=[5560001,12910009,15640003,23870002]
colors=['k','r','g','b']
label=['07/10/2012','07/15/2014','04/14/2015','07/17/2017']
for ichip,chip in enumerate(['a','b','c']) :
    for idark,dark in enumerate(darks) :
        a=fits.open('apDarkRate-{:s}-{:08d}.fits'.format(chip,dark))[0].data
        ax[ichip].hist(a.flatten(),bins=np.arange(-0.5,10,0.1),histtype='step',log=True,color=colors[idark])
        ax[ichip].set_xlabel('DN per 10.6s readout')
        if ichip == 0 : ax[ichip].set_ylabel('Number of pixels')
        ax[ichip].text(0.8,0.8-idark*0.05,label[idark],ha='right',transform=ax[ichip].transAxes,color=colors[idark])
    ax[ichip].text(0.5,0.9,'APOGEE-N chip '+chip,ha='center',transform=ax[ichip].transAxes)
    plt.draw()


