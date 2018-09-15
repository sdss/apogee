from tools import html
from tools import plots
from astropy.io import fits
import pdb
import numpy as np
import matplotlib.pyplot as plt

# routine for comparing output of syntheses with different dw

def wave_fits(hdr) :
    """ Returns wavelength array given FITS header
    """
    wave = hdr['CRVAL1']+np.arange(hdr['NAXIS1'])*hdr['CDELT1']
    if hdr['LOGW'] == 0 :
      return wave
    elif hdr['LOGW'] == 1 :
      return 10.**wave
    elif hdr['LOGW'] == 2 :
      return np.exp(wave)

# list of file name suffixes and microturbulences
#suffix='cn'
#vplist=['03','06','12','24','48']
suffix='c'
vplist=['05','10','20','40','80']
y=[]
ylab=[]

# loop through different dw
for t1 in ['01','02','03','05','10'] :
  x=[]
  xlab=[]
  for t2 in ['01','02','03','05','10'] :
    # loop through different dw for comparison
    fig,ax=plots.multi(1,5,hspace=0.001) 
    for iplt,vp in enumerate(vplist) :
      a=fits.open('test'+t1+suffix+'_lsfcombo5/new_ap00cp00np00vp'+vp+'.fits')[0]
      b=fits.open('test'+t2+suffix+'_lsfcombo5/new_ap00cp00np00vp'+vp+'.fits')[0]

      # loop over 3 Teff, 3 logg, 3 [M/H]
      for i in range(3) :
        for j in range(3) :
          for k in range(3) :
            if k == 0 : color = 'r'
            elif k == 1 : color = 'g'
            elif k == 2 : color = 'b'
            plots.plotl(ax[iplt],wave_fits(a.header),a.data[i,j,k,:]/b.data[i,j,k,:],yr=[0.95,1.05],color=color)
            ax[iplt].text(0.05,0.9,r'$v_{micro}= $'+'{:3.1f}'.format(float(vp)/10.),transform=ax[iplt].transAxes,va='top')
      fig.savefig('dw/'+t1+'_'+t2+suffix+'.png')
      plt.close(fig)
    x.append('dw/'+t1+'_'+t2+suffix+'.png')
    xlab.append('dw={:4.3f}'.format(float(t2)/100.))
  y.append(x)
  ylab.append(t1)

# create output HTML table with figures and links
html.htmltab(y,file='dw'+suffix+'.html',xtitle=xlab,ytitle=ylab)

