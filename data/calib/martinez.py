from astropy.io import ascii, fits
from tools import match,plots
import os

a=ascii.read('martinez.txt')
b=ascii.read('koi.txt')
c=ascii.read('../apokasc/id_xmatch_apokasc.csv')

# get KIC of Martinez objects
i1,i2=match.match(a['KOI'],b['col2'])
kic=b['col1'][i2]
mart=a[i1]

# get 2MASS 
j1,j2=match.match(kic,c['KEPLER_INT'])
twomass=c['2MASS_ID'][j2]
mart=mart[j1]


allstar=fits.open(os.environ['APOGEE_ASPCAP']+'/dr17/synspec/allStar-dr17-synspec.fits')[1].data
i,j=match.match(allstar['APOGEE_ID'],twomass)

fig,ax=plots.multi(1,1)
plots.plotc(ax,allstar['FPARAM'][i,3],allstar['FPARAM'][i,0]-mart['Teff'][j],allstar['FPARAM'][i,0],
            xt='[M/H]',yt='ASPCAP - martinez',zt='Teff',colorbar=True)
fig.savefig('martinez1.png')
fig,ax=plots.multi(1,1)
plots.plotc(ax,allstar['FPARAM'][i,0],allstar['FPARAM'][i,1],allstar['FPARAM'][i,0]-mart['Teff'][j],
    xt='Teff',yt='log g',zt='ASPCAP-Martinez',zr=[-200,200],xr=[8000,3000],yr=[6,-1],colorbar=True)
fig.savefig('martinez2.png')



