from astropy.table import Table, Column, vstack
from astropy.io import ascii, fits
from tools import match
import numpy as np
import pdb

xmatch=ascii.read('../apokasc/id_xmatch_apokasc.csv')

# Berger et al.
a=ascii.read('GKSPCPapTable2_Final.csv')
i1,i2=match.match(a['KIC'],xmatch['KEPLER_INT'])
gd=np.where((a['iso_teff_err1'][i1]<100)&(a['iso_teff_err2'][i1]<100))[0]
print('number of berger stars: ', len(gd))
berger_tab=Table()
berger_tab.add_columns([Column(xmatch['2MASS_ID'][i2[gd]],name='APOGEE_ID'),
                Column(a['iso_teff'][i1[gd]],name='TEFF',dtype=np.float32),
                Column(a['iso_teff_err1'][i1[gd]],name='TEFF_ERR',dtype=np.float32),
                Column(np.zeros(len(i2[gd]))+1,name='SOURCE',dtype=np.int)]
                )

berger_tab.write('teff_berger.fits')
