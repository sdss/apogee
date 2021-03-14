from astropy.table import Table, Column, vstack
from astropy.io import ascii, fits
from tools import match
import numpy as np
import pdb

xmatch=ascii.read('../apokasc/id_xmatch_apokasc.csv')

# APOKASC3 gold sample from Marc P
a=ascii.read('../apokasc/APOKASC_cal')
i1,i2=match.match(a['KIC'],xmatch['KEPLER_INT'])
print('number of APOKASC stars: ', len(i1))

pdb.set_trace()
col=Column(xmatch['2MASS_ID'][i2],name='APOGEE_ID')
apokasc_tab=Table()
apokasc_tab.add_columns([Column(xmatch['2MASS_ID'][i2],name='APOGEE_ID'),
                Column(a['LOGG_S'][i1],name='LOGG',dtype=np.float32),
                Column(a['SIG_LOGG_S'][i1],name='LOGG_ERR',dtype=np.float32),
                Column(np.zeros(len(i2)),name='SOURCE',dtype=np.int)]
                )

a['APOGEE_ID'] = col


# Berger et al.
a=ascii.read('GKSPCPapTable2_Final.csv')
i1,i2=match.match(a['KIC'],xmatch['KEPLER_INT'])
gd=np.where((a['iso_logg_err1'][i1]<0.03)&(a['iso_logg_err2'][i1]<0.03))[0]
print('number of berger stars: ', len(gd))
berger_tab=Table()
berger_tab.add_columns([Column(xmatch['2MASS_ID'][i2[gd]],name='APOGEE_ID'),
                Column(a['iso_logg'][i1[gd]],name='LOGG',dtype=np.float32),
                Column(a['iso_logg_err1'][i1[gd]],name='LOGG_ERR',dtype=np.float32),
                Column(np.zeros(len(i2[gd]))+1,name='SOURCE',dtype=np.int)]
                )

#APOKASC LOGG_DW
a=fits.open('../apokasc/APOKASC_cat_v6.6.5.fits.gz')[1].data
gd=np.where((a['LOGG_DW'] > -1)& (a['LOGG_DW'] < 6) )[0]
print('number of LOGG_DW stars: ', len(gd))
logg_dw=Table()
logg_dw.add_columns([Column(a['2MASS_ID'][gd],name='APOGEE_ID'),
                Column(a['LOGG_DW'][gd],name='LOGG',dtype=np.float32),
                Column(a['LOGG_DW_SYSERR'][gd],name='LOGG_ERR',dtype=np.float32),
                Column(np.zeros(len(gd))+2,name='SOURCE',dtype=np.int)]
                )

# cluster physical gravities
a=ascii.read('clusters_physical_logg.txt')
physical_tab=Table()
physical_tab.add_columns([Column(a['col1'],name='APOGEE_ID'),
                Column(a['col2'],name='LOGG',dtype=np.float32),
                Column(np.zeros(len(a))-1.,name='LOGG_ER',dtype=np.float32),
                Column(np.zeros(len(a))+3,name='SOURCE',dtype=np.int)]
                )

tab=vstack([apokasc_tab,berger_tab,logg_dw,physical_tab])
tab.write('logg_cal.fits')
