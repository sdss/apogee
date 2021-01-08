from astropy.table import Table, Column
import numpy as np

a=Table.read('logg_cal.csv')
col = Column(np.full([len(a)],''),name='CONS_EVSTATES',dtype='S8')
a.add_column(col)
j=np.where(a['ES_AST']==2)[0]
a['CONS_EVSTATES'][j] = 'RC'
j=np.where(a['ES_AST']==1)[0]
a['CONS_EVSTATES'][j] = 'RGB'

col = Column(np.full([len(a)],np.nan),name='APOKASC2_LOGG',dtype=np.float32)
a.add_column(col)
a['APOKASC2_LOGG'] = a['LOGG_SEIS']
a['LOGG_DW'] = a['LOGG_SEIS']

a.write('logg_cal.fits',overwrite=True)


