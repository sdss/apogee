import numpy as np
from astropy.table import Column
from astropy.io import fits
from apogee.utils import bitmask
from sdss import yanny
import os
import pdb

plans = None

def add_design(tab) :
    """ Add design information for objects
    """
    global plans

    for col in ['MIN_H','MAX_H','MIN_JK','MAX_JK'] :
        try : tab.add_column(Column(name=col,dtype=np.float32,length=len(tab)) )
        except : pass
        if 'MIN' in col : tab[col] = -9999.99
        else : tab[col] = 9999.99

    if tab['TELESCOPE'][0] == 'apo1m' : return

    if plans is None : plans = yanny.yanny(os.environ['PLATELIST_DIR']+'/platePlans.par')['PLATEPLANS']
    apogee_design = fits.open(os.environ['APOGEE_TARGET']+'/apogeeDesign.fits')[1].data
    apogee2_design = fits.open(os.environ['APOGEE_TARGET']+'/apogee2Design.fits')[1].data

    plate = int(tab['PLATE'][0]) 
    iplan = np.where(np.array(plans['plateid']) == plate)[0]
    if len(iplan) != 1 : print('zero or multiple plates found in platePlans')
    else : iplan=iplan[0]

    if 'apogee2' in tab['SURVEY'][0] :
        idesign=np.where(apogee2_design['design_id'] == plans['designid'][iplan])[0][0]
        min_h = apogee2_design['cohort_min_h'][idesign]
        max_h = apogee2_design['cohort_max_h'][idesign]
        targflag = bitmask.Apogee2Target1()
        j=np.where(tab['APOGEE_TARGET1']&targflag.getval('APOGEE2_SHORT'))[0]
        tab['MIN_H'][j] = min_h[0]
        tab['MAX_H'][j] = max_h[0]
        j=np.where(tab['APOGEE_TARGET1']&targflag.getval('APOGEE2_MEDIUM'))[0]
        tab['MIN_H'][j] = min_h[1]
        tab['MAX_H'][j] = max_h[1]
        j=np.where(tab['APOGEE_TARGET1']&targflag.getval('APOGEE2_LONG'))[0]
        tab['MIN_H'][j] = min_h[2]
        tab['MAX_H'][j] = max_h[2]
        j=np.where(tab['APOGEE_TARGET1']&targflag.getval('APOGEE2_ONEBIN_GT_0_5'))[0]
        tab['MIN_JK'][j] = 0.5
        j=np.where(tab['APOGEE_TARGET1']&targflag.getval('APOGEE2_TWOBIN_0_5_TO_0_8'))[0]
        tab['MIN_JK'][j] = 0.5
        tab['MAX_JK'][j] = 0.8
        j=np.where(tab['APOGEE_TARGET1']&targflag.getval('APOGEE2_TWOBIN_GT_0_8'))[0]
        tab['MIN_JK'][j] = 0.8
        j=np.where(tab['APOGEE_TARGET1']&targflag.getval('APOGEE2_ONEBIN_GT_0_3'))[0]
        tab['MIN_JK'][j] = 0.3
              
    else :
        idesign=np.where(apogee_design['design_id'] == plans['designid'][iplan])[0][0]
        min_h = [apogee_design['short_cohort_min_h'][idesign],
                 apogee_design['medium_cohort_min_h'][idesign],
                 apogee_design['long_cohort_min_h'][idesign]]
        max_h = [apogee_design['short_cohort_max_h'][idesign],
                 apogee_design['medium_cohort_max_h'][idesign],
                 apogee_design['long_cohort_max_h'][idesign]]
        min_jk = apogee_design['DEREDDENED_MIN_J_KS_COLOR'][idesign]
        targflag = bitmask.ApogeeTarget1()
        j=np.where(tab['APOGEE_TARGET1']&targflag.getval('APOGEE_SHORT'))[0]
        tab['MIN_H'][j] = min_h[0]
        tab['MAX_H'][j] = max_h[0]
        tab['MIN_JK'][j] = min_jk
        j=np.where(tab['APOGEE_TARGET1']&targflag.getval('APOGEE_INTERMEDIATE'))[0]
        tab['MIN_H'][j] = min_h[1]
        tab['MAX_H'][j] = max_h[1]
        tab['MIN_JK'][j] = min_jk
        j=np.where(tab['APOGEE_TARGET1']&targflag.getval('APOGEE_LONG'))[0]
        tab['MIN_H'][j] = min_h[2]
        tab['MAX_H'][j] = max_h[2]
        tab['MIN_JK'][j] = min_jk
