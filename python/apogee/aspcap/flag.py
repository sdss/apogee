import numpy as np

from apogee.utils import bitmask
from apogee.aspcap import aspcap
from apogee.aspcap import teff

def aspcapflag(aspcapfield) :
    """ Set bits in ASPCAPFLAG
    """

    parambitmask=bitmask.ParamBitMask()
    aspcapbitmask=bitmask.AspcapBitMask()

    gd=np.where((aspcapfield['ASPCAPFLAG'] & aspcapbitmask.getval('NO_ASPCAP_RESULT') ==0) &
                (aspcapfield['ASPCAPFLAG'] & aspcapbitmask.getval('NO_GRID') == 0))  [0]

    # set ASPCAPFLAG bits for grid edge with final adopted grid
    for iparam,flagname in enumerate(aspcap.params()[2]) :
        j=np.where(aspcapfield['PARAMFLAG'][gd,iparam] & parambitmask.getval('GRIDEDGE_BAD') )[0]
        aspcapfield['ASPCAPFLAG'][gd[j]] |= aspcapbitmask.getval(flagname+'_BAD')
        print(flagname+'_BAD',len(j))
        j=np.where(aspcapfield['PARAMFLAG'][gd,iparam] & parambitmask.getval('GRIDEDGE_WARN') )[0]
        aspcapfield['ASPCAPFLAG'][gd[j]] |= aspcapbitmask.getval(flagname+'_WARN') 
        print(flagname+'_WARN',len(j))
    
    # chi**2 
    tmp= aspcapfield['ASPCAP_CHI2']/(aspcapfield['SNR']/100.)**2
    j=np.where(tmp[gd] > 50)[0]
    aspcapfield['ASPCAPFLAG'][gd[j]] |= aspcapbitmask.getval('CHI2_BAD')
    print('CHI2_BAD',len(j))
    j=np.where((tmp[gd] > 30) & (aspcapfield['ASPCAPFLAG'][gd]&aspcapbitmask.getval('CHI2_BAD')==0) )[0]
    aspcapfield['ASPCAPFLAG'][gd[j]] |= aspcapbitmask.getval('CHI2_WARN')
    print('CHI2_WARN',len(j))

    # rotation in giant grids
    j=np.where( ( (np.core.defchararray.find(aspcapfield['ASPCAP_CLASS'][gd].astype(str),'GKg') >=0)  |
                  (np.core.defchararray.find(aspcapfield['ASPCAP_CLASS'][gd].astype(str),'Mg') >=0) ) &
                (aspcapfield['RV_CCFWHM'][gd]/aspcapfield['RV_AUTOFWHM'][gd] > 2.0 ) ) [0]
    aspcapfield['ASPCAPFLAG'][gd[j]] |= aspcapbitmask.getval('ROTATION_BAD')
    print('ROTATION_BAD',len(j))
    j=np.where( ( (np.core.defchararray.find(aspcapfield['ASPCAP_CLASS'][gd].astype(str),'GKg') >=0)  |
                  (np.core.defchararray.find(aspcapfield['ASPCAP_CLASS'][gd].astype(str),'Mg') >=0) ) &
                (aspcapfield['RV_CCFWHM'][gd]/aspcapfield['RV_AUTOFWHM'][gd] > 1.5 ) &
                (aspcapfield['ASPCAPFLAG'][gd]&aspcapbitmask.getval('ROTATION_BAD')==0) ) [0]
    aspcapfield['ASPCAPFLAG'][gd[j]] |= aspcapbitmask.getval('ROTATION_WARN')
    print('ROTATION_WARN',len(j))

    # S/N
    j=np.where(aspcapfield['SNR'][gd] < 30)[0]  
    aspcapfield['ASPCAPFLAG'][gd[j]] |= aspcapbitmask.getval('SN_BAD')
    print('SNR_BAD',len(j))
    j=np.where((aspcapfield['SNR'][gd] < 70) & ( aspcapfield['ASPCAPFLAG'][gd]&aspcapbitmask.getval('SN_BAD')) )[0]  
    aspcapfield['ASPCAPFLAG'][gd[j]] |= aspcapbitmask.getval('SN_WARN')
    print('SNR_WARN',len(j))

    #color-Teff
    jk0=(aspcapfield['J']-aspcapfield['K'])-1.5*np.clip(aspcapfield['AK_TARG'],0.,None)
    j = np.where( (aspcapfield['H'][gd] < 90) & (aspcapfield['AK_TARG'][gd] != 0.) &
                  (np.abs(aspcapfield['FPARAM'][gd,0]-teff.cte_ghb(jk0[gd],aspcapfield['FPARAM'][gd,3])[0]) > 1000) )[0]
    aspcapfield['ASPCAPFLAG'][gd[j]] |= aspcapbitmask.getval('COLORTE_BAD')
    print('COLORTE_BAD',len(j))
    j = np.where( (aspcapfield['H'][gd] < 90) & (aspcapfield['AK_TARG'][gd] != 0.) &
                  (np.abs(aspcapfield['FPARAM'][gd,0]-teff.cte_ghb(jk0[gd],aspcapfield['FPARAM'][gd,3])[0]) > 500) &
                  (aspcapfield['ASPCAPFLAG'][gd] & aspcapbitmask.getval('COLORTE_BAD')==0) )[0]
    aspcapfield['ASPCAPFLAG'][gd[j]] |= aspcapbitmask.getval('COLORTE_WARN')
    print('COLORTE_WARN',len(j))

    # bad targets
    j = np.where( (np.core.defchararray.find(aspcapfield['TARGFLAGS'][gd].astype(str),'EXTENDED') >=0) |
                  (np.core.defchararray.find(aspcapfield['TARGFLAGS'][gd].astype(str),'EMBEDDED') >=1) |
                  (np.core.defchararray.find(aspcapfield['TARGFLAGS'][gd].astype(str),'APOGEE2_M31') >=1) |
                  (np.core.defchararray.find(aspcapfield['TARGFLAGS'][gd].astype(str),'APOGEE2_M33') >=1) |
                  (np.core.defchararray.find(aspcapfield['TARGFLAGS'][gd].astype(str),'APOGEE2_QSO') >=1) ) [0]
    aspcapfield['ASPCAPFLAG'][gd[j]] |= aspcapbitmask.getval('PROBLEM_TARGET')
    print('PROBLEM_TARGET',len(j))

    # star level flag
    j = np.where( aspcapfield['ASPCAPFLAG']&aspcapbitmask.badval() > 0)[0]
    aspcapfield['ASPCAPFLAG'][j] |= aspcapbitmask.getval('STAR_BAD')
    print('STAR_BAD',len(j))
    j = np.where( (aspcapfield['ASPCAPFLAG']&aspcapbitmask.warnval() > 0) &
                  (aspcapfield['ASPCAPFLAG']&aspcapbitmask.getval('STAR_BAD') == 0) )[0]
    aspcapfield['ASPCAPFLAG'][j] |= aspcapbitmask.getval('STAR_WARN')
    print('STAR_WARN',len(j))

    #populate character ASPCAPFLAGS
    maxlen=0
    for istar in range(len(aspcapfield)) :
        # set ASPCAPFLAGS character string
        aspcapflags=aspcapbitmask.getname(aspcapfield['ASPCAPFLAG'][istar])
        if len(aspcapflags) > maxlen: 
            maxlen=len(aspcapflags)
            print(aspcapflags)
        aspcapfield['ASPCAPFLAGS'][istar] = aspcapflags
    print('maxlen: ', maxlen)

