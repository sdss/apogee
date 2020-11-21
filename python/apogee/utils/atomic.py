import numpy as np

def periodic(n) :
    """ Routine to get element name / atomic number conversion
    """
    elem=np.array(['H','He','Li','Be','B','C','N','O','F','Ne',
                   'Na','Mg','Al','Si','P','S','Cl','Ar',
                   'K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
                   'Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe',
                   'Cs','Ba','La','Ce','Pr','Nd'])
    if isinstance(n,str) :
        j=np.where(elem == n)[0]
        return j+1
    else :
        if n == 0 :
            return ''    
        else :
            return elem[n-1]

def solar(el=None) :
    """ Return solar abundances
    """
    sunabund_2007  =  np.array([
       12.00, 10.93,  1.05,  1.38,  2.70,  8.39,  7.78,  8.66,  4.56,   #  1 -  9
        7.84,  6.17,  7.53,  6.37,  7.51,  5.36,  7.14,  5.50,  6.18,   # 10 - 18
        5.08,  6.31,  3.17,  4.90,  4.00,  5.64,  5.39,  7.45,  4.92,   # 19 - 27
        6.23,  4.21,  4.60,  2.88,  3.58,  2.29,  3.33,  2.56,  3.25,   # 28 - 36
        2.60,  2.92,  2.21,  2.58,  1.42,  1.92, -99.0,  1.84,  1.12,   # 37 - 45
        1.66,  0.94,  1.77,  1.60,  2.00,  1.00,  2.19,  1.51,  2.24,   # 46 - 54
        1.07,  2.17,  1.13,  1.70,  0.58,  1.45, -99.0,  1.00,  0.52,   # 55 - 63
        1.11,  0.28,  1.14,  0.51,  0.93,  0.00,  1.08,  0.06,  0.88,   # 64 - 72
       -0.17,  1.11,  0.23,  1.25,  1.38,  1.64,  1.01,  1.13,  0.90,   # 73 - 81
        2.00,  0.65, -99.0, -99.0, -99.0, -99.0, -99.0, -99.0,  0.06,   # 82 - 90
       -99.0, -0.52, -99.0, -99.0, -99.0, -99.0, -99.0, -99.0, -99.0 ]) # 91 - 99
    if el is None : return sunabund_2007
    n = periodic(el)
    return sunabund_2007[n-1]

def rydberg(n1,n2) :
    """ Rydberg formula to give H wavelengths for transitions between levels n1 and n2
    """
    r=1.0973731568e7
    me=9.109382e-31
    mprot=1.672621e-27
    rm=r/(1+me/mprot)

    l=rm*(1./n1**2-1./n2**2)
    w=1./l*1.e10
    return w

def hlines(plot=None,yloc=0.,n1=4,n2=range(11,22)) :
    """ Return approximate (Rydberg) location of H lines, defaulting to lines in APOGEE spectra
    """
    h=[]
    for n in n2 :
        h.append(rydberg(n1,n))
    h=np.array(h)
    if plot is not None :
        plots.plotp(plot,h,h*0.+yloc)
    return h
