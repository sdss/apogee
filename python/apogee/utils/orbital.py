import numpy 
import galpy
from galpy.orbit import Orbit
from galpy.potential import MWPotential2014

def parameters(tab,dist='gaia') :
    """ Add orbital parameters from galpy quick computation
        needs distances
    """

    for tag in ['ECC','R_PERI','R_AP','Z_MAX'] :
        try: tab.remove_column(tag)
        except KeyError: pass
        tab[tag] = numpy.nan

    if dist == 'gaia' :
        dist=tab['GAIAEDR3_R_MED_GEO']
    else :
        dist=tab['diso'][:,1]

    n = len(tab)
    # do 10000 at a time
    for i1 in numpy.arange(1,n+10000,10000) :
        i2=numpy.min([i1+10000,len(tab)])
        print(i1,i2)
        o=Orbit(numpy.dstack([tab['RA'][i1:i2],tab['DEC'][i1:i2],
                dist[i1:i2]/1000.,tab['GAIAEDR3_PMRA'][i1:i2],tab['GAIAEDR3_PMDEC'][i1:i2],
                tab['VHELIO_AVG'][i1:i2]])[0],ro=8.,vo=220.,radec=True)
        tab['ECC'][i1:i2] =  o.e(analytic=True,pot=MWPotential2014)
        tab['R_PERI'][i1:i2] = o.rperi(analytic=True,pot=MWPotential2014)
        tab['R_AP'][i1:i2] = o.rap(analytic=True,pot=MWPotential2014)
        tab['Z_MAX'][i1:i2] = o.zmax(analytic=True,pot=MWPotential2014)

    return tab

