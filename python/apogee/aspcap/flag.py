from apogee.utils import bitmask
from apogee.aspcap import aspcap

def aspcap(aspcapfield) :
    """ Set bits in ASPCAPFLAG
    """

    # set ASPCAPFLAG bits for grid edge with final adopted grid
    parambitmask=bitmask.ParamBitMask()
    for istar,star in enumerate(aspcapfield) :
        for iparam,flagname in enumerate(aspcap.params()[2]) :
            if aspcapfield['PARAMFLAG'][istar,iparam] & parambitmask.getval('GRIDEDGE_BAD') :
                aspcapfield['ASPCAPFLAG'][istar] |= aspcapmask.getval(flagname+'_BAD')
            if aspcapfield['PARAMFLAG'][istar,iparam] & parambitmask.getval('GRIDEDGE_WARN') :
                aspcapfield['ASPCAPFLAG'][istar] |= aspcapmask.getval(flagname+'_WARN')

        # set ASPCAPFLAGS character string
        aspcapfield['ASPCAPFLAGS'][istar] = aspcapmask.getname(aspcapfield['ASPCAPFLAG'][istar])

