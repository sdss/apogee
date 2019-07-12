# encoding: utf-8
#
# @Author: Jon Holtzman
# @Date: February 2019
# @Filename: bc.py
# @License: BSD 3-Clause
# @Copyright: Jon Holtzman

# routines for APOGEE wavelength calibration

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import argparse
import barycorrpy
import os
import pdb

def getbc(ra,dec,jd,obs='APO') :
    """ Get barycentric correction using barycorrpy from Eastman & Wright,
        as implemented by Kanodia & Wright, RNAAS 2, 1
    """
    if obs == 'APO' :
        longitude = 360. - (105. + 49./60. + 13/3600.)
        latitude = 32. + 46/60. + 49./3600.
        altitude = 2788.
    elif obs == 'LCO' :
        longitude = 360. - (70 + 41/60. + 33.36/3600.)
        latitude = -1 * (29 + 0./60. + 52.56/3600.)
        altitude = 2380.
    else :
        print('Unknown observatory')
        pdb.set_trace()

    out=barycorrpy.get_BC_vel(JDUTC=jd,ra=ra,dec=dec, 
                              longi=longitude,lat=latitude,alt=altitude)
    return out[0]

def main(args) :

    parser = argparse.ArgumentParser(
        description='Gets barycentric correction')

    parser = argparse.ArgumentParser(description="get BC",
                                     prog=os.path.basename(args[0]),
                                     usage="bc ra dec jd obs")
    #parser.add_argument("ra", type=float, help="RA")
    #parser.add_argument("dec", type=float, help="DEC")
    #parser.add_argument("jd", type=float, help="JD")
    #parser.add_argument("obs", type=str, help="Observatory")
    #out=getbc(args.ra,args.dec,args.jd,obs=args.obs)[0]/1000.
    parser.add_argument("file", type=str, help="input file")
    parser.add_argument("--out", type=str, help="Output file")
    args = parser.parse_args()
    fp = open(args.file,'r')
    if args.out is not None : fout=open(args.out,'w')
    for line in fp:
        ra,dec,jd,obs = line.split()
        out=getbc(float(ra),float(dec),float(jd),obs=obs)[0]/1000.
        if args.out is None :
            print(out)
        else :
            fout.write("{:12.6f}\n".format(out))
    if args.out is not None : fout.close()

