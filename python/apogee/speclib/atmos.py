# encoding: utf-8
#
# @Author: Jon Holtzman
# @Date: March 2018
# @Filename: atmos.py
# @License: BSD 3-Clause
# @Copyright: Jon Holtzman

#Routines for working with model atmospheres, including
#finding and filling "holes" in a rectangular atmospheres grid.
#Main routine is fill_holes(). Takes as input the grid parameters (start, end, step),
#looks for existing atmospheres file. If file not found, adopt "nearest" atmosphere
#with algorithm as given in find_filler(), and copy this into the missing
#file location. For Kurucz models, need to then "fix" this file for appropriate
#parameters; MARCS models are left as is.
#
#Also outputs holefile that gives location of holes and how they were filled.
#
#J. Holtzman, working off of Neville Shane's IDL scripts to do some of these tasks.

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import os
import numpy as np
import sys
import subprocess
import shutil
import argparse
import pdb
from astropy.io import fits

def cval(x) :
    """ routine to convert value to "Kurucz-style" string, i.e. mXX or pXX

    Args:
        x (float) : value to put into output format
    Returns :
        (str) : string with output
    """
    if x < -0.000001 :
      prefix = 'm'
    else :
      prefix = 'p'
    return prefix+'{:02d}'.format(int(round(abs(x)*10.)))

def filename(teff,logg,z,c,a,model='Unknown',vers='x3') :
    """Function to return atmosphere file name for either MARCS-style or KURUCZ-style 

    Args:
        teff (int) : effective temperature
        logg (float) : surface gravity
        z (float) : [M/H]
        c (float) : [C/M]
        model (str) : atmosphere type (default='Unknown' --> kurucz)
        vers (str) : for MARCS, "version number" for filename (default='x3')
    """
    if model.upper() == 'MARCS' :
        # MARCS file names
        if logg <= 3.001 :
            prefix = 's'
            mass = 1.
            t = 2
        else :
            prefix = 'p'
            mass = 0.
            t = 1
        dir = 'mod_z{:+.2f}'.format(z)
        file = (dir+'/'+prefix+'{:04d}_g{:+.1f}_m{:.1f}_t{:02d}_'+vers+'_z{:+.2f}_a{:+.2f}_c{:+.2f}_n+0.00_o{:+.2f}_r+0.00_s+0.00').format(int(teff),logg,mass,t,z,a,c,a) + '.mod'
    else :
        # KURUCZ file names
        dir = 'm'+cval(z)+'c'+cval(c)+'o'+cval(a)
        file = dir + '/' + 'am'+cval(z)+'c'+cval(c)+'o'+cval(a)+'t'+str(int(teff))+'g'+'{:02d}'.format(int(logg*10.))+'v20.mod'

    return file

def find_filler(all,i,model='Unknown') :
    """ Find "closest" existing model for an atmosphere hole using "distance" scheme suggested largely by Matt Shetrone 

    Args :
        all (recarray) : input list of atmospheres
        i (int) : index of hole
        model (str) : type of atmopshere (default='Unknown' --> kurucz)
    Returns :
        fill : index of "closest" atmosphere
        dist : "distance" to closest atmosphere
    """

    # calculate distance of all models from the hole
    dist = 0.7*abs(all['z'][i]-all['z'])/1.00 + 0.4*abs(all['a'][i]-all['a'])/1.00 + 0.17*abs(all['c'][i]-all['c'])/1.00 + 0.62*abs(all['teff'][i]-all['teff'])/100. + 1.5*abs(all['logg'][i]-all['logg'])
    # penalize models where [c/m]-[alpha/m] have opposite sign
    if all['c'][i] - all['a'][i] <= 0. :
        bd = np.where((all['metric'] == 0.) & ((all['c']-all['a']) > 0.))
    else :
        bd = np.where((all['metric'] == 0.) & ((all['c']-all['a']) <= 0.))
    dist[bd] +=  4*0.4*abs(all['a'][i]-all['a'][bd])/1.00 + 4*0.17*abs(all['c'][i]-all['c'][bd])/1.00

    # search only through computed models and only where alpha and m go in opposite directions
    #gd = np.where((all['metric'] == 0.) & ((all['z'][i]-all['z'])*(all['a'][i]-all['a']) <= 0.))
    # search only through computed models 
    gd = np.where(all['metric'] == 0.)
    fill = gd[0][np.argmin(dist[gd])]
    return fill, dist[fill]

def writetmp(tmp,file) :
    """ Auxiliary routine to open file, write string, and close file

    Args:
        temp (str) : string to write
        file (str) : file to write to
    """

    fp = open(file,'w')
    fp.write(tmp)
    fp.close()

def update_header_metallicity(file,z) :
    """ update metallicity in Kurucz atmosphere file 

    Args :
        file (str) : name of file to update
        z (float) : new metallicity
    """
    out=subprocess.check_output(['grep','TITLE',file])
    old_m_h = out.split()[1].strip('[]')
    subprocess.check_output(['sed','-i',
       's/\['+old_m_h+'\]/\['+'{:.2f}'.format(z)+'\]/g',file])

    #update ABUNDANCE SCALE
    hole_scale = 10**z
    if hole_scale < 10 :
       str = '{:7.5f}'.format(10**z)
    else :
       str = '{:8.5f}'.format(10**z)
    tmp=subprocess.check_output(['gawk',
       '/ABUNDANCE SCALE/{gsub ($3,'+str+')};{print}',file])
    writetmp(tmp,file)


def update_abundance_change(linestring,column,delta,file) :
    """ substitute value in given column of line with given string 
    """
    ac_old = subprocess.check_output(['gawk','/'+linestring+'/{print $'+column+'}', file])
    str = '{:.2f}'.format(float(ac_old)+delta)
    tmp = subprocess.check_output(['gawk','/'+linestring+'/{gsub ($'+column+',"'+str+'")};{print}',file])
    writetmp(tmp,file)

def update_header_carbon(file,delta_c) :
    """Update carbon abundances for Kurucz file
    """
    update_abundance_change('ABUNDANCE CHANGE  3','10',delta_c,file)

def update_header_alpha(file,delta_a) :
    """ Update alpha element abundances for Kurucz file
    """
    update_abundance_change('ABUNDANCE CHANGE  3','14',delta_a,file)
    update_abundance_change('ABUNDANCE CHANGE  9','10',delta_a,file)
    update_abundance_change('ABUNDANCE CHANGE  9','14',delta_a,file)
    update_abundance_change('ABUNDANCE CHANGE 15','6',delta_a,file)
    update_abundance_change('ABUNDANCE CHANGE 15','14',delta_a,file)
    update_abundance_change('ABUNDANCE CHANGE 21','6',delta_a,file)

def update_header_h_he(file) :
    """ Update H and He abundances in Kurucz atmosphere file
    """
    h_old = subprocess.check_output(['gawk','/ABUNDANCE SCALE/{print $7}',file])
    he_old = subprocess.check_output(['gawk','/ABUNDANCE SCALE/{print $9}',file])
    abun_scale = subprocess.check_output(['gawk','/ABUNDANCE SCALE/{print $3}',file])

    #calculate new values
    #sum up all 10^abundance changes above He
    sum_ac=subprocess.check_output(['gawk', 'NR==6,NR==22{for (i=4; i<=14;i+=2) if (10^$i < 1) sum+=10^$i};END {print sum}',file])
    sum_el = float(abun_scale) * float(sum_ac)

    h_new = (1 - sum_el)/(1 + float(he_old)/float(h_old))
    he_new = '{:7.5f}'.format(h_new * float(he_old)/float(h_old))
    h_new = '{:7.5f}'.format((1 - sum_el)/(1 + float(he_old)/float(h_old)))

    #update header with new values
    tmp = subprocess.check_output(['gawk','/ABUNDANCE SCALE/{gsub ($7,"'+h_new+'")};{print}',file])
    writetmp(tmp,file)
    tmp = subprocess.check_output(['gawk','/ABUNDANCE SCALE/{gsub ($9,"'+he_new+'")};{print}',file])
    writetmp(tmp,file)

def pars(s) :
    """ routine to parse command-line triples
    """
    try :
        x, y, z = s.split(',')
        return x, y, z
    except :
        raise argparse.ArgumentTypeError("range must be specified as number,start,delta")

def fill_holes(argv) :
    """ Main routine to create "holefile" for grid as specified on command line.
    Creates links for the missing atmospheres, and outputs FITS and ASCII file
    giving the grid and hole filling information.
    """

    # set up arguments, all required!
    parser = argparse.ArgumentParser()
    parser.add_argument("--z",required=True,type=pars)
    parser.add_argument("--teff",required=True,type=pars)
    parser.add_argument("--logg",required=True,type=pars)
    parser.add_argument("--alpha",required=True,type=pars)
    parser.add_argument("--carbon",required=True,type=pars)
    parser.add_argument("--model",required=True)
    parser.add_argument("--dir",required=True)
    parser.add_argument("--fits",required=False)
    args=parser.parse_args()

    # transfer argument names
    model = args.model
    dir = args.dir + '/'

    # set up the grid value arrays based on command-line input
    all_z = float(args.z[1]) + np.arange(int(args.z[0]))*float(args.z[2])
    all_teff = int(args.teff[1]) + np.arange(int(args.teff[0]))*int(args.teff[2])
    all_logg = float(args.logg[1]) + np.arange(int(args.logg[0]))*float(args.logg[2])
    all_a = float(args.alpha[1]) + np.arange(int(args.alpha[0]))*float(args.alpha[2])
    all_c = float(args.carbon[1]) + np.arange(int(args.carbon[0]))*float(args.carbon[2])
    ntot = int(args.z[0]) * int(args.teff[0]) * int(args.logg[0]) * int(args.alpha[0]) * int(args.carbon[0])
    all = np.zeros(ntot, dtype = {'names':['teff','logg','z','c','a','metric'], 
       'formats':['i4','f4','f4','f4','f4','f4']} )

    print('# [M/H] [C/M] [A/M] Teff logg metric Delta(M) Delta(C) Delta(A) Delta(Teff) Delta(logg) (SPECLIB_VERS: '+os.environ['APOGEE_VER']+')')

    # loop over grid finding computed models and filling up array of all models
    nmissing = 0
    i=0
    for teff in all_teff :
       for z in all_z :
          for c in all_c :
             for a in all_a :
                for logg in all_logg :

                    file = filename(teff,logg,z,c,a,model=model)
                    all['teff'][i] = teff
                    all['logg'][i] = logg
                    all['z'][i] = z
                    all['c'][i] = c
                    all['a'][i] = a
                    all['metric'][i] = -1.
                    # see if model exists 
                    if os.path.exists(dir+file) :
                      # set "metric" to zero
                      all['metric'][i] = 0.
                    i += 1

    # now loop through and try to fill the missing models
    for i in range(ntot) :
        file =filename(all['teff'][i],all['logg'][i],all['z'][i],all['c'][i],all['a'][i],model=model)
        if all['metric'][i] < 0 :
            fill, dist = find_filler(all,i,model=model) 
            all['metric'][i] = dist
            file =file+'.filled'
            fillfile = filename(all['teff'][fill],all['logg'][fill],all['z'][fill],all['c'][fill],all['a'][fill],model=model)
            try :
                os.mkdir(os.path.dirname(dir+file))
            except OSError :
                pass
            #shutil.copyfile(dir+fillfile,dir+file)
            try :
                os.remove(dir+file)
            except OSError :
                pass
            os.symlink('../'+dir+fillfile,dir+file)

            # fix up headers of filled models for new abundances
            if model.upper() == 'KURUCZ' :
                #If metallicity has changed, need to update header
                if all['z'][i] != all['z'][fill] :
                  update_header_metallicity(dir+file,all['z'][i])
                #If [C/M] has changed, need to update header
                if all['c'][i] != all['c'][fill] :
                  update_header_carbon(dir+file,all['c'][i]-all['c'][fill])
                #If [A/M] has changed, need to update header
                if all['a'][i] != all['a'][fill] : 
                  update_header_alpha(dir+file,all['a'][i]-all['a'][fill])
                #Update ABUNDANCE CHANGE 1 and 2 in header
                update_header_h_he(dir+file)
        else :
            fill = i

        # summary output for holefile
        print ('{:+.2f} {:+.2f} {:+.2f} {:04d} {:+.1f} {:.2f}'.format(all['z'][i],all['c'][i],all['a'][i],all['teff'][i],all['logg'][i],all['metric'][i])+ 
               '  {:+.2f} {:+.2f} {:+.2f} {:4d} {:+.1f}'.format(all['z'][i]-all['z'][fill],all['c'][i]-all['c'][fill],all['a'][i]-all['a'][fill],all['teff'][i]-all['teff'][fill],all['logg'][i]-all['logg'][fill]) )

    # FITS file output
    out=np.zeros([len(all_c),len(all_a),len(all_z),len(all_logg),len(all_teff)],dtype=np.float32)
    i=0
    for teff in range(len(all_teff)) :
       for z in range(len(all_z)) :
          for c in range(len(all_c)) :
             for a in range(len(all_a)) :
                 for logg in range(len(all_logg)) :
                     out[c,a,z,logg,teff] = all['metric'][i]
                     i+=1
    if args.fits is not None :
        hd=fits.PrimaryHDU(out)
        hd.header['CTYPE1']='TEFF'
        hd.header['CRVAL1']=float(args.teff[1])
        hd.header['CDELT1']=float(args.teff[2])
        hd.header['CTYPE2']='LOGG'
        hd.header['CRVAL2']=float(args.logg[1])
        hd.header['CDELT2']=float(args.logg[2])
        hd.header['CTYPE3']='METALS'
        hd.header['CRVAL3']=float(args.z[1])
        hd.header['CDELT3']=float(args.z[2])
        hd.header['CTYPE4']='O Mg Si S Ca Ti'
        hd.header['CRVAL4']=float(args.alpha[1])
        hd.header['CDELT4']=float(args.alpha[2])
        hd.header['CTYPE5']='C'
        hd.header['CRVAL5']=float(args.carbon[1])
        hd.header['CDELT5']=float(args.carbon[2])
        hd.writeto(args.fits,overwrite=True)


if __name__ == "__main__" :
    fill_holes(sys.argv)
