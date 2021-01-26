# encoding: utf-8
#
# @Author: Jon Holtzman
# @Date: March 2018
# @Filename: synth.py
# @License: BSD 3-Clause
# @Copyright: Jon Holtzman

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import numpy as np
import pdb
from tools import match
from apogee.aspcap import aspcap
from apogee.utils import bitmask
import glob
try: import corner
except: pass
import os
import subprocess
import tempfile

def writespec(name,data) :
    """ Writes FERRE 'spectrum' file with input data, one line per star
    """
    f=open(name,'w')
    for spec in data :
        for pix in np.arange(spec.shape[0]) :
            f.write('{:14.6e}'.format(spec[pix]))
        f.write('\n')
    f.close()
    return 


def writenml(outfile,file,libhead,ncpus=2,nruns=1,inter=3,pca=1,errbar=1,indi=None,indv=None,filterfile=None,f_format=1,f_access=0,f_sort=None,
               init=None,indini=None,renorm=None,obscont=0,rejectcont=0,algor=1,nov=None,stopcr=None,ttie=None) :
    """ Writes FERRE control file
    """
    f=open(outfile,'w')
    f.write(' &LISTA\n')
    ndim=libhead['N_OF_DIM']
    f.write(' NDIM = {:d}\n'.format(ndim))
    if indi is not None : f.write(' INDI = '+np.array2string(np.array(indi)).strip('[]')+'\n')
    if nov is None : nov=ndim
    if indv is None : 
        f.write((' NOV = {:2d}\n').format(nov))
        f.write(' INDV = '+np.array2string(np.arange(1,ndim+1)).strip('[]')+'\n')
    else :
        nov = len(indv)
        f.write((' NOV = {:2d}\n').format(nov))
        f.write(' INDV = '+np.array2string(np.array(indv)).strip('[]')+'\n')
    f.write(" SYNTHFILE(1) = '"+libhead['FILE']+"'\n")
    if filterfile is not None : f.write(" FILTERFILE = '"+filterfile+"'\n")
    f.write(" PFILE = '"+file+".ipf'\n")
    f.write(" OFFILE = '"+file+".mdl'\n")
    if nov > 0 :
        f.write(" ERFILE = '"+file+".err'\n")
        f.write(" OPFILE = '"+file+".spm'\n")
    if renorm is not None :
        f.write(" FFILE = '"+file+".obs'\n")
        f.write(' CONT = 1\n')
        f.write(' NCONT = {:d}\n'.format(renorm))
        f.write(' OBSCONT = {:d}\n'.format(obscont))
        f.write(' REJECTCONT = {:f}\n'.format(rejectcont))
        f.write(" SFFILE = '"+file+".frd'\n")
    elif nov > 0 :
        f.write(" FFILE = '"+file+".frd'\n")
    if nov > 0 :
        f.write(' ERRBAR = {:d}\n'.format(errbar))
        if init is not None: f.write(' INIT = {:d}\n'.format(init))
        f.write(' ALGOR = {:2d}\n'.format(algor))
        if indini is not None :
            f.write(' INDINI = '+np.array2string(np.array(indini)).strip('[]')+'\n')
            nruns = 1
            for term in indini : nruns = nruns * term
        f.write(' NRUNS = {:2d}\n'.format(nruns))
    if stopcr is not None : f.write(' STOPCR = {:f}\n'.format(stopcr))
    if ttie is not None :
        j=np.where(np.array(ttie) >0)[0]
        f.write(' NTIE = {:2d}\n'.format(len(j)))
        f.write(' TYPETIE = 1\n')
        for i,tie in enumerate(np.array(ttie)[j]) :
            f.write(' INDTIE({:d}) = {:d}\n'.format(i+1,tie))
            f.write(' TTIE0({:d}) = 0.\n'.format(i+1,tie))
            f.write(' TTIE({:d},{:d}) = -1.\n'.format(i+1,indv[0],tie))

    f.write(' NTHREADS = {:2d}\n'.format(ncpus))
    f.write(' COVPRINT = 1\n')
    f.write(' PCAPROJECT = 0\n')
    f.write(' PCACHI = 0\n')
    f.write(' INTER = {:d}\n'.format(inter))
    f.write(' F_FORMAT = {:d}\n'.format(f_format))
    f.write(' F_ACCESS = {:d}\n'.format(f_access))
    if f_sort is not None : f.write(' F_SORT = {:d}\n'.format(f_sort))
    f.write(' /\n')
    f.close()

def clip(x,lim,eps=None) :
    """ Utility routine to clip values within limits, and move slightly off edges if requested
    """
    # set negative zero to zero
    if np.isclose(x,0.) : x=0.
    # clip to limits
    tmp=np.max([lim[0],np.min([lim[1],x])])
    # move off limit if requested
    if eps is not None :
        if np.isclose(tmp,lim[0],atol=0.001) : tmp+=eps
        if np.isclose(tmp,lim[1],atol=0.001) : tmp-=eps
    return tmp


def writeipf(name,libfile,stars,param=None) :
    """ Writes FERRE input file
    """
    # get library headers and load wavelength array
    libhead0, libhead=rdlibhead(libfile)

    params=aspcap.params()[0]
    nparams=libhead0['N_OF_DIM']
    # get the index numbers in input parameter array for correct library parameter order
    index=np.zeros(nparams,dtype=int)
    for i in range(nparams) :
        index[i] = np.where(params == libhead0['LABEL'][i].decode())[0]
    # if input parameters aren't specified, use zeros
    if param is None :
        param=np.zeros([len(stars),len(params)])
    # write the IPF file
    f=open(name+'.ipf','w')
    for i,star in enumerate(stars) :
        f.write('{:<40s}'.format(star))
        for ipar in range(nparams) : 
            lims = [libhead0['LLIMITS'][ipar],libhead0['LLIMITS'][ipar]+libhead0['STEPS'][ipar]*(libhead0['N_P'][ipar]-1)]
            f.write('{:12.3f}'.format(clip(param[i][index[ipar]],lims,eps=0.001)))
        f.write('\n')
    f.close()
        
def readmcmc(name,libfile,nov=None,burn=500) :
    libhead0, libhead=rdlibhead(libfile)

    files=glob.glob(name+'.chain*.dat')
    print(files)
    alldat=[]
    for file in files :
        a=np.loadtxt(file,skiprows=1)
        alldat.extend(a[burn:,:])
    alldat=np.array(alldat)
    if nov is None : nov = np.arange(libhead0['N_OF_DIM'])+1
    pdb.set_trace()
    # transform from normalized to true parameters
    labels=[]
    for i,ipar in enumerate(nov) :
        alldat[:,i+2] = libhead0['LLIMITS'][ipar-1] + alldat[:,i+2]*libhead0['STEPS'][ipar-1]*(libhead0['N_P'][ipar-1]-1)
        labels.append(libhead0['LABEL'][ipar-1])
    corner.corner(alldat[:,2:],labels=labels,show_titles=True)

def read(name,libfile) :

    """ Read all of the FERRE files associated with a FERRE run
    """
    # get library headers and load wavelength array
    libhead0, libhead=rdlibhead(libfile)
    wave=[]
    for ichip in range(len(libhead)) :
        wave.extend(libhead[ichip]['WAVE'][0] + np.arange(libhead[ichip]['NPIX'])*libhead[ichip]['WAVE'][1])
    wave=10.**np.array(wave)
    nwave=len(wave)
    nparam=libhead0['N_OF_DIM']

    # input ipf and output spm files, match object names
    ipfobj,ipf=readferredata(name+'.ipf')
    spmobj,spm=readferredata(name+'.spm')
    i1,i2=match.match(ipfobj,spmobj)
    nobj=len(ipfobj)
    param=spm[i2,0:nparam]
    paramerr=spm[i2,nparam:2*nparam]
    chi2=10.**spm[i2,2*nparam+2]
    covar=spm[i2,2*nparam+3:2*nparam+3+nparam**2]
    covar=np.reshape(covar,(nobj,nparam,nparam))

    # load param array
    params,tagnames,flagnames=aspcap.params()
    ntotparams=len(params)
    index=np.zeros(ntotparams,dtype=int)
    a=np.zeros(nobj, dtype=[('APOGEE_ID','S100'),
                              ('FPARAM','f4',(ntotparams)),
                              ('FPARAM_COV','f4',(ntotparams,ntotparams)),
                              ('ASPCAP_CHI2','f4'),
                              ('PARAMFLAG','i4',(ntotparams)),
                              ('ASPCAPFLAG','i4')])
    a['APOGEE_ID']=ipfobj[i1]
    a['ASPCAP_CHI2']=chi2

    parammask=bitmask.ParamBitMask()
    aspcapmask=bitmask.AspcapBitMask()
    for i in range(nparam) :
        # load FPARAM and FPARAM_COV into the right order of parameters
        pname = libhead0['LABEL'][i].decode()
        index = np.where(params == pname)[0][0]
        val = param[:,i]
        a['FPARAM'][:,index] = val
        for j in range(nparam) :
            jindex = np.where(params == libhead0['LABEL'][j].decode())[0][0]
            a['FPARAM_COV'][:,index,jindex]=covar[:,i,j]

        # check for grid edge and flag
        bad = np.where((val < libhead0['LLIMITS'][i]+libhead0['STEPS'][i]/8.) )[0]
        if pname != 'N' and pname !='LOG10VDOP' and pname != 'LGVSINI' :
            a['PARAMFLAG'][bad,index] |= parammask.getval('GRIDEDGE_BAD')
            a['ASPCAPFLAG'][bad] |= aspcapmask.getval(flagnames[index]+'_BAD')
        bad = np.where((val > libhead0['LLIMITS'][i]+libhead0['STEPS'][i]*(libhead0['N_P'][i]-1-1/8.)) )[0]
        if pname != 'N' :
            a['PARAMFLAG'][bad,index] |= parammask.getval('GRIDEDGE_BAD')
            a['ASPCAPFLAG'][bad] |= aspcapmask.getval(flagnames[index]+'_BAD')
        warn = np.where((a['PARAMFLAG'][:,index]&parammask.getval('GRIDEDGE_BAD') ==0) &
                        (val < libhead0['LLIMITS'][i]+libhead0['STEPS'][i]) |
                        (val > libhead0['LLIMITS'][i]+libhead0['STEPS'][i]*(libhead0['N_P'][i]-1-1.)) )[0]
        a['PARAMFLAG'][warn,index] |= parammask.getval('GRIDEDGE_WARN')
        a['ASPCAPFLAG'][warn] |= aspcapmask.getval(flagnames[index]+'_WARN')
        bad = np.where(val < -999)[0]
        a['PARAMFLAG'][bad,index] |= parammask.getval('FERRE_FAIL')
        a['ASPCAPFLAG'][bad] |= aspcapmask.getval(flagnames[index]+'_BAD')
        a['ASPCAPFLAG'][bad] |= aspcapmask.getval('FERRE_FAIL')


    # put spectral data into a structured array
    form='{:d}f4'.format(nwave)
    sform='{:d}f4'.format(spm.shape[1])
    out=np.empty(nobj, dtype=[('obj','S24'),('spm',sform),('obs',form),('frd',form),('err',form),('mdl',form),('chi2',form)])
    out['obj']=ipfobj
    out['spm'][i1,:]=spm[i2,:]
    out['obs'][i1,:]=readspec(name+'.obs')[i1,:]
    out['frd'][i1,:]=readspec(name+'.frd')[i2,:]
    out['err'][i1,:]=readspec(name+'.err')[i1,:]
    out['mdl'][i1,:]=readspec(name+'.mdl')[i2,:]
    out['chi2'][i1,:]=(out['obs']-out['mdl'])**2/out['err']**2

    return a,out,wave

def readspec(name) :
    """ Read a single file with FERRE-format spectra, and return as 2D array [nspec,nwave]
    """
    f=open(name)
    data=[]
    for line in f :
        spec=np.array(line.split())
        spec=spec.astype(float)
        data.append(spec)
    return np.array(data)

def readmask(name) :
    """ Read a single FERRE mask file
    """
    f=open(name)
    data=[]
    for line in f :
        data.append(float(line))
    return np.array(data)

def writemask(name, mask) :
    '''
    Write a single FERRE ask file
    '''
    f=open(name,'w')
    for m in mask:
        f.write('{:8.3f}\n'.format(m))
    f.close()

def readferredata(name) :
    """ Read a single file with FERRE-format data, and return as 2D array [nspec,nwave]
    """
    f=open(name)
    alldata=[]
    allobj=[]
    for line in f :
        a=line.split()
        obj=np.array(a[0])
        data=np.array(a[1:len(a)]).astype(float)
        alldata.append(data)
        allobj.append(obj)
    return np.array(allobj),np.array(alldata)


def rdsinglehead(f) :
    """ Read a FERRE library header into a dictionary
    """
    dict={}
    for line in f:  
      # read until we hit a '/' in first character
      if line[1] != '/' and line[0] != '/' : 
        words=line.split()
        nwords=len(words)
        card=words[0].upper()
        if card[-1]=='=' : 
            # adjust if = is not separated from card name
            card = card[0:-1]
            nwords+=1
        if nwords>2 :
            # we have a value, set it to characters to right of = sign, not including newling
            value=line[line.find('=')+1:len(line)-1]
            if value.find("'") >= 0 :
                # we have a string, strip it
                val=value.replace("'"," ").strip()
                if card.find('LABEL') >= 0 :
                    ilab= int(card[card.find('(')+1:card.find(')')])
                    label[ilab-1] = val
            else :
                # we have numbers, put into array of correct type
                vals=value.split()
                n=len(vals)
                try:
                    val=np.array(vals).astype(int)
                except :
                    val=np.array(vals).astype(float)
                # if we just have one variable, we don't want an array
                if n == 1: val=val[0]
            if card == 'N_OF_DIM' : label=np.zeros(val,dtype='S24')
            if card.find('LABEL') < 0 :
                # add card, value to dictionary if not a LABEL() card
                dict[card]=val
      else:
          break
    try:
        dict['LABEL'] = label
    except:
        pass
    return dict

def rdlibhead(name) :
    """ Read a full FERRE library header with multi-extensions

    Returns:
       libstr0, libstr : first header, then list of extension headers; headers returned as dictionaries
    """
    try:
        f=open(name,'r')
    except:
        print('cannot open file',name)
        return
    libstr0=rdsinglehead(f)
    libstr0['FILE'] = name
    try:
        multi=libstr0['MULTI']
        libstr=[]
        for imulti in range(multi) :
            libstr.append(rdsinglehead(f))
        return libstr0,libstr
    except:
        return libstr0

def elemmask(libfile,wspec,spec,dw=0.1,thresh=0.,out=None) :
    """ writes FERRE binary (0/1) mask file for wavelengths from input libfile, where input spectrum > thresh within wavelength dw
    """
    libhead0,libhead=rdlibhead(libfile)
    mask = []
    for chip in range(libhead0['MULTI']) :
        wave=libhead[chip]['WAVE'][0]+np.arange(libhead[chip]['NPIX'])*libhead[chip]['WAVE'][1]
        if libhead[chip]['LOGW'] == 1 : wave = 10.**wave
        for w in wave :
            j=np.where(abs(wspec-w) < dw)[0]
            if len(j) > 0 :
                if spec[j].max() > thresh :
                    mask.append(0.)
                else :
                    mask.append(1.)
    if out is not None :
        writemask(out,mask)

    return np.array(mask)
            

def wrhead(planstr,file,npca=None,npix=None,wchip=None,cont=None) :
    """ Write header of library file given plan structure
    """
    ndim=0
    name=[]
    llimits=[]
    steps=[]
    n=[]
    idim = 0
    for dim in ['oa','vt','cm','nm','am','rot','mh','logg','teff'] :
        if int(planstr['n'+dim]) > 1 : 
            ndim+=1
            n.append(int(planstr['n'+dim]))
            llimits.append(float(planstr[dim+'0']))
            steps.append(float(planstr['d'+dim]))
            if dim == 'oa' : name.append('O')
            if dim == 'vt' : name.append('LOG10VDOP')
            if dim == 'cm' : name.append('C')
            if dim == 'nm' : name.append('N')
            if dim == 'am' : name.append('O Mg Si S Ca Ti')
            if dim == 'rot' : name.append('LGVSINI')
            if dim == 'mh' : name.append('METALS')
            if dim == 'logg' : name.append('LOGG')
            if dim == 'teff' : name.append('TEFF')
    fp=open(file,'w')
    fp.write(" &SYNTH\n")
    fp.write(" ID = '"+file+"'\n")
    if wchip is not None : 
        nchips = len(wchip)
        fp.write(" MULTI = {:d}\n".format(nchips))
    if npca is not None:
        fp.write(" NPCA = ".format(ndim))
        for item in npca : fp.write('{:4d}'.format(item))
        fp.write("\n")
    if npix is not None: fp.write(" NPIX = {:d}\n".format(npix))
    fp.write(" N_OF_DIM = {:d}\n".format(ndim))
    fp.write(" N_P = ")
    for item in n : fp.write('{:4d}'.format(item))
    for i,item in enumerate(name) : fp.write("\n LABEL({:1d}) = '{:s}'".format(i+1,item))
    fp.write("\n LLIMITS = ")
    for item in llimits : fp.write('{:11.5f}'.format(item))
    fp.write("\n STEPS = ")
    for item in steps : fp.write('{:11.5f}'.format(item))
    fp.write("\n /\n")
    if wchip is not None :
        for chip,c in zip(wchip,cont) :
            fp.write(" &SYNTH\n")
            fp.write(" ID = '"+file+"'\n")
            fp.write(" N_OF_DIM = {:d}\n".format(ndim))
            fp.write(" N_P = ")
            for item in n : fp.write('{:4d}'.format(item))
            for i,item in enumerate(name) : fp.write("\n LABEL({:1d}) = '{:s}'".format(i+1,item))
            fp.write("\n LLIMITS = ")
            for item in llimits : fp.write('{:11.5f}'.format(item))
            fp.write("\n STEPS = ")
            for item in steps : fp.write('{:11.5f}'.format(item))
            fp.write("\n")
            fp.write(" NPIX = {:d}\n".format(chip[0]))
            fp.write(" WAVE = {:16.7f}{:16.7e}\n".format(chip[1],chip[2]))
            fp.write(" LOGW = {:d}\n".format(1))
            fp.write(" VACUUM = {:d}\n".format(1))
            fp.write(" RESOLUTION = {:d}\n".format(22500))
            fp.write(" CONTINUUM = {:6d}{:6d}{:7.2f}{:7.2f}\n".format(int(round(c[0])),int(round(c[1])),c[2],c[3]))
            fp.write(" /\n")
    fp.close()

def interp(libfile,params,renorm=4,obscont=0,rejectcont=0.3) :
    """ Use FERRE to get an interpolated spectrum
    """
    tmpname = tempfile.NamedTemporaryFile().name
    libhead0,libhead=rdlibhead(libfile)
    writenml(tmpname+'.nml',tmpname,libhead0,ncpus=1,nruns=1,inter=3,pca=1,errbar=1,indi=None,indv=None,filterfile=None,f_format=1,f_access=1,f_sort=None,
               init=None,indini=None,renorm=renorm,obscont=obscont,rejectcont=rejectcont,algor=1,nov=0,stopcr=None,ttie=None) 
    writeipf(tmpname,libfile,['dummy'],param=[params]) 
    ret = subprocess.call(['ferre.x',tmpname+'.nml'],shell=False)
    mdl = readspec(tmpname+'.mdl')
    for ext in ['nml','ipf','mdl'] : os.remove(tmpname+'.'+ext)

    return mdl
