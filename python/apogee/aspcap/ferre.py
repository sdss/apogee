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
#from apogee.tools import match
from apogee.aspcap import aspcap

def writespec(name,data) :
    """ Writes FERRE 'spectrum' file with input data, one line per star
    """
    f=open(name,'w')
    for spec in data :
        for pix in np.arange(spec.shape[0]) :
            f.write('{:12.6f}'.format(spec[pix]))
        f.write('\n')
    f.close()
    return 


def writenml(outfile,file,libhead,ncpus=2,nruns=1,interord=3,direct=1,pca=1,errbar=1,indi=None,indv=None,filterfile=None,f_format=1,f_access=0,
               init=None,indini=None,renorm=None) :
    """ Writes FERRE control file
    """

    f=open(outfile,'w')
    f.write(' &LISTA\n')
    ndim=libhead['N_OF_DIM']
    f.write(' NDIM = {:d}\n'.format(ndim))
    if indi is not None : f.write((' INDI = '+'{:2d}'*ndim+'\n').format(indi))
    if indv is not None : 
        f.write((' NOV = {:2d}\n').format(len(indv)))
        if len(indv) > 0 : f.write(' INDV = '+np.array2string(np.array(indv)).strip('[]')+'\n')
    else :
        f.write((' NOV = {:2d}\n').format(ndim))
        f.write(' INDV = '+np.array2string(np.arange(1,ndim+1)).strip('[]')+'\n')
    f.write(" SYNTHFILE(1) = '"+libhead['FILE']+"'\n")
    if filterfile is not None : f.write(' FILTERFILE = '+filterfile+'\n')
    f.write(" PFILE = '"+file+".ipf'\n")
    f.write(" ERFILE = '"+file+".err'\n")
    f.write(" OPFILE = '"+file+".spm'\n")
    f.write(" OFFILE = '"+file+".mdl'\n")
    f.write(' ERRBAR = {:d}\n'.format(errbar))
    if renorm is not None :
        f.write(' CONT = 1\n')
        f.write(' NCONT = {:d}\n'.format(renorm))
        f.write(' OBSCONT = 0\n')
        f.write(" FFILE = '"+file+".obs'\n")
        f.write(" SFFILE = '"+file+".frd'\n")
    else :
        f.write(" FFILE = '"+file+".frd'\n")
    if init is not None: f.write(' INIT = {:d}\n'.format(init))
    if indini is not None :
        f.write(' INDINI = '+np.array2string(np.array(indini)).strip('[]')+'\n')
        nruns = 1
        for term in indini : nruns = nruns * term
    f.write(' NRUNS = {:2d}\n'.format(nruns))
    f.write(' NTHREADS = {:2d}\n'.format(ncpus))
    f.write(' PCAPROJECT = 0\n')
    f.write(' PCACHI = 0\n')
    f.write(' INTER = {:d}\n'.format(interord))
    f.write(' F_FORMAT = {:d}\n'.format(f_format))
    f.write(' F_ACCESS = {:d}\n'.format(f_access))
    f.write(' /\n')
    f.close()

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
        index[i] = np.where(params == libhead0['LABEL'][i])[0]
    # if input parameters aren't specified, use zeros
    if param is None :
        param=np.zeros([len(stars),len(params)])
    # write the IPF file
    f=open(name+'.ipf','w')
    for i,star in enumerate(stars) :
        f.write('{:<40s}'.format(star))
        for ipar in range(nparams) : 
            f.write('{:12.3f}'.format(param[i][index[ipar]]))
        f.write('\n')
    f.close()
        

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
    param=spm[:,0:nparam]
    paramerr=spm[:,nparam:2*nparam]
    chi2=10.**spm[:,2*nparam+2]
    covar=spm[:,2*nparam+3:2*nparam+3+nparam**2]
    covar=np.reshape(covar,(nobj,nparam,nparam))

    # load param array
    params=aspcap.params()[0]
    ntotparams=len(params)
    index=np.zeros(ntotparams,dtype=int)
    a=np.zeros(nobj, dtype=[('APOGEE_ID','S100'),
                              ('FPARAM','f4',(ntotparams)),
                              ('FPARAM_COV','f4',(ntotparams,ntotparams)),
                              ('PARAM_CHI2','f4')])
    a['APOGEE_ID']=ipfobj
    for i in range(ntotparams) :
        try :
            index[i] = np.where(libhead0['LABEL'] == params[i])[0]
            a['FPARAM'][:,i] = param[:,index[i]]
            for j in range(ntotparams) :
                a['FPARAM_COV'][:,i,j]=covar[:,index[i],index[j]]
        except :
            index[i] = -1
    a['PARAM_CHI2']=chi2

    # put it all into a structured array
    form='{:d}f4'.format(nwave)
    sform='{:d}f4'.format(spm.shape[1])
    out=np.empty(nobj, dtype=[('obj','S24'),('spm',sform),('obs',form),('err',form),('mdl',form),('chi2',form)])
    out['obj']=ipfobj
    out['spm'][i1,:]=spm[i2,:]
    out['obs']=readspec(name+'.frd')[i2,:]
    out['err']=readspec(name+'.err')[i2,:]
    out['mdl']=readspec(name+'.mdl')[i2,:]
    out['chi2']=(out['obs']-out['mdl'])**2/out['err']**2

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
      if line[1] != '/' : 
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
            

#def plotspec(w,spec,n=0) :
#    fig,ax=plots.multi(1,2)
#    plots.plotl(ax[0],w,spec['obs'][n,:],yr=[0,1.3]) 
#    plots.plotl(ax[0],w,spec['err'][n,:]) 
#    plots.plotl(ax[0],w,spec['mdl'][n,:]) 
#    plots.plotl(ax[1],w,spec['chi2'][n,:])



