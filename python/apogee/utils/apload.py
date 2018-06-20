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

from astropy.io import fits
import os
from sdss_access.path import path
from sdss_access.sync.http import HttpAccess
import pdb
import sys

apred = 'r6'
apstar = 'stars'
aspcap = 'l30e'
results = 'l30e.2'
instrument = 'apogee-n'

if os.getenv('APOGEE_REDUX') is None :
    print('you must set environment variable APOGEE_REDUX!')
print('\nusing defaults: apred = ', apred, 'apstar = ', apstar, 'aspcap = ', aspcap, 'results = ',results, 'instrument = ',instrument)
print("you can change using, e.g. apload.apred='r6'\n")

def dr10() :
    global apred, aspcap, results
    apred='r3'
    apstar='s3'
    aspcap='v304'
    results='v304'
    print(apred)

def dr12() :
    global apred, aspcap, results
    apred='r5'
    aspcap='l25_6d'
    results='v603'
    print(apred)

def dr13() :
    global apred, aspcap, results
    apred='r6'
    aspcap='l30e'
    results='l30e.2'

def dr14() :
    global apred, aspcap, results
    apred='r8'
    aspcap='l31c'
    results='l31c.2'

def printerror() :
    print('cannot find file: do you have correct version? permission? wget authentication?')

def allStar(hdu=None) :
    ''' Read allStar file (downloading if necesssary)'''
    try :
        file = allfile('allStar',results=results,apred=apred,apstar=apstar,aspcap=aspcap)
        return _readhdu(file,hdu=hdu)
    except :
        printerror()

def allVisit(hdu=None) :
    ''' Read allVisit file (downloading if necesssary)'''
    try :
        file = allfile('allVisit',results=results,apred=apred,apstar=apstar,aspcap=aspcap)
        return _readhdu(file,hdu=hdu)
    except :
        printerror()

def allPlates(hdu=None) :
    ''' Read allPlates file (downloading if necesssary)'''
    try :
        file = allfile('allPlates',results=results,apred=apred,apstar=apstar,aspcap=aspcap)
        return _readhdu(file,hdu=hdu)
    except :
        printerror()

def allExp(hdu=None) :
    ''' Read allExp file (downloading if necesssary)'''
    try :
        file = allfile('allExp',results=results,apred=apred)
        return _readhdu(file,hdu=hdu)
    except :
        printerror()

def allSci(hdu=None) :
    ''' Read allSci file (downloading if necesssary)'''
    try :
        file = allfile('allSci',results=results,apred=apred)
        return _readhdu(file,hdu=hdu)
    except :
        printerror()

def allCal(hdu=None) :
    ''' Read allCal file (downloading if necesssary)'''
    try :
        file = allfile('allCal',results=results,apred=apred)
        return _readhdu(file,hdu=hdu)
    except :
        printerror()

def apFlat(*args,**kwargs) :
    """
    NAME: apload.apFlat
    PURPOSE:  read apFlat file (downloading if necessary)
    USAGE:  ret = apload.apFlat(imagenumber[,hdu=N,tuple=True])
    RETURNS: if hdu==None : dictionary of ImageHDUs (all extensions) 
                            for chips 'a', 'b', 'c'
             if hdu=N : returns dictionaries (data, header) for specified HDU
             if tuple=True : returns tuples rather than dictionaries
    """
    if len(args) != 1 :
        print('Usage: apFlat(imagenumber)')
    else :
        try :
            file = allfile(
               'Flat',num=args[0],mjd=cmjd(args[0]),chips=True,
               apred=apred,apstar=apstar,aspcap=aspcap,results=results,dr='collab')
            return _readchip(file,'Flat',**kwargs)
        except :
            printerror()

def apWave(*args,**kwargs) :
    """
    NAME: apload.apWave
    PURPOSE:  read apWave file (downloading if necessary)
    USAGE:  ret = apload.apWave(imagenumber[,hdu=N,tuple=True])
    RETURNS: if hdu==None : dictionary of ImageHDUs (all extensions) 
                            for chips 'a', 'b', 'c'
             if hdu=N : returns dictionaries (data, header) for specified HDU
             if tuple=True : returns tuples rather than dictionaries
    """
    if len(args) != 1 :
        print('Usage: apWave(imagenumber)')
    else :
        try :
            file = allfile(
               'Wave',num=args[0],mjd=cmjd(args[0]),chips=True,
               apred=apred,apstar=apstar,aspcap=aspcap,results=results,dr='collab')
            return _readchip(file,'Wave',**kwargs)
        except :
            printerror()

def apLSF(*args,**kwargs) :
    """
    NAME: apload.apLSF
    PURPOSE:  read apLSF file (downloading if necessary)
    USAGE:  ret = apload.apLSF(imagenumber[,hdu=N,tuple=True])
    RETURNS: if hdu==None : dictionary of ImageHDUs (all extensions) 
                            for chips 'a', 'b', 'c'
             if hdu=N : returns dictionaries (data, header) for specified HDU
             if tuple=True : returns tuples rather than dictionaries
    """
    if len(args) != 1 :
        print('Usage: apLSF(imagenumber)')
    else :
        try :
            file = allfile(
               'LSF',num=args[0],mjd=cmjd(args[0]),chips=True,
               apred=apred,apstar=apstar,aspcap=aspcap,results=results,dr='collab')
            return _readchip(file,'LSF',**kwargs)
        except :
            printerror()

def apPSF(*args,**kwargs) :
    """
    NAME: apload.apPSF
    PURPOSE:  read apPSF file (downloading if necessary)
    USAGE:  ret = apload.apPSF(imagenumber[,hdu=N,tuple=True])
    RETURNS: if hdu==None : dictionary of ImageHDUs (all extensions) 
                            for chips 'a', 'b', 'c'
             if hdu=N : returns dictionaries (data, header) for specified HDU
             if tuple=True : returns tuples rather than dictionaries
    """
    if len(args) != 1 :
        print('Usage: apPSF(imagenumber)')
    else :
        try :
            file = allfile(
               'PSF',num=args[0],mjd=cmjd(args[0]),chips=True,
               apred=apred,apstar=apstar,aspcap=aspcap,results=results,dr='collab')
            return _readchip(file,'PSF',**kwargs)
        except :
            printerror()

def apEPSF(*args,**kwargs) :
    """
    NAME: apload.apEPSF
    PURPOSE:  read apEPSF file (downloading if necessary)
    USAGE:  ret = apload.apEPSF(imagenumber[,hdu=N,tuple=True])
    RETURNS: if hdu==None : dictionary of ImageHDUs (all extensions) 
                            for chips 'a', 'b', 'c'
             if hdu=N : returns dictionaries (data, header) for specified HDU
             if tuple=True : returns tuples rather than dictionaries
    """
    if len(args) != 1 :
        print('Usage: apEPSF(imagenumber)')
    else :
        try :
            file = allfile(
               'EPSF',num=args[0],mjd=cmjd(args[0]),chips=True,
               apred=apred,apstar=apstar,aspcap=aspcap,results=results,dr='collab')
            return _readchip(file,'EPSF',**kwargs)
        except :
            printerror()

def ap1D(*args,**kwargs) :
    """
    NAME: apload.ap1D
    PURPOSE:  read ap1D file (downloading if necessary)
    USAGE:  ret = apload.ap1D(imagenumber[,hdu=N,tuple=True])
    RETURNS: if hdu==None : dictionary of ImageHDUs (all extensions) 
                            for chips 'a', 'b', 'c'
             if hdu=N : returns dictionaries (data, header) for specified HDU
             if tuple=True : returns tuples rather than dictionaries
    """
    if len(args) != 1 :
        print('Usage: ap1D(imagenumber)')
    else :
        try :
            file = allfile(
               '1D',num=args[0],mjd=cmjd(args[0]),chips=True,
               apred=apred,apstar=apstar,aspcap=aspcap,results=results,dr='collab')
            return _readchip(file,'1D',**kwargs)
        except :
            printerror()


def ap2D(*args,**kwargs) :
    """
    NAME: apload.ap2D
    PURPOSE:  read ap2D file (downloading if necessary)
    USAGE:  ret = apload.ap2D(imagenumber)
    RETURNS: if hdu==None : dictionary of ImageHDUs (all extensions) 
                            for chips 'a', 'b', 'c'
             if hdu=N : returns dictionaries (data, header) for specified HDU
             if tuple=True : returns tuples rather than dictionaries
    """
    fz=''
    for key in kwargs : 
        if key == 'fz' : fz='fz'
    if len(args) != 1 :
        print('Usage: ap2D(imagenumber)')
    else :
        try :
            file = allfile(
               '2D'+fz,num=args[0],mjd=cmjd(args[0]),chips=True,
               apred=apred,apstar=apstar,aspcap=aspcap,results=results,dr='collab')
            print('file: ', file)
            return _readchip(file,'2D',**kwargs)
        except :
            printerror()

def ap2Dmodel(*args,**kwargs) :
    """
    NAME: apload.ap2Dmodel
    PURPOSE:  read ap2Dmodel file (downloading if necessary)
    USAGE:  ret = apload.ap2Dmodel(imagenumber)
    RETURNS: if hdu==None : dictionary of ImageHDUs (all extensions) 
                            for chips 'a', 'b', 'c'
             if hdu=N : returns dictionaries (data, header) for specified HDU
             if tuple=True : returns tuples rather than dictionaries
    """
    if len(args) != 1 :
        print('Usage: ap2Dmodel(imagenumber)')
    else :
        try :
            file = allfile(
               '2Dmodel',num=args[0],mjd=cmjd(args[0]),chips=True,
               apred=apred,apstar=apstar,aspcap=aspcap,results=results,dr='collab')
            return _readchip(file,'2Dmodel',**kwargs)
        except :
            printerror()

def apCframe(*args, **kwargs) :
    """
    NAME: apload.apCframe
    PURPOSE:  read apCframe file (downloading if necessary)
    USAGE:  ret = apload.apCframe(plate,mjd,imagenumber[,hdu=N,tuple=True])
    RETURNS: if hdu==None : dictionary of ImageHDUs (all extensions) 
                            for chips 'a', 'b', 'c'
             if hdu=N : returns dictionaries (data, header) for specified HDU
             if tuple=True : returns tuples rather than dictionaries
    """
    if len(args) != 4 :
        print('Usage: apCframe(field,plate,mjd,imagenumber)')
    else :
        try :
            file = allfile(
               'Cframe',field=args[0],plate=args[1],mjd=args[2],num=args[3],chips=True,
               apred=apred,apstar=apstar,aspcap=aspcap,results=results,dr='collab')
            return _readchip(file,'Cframe',**kwargs)
        except :
            printerror()

def apPlate(*args, **kwargs) :
    """
    NAME: apload.apPlate
    PURPOSE:  read apPlate file (downloading if necessary)
    USAGE:  ret = apload.ap2D(plate,mjd[,hdu=N,tuple=True])
    RETURNS: if hdu==None : dictionary of ImageHDUs (all extensions) 
                            for chips 'a', 'b', 'c'
             if hdu=N : returns dictionaries (data, header) for specified HDU
             if tuple=True : returns tuples rather than dictionaries
    """
    if len(args) != 2 :
        print('Usage: apPlate(plate,mjd)')
    else :
        try :
            file = allfile(
               'Plate',plate=args[0],mjd=args[1],chips=True,telescope='apo25m',
               apred=apred,apstar=apstar,aspcap=aspcap,results=results,dr='collab')
            return _readchip(file,'Plate',**kwargs)
        except :
            printerror()

def apVisit(*args, **kwargs) :
    """
    NAME: apload.apVisit
    PURPOSE:  read apVisit file (downloading if necessary)
    USAGE:  ret = apload.apVisit(plate,mjd,fiber,[hdu=N])
    RETURNS: if hdu==None : ImageHDUs (all extensions)
             if hdu=N : returns (data, header) for specified HDU
    """
    if len(args) != 3 :
        print('Usage: apVisit(plate,mjd,fiber)')
    else :
        try :
            file = allfile(
               'Visit',plate=args[0],mjd=args[1],fiber=args[2],
               apred=apred,apstar=apstar,aspcap=aspcap,results=results,dr='collab')
            return _readhdu(file,**kwargs)
        except :
            printerror()

def apVisit1m(*args, **kwargs) :
    """
    NAME: apload.apVisit
    PURPOSE:  read apVisit file (downloading if necessary)
    USAGE:  ret = apload.apVisit(program,mjd,object,[hdu=N])
    RETURNS: if hdu==None : ImageHDUs (all extensions)
             if hdu=N : returns (data, header) for specified HDU
    """
    if len(args) != 3 :
        print('Usage: apVisit1m(program,mjd,object)')
    else :
        try :
            file = allfile(
               'Visit1m',plate=args[0],mjd=args[1],obj=args[2],telescope='apo1m',
               apred=apred,apstar=apstar,aspcap=aspcap,results=results,dr='collab')
            return _readhdu(file,**kwargs)
        except :
            printerror()

def apVisitSum(*args, **kwargs) :
    """
    NAME: apload.apVisitSum
    PURPOSE:  read apVisitSum file (downloading if necessary)
    USAGE:  ret = apload.apVisitSum(location,plate,mjd)
    RETURNS: if hdu==None : ImageHDUs (all extensions)
             if hdu=N : returns (data, header) for specified HDU
    """
    if len(args) != 3 :
        print('Usage: apVisitSum(location,plate,mjd)')
    else :
        try :
            file = allfile(
               'VisitSum',location=args[0],plate=args[1],mjd=args[2],
               apred=apred,apstar=apstar,aspcap=aspcap,results=results,dr='collab')
            return _readhdu(file,**kwargs)
        except :
            printerror()

def apStar(*args, **kwargs) :
    """
    NAME: apload.apStar
    PURPOSE:  read apStar file (downloading if necessary)
    USAGE:  ret = apload.apStar(location,object)
    RETURNS: if hdu==None : ImageHDUs (all extensions)
             if hdu=N : returns (data, header) for specified HDU
    """
    if len(args) != 2 :
        print('Usage: apStar(location,object)')
    else :
        try :
            file = allfile(
               'Star',location=args[0],obj=args[1],
               apred=apred,apstar=apstar,aspcap=aspcap,results=results,dr='collab')
            return _readhdu(file,**kwargs)
        except :
            printerror()

def apStar1m(*args, **kwargs) :
    """
    NAME: apload.apStar1m
    PURPOSE:  read apStar file (downloading if necessary)
    USAGE:  ret = apload.apStar1m(location,object)
    RETURNS: if hdu==None : ImageHDUs (all extensions)
             if hdu=N : returns (data, header) for specified HDU
    """
    if len(args) != 2 :
        print('Usage: apStar(location,object)')
    else :
        try :
            file = allfile(
               'Star',location=args[0],obj=args[1],telescope='apo1m',
               apred=apred,apstar=apstar,aspcap=aspcap,results=results,dr='collab')
            return _readhdu(file,**kwargs)
        except :
            printerror()

def aspcapStar(*args, **kwargs) :
    """
    NAME: apload.aspcapStar
    PURPOSE:  read aspcapStar file (downloading if necessary)
    USAGE:  ret = apload.aspcapStar(location,object)
    RETURNS: if hdu==None : ImageHDUs (all extensions)
             if hdu=N : returns (data, header) for specified HDU
    """
    if len(args) != 2 :
        print('Usage: aspcapStar(location,object)')
    else :
        try :
            file = allfile(
               'aspcapStar',location=args[0],obj=args[1],
               apred=apred,apstar=apstar,aspcap=aspcap,results=results,dr='collab')
            return _readhdu(file,**kwargs)
        except :
            printerror()

def apField(*args, **kwargs) :
    """
    NAME: apload.apField
    PURPOSE:  read apField file (downloading if necessary)
    USAGE:  ret = apload.apField(location)
    RETURNS: if hdu==None : ImageHDUs (all extensions)
             if hdu=N : returns (data, header) for specified HDU
    """
    if len(args) != 1 :
        print('Usage: apField(location)')
    else :
        try :
            file = allfile(
               'Field',location=args[0],
               apred=apred,apstar=apstar,aspcap=aspcap,results=results,dr='collab')
            return _readhdu(file,**kwargs)
        except :
            printerror()

def aspcapField(*args, **kwargs) :
    """
    NAME: apload.aspcapField
    PURPOSE:  read aspcapField file (downloading if necessary)
    USAGE:  ret = apload.aspcapField(location)
    RETURNS: if hdu==None : ImageHDUs (all extensions)
             if hdu=N : returns (data, header) for specified HDU
    """
    if len(args) != 1 :
        print('Usage: aspcapField(location)')
    else :
        try :
            file = allfile(
               'aspcapField',location=args[0],
               apred=apred,apstar=apstar,aspcap=aspcap,results=results,dr='collab')
            return _readhdu(file,**kwargs)
        except :
            printerror()

def cmjd(frame) :
    """ Get chracter MJD from frame number """
    num = (frame - frame%10000 ) / 10000
    return('{:05d}'.format(int(num)+55562) )

def _readchip(file,root,hdu=None,tuple=None,fz=None) :
    """ low level routine to read set of 3 chip files and return data as requested"""
    print('Reading from file: ', file)
    try:
        print (file.replace(root,root+'-a'))
        a=fits.open(file.replace(root,root+'-a'))
        print (file.replace(root,root+'-b'))
        b=fits.open(file.replace(root,root+'-b'))
        print (file.replace(root,root+'-c'))
        c=fits.open(file.replace(root,root+'-c'))
    except:
        print("Can't open file: ", file)
        return(0)

    if hdu is None :
        if tuple :
           print('file: ', file,' read into tuple')
           return a,b,c 
        else :
           print('file: ', file,' read into dictionary with entries a, b, c')
           return {'a' : a, 'b' : b, 'c' : c, 'filename' : file}
    else :
        a[hdu].header.set('filename',os.path.basename(file.replace(root,root+'-a')))
        b[hdu].header.set('filename',os.path.basename(file.replace(root,root+'-b')))
        c[hdu].header.set('filename',os.path.basename(file.replace(root,root+'-c')))
        if tuple :
            data =( a[hdu].data, b[hdu].data, c[hdu].data)
            header =( a[hdu].header, b[hdu].header, c[hdu].header)
        else :
            data ={'a' : a[hdu].data, 'b' : b[hdu].data, 'c': c[hdu].data}
            header = {'a' : a[hdu].header, 'b' : b[hdu].header, 'c': c[hdu].header}
        a.close()
        b.close()
        c.close()
        return data, header

def _readhdu(file,hdu=None,verbose=False) :
    '''
    internal routine for reading all HDU or specified HDU and returning data and header
    '''
    if verbose :
        print('Reading from file: ', file)
    if hdu is None :
        fits.open(file)
        return fits.open(file)
    else :
        hd = fits.open(file)
        data = hd[hdu].data 
        header = hd[hdu].header
        hd.close()
        return data, header

def allfile(root,dr=None,apred=None,apstar=None,aspcap=None,results=None,location=None,obj=None,plate=None,mjd=None,num=None,telescope='apo25m',fiber=None,chips=False,field=None) :
    '''
    Uses sdss_access to create filenames and download files if necessary
    '''
    print('allfile...')
    sdss_path=path.Path()
    http_access=HttpAccess(verbose=True)
    print('http_access.remote..')
    http_access.remote()
    sys.stdout.flush()
    if instrument == 'apogee-n' :
        if root == 'R' :
            prefix='ap'
        else :
            prefix='ap'
    else :
        prefix='as'
        telescope='lco25m'

    if 'all' in root or 'aspcap' in root or 'cannon' in root :
        sdssroot = root 
    else :
        sdssroot = 'ap'+root

    if chips == False :
        # First make sure the file doesn't exist locally
        filePath = sdss_path.full(sdssroot,apred=apred,apstar=apstar,aspcap=aspcap,results=results,
                location=location,obj=obj,plate=plate,mjd=mjd,num=num,telescope=telescope,fiber=fiber,prefix=prefix,instrument=instrument)
        print('filePath',filePath)
        if os.path.exists(filePath) is False: 
            downloadPath = sdss_path.url(sdssroot,apred=apred,apstar=apstar,aspcap=aspcap,results=results,
                location=location,obj=obj,plate=plate,mjd=mjd,num=num,telescope=telescope,fiber=fiber,prefix=prefix,instrument=instrument)
            http_access.get(sdssroot,apred=apred,apstar=apstar,aspcap=aspcap,results=results,
                location=location,obj=obj,plate=plate,mjd=mjd,num=num,telescope=telescope,fiber=fiber,prefix=prefix,instrument=instrument)
        return filePath
    else :
        for chip in ['a','b','c'] :
            print(chip,root,num,mjd,prefix)
            filePath = sdss_path.full(sdssroot,apred=apred,apstar=apstar,aspcap=aspcap,results=results,
                location=location,obj=obj,plate=plate,mjd=mjd,num=num,telescope=telescope,fiber=fiber,
                chip=chip,prefix=prefix,instrument=instrument)
            print('filePath: ', filePath)
            if os.path.exists(filePath) is False: 
                http_access.get(sdssroot,apred=apred,apstar=apstar,aspcap=aspcap,results=results,
                    location=location,obj=obj,plate=plate,mjd=mjd,num=num,telescope=telescope,fiber=fiber,
                    chip=chip,prefix=prefix,instrument=instrument)
        return filePath.replace('-c','')

