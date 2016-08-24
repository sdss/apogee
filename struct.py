"""
Utilities for numpy structured arrays
"""

import numpy as np
import glob
import sys
import pdb
from astropy.io import fits

def pformat(file,val,iformat,fformat,sformat) :
    """ Utility routine for printing a value """
    #print('type: ', val, type(val))
    if isinstance(val,np.ndarray) :
        for v in val :
            #call pformat in case of Nd array
            pformat(file,v,iformat,fformat,sformat)
    elif isinstance(val,(float,np.float32)) :
        file.write(fformat.format(val))
    elif isinstance(val,(int,np.int32,np.int16)) :
        file.write(iformat.format(val))
    else :
        file.write(sformat.format(str(val)))

def list(s,cols=None, cond=None, ind=None, table=None, iformat='{:6d}',fformat='{:8.2f}', sformat='{:<12s}',file=None) :
    """
    List columns of requested row

    Args:
      s (numpy structured array)  : array to show

    Keyword args:
      cols : list of columns to list
      cond : tuple of (column, value); rows where column==value will be listed
      ind : list of index(es) to print
      table : use table format

    Returns:
    """
    if file is None :
       f = sys.stdout
    else :
       f = open(file,'w')

    # Use input columns if given, else all columns
    if cols is None :
        cols = s.names

    # Use condition if given, else specified index (default ind=0)
    if cond is not None :
        inds = np.where(s[cond[0]] == cond[1])[0]
    elif ind is not None :
        try :
            test=len(ind)
        except :
            ind=[ind]
        inds = np.array(ind)
    else :
        inds = np.arange([0])
   
    # if not specified, use table format for multiple entries
    if table is None :
        if len(inds) > 1 :
            table = True
        else :
            table = False

    # in table format, print column names 
    if table :
        for col in cols :
            f.write(sformat.format(col))
        f.write('\n')

    # print
    for i in inds :
        for col in cols :
            if not table :
                f.write(sformat.format(col))
            pformat(f,s[i][col],iformat,fformat,sformat)
            if not table : 
                f.write('\n')
        if table :
            f.write('\n')

def add_cols(a,b):
    """ 
    Add empty columns from recarray b to existing columns from a,
    return new recarray
    """

    # need to handle array elements properly
    newdtype = []
    names = a.dtype.names+b.dtype.names
    descrs = a.dtype.descr+b.dtype.descr
    for i in range(len(descrs)) :
        name= names[i]
        desc= descrs[i]
        if i < len(a.dtype.names) :
            shape= a[name][0].shape
        else :
            shape= b[name][0].shape
        if len(desc) > 2 :
            newdtype.append((desc[0],desc[1],shape))
        else :
            newdtype.append((desc[0],desc[1]))
    # create new array
    newrecarray = np.empty(len(a), dtype = newdtype)
    # fill in all of the old columns
    #print('copying...')
    for name in a.dtype.names:
         #print(name)
         newrecarray[name] = a[name]
    return newrecarray


def append(a,b) :
    '''
    Append two structured arrays, checking for same fields, and increasing size of character fields
    as needed
    '''

    dt_a=a.dtype.descr
    dt_b=b.dtype.descr
    if len(dt_a) != len(dt_b) :
        print("structures don't have same number of fields")

    dt=dt_a
    for i in range(len(dt_a)) :
        if dt_a[i][0] != dt_b[i][0] :
            print("fields don't match",i,dt_a[i],dt_b[i])
        elif dt_a[i][1].find('S') >= 0 :
            j=dt_a[i][1].find('S')
            n=len(dt_a[i][1])
            s_a=int(dt_a[i][1][j+1:n])
            n=len(dt_b[i][1])
            s_b=int(dt_b[i][1][j+1:n])
            dt[i]=(dt_a[i][0],dt_a[i][1][0:j+1]+'{:<d}'.format(max([s_a,s_b])))
            #print(dt_a[i][0],dt_a[i][1],dt_b[i][0],dt_b[i][1],dt[i][1],s_a,s_b)
    dt=np.dtype(dt)
    return np.append(a.astype(dt),b.astype(dt))

def concat(files,hdu=1) :
    '''
    Create concatenation of structures from an input list of files files; structures must have identical tags

    Args:
        files : single file name(str) or list of files, can include wildcards (expanded using glob)
  
    Keyword args:
        hdu=  : specifies which HDU to read/concatenation (default=1)

    Returns:
        structure with concatenated records
    '''
    if type(files) == str:
        files=[files]
    allfiles=[]
    for file in files :
        allfiles.extend(glob.glob(file))
    if len(allfiles) == 0 :
        print('no files found!',file)
        return

    for file in allfiles :
        print(file)
        a=fits.open(file)[hdu].data
        if file == allfiles[0] :
            all=a
        else :
            all=append(all,a)
        print len(all), len(a)
    return all


def wrfits(a,file) :
    '''
    write input FITS structure to file
    '''
    tab=fits.BinTableHDU.from_columns(a)
    tab.writeto(file,clobber=True)

