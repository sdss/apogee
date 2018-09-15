import os
import numpy as np
from astropy.io import ascii
import astropy.table as table

def read(infile,columns=None,agerange=[0,20.]) :
    """
    Routine to read a Padova isochrone file and return
    a numpy structured array with the contents

    Args: 
        file name : input data file name
        columns=[list] : list of columns to extract
        age = age : single age to extract

    Returns: 
        structured array with isochrone data
    """
    if os.getenv('ISOCHRONE_DIR') != "" and infile.find('/') < 0 :
        data=ascii.read(os.getenv('ISOCHRONE_DIR')+'/'+infile,
             names=['z','age','mini','mact','logl','logte','logg',
                    'mbol','u','b','v','r','i','j','h','k','intimf','stage'])
    else :
        data=ascii.read(infile,
             names=['z','age','mini','mact','logl','logte','logg',
                    'mbol','u','b','v','r','i','j','h','k','intimf','stage'])

    # add some derived columns
    data.add_column(table.column.Column(name='feh',data=np.log10(data['z']/0.0152)))
    data.add_column(table.column.Column(name='teff',data=10.**(data['logte'])))
    data.add_column(table.column.Column(name='jk',data=data['j']-data['k']))

    # select ages within specified age range
    gd = np.where((data['age'] >=agerange[0]) & (data['age'] <= agerange[1]))
    data=data[gd]

    # default columns
    if columns is None:
        # can set default columns here, or keep all quantitites
        #data.keep_columns(['Z','age','logte','logl','intimf','stage'])
        pass
    # option to extract specified columns
    else:
        data.keep_columns(columns)

    return data


