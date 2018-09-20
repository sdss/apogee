#!/usr/bin/env python
"""
This script will transform a plPlugMapM file to a plPlugMapA file for
APOGEE, adding the 2-Mass JHK magnitudes to the table.

Usage: Specify the plateId

    plugmapm2a.py -p 4385


Script prerequisites:
---------------------
setup sdss_python_module
"""

import os
import sys
import time
import types
import string

try:
    from sdss.internal.database.connections.APODatabaseUserLocalConnection import db # access to engine, metadata, Session
except ImportError:
    print 'Error on import - did you "setup sdss_python_module" before running this script??\n'
    sys.exit(1)

try:
    print "Importing PlateDB"
    from sdss.internal.database.apo.platedb.ModelClasses import *
except ImportError:
    print 'Could not create ModelClasses - did you "setup sdss_python_module" before running this script??\n'
    try:
        db
    except:
        db.engine.dispose() # avoid "unexpected EOF on client connection" on db server
    sys.exit(1)

try:
  from sdss.apogee.makeApogeePlugMap import *
except ImportError:
    print 'Error on import (sdss.apogee.makeApogeePlugMap) - did you "setup sdss_python_module" before running this script??\n'
    sys.exit(1)

import sqlalchemy
from sqlalchemy import not_, and_, or_
#from sdss.utilities import yanny

import pyfits
import traceback

def main(argv=None):
    from optparse import OptionParser
    if argv is None: argv = sys.argv[1:]
 
    plateId = None
    usage_text = '%s plateId ' % sys.argv[0]
    description_text = "Append the APOG/EE specifc information to a plPlugMapM file and writes a new " + \
                "plPlugMapA file."
    
    parser = OptionParser(usage=usage_text, description=description_text)
    
    parser.add_option("-o", "--overwrite",
                      action="store_true", # for boolean options
                      dest="overwrite", # this will be the variable name
                      default=False,
                      help="overwrite any existing plPlugMapA file")
    
    parser.add_option("-p", "--plateid",
                      dest="plateId",    # replaces "action", no default
                      help="plateId to process")
    
    global options
    (options, args) = parser.parse_args()
    
    if (len(sys.argv) > 1):
        plateId = sys.argv[1]
    else:
        # get the info from the command line arguments
        if options.plateId:
            plateId = int(options.plateId)
    
    # make sure we have all the parameters to run
    if (plateId == None):
        print
        print "Please specify the plateId to process."
        print
        print "Enter '%s --help' for more details." % sys.argv[0]
        print
        # db.engine.dispose() # avoid "unexpected EOF on client connection" on db server
        sys.exit()
    
    mysession = db.Session()
    
    # look for all the matching entries
    pm = mysession.query(PlPlugMapM).join(Plugging).join(Plate).filter(Plate.plate_id==plateId).\
            order_by(PlPlugMapM.fscan_mjd.desc()).order_by(PlPlugMapM.fscan_id.desc())
    
    if pm.count() == 0:
        print 'No entries for plate ',plateId
        sys.exit()
    else:
        for i in range(pm.count()):
            print '%d   pm.filename=%s     pm.fscan_mjd=%d    pm.fscan_id=%d' % \
                (i+1,pm[i].filename,pm[i].fscan_mjd,pm[i].fscan_id)
    
    id=''
    id = raw_input("Select file to modify [1]: ")
    if id == '':
        id=0
    else:
        id = int(id)-1
    
    if id+1 > pm.count():
        print 'wrong pm selected'
        sys.exit()
    
    pm=pm[id]
    
    fname = pm.filename
    p = fname.find('MapM')
    fname  = os.path.join('/data-ql/plugmaps/',fname[0:p+3]+'A'+fname[p+4:])
    makeApogeePlugMap(mysession, pm, fname)
    
    print 'wrote ',fname

if __name__ == '__main__':
    main()
