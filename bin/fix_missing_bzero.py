#!/usr/bin/env python
"""
To fix the files in MJD 56532-56540 where the bzero field wasn't written
during the annotation, due to the pyfits 2.4->3.1 change.
"""
import os
import sys
import getpass
import socket
import glob
import pyfits

suffix = 'badheaders'

def get_yes(text):
    getOk = raw_input(text)
    if not getOk.lower().startswith("y"):
        return False
    else:
        return True

def fix_header(hdulist):
    """
    Replace the missing BZERO/BSCALE values in the header,
    in the correct place.
    """
    bscale = int(hdulist[0].header.get('BSCALE',1))
    bzero = int(hdulist[0].header.get('BZERO',32768))
    del hdulist[0].header['BSCALE']
    del hdulist[0].header['BZERO']
    hdulist[0].header.update('BSCALE',bscale,after='GCOUNT')
    hdulist[0].header.update('BZERO',bzero,after='BSCALE')
    return hdulist

def fix_file(oldfile,newfile,force=False):
    """Fix the headers of oldfile, writing them to newfile."""
    pre = '' if force else '(NOT) '
    print "%sFixing header: %s -> %s"%(pre,oldfile,newfile)
    if not force: oldfile = newfile
    if force:
        hdulist = pyfits.open(oldfile,uint16=True,do_not_scale_image_data=True,checksum=True)
        hdulist = fix_header(hdulist)
        hdulist.writeto(newfile,clobber=False,output_verify='warn',checksum=True)
        os.chmod(newfile,0o444)

def move_file(file,force=False):
    """Move file to FILEDIR/badheaders/."""
    pre = '' if force else '(NOT) '
    newdir = os.path.join(os.path.split(file)[0],suffix)
    if not os.path.exists(newdir):
        os.mkdir(newdir)
    newfile = os.path.join(newdir,os.path.split(file)[-1])
    print "%sMoving: %s -> %s"%(pre,file,newfile)
    if force:
        os.rename(file,newfile)
        os.chmod(newfile,0o444)
    return newfile

def file_is_ok(file):
    """Return True if file already has BZERO in its header."""
    if pyfits.getheader(file).get('BZERO',None):
        print "%s already has BZERO. Not processing."%file
        return True
    else:
        return False

def fix_files(files,force=False):
    """Fix the files, asking about each one unless force is set."""
    for file in files:
        if file_is_ok(file):
            continue
        if not force:
            forceThis = get_yes("Process %s y/[n]? "%file)
        else:
            forceThis = True
        moved = move_file(file,forceThis)
        fix_file(moved,file,forceThis)

def main(argv=None):
    from optparse import OptionParser
    if argv is None: argv = sys.argv[1:]

    usage = "%prog [OPTIONS] DIR1 [DIR2 [DIR3 [...]]]"
    usage += "\n\nRepair BZERO/BSCALE headers in the apRaw files in DIR, after"
    usage += "\nmoving them to a backup directory."
    usage += "\nRequest permission to process each file, unless -f is specified."
    usage += "\n\nDIR example (on apogee-ql): /data-ql/data/56535"
    parser = OptionParser(usage)
    parser.add_option('--yes',dest='yes',action='store_true',default=False,
                      help='Process all the files in a directory without asking (%default).')
    (opts,args) = parser.parse_args(args=argv)
    
    #if opts.yes and getpass.getuser() != 'sdss3' and 'apogee-ql' not in socket.gethostname():
    #    print "If 'yes', we must be run as sdss3@apogee-ql!"
    #    sys.exit(-2)
    
    if len(args) == 0:
        print "Need at least one directory as an argument."
        sys.exit(-1)
    for directory in args:
        print "Processing:",directory
        files = sorted(glob.glob(os.path.join(directory,'apRaw*.fits')))
        fix_files(files,opts.yes)

if __name__ == "__main__":
    main()
