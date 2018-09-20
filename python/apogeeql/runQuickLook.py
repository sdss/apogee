#!/usr/bin/env python
#Author: Neville Shane, University of Virginia, nss5b@virginia.edu
"""
This script will run the apogee quicklook process manually
If an output file is specified, all quicklook data, with the exception of large arrays,
will be written to that file in the form of CSV tables.

Usage: Specify the MJD, and optionally plate id (-p plateId) or exposure number (-e exposureId).  
Use the -d (--dbinsert) flag to specify a DB insert to QuicklookDB.
Use the -o (--outputfile) flag to specify the name of the file to which you wish to write the quicklook output.
Use the -r (--quickred) flag if you also wish to run Quickred.
Use the -s (--qr_outfile) flag to specify the name of the file to which you wish to write the quickred output.

	runQuickLook mjd -p plateId -e exposureId -o output_file -d -r -p quickred_output_file

Script prerequisites:
----------------------
setup apogeereduce
setup apgquicklook
setup sdss_python_module
"""

from optparse import OptionParser
import sys, os, tempfile, glob, shutil, subprocess
import logging
import pyfits
import traceback
import sqlalchemy

try:
    from sdss.internal.database.connections.APODatabaseAdminLocalConnection import db # access to engine, metadata, Session
#	from sdss.internal.database.connections.APODatabaseDevAdminLocalConnection import db as db_dev #dev db for testing
except ImportError:
    print 'Error on import - did you "setup sdss_python_module" before running this script??\n'
    sys.exit(1)

try:
    from sdss.internal.database.apo.platedb.ModelClasses import *
    from sdss.apogee.addExposure import *
    from sdss.apogee.makeApogeePlugMap import *
except ImportError:
    print 'Error on import - did you "setup sdss_python_module" before running this script??\n'
    try:
        db
    except:
        db.engine.dispose() # avoid "unexpected EOF on client connection" on db server
    sys.exit(1)

#initial setup
mjd = None
plate_id = None
exposure_id = None
output_file = None
qr_output_file = None
prevPlateId = None
prevCartId = None
prevPlugging = None
prevPointing = None
prevScanId = None
prevScanMJD = None
surveyLabel='APOGEE-2'
startOfSurvey = 55562

# Handle inputs
usage_text = '%s <mjd to process>' % sys.argv[0]
description_text = "This script will run the apogee quicklook process for the given MJD/plate/exposure."

parser = OptionParser(usage=usage_text, description=description_text)
    
parser.add_option("-p", "--plateid",
                      dest="plate_id",    # replaces "action", no default
                      help="plateId to process") 
parser.add_option("-e", "--exposureid",
                      dest="exposure_id",    # replaces "action", no default
                      help="exposureId to process") 

parser.add_option("-o", "--outputfile",
						dest="output_file",
						help="write quicklook output to this file")

parser.add_option("-d", "--dbinsert",
						dest="dbinsert",
						action="store_true",
						default=False,
						help="write to database if true")

parser.add_option("-r", "--quickred",
						dest="quickred",
						action="store_true",
						default=False,
						help="also run quickred")

parser.add_option("-s", "--qroutfile",
						dest="qr_output_file",
						help="write quickred output to this file")

global options
(options, args) = parser.parse_args()


# get the info from the command line arguments
if (len(sys.argv) > 1):
	mjd = sys.argv[1]
if options.plate_id:
    plate_id = int(options.plate_id)
if options.exposure_id:
    exposure_id = int(options.exposure_id)
if options.output_file:
	output_file = options.output_file
dbinsert = options.dbinsert
runquickred = options.quickred
if options.qr_output_file:
	qr_output_file = options.qr_output_file


if (mjd == None):
    print
    print "Please specify the mjd to process."
    print
    print "Enter '%s --help' for more details." % sys.argv[0]
    print
    sys.exit()

try:
	data_dir = os.environ["APQLDATA_DIR"]
except:
	raise RuntimeError("Failed: APQLDATA_DIR is not defined")

try:
	spectro_dir = os.environ["APQLSPECTRO_DIR"]
except:
	raise RuntimeError("Failed: APQLSPECTRO_DIR is not defined")

try:
	archive_dir = os.environ["APQLARCHIVE_DIR"]
except:
	raise RuntimeError("Failed: APQLARCHIVE_DIR is not defined")

try:
	quickred_dir = os.environ["APQLQUICKRED_DIR"]
except:
	raise RuntimeError("Failed: APQLQUICKRED_DIR is not defined")

plugmap_dir = '/data-ql/plugmaps/'
ics_dir = '/data-ics/'
raw_dir = os.path.join(data_dir,mjd)

mysession = db.Session()

#Verify that all the UTR files were copied from the ICS (won't do anything if there are no longer files
#in the ICS directory)
dayOfSurvey = str(int(mjd) - startOfSurvey)
indir = os.path.join(ics_dir,dayOfSurvey)

if exposure_id is not None:
	lst = glob.glob(os.path.join(indir,'apRaw-'+str(exposure_id)+'*.fits'))
else:
	lst = glob.glob(os.path.join(indir,'apRaw*.fits'))
lst.sort()
count=0
for infile in lst:
	# check that the file exists in the outdir
	outfile = os.path.join(raw_dir,os.path.basename(infile))
	if not os.path.exists(outfile):
		count+=1
		#just copy the file (no appending of fits keywords)
		shutil.copy(infile,outfile)
if count > 0:
	print '%d missing UTR files were copied from ICS directory' % (count)


#open a temporary file that will contain a list fo exposures, etc, to send to apql_wrapper_manual.pro
f_ql = tempfile.NamedTemporaryFile(delete=False)
listfile = f_ql.name
#if we are also running quickred, open a file taht will contain a list to send to apqr_wrapper_manual.pro
if runquickred:
	f_qr =  tempfile.NamedTemporaryFile(delete=False)
	qrlistfile = f_qr.name

# get list of exposures for the mjd
if exposure_id is not None:
	rawfiles = glob.glob(os.path.join(raw_dir,'apRaw-'+str(exposure_id)+'*.fits'))
	exposures = sorted(set([os.path.basename(fn).split('-')[1] for fn in rawfiles]))
	if len(exposures) == 0:
		raise RuntimeError("Exposure %s not found for MJD %s" % (exposure_id, mjd))
else:
	rawfiles = glob.glob(os.path.join(raw_dir,'apRaw*.fits'))
	exposures = sorted(set([os.path.basename(fn).split('-')[1] for fn in rawfiles]))
	if len(exposures) == 0:
		raise RuntimeError("No exposures found for MJD %s" % (mjd))

count_exp = 0
# for each exposure
for exp in exposures:

	#get raw file and read cartidge and pointing from header
	exp_files = sorted(glob.glob(os.path.join(raw_dir,'apRaw-'+str(exp)+'*.fits')))
	hdulist=pyfits.open(exp_files[0])
	#check that the raw file contains PLATEID in the header so taht we can find plugmap
	if 'PLATEID' in hdulist[0].header:
		plateId = hdulist[0].header['PLATEID']

		#if user has specified a plate and this file is for a different plate, then skip
		if plate_id is not None and plateId != plate_id:
			continue
		
		survey=mysession.query(Survey).join(PlateToSurvey).join(Plate).filter(Plate.plate_id==plateId)
		if survey.count() > 0:
			if survey[0].label.upper().find("APOGEE") == -1 and survey[0].label.upper().find("MANGA") == -1:
				 # not an apogee or marvels plate - just skip
				 continue

		print 'Processing exposure ',exp

		cartId = hdulist[0].header['CARTID']
		pointing = hdulist[0].header['POINTING']
		
		# starttime is MJD in seconds
		startTime = int(mjd)*24.0*3600.0
		time_string = hdulist[0].header['DATE-OBS']
		p = time_string.find(':')
		if p > 0:
			hours = float(time_string[p-2:p])
			minutes = float(time_string[p+1:p+3])
			seconds = float(time_string[p+4:])
			startTime = startTime + seconds + (minutes + hours*60.0) * 60.0
		
		expTime = hdulist[0].header.get('EXPTIME')
		expType = hdulist[0].header.get('IMAGETYP')
	#if the file has been copied straight from ICS without appending the extra fits keywords
	#ask the user to manually enter the plugmap filename
	else:
		print 'File %s does not contain PLATEID in header' %(exp_files[0])
		plateId = 999
		cartId = 999
		pointing = 999
		found = 0
		skip = 0
		while found == 0:
			if prevPlugging is not None:
				ans = raw_input('Use (p)revious plugmap file, %s, enter new plugmap_m filename, or (s)kip: ' %(prevPlugging))
			else:
				ans = raw_input('Enter plugmap_m filename, or (s)kip: ')
			if ans == 's' or ans == 'S':
				skip = 1
				found = 1
			elif ans == 'p' or ans == 'P':
				plateId = prevPlateId
				cartId = prevCartId
				pointing = prevPointing
				found = 1
			else:
				pm_ans = mysession.query(PlPlugMapM).filter(PlPlugMapM.filename==ans)
				if pm_ans.count() == 0:
					 print "No plugmap found with name %s" % (ans)
				else:
					pm_ans = pm_ans[0]
					found = 1
		if skip == 1:
			continue


		


		# starttime is MJD in seconds
		startTime = int(mjd)*24.0*3600.0
		time_string = hdulist[0].header['DATE-OBS']
		p = time_string.find(':')
		if p > 0:
			hours = float(time_string[p-2:p])
			minutes = float(time_string[p+1:p+3])
			seconds = float(time_string[p+4:])
			startTime = startTime + seconds + (minutes + hours*60.0) * 60.0
		
		expTime = hdulist[0].header.get('EXPTIME')
		expType = hdulist[0].header.get('IMAGETYP')

	count_exp += 1
	# get plugmap
	if plateId != prevPlateId or cartId != prevCartId or pointing != prevPointing: 
		if plateId is 999:
			#use plugmap entered manually
			pm = pm_ans
			plateId = plate_id
		else:
			pm = mysession.query(PlPlugMapM).join(Plugging,Plate,Cartridge).\
				filter(Plate.plate_id==plateId).\
				filter(Cartridge.number==cartId).\
				filter(PlPlugMapM.pointing_name==pointing).\
				order_by(PlPlugMapM.fscan_mjd.desc()).order_by(PlPlugMapM.fscan_id.desc())

			if pm.count() == 0:
				raise RuntimeError("No plugmap found for plate %d cartidge %d plugging %s" % (plateId, cartId, pointing))
			else:
				pm = pm[0]		

		# create apogee plugmap file
		fname = pm.filename
		p = fname.find('MapM')
		fname_a  = os.path.join(plugmap_dir,fname[0:p+3]+'A'+fname[p+4:])
		print 'Creating APOGEE plugmap file', fname_a
		makeApogeePlugMap(mysession,pm,fname_a)
		
		prevPlateId = plateId
		prevCartId = cartId
		prevPlugging = fname
		prevPointing = pointing
		prevScanId = pm.fscan_id
		prevScanMJD = pm.fscan_mjd
		
		# write a copy to the archive directory
		# define the current mjd archive directory to store the plPlugMapA file
		arch_dir = os.path.join(archive_dir, mjd)
		if not os.path.isdir(arch_dir):
			os.mkdir(arch_dir, 0o0775)

		res=os.path.split(fname_a)
		archivefile = os.path.join(arch_dir,res[1])
		print 'Archiving APOGEE plugmap file to ',archivefile
		shutil.copyfile(fname_a,archivefile)

	# see if exposure entry already exists in DB  
	exposure = mysession.query(Exposure).filter(Exposure.exposure_no==exp)
	if exposure.count() != 0:
		exp_pk = exposure[0].pk
	else:
		if dbinsert:
			# create exposure entry in DB
			survey = mysession.query(Survey).filter(Survey.label==surveyLabel)
			survey = survey[0]
			# check if exposure already exists in DB, otherwise create it
			exp_obj = mysession.query(Exposure).filter(Exposure.survey_pk==survey.pk).filter(Exposure.exposure_no==int(exp))
			if exp_obj.count() == 1:
				exp_pk=exp_obj[0].pk
			elif exp_obj.count() == 0: 	
				exp_pk = addExposure(mysession, prevScanId, prevScanMJD, prevPlateId, mjd, int(exp), surveyLabel, startTime, expTime, expType, 'Manual apogeeQL')
			else:
				raise RuntimeError("ERROR: Multiple exposures already exist for exposure number %d, survey %s" \
								 % (exposureNo, survey.label))
		else:
			exp_pk = 'None'

	plugfile = fname_a

	# write mjd, exp, plate_id and plug file to tmp file to send to apql_wrapper_manual 
	f_ql.write('%s, %s, %s, %s\n' %(mjd, exp, plateId, plugfile))
	# if we are also running quickred, write mjd, exp, exp_pk, and plug file to tmp file to send to apqr_wrapper_manual
	if runquickred:
		f_qr.write('%s, %s, %s, %s\n' %(mjd, exp, exp_pk, plugfile))

f_ql.close()
if runquickred: 
	f_qr.close()

if count_exp > 0 :
	# run apql_wrapper_manual IDL code
	if output_file is not None:
		ql_cmd = 'idl -e "apql_wrapper_manual,\'%s\',no_dbinsert=%i,outfile=\'%s\',data_dir=\'%s\',spectro_dir=\'%s\'"' \
			 % (listfile, not dbinsert, output_file, data_dir, spectro_dir)
	else:
		ql_cmd = 'idl -e "apql_wrapper_manual,\'%s\',no_dbinsert=%i,data_dir=\'%s\',spectro_dir=\'%s\'"'  \
			 % (listfile, not dbinsert, data_dir, spectro_dir)
	print ql_cmd

	ql_process = subprocess.Popen(ql_cmd, stderr=subprocess.PIPE, shell=True)
	output=ql_process.communicate()[0] 
	if output is not None:
		print output

	# if quickred requested, also run apqr_wrapper_manual IDL code
	if runquickred:
		print 'RUNNING QUICKRED'
		if qr_output_file is not None:
			qr_cmd = 'idl -e "apqr_wrapper_manual,\'%s\',no_dbinsert=%i,outfile=\'%s\',data_dir=\'%s\',spectro_dir=\'%s\',archive_dir=\'%s\',quickred_dir=\'%s\'"' \
				% (qrlistfile, not dbinsert, qr_output_file, data_dir, spectro_dir, archive_dir, quickred_dir)
		else:
			qr_cmd = 'idl -e "apqr_wrapper_manual,\'%s\',no_dbinsert=%i,data_dir=\'%s\',spectro_dir=\'%s\',archive_dir=\'%s\',quickred_dir=\'%s\'"' \
				% (qrlistfile, not dbinsert, data_dir, spectro_dir, archive_dir, quickred_dir)
		print qr_cmd
		
		qr_process = subprocess.Popen(qr_cmd, stderr=subprocess.PIPE, shell=True)		
		output=qr_process.communicate()[0] 
		if output is not None:
			print output


else:
	raise RuntimeError("No APOGEE or MANGA exposures found for plate %s on MJD %s" % (plate_id, mjd))

# delete tmp list files
os.remove(listfile)
if runquickred:
	os.remove(qrlistfile)
# close db connection
db.engine.dispose()
