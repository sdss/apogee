#!/usr/bin/env python

from twisted.internet.protocol import Protocol, Factory, ClientFactory
from twisted.protocols.basic import LineReceiver
from twisted.internet import reactor
from twisted.internet.protocol import Protocol, Factory

from sdss.internal.database.connections.APODatabaseAdminLocalConnection import db # access to engine, metadata, Session
from sdss.internal.database.apo.platedb.ModelClasses import *
from sdss.apogee.addExposure import addExposure

import sqlalchemy

import opscore.actor.model
import opscore.actor.keyvar

import actorcore.Actor
import actorcore.CmdrConnection as actorCmdrConnection
import actorcore.utility.fits as actorFits

import pyfits, warnings
import actorkeys
import traceback

#
# Import sdss3logging before logging if you want to use it
#
import logging
import os, signal, subprocess, tempfile, shutil, glob
import time
import types
from sdss.utilities import yanny
import RO.Astro.Tm.MJDFromPyTuple as astroMJD
from astropy.time import Time


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class QuickLookLineServer(LineReceiver):
    def connectionMade(self):
        # the IDL end-of-line string is linefeed only (no carriage-return) 
        # we need to change the delimiter to be able to detect a full line
        # and have lineReceived method be called (otherwise it just sits there)
        self.delimiter = '\n'
        self.peer = self.transport.getPeer()
        self.factory.qlActor.qlSources.append(self)
        logging.info("Connection from %s %s" % (self.peer.host, self.peer.port))
        print "Connection from ", self.peer.host, self.peer.port
        # ping the quicklook
        if self.factory.qlActor.ql_pid > 0:
           for s in self.factory.qlActor.qlSources:
              s.sendLine('PING')
              s.sendLine('STARTING')

    def lineReceived(self, line):
        if line.upper()=='PONG':
            # normal response from aliveness test
            self.factory.qlActor.watchDogStatus = True
        elif line.upper()=='QUIT':
            # request to disconnect
            self.transport.loseConnection()
        elif line.upper()=='STARTED':
            # we got the initial response from the apql_wrapper
            # send the PointingInfo stuff
            if self.factory.qlActor.prevPlate != -1:
               for s in self.factory.qlActor.qlSources:
                  s.sendLine('plugMapInfo=%s,%s,%s,%s' % (self.factory.qlActor.prevPlate, \
                       self.factory.qlActor.prevScanMJD, self.factory.qlActor.prevScanId, \
                       self.factory.qlActor.plugFname))
        elif line == "callback":
            logging.info("preparing callback")
            reactor.callLater(5.0,self.sendcomment)
        else:
            # assume the messages are properly formatted to pass along
            # print 'Received from apql_wrapper.pro: ',line
            self.factory.qlActor.bcast.finish(line)

    def connectionLost(self, reason):
        logging.info("Disconnected from %s %s"  % (self.peer.port, reason.value))

    def sendcomment(self):
        logging.info("in the callback routine")


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class QLFactory(ClientFactory):
    protocol = QuickLookLineServer
    def __init__(self, qlActor): 
        self.qlActor=qlActor

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class QuickRedLineServer(LineReceiver):
    def connectionMade(self):
        # the IDL end-of-line string is linefeed only (no carriage-return) 
        # we need to change the delimiter to be able to detect a full line
        # and have lineReceived method be called (otherwise it just sits there)
        self.delimiter = '\n'
        self.peer = self.transport.getPeer()
        self.factory.qrActor.qrSources.append(self)
        logging.info("apqr_wrapper Connection from %s %s" % (self.peer.host, self.peer.port))
        print "apqr_wrapper -> Connection from ", self.peer.host, self.peer.port
        # ping the quicklook
        if self.factory.qrActor.ql_pid > 0:
           for s in self.factory.qrActor.qrSources:
              s.sendLine('PING')
              s.sendLine('STARTING')

    def lineReceived(self, line):
        if line.upper()=='PONG':
            # normal response from aliveness test
            self.factory.qrActor.watchDogStatus = True
        if line=='quit':
            # request to disconnect
            self.transport.loseConnection()
        elif line.upper()=='STARTED':
            # we got the initial response from the apql_wrapper
            # send the PointingInfo stuff
            if self.factory.qrActor.prevPlate != -1:
               for s in self.factory.qrActor.qrSources:
                  s.sendLine('plugMapInfo=%s,%s,%s,%s' % (self.factory.qrActor.prevPlate, \
                       self.factory.qrActor.prevScanMJD, self.factory.qrActor.prevScanId, \
                       self.factory.qrActor.plugFname))
        elif line == "callback":
            logging.info("preparing callback")
            reactor.callLater(5.0,self.sendcomment)
        else:
            # assume the messages are properly formatted to pass along
            # print 'Received from apqr_wrapper.pro: ',line
            self.factory.qrActor.bcast.finish(line)

    def connectionLost(self, reason):
        logging.info("Disconnected from %s %s"  % (self.peer.port, reason.value))

    def sendcomment(self):
        logging.info("in the callback routine")


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class QRFactory(ClientFactory):
    protocol = QuickRedLineServer
    def __init__(self, qrActor): 
        self.qrActor=qrActor

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class Apogeeql(actorcore.Actor.Actor):

   models = {}
   ql_pid = 0
   qlSources=[]
   qr_pid = 0
   qrSources=[]
   prevCartridge=-1
   prevPlate=-1
   prevScanId = -1
   prevScanMJD = -1
   prevPointing='A'
   pluggingPk = 0
   apogeeSurveyPk=1   # use the known APOGEE survey value in the db as default
   plugFname=''
   inst = ''
   startExp = False
   endExp = False
   expState=''
   expType=''
   numReadsCommanded=0
   actor=''
   obs_pk = 0
   exp_pk = 0
   cdr_dir = ''
   summary_dir = ''
   frameid = ''
   ditherPos = 0.0
   namedDitherPos = ''

   # setup the variables for the watchdog to the IDL code
   watchDogStatus=True  # if True, the idl code is stil responding
   watchDogTimer=2.0    # should return a reply within this time if still alive

   def __init__(self, name, productName=None, configFile=None, debugLevel=20):
      actorcore.Actor.Actor.__init__(self, name, productName=productName, configFile=configFile)

      Apogeeql.actor=self
      self.headURL = '$HeadURL: https://svn.sdss.org/repo/operations/general/actors/apogeeql/trunk/python/apogeeql/apogeeql_main.py $'

      self.logger.setLevel(debugLevel)
      self.logger.propagate = True
      self.logger.info('Starting the Apogeeql Actor ...')

      # only the ics_datadir is defined in the cfg file - expect the datadir to be an
      # environment variable set by the apgquicklook setup 
      # (don't want to keep it in 2 different places)
      try:
         self.datadir = os.environ["APQLDATA_DIR"]
      except:
         self.logger.error("Failed: APQLDATA_DIR is not defined")
         traceback.print_exc()

      # environment variable set by the apgquicklook setup (don't want to keep it in 2 different places)
      try:
         self.archive_dir = os.environ["APQLARCHIVE_DIR"]
      except:
         self.logger.error("Failed: APQLARCHIVE_DIR is not defined")
         traceback.print_exc()

      # MJD value of start of survey (Jan 1 2011) for filenames (55562)
      self.startOfSurvey = self.config.get('apogeeql','startOfSurvey')
      self.snrAxisRange = [self.config.get('apogeeql','snrAxisMin'), self.config.get('apogeeql','snrAxisMax')]
      self.rootURL = self.config.get('apogeeql','rootURL')
      self.plugmap_dir = self.config.get('apogeeql','plugmap_dir')
      self.ics_datadir = self.config.get('apogeeql','ics_datadir')
      self.cdr_dir = self.config.get('apogeeql','cdr_dir')
      self.summary_dir = self.config.get('apogeeql','summary_dir')
      self.delFile = self.config.get('apogeeql','delFile')
      self.qlPort = self.config.getint('apogeeql', 'qlPort') 
      self.qlHost = self.config.get('apogeeql', 'qlHost') 
      self.qrPort = self.config.getint('apogeeql', 'qrPort') 
      self.qrHost = self.config.get('apogeeql', 'qrHost') 
      self.criticalDiskSpace = self.config.get('apogeeql', 'criticalDiskSpace') 
      self.seriousDiskSpace = self.config.get('apogeeql', 'seriousDiskSpace') 
      self.warningDiskSpace = self.config.get('apogeeql', 'warningDiskSpace')

      #
      # Explicitly load other actor models. We usually need these for FITS headers.
      #
      self.models = {}
      # for actor in ["mcp", "guider", "platedb", "tcc", "apo", "apogeetest"]:
      for actor in ["mcp", "guider", "platedb", "tcc", "apogee", "apogeecal", "sop"]:
         self.models[actor] = opscore.actor.model.Model(actor)

      #
      # register the keywords that we want to pay attention to
      #
      self.models["tcc"].keyVarDict["inst"].addCallback(self.TCCInstCB, callNow=False)
      self.models["platedb"].keyVarDict["pointingInfo"].addCallback(self.PointingInfoCB, callNow=True)
      self.models["apogee"].keyVarDict["exposureState"].addCallback(self.ExposureStateCB, callNow=False)
      self.models["apogee"].keyVarDict["exposureWroteFile"].addCallback(self.exposureWroteFileCB, callNow=False)
      self.models["apogee"].keyVarDict["exposureWroteSummary"].addCallback(self.exposureWroteSummaryCB, callNow=False)
      self.models["apogee"].keyVarDict["ditherPosition"].addCallback(self.ditherPositionCB, callNow=False)

      #
      # Connect to the platedb
      #
      self.mysession = db.Session()


   @staticmethod
   def TCCInstCB(keyVar):
      '''callback routine for tcc.inst'''

      # print 'TCCInstCB keyVar=',keyVar
      Apogeeql.inst = keyVar[0]

   @staticmethod
   def PointingInfoCB(keyVar):
      '''callback routine for platedb.pointingInfo'''
      # plate_id, cartridge_id, pointing_id, boresight_ra, boresight_dec, hour_angle, temperature, wavelength
      # print "PointingInfoCB=",keyVar

      # if pointing is None than just skip this
      if keyVar[1] == None:
         return

      plate = int(keyVar[0])
      cartridge = int(keyVar[1])
      pointing = str(keyVar[2])

      """
      MODIFIED TO TESTING QUICKLOOK ON TEST DATABASE WITH SIMULATED DATA AND FAKE ICS

      cartridge = 3
      plate = 4918
      pointing = 'A'

      cartridge = 1
      plate = 4929
      pointing = 'A'
      """

      # print Apogeeql.actor.models['platedb'].keyVarDict['activePlugging']

      # find the platedb.survey.pk corresponding to APOGEE (-2)
      survey = Apogeeql.actor.mysession.query(Survey).filter(Survey.label=='APOGEE-2')
      if survey.count() > 0:
         Apogeeql.actor.apogeeSurveyPk = survey[0].pk

      if plate != Apogeeql.prevPlate or cartridge != Apogeeql.prevCartridge or pointing != Apogeeql.prevPointing:
         # we need to ignore all plates that are not for APOGEE or MANGA
         survey=Apogeeql.actor.mysession.query(Survey).join(PlateToSurvey).join(Plate).filter(Plate.plate_id==plate)
         if survey.count() > 0:
              if survey[0].label.upper().find("APOGEE") == -1 and survey[0].label.upper().find("MANGA") == -1:
                  # not an apogee or marvels plate - just skip
                  return

         # we need to extract and pass a new plugmap to IDL QuickLook
         pm = Apogeeql.actor.getPlPlugMapM(Apogeeql.actor.mysession, cartridge, plate, pointing)

         # routine returns a yanny par file
         # replace plPlugMapM-xxxx by plPlugMapA-xxxx
         # str() to convert from unicode, as twisted can't take it.
         fname = str(pm.filename)
         p = fname.find('MapM')
         fname  = os.path.join(Apogeeql.actor.plugmap_dir,fname[0:p+3]+'A'+fname[p+4:])

         # print 'fname=',fname
         Apogeeql.actor.makeApogeePlugMap(pm, fname)

         # pass the info to IDL QL
         for s in Apogeeql.actor.qlSources:
             s.sendLine('plugMapInfo=%s,%s,%s,%s' % (plate, pm.fscan_mjd, pm.fscan_id, fname))
         for s in Apogeeql.actor.qrSources:
             s.sendLine('plugMapInfo=%s,%s,%s,%s' % (plate, pm.fscan_mjd, pm.fscan_id, fname))

         # print 'plugMapFilename=%s' % (fname)
         Apogeeql.prevPointing = pointing
         Apogeeql.prevPlate = plate
         Apogeeql.prevCartridge = cartridge
         Apogeeql.prevScanId = pm.fscan_id
         Apogeeql.prevScanMJD = pm.fscan_mjd
         # the plugging_pk is needed to find the right observation_pk
         Apogeeql.pluggingPk = pm.plugging_pk
         Apogeeql.plugFname = fname


   @staticmethod
   def ExposureStateCB(keyVar):
      '''callback routine for apogeeICC.exposureState '''
      # exposureState expState expType nReads expName
      #
      # exposureState EXPOSING SCIENCE 50 apRaw-0054003      -> mark start of exposure
      # utrReadState apRaw-0054003 READING 1 50              -> mark first UTR read status
      # utrReadState apRaw-0054003 SAVING 1 50              
      # utrReadState apRaw-0054003 DONE 1 50              
      # exposureWroteFile apRaw-0054003-001.fits             -> mark first UTR read
      # exposureWroteFile apRaw-0054003-002.fits             -> mark first UTR read
      # exposureWroteSummary apRaw-0054003.fits              -> mark first CDS read
      # ...
      # exposureWroteFile apRaw-0054003-050.fits             -> mark last requested UTR read
      # exposureWroteSummary apRaw-0054003.fits              -> mark first UTR read
      # exposureState DONE SCIENCE 50 apRaw-0054003          -> mark end of exposure
      # test existance of passed variable
      # print 'ExposureStateCB keyVar=',keyVar
      if not keyVar.isGenuine:
         return

      Apogeeql.expState=keyVar[0]
      if Apogeeql.expState.upper() == 'EXPOSING':
         # check the communication to IDL quicklook
         # Apogeeql.startQuickLook(Apogeeql.actor)
         Apogeeql.startExp = True
         Apogeeql.startExp = True
         Apogeeql.endExp = False
         Apogeeql.expType = keyVar[1].upper()
         Apogeeql.numReadsCommanded = int(keyVar[2])
         # we may have to do something special for QL at the start of a new exposure
      elif Apogeeql.expState.upper() in ['DONE', 'STOPPED', 'FAILED']:
         # ignore if we weren't actually exposing
         if Apogeeql.startExp == True:
            Apogeeql.startExp = False
            Apogeeql.endExp = True
            Apogeeql.expType = keyVar[1].upper()
            filebase = keyVar[3]
            # make sure we have all the UTR files before bundling
            Apogeeql.completeUTR(Apogeeql.actor,filebase)
            Apogeeql.numReadsCommanded = 0
            res = keyVar[3].split('-')
            Apogeeql.frameid = res[1][:8]
            mjd5 = int(Apogeeql.frameid[:4]) + int(Apogeeql.actor.startOfSurvey)
            # do something for the quickreduction at the end of an exposure
            for s in Apogeeql.qlSources:
               s.sendLine('UTR=DONE')
            for s in Apogeeql.qrSources:
               s.sendLine('UTR=DONE,%s,%d,%s' % (Apogeeql.frameid, mjd5, Apogeeql.exp_pk))
      elif Apogeeql.expState.upper() != 'STOPPING':
         # when a stop was requested, a couple of images will still be coming in
         # only change the startExp in other cases
         Apogeeql.startExp = False


   @staticmethod
   def exposureWroteFileCB(keyVar):
      '''callback routine for apogeeICC.exposureState '''
      # exposureWroteFile apRaw-0054003-001.fits   -> mark first UTR read
      # exposureWroteFile apRaw-0054003-002.fits   -> mark first UTR read
      # exposureWroteSummary apRaw-0054002.fits    -> mark first CDS read
      # ...
      # exposureWroteFile apRaw-0054003-050.fits   -> mark last requested UTR read
      # exposureWroteSummary apRaw-0054003.fits    -> mark last CDS read
      # exposureState DONE SCIENCE 50 apRaw-0054003          -> mark end of exposure
      # test existance of passed variable

      # if we're not actually exposing than skip (we get these messages from the hub when
      # starting the apogeeql if apogee ics is already started)
      # print "exposureWroteFileCB=",keyVar
      if not keyVar.isGenuine:
         return

      filename=keyVar[0]
      if filename == None:
         return

      if len(filename) == 0:
         print "exposureWroteFileCB  -> Null filename received"
         return

      # create a new FITS file by appending the telescope fits keywords
      newfilename, starttime, exptime = Apogeeql.appendFitsKeywords(Apogeeql.actor,filename)

      #Don't create a new exposure if the exposure is not an APOGEE or MANGA object
      if Apogeeql.prevPlate == -1:
         print "exposureWroteFileCB  -> Not an APOGEE/MANGA Object"
         return

      # get the mjd from the filename
      res=filename.split('-')
      try:
         mjd = int(res[1][:4]) + int(Apogeeql.actor.startOfSurvey)
         readnum=int(res[2].split('.')[0])
         expnum=int(res[1])
      except:
         raise RuntimeError( "The filename doesn't match expected format (%s)" % (filename)) 

      if readnum == 1 or Apogeeql.exp_pk == 0:

         #Create new exposure object
         #currently hard-coded for survey=APOGEE-2
         surveyLabel='APOGEE-2'
         try:
             Apogeeql.exp_pk = addExposure(Apogeeql.actor.mysession, Apogeeql.prevScanId, Apogeeql.prevScanMJD, 
                                           Apogeeql.prevPlate, mjd, expnum, surveyLabel, starttime, exptime, Apogeeql.expType, 'apogeeQL')
         except RuntimeError as e:
             Apogeeql.actor.logger.error('Failed in call addExposure for exposureNo %d' %expnum)
             Apogeeql.actor.logger.error('Exception: %s'%e)
             raise RuntimeError('Failed in call addExposure for exposureNo %d' %expnum +'\n'+str(e))


      for s in Apogeeql.qlSources:
         s.sendLine('UTR=%s,%d,%d,%d' % (newfilename, Apogeeql.exp_pk, readnum, Apogeeql.numReadsCommanded))

   @staticmethod
   def exposureWroteSummaryCB(keyVar):

      '''callback routine for apogeeICC.exposureWroteSummary '''
      # exposureWroteFile apRaw-0054003-001.fits   -> mark first UTR read
      # exposureWroteFile apRaw-0054003-002.fits   -> mark first UTR read
      # exposureWroteSummary apRaw-0054003.fits    -> mark first CDS read
      # ...
      # exposureWroteFile apRaw-0054003-050.fits   -> mark last requested UTR read
      # exposureWroteSummary apRaw-0054003.fits    -> mark first CDS read
      # exposureState DONE SCIENCE 50 apRaw-0054003          -> mark end of exposure
      # test existance of passed variable
      # print "exposureWroteSummaryCB=",keyVar
      if not keyVar.isGenuine:
         return

      filename=keyVar[0]
      res=(filename.split('-'))[1].split('.')
      dayOfSurvey = res[0][:4]
      mjd = int(dayOfSurvey) + int(Apogeeql.actor.startOfSurvey)
      indir=Apogeeql.actor.summary_dir
      outdir = Apogeeql.actor.cdr_dir

      try:
         indir  = os.path.join(indir,dayOfSurvey)
         infile= os.path.join(indir,filename)
         outdir = os.path.join(outdir,str(mjd))
         outfile= os.path.join(outdir,filename)
         if not os.path.isdir(outdir):
            os.mkdir(outdir)
            # print 'Directory created at: ' + dest
         # t0=time.time()
         shutil.copy2(infile,outfile)
         # print "shutil.copy2 took %f seconds" % (time.time()-t0)
      except:
         raise RuntimeError( "Failed to copy the summary file (%s)" % (filename)) 

   @staticmethod
   def ditherPositionCB(keyVar):

      '''callback routine for apogeeICC.ditherPosition '''
      # ditherPosition=13.9977,A
      # test existance of passed variable
      # print "ditherPositionCB=",keyVar
      if not keyVar.isGenuine:
         return

      # save the current dither pixel and named position
      Apogeeql.ditherPos = float(keyVar[0])
      Apogeeql.namedDitherPos = keyVar[1]
      for s in Apogeeql.actor.qlSources:
          s.sendLine('ditherPosition=%f,%s' % (Apogeeql.ditherPos, Apogeeql.namedDitherPos))


   def connectQuickLook(self):
      '''open a socket through twisted to send/receive information to/from apogee_IDL'''
      # get the port from the configuratio file 
      reactor.listenTCP(self.qlPort, QLFactory(self))

   def startQuickLook(self):
      '''Open a twisted reactor to communicate with IDL socket'''
      #
      # check if an apogeeql_IDL process is already running before starting a new one
      if self.ql_pid > 0:
         self.stopQuickLook()

      # spawn the apogeeql IDL process and don't wait for its completion
      try:
         # get the string corresponding to the command to start the IDL quicklook process
         qlCommand = self.config.get('apogeeql','qlCommandName')
         qlCommand = qlCommand.strip('"')
         # this adds the arguments to the IDL command line
         qlCommand += " -args %s %s" % (self.qlHost, self.qlPort)
         # Popen does NOWAIT by default
         ql_process = subprocess.Popen(qlCommand.split(), stderr=subprocess.STDOUT)
         self.ql_pid = ql_process.pid
      except:
         self.logger.error("Failed to start the apogeeql_IDL process")
         traceback.print_exc()


   def stopQuickLook(self):
      '''If a quickLook IDL process already exists - just kill it (for now)'''
      # first send a QUIT command to the program
      for s in self.qlSources:
         s.sendLine('QUIT')

      if self.ql_pid == 0:
         p1=subprocess.Popen(['/bin/ps'],stdout=subprocess.PIPE)
         # apql_wrapper is the name of the program ran when starting IDL (in apogeeql.cfg)
         processName = (self.config.get('apogeeql','qlCommandName')).split()
         if '-e' in processName:
            # look at the name of the program called when IDL is started
            pos = processName[processName.index('-e')+1]
         else:
            pos = 0
         args = ['grep',processName[pos]]
         p2=subprocess.Popen(args,stdin=p1.stdout,stdout=subprocess.PIPE)
         output=p2.communicate()[0]
         p2.kill()
         p1.kill()
         if len(output) > 0:
            # process exists -> kill it
            self.ql_pid = output.split()[0]

      os.kill(self.ql_pid,signal.SIGKILL)
      killedpid, stat = os.waitpid(self.ql_pid, os.WNOHANG)
      if killedpid == 0:
         # failed in killing old IDL process
         # print error message here
         self.logger.warn("Unable to kill existing apogeeql_IDL process %s" % self.ql_pid)
      else:
         self.ql_pid = 0

   def connectQuickReduce(self):
      '''open a socket through twisted to send/receive information to/from apqr_wrapper IDL'''
      # get the port from the configuratio file 
      reactor.listenTCP(self.qrPort, QRFactory(self))

   def startQuickReduce(self):
      '''Open a twisted reactor to communicate with IDL socket'''
      #
      # check if an apogeeql_IDL process is already running before starting a new one
      if self.qr_pid > 0:
         self.stopQuickReduce()

      # spawn the apgreduce IDL process and don't wait for its completion
      try:
         # get the string corresponding to the command to start the IDL quicklook process
         qrCommand = self.config.get('apogeeql','qrCommandName')
         qrCommand = qrCommand.strip('"')
         # this adds the arguments to the IDL command line
         qrCommand += " -args %s %s" % (self.qrHost, self.qrPort)
         # Popen does NOWAIT by default
         qr_process = subprocess.Popen(qrCommand.split(), stderr=subprocess.STDOUT)
         self.qr_pid = qr_process.pid
      except:
         self.logger.error("Failed to start the apogeeqr_IDL process")
         traceback.print_exc()


   def stopQuickReduce(self):
      '''If a quickReduce IDL process already exists - just kill it (for now)'''
      # first send a QUIT command to the program
      for s in self.qrSources:
         s.sendLine('QUIT')

      if self.qr_pid == 0:
         p1=subprocess.Popen(['/bin/ps'],stdout=subprocess.PIPE)
         # apqr_wrapper is the name of the program ran when starting IDL (in apogeeqr.cfg)
         processName = (self.config.get('apogeeqr','qlCommandName')).split()
         if '-e' in processName:
            # look at the name of the program called when IDL is started
            pos = processName[processName.index('-e')+1]
         else:
            pos = 0
         args = ['grep',processName[pos]]
         p2=subprocess.Popen(args,stdin=p1.stdout,stdout=subprocess.PIPE)
         output=p2.communicate()[0]
         p2.kill()
         p1.kill()
         if len(output) > 0:
            # process exists -> kill it
            self.qr_pid = output.split()[0]

      os.kill(self.qr_pid,signal.SIGKILL)
      killedpid, stat = os.waitpid(self.qr_pid, os.WNOHANG)
      if killedpid == 0:
         # failed in killing old IDL process
         # print error message here
         self.logger.warn("Unable to kill existing apogeeql_IDL process %s" % self.qr_pid)
      else:
         self.qr_pid = 0

   def periodicStatus(self):
      '''Run some command periodically'''
      self.callCommand('update')
      reactor.callLater(int(self.config.get(self.name, 'updateInterval')), self.periodicStatus)

   def periodicDisksStatus(self):
      '''Check on the disk free space '''
      self.callCommand('checkdisks')
      reactor.callLater(int(self.config.get(self.name, 'diskAlarmInterval')), self.periodicDisksStatus)

   def sendAliveTest(self):
      '''Run some command periodically'''
      self.watchDogStatus = False
      s.sendLine('PING')
      reactor.callLater(int(self.watchDogTimer), self.isApqlAlive)

   def isApqlAlive(self):
      '''Run some command periodically'''
      if not self.watchDogStatus:
         # APQL is not reponding - > restart it
         self.logger.warn("APOGEEQL -> IDL_QL process not responding ... restarting ...")
         self.startQuickLook()

   def connectionMade(self):
      '''Runs this after connection is made to the hub'''
      #
      # Schedule an update.
      #
      # reactor.callLater(3, self.periodicStatus)

   def appendFitsKeywords(self, filename):
      '''make a copy of the input FITS file with added keywords'''

      # we need to form the paths where the file can be found and written
      # expecting something like: apRaw-DDDDXXXX-RRR.fits
      # where:
      # DDDD is a 4 digit day number starting with 0000 for Jan 1, 2011 which is MJD=55562
      # XXXX is a 4 digit exposure number of the day, starting with 0001 
      # RRR is a 3 digit read number starting with 001 for the first read
      #
      # We are using the SDSS version of MJD which is MJD+0.3 days, which
      # rolls over at 9:48 am MST.
      #
      # The directories from the ICS are named from the day-of-survey 4 digits
      # The directories for the quicklook and archive are by the 5 digits MJD
      #
      res=filename.split('-')
      indir=self.ics_datadir
      outdir = self.datadir

      try:
         indir  = os.path.join(self.ics_datadir,res[1][:4])
         mjd = int(res[1][:4])+int(self.startOfSurvey)
         outdir = os.path.join(self.datadir,str(mjd))
      except:
         raise RuntimeError( "The filename doesn't match expected format (%s)" % (filename)) 

      # create directory if it doesn't exist
      if not os.path.isdir(outdir):
         os.mkdir(outdir)
      outFile = os.path.join(outdir, filename)
      filename = os.path.join(indir, filename)

      # Since the Java nom.tam.fits library used by the ICS has an incompatible
      # checksum calculation with pyfits, we're rolling our own here

      # first extract the value of the checksum from the fits header (pyfits.getval removes
      # any checksum or datasum keywords)
      f=open(filename,'rb')
      checksum = None
      # only read the first 72 lines (which should be the whole header plus padding)
      for p in range(72):
          line = f.read(80)
          if line[0:8] == 'END     ':
              break
          if line[0:8] == 'CHECKSUM':
              checksum = line[11:27]
              cs_comment = line[33:80]

      f.close()
      # don't touch the data, which is supposed to be uint16s.
      hdulist = pyfits.open(filename, do_not_scale_image_data=True, uint16=True)
      if checksum != None:
          # validate the value of the checksum found (corresponding to DATASUM in pyfits)
          # calulate the datasum
          ds = hdulist[0]._calculate_datasum('standard')

          # add a new CHECKSUM line to the header (pyfits.open removes it) with same comment
          hdulist[0].header.update("CHECKSUM",'0'*16, cs_comment)

          # calulate a new checksum
          cs=hdulist[0]._calculate_checksum(ds,'standard')
          if cs != checksum:
              self.logger.error("CHECKSUM Failed for file %s" % (filename))
          else:
              # self.logger.info("CHECKSUM checked ok")
              f = open(os.path.join(outdir,Apogeeql.actor.delFile),'a')
              f.write(filename+'\n')
              f.close()

      # force these to be ints:
      # As of August 2013, the ICS writes them both as floats, but the
      # FITS standard wants them to be ints.
      bscale = int(hdulist[0].header.get('BSCALE',1))
      bzero = int(hdulist[0].header.get('BZERO',32768))
      del hdulist[0].header['BSCALE']
      del hdulist[0].header['BZERO']
      hdulist[0].header.update('BSCALE',bscale,after='GCOUNT')
      hdulist[0].header.update('BZERO',bzero,after='BSCALE')
      
      hdulist[0].header.update('TELESCOP' , 'SDSS 2-5m')
      hdulist[0].header.update('FILENAME' ,outFile)
      hdulist[0].header.update('EXPTYPE' ,self.expType)
      
      # get the calibration box status
      lampqrtz, lampune, lampthar, lampshtr, lampcntl = self.getCalibBoxStatus()
      hdulist[0].header.update('LAMPQRTZ',lampqrtz, 'CalBox Quartz Lamp Status')
      hdulist[0].header.update('LAMPUNE',lampune, 'CalBox UNe Lamp Status')
      hdulist[0].header.update('LAMPTHAR',lampthar, 'CalBox ThArNe Lamp Status')
      hdulist[0].header.update('LAMPSHTR',lampshtr, 'CalBox Shutter Lamp Status')
      hdulist[0].header.update('LAMPCNTL',lampcntl, 'CalBox Controller Status')

      """
      # guider i refractionCorrection=1.00000
      refraction = Apogeeql.actor.models['guider'].keyVarDict['refractionCorrection'][0]
      refraction = numpy.nan_to_num(refraction)
      hdulist[0].header.update('REFRACOR',refraction, 'guider refractionCorrection')
      """

      # guider i seeing=2.09945
      seeing = Apogeeql.actor.models['guider'].keyVarDict['seeing'][0]
      try:
          seeing=float(seeing)
          if seeing == float('NaN'):
              seeing=0.0
      except:
          seeing=0.0
      hdulist[0].header.update('SEEING',seeing, 'RMS seeing from guide fibers')

      # starttime is MJD in seconds      
      time_string = hdulist[0].header['DATE-OBS']
      starttime = Time(time_string,format='isot',scale='utc').mjd * 86400.0

      exptime = hdulist[0].header['exptime']

      cards=[]
      cards.extend(actorFits.mcpCards(self.models, cmd=self.bcast))
      cards.extend(actorFits.tccCards(self.models, cmd=self.bcast))
      cards.extend(actorFits.plateCards(self.models, cmd=self.bcast))

      for name, val, comment in cards:
          try:
              hdulist[0].header.update(name, val, comment)
          except:
              self.logger.warn('text="failed to add card: %s=%s (%s)"' % (name, val, comment))

      # repair (if possible) any problems with the header (mainly OBSCMNT too long)
      hdulist.verify('fix')
      hdulist.writeto(outFile, clobber=True, output_verify='warn', checksum=True)
      os.chmod(outFile,0o444) # all read only
      return outFile, starttime, exptime


   def getPlPlugMapM(self, session, cartridgeId, plateId, pointingName):
       """Return the plPlugMapM given a plateId and pointingName"""

       try:
           pm = session.query(PlPlugMapM).join(Plugging,Plate,Cartridge,ActivePlugging).\
                   filter(Plate.plate_id==plateId).\
                   filter(Cartridge.number==cartridgeId).\
                   filter(PlPlugMapM.pointing_name==pointingName).order_by(PlPlugMapM.fscan_mjd.desc()).\
                   order_by(PlPlugMapM.fscan_id.desc()).one()
       except sqlalchemy.orm.exc.NoResultFound:
           raise RuntimeError("NO plugmap from for plate %d" % (plateId))
       except sqlalchemy.orm.exc.MultipleResultsFound:
           # raise RuntimeError, ("More than one plugmap from for plate %d" % (plateId))
           # use thae last entry hoping all is well
           pm = session.query(PlPlugMapM).join(Plugging,Plate,Cartridge,ActivePlugging).\
                   filter(Plate.plate_id==plateId).\
                   filter(Cartridge.number==cartridgeId).\
                   filter(PlPlugMapM.pointing_name==pointingName).order_by(PlPlugMapM.fscan_mjd.desc()).\
                   order_by(PlPlugMapM.fscan_id.desc())
           pm=pm[0]

       return pm

   def makeApogeePlugMap(self, plugmap, newfilename):
       """Return the plPlugMapM given a plateId and pointingName"""

       from sqlalchemy import and_

       # replace the \r with \n to get yanny to parse the text properly
       import re
       data = re.sub("\r", "\n", plugmap.file)

       # append to the standard plPlugMap to add 2mass_style and J, H, Ks mags
       par = yanny.yanny()
       par._contents = data
       par._parse()
       p0=par
       # update the definition of PLUGMAPOBJ
       pos=0
       for t in p0.tables():
          if t=='PLUGMAPOBJ':
             p0['symbols']['struct'][pos] = (p0['symbols']['struct'][pos]).replace('secTarget;','secTarget;\n char tmass_style[30];')
             break
          else:
             pos+=1
       
       p0['symbols']['PLUGMAPOBJ'].append('tmass_style')
       p0['PLUGMAPOBJ']['tmass_style']=[]
       for i in range(p0.size('PLUGMAPOBJ')):
          p0['PLUGMAPOBJ']['tmass_style'].append('-')

       # get the needed information from the plate_hole 
       ph = self.mysession.query(Fiber).join(PlateHole).\
             filter(Fiber.pl_plugmap_m_pk==plugmap.pk).order_by(Fiber.fiber_id).\
             values('fiber_id','tmass_j','tmass_h','tmass_k','apogee_target1','apogee_target2')
                         
       # we'll use the target1 and target2 to define the type of target
       # these are 32 bits each with each bit indicating a type
       skymask = 16
       hotmask = 512
       extmask = 1024
       starmask = skymask | hotmask

       # loop through the list and update the PLUGMAPOBJ
       tmass_style = 'Unknown'
       for fid, j_mag, h_mag, k_mag, t1, t2 in ph:
           count = p0['PLUGMAPOBJ']['fiberId'].count(fid)
           if count >= 1:
              # we have more than one entry for this fiberId -> get the APOGEE 
              ind = -1
              for i in range(count):
                  pos = p0['PLUGMAPOBJ']['fiberId'][ind+1:].index(fid)
                  ind = pos+ind+1
                  if p0['PLUGMAPOBJ']['spectrographId'][ind] == 2:
                      break

                 # print "fid=%d    t1=%d   t2=%d" % (fid,t1,t2)
                 # only modify the fibers for APOGEE (2) that are not sky fibers
              if p0['PLUGMAPOBJ']['spectrographId'][ind] == 2 and p0['PLUGMAPOBJ']['objType'][ind] != 'SKY':
                  if not (j_mag and h_mag and k_mag):
                      #cmd.warn('text="some IR mags are bad: j=%s h=%s k=%s"' % (j_mag, h_mag, k_mag))
                      logging.warn('text="some IR mags are bad: j=%s h=%s k=%s"' % (j_mag, h_mag, k_mag))
                  p0['PLUGMAPOBJ']['mag'][ind][0] = j_mag if j_mag else 0.0
                  p0['PLUGMAPOBJ']['mag'][ind][1] = h_mag if h_mag else 0.0
                  p0['PLUGMAPOBJ']['mag'][ind][2] = k_mag if k_mag else 0.0
                  p0['PLUGMAPOBJ']['tmass_style'][ind] = tmass_style
                  if (t2 & skymask) > 0:
                     p0['PLUGMAPOBJ']['objType'][ind] = 'SKY'
                  elif (t2 & hotmask) > 0:
                     p0['PLUGMAPOBJ']['objType'][ind] = 'HOT_STD'
                  elif (t1 & extmask) > 0:
                     p0['PLUGMAPOBJ']['objType'][ind] = 'EXTOBJ'
                  elif (t2 & starmask) == 0 and (t1 & extmask) ==0:
                    p0['PLUGMAPOBJ']['objType'][ind] = 'STAR'
        
       # delete file if it already exists
       if os.path.isfile(newfilename):
          os.remove(newfilename)

       p0.write(newfilename)

       # write a copy to the archive directory
       # define the current mjd archive directory to store the plPlugMapA file
       mjd = astroMJD.mjdFromPyTuple(time.gmtime())
       fmjd = str(int(mjd + 0.3))
       arch_dir = os.path.join(self.archive_dir, fmjd)
       if not os.path.isdir(arch_dir):
           os.mkdir(arch_dir, 0o0775)

       res=os.path.split(newfilename)
       archivefile = os.path.join(arch_dir,res[1])
       # delete file if it already exists
       if os.path.isfile(archivefile):
          os.remove(archivefile)
       p0.write(archivefile)
 
       return 

   def getCalibBoxStatus(self):
       """Insert a new row in the platedb.exposure table """

       # get the lamp status from the actor
       lampqrtz = Apogeeql.actor.models['apogeecal'].keyVarDict['calSourceStatus'][0]
       lampune  = Apogeeql.actor.models['apogeecal'].keyVarDict['calSourceStatus'][1]
       lampthar = Apogeeql.actor.models['apogeecal'].keyVarDict['calSourceStatus'][2]
       lampshtr = Apogeeql.actor.models['apogeecal'].keyVarDict['calShutter'][0]
       lampcntl = Apogeeql.actor.models['apogeecal'].keyVarDict['calBoxController'][0]

       return lampqrtz, lampune, lampthar, lampshtr, lampcntl

   def completeUTR(self,filebase=None):
       """Verifies that all of the UTR files where copied from the ICS"""

       # expecting something like: apRaw-DDDDXXXX
       # make sure a filebase was passed
       if not filebase: 
           return

       res=filebase.split('-')
       try:
          indir  = os.path.join(self.ics_datadir,res[1][:4])
          mjd = int(res[1][:4])+int(self.startOfSurvey)
          outdir = os.path.join(self.datadir,str(mjd))
       except:
          raise RuntimeError( "The filename doesn't match expected format (%s)" % (filename)) 
 
       lst = glob.glob(os.path.join(indir,filebase+'*.fits'))
       lst.sort()
       count=0
       for infile in lst:
           # check that the file exists in the outdir
           outfile = os.path.join(outdir,os.path.basename(infile))
           if not os.path.exists(outfile):
               # need to annotate and copy the file
               count+=1
               # should we try this or just make a copy without the annotation?
               newfilename, starttime, exptime = self.appendFitsKeywords(os.path.basename(infile))
               if not os.path.exists(outfile):
                   # copy the file if appendFitsKeywords did not work
                   shutil.copy(infile,outdir)
 
       if count > 0:
           self.bcast.warn('text="%s had %d missing UTR"' % (filebase,count))
           self.logger.info('APOGEEQL -> had %d missing UTR' % (count))
       return



#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def kill_handler(signum, frame):
   # send the quit message to the idl routines
   print 'Signal handler called with signal', signum
   logging.info("Stopping the IDL processes ... ")
   # stopping apqr_wrapper

   p1=subprocess.Popen(['/bin/ps'],stdout=subprocess.PIPE)
   # apqr_wrapper is the name of the program ran when starting IDL (in apogeeqr.cfg)
   processName = (self.config.get('apogeeql','qrCommandName')).split()
   if '-e' in processName:
       # look at the name of the program called when IDL is started
       pos = processName[processName.index('-e')+1]
   else:
       pos = 0
   args = ['grep',processName[pos]]
   p2=subprocess.Popen(args,stdin=p1.stdout,stdout=subprocess.PIPE)
   output=p2.communicate()[0]
   p2.kill()
   p1.kill()
   if len(output) > 0:
       # process exists -> kill it
       qr_pid = output.split()[0]
       os.kill(qr_pid,signal.SIGKILL)
       killedpid, stat = os.waitpid(qr_pid, os.WNOHANG)
       if killedpid == 0:
          # failed in killing old IDL process
          # print error message here
          logging.warn("Unable to kill existing apogeeql_IDL process %s" % qr_pid)


   p1=subprocess.Popen(['/bin/ps'],stdout=subprocess.PIPE)
   # apqr_wrapper is the name of the program ran when starting IDL (in apogeeqr.cfg)
   processName = (self.config.get('apogeeql','qlCommandName')).split()
   if '-e' in processName:
       # look at the name of the program called when IDL is started
       pos = processName[processName.index('-e')+1]
   else:
       pos = 0
   args = ['grep',processName[pos]]
   p2=subprocess.Popen(args,stdin=p1.stdout,stdout=subprocess.PIPE)
   output=p2.communicate()[0]
   p2.kill()
   p1.kill()
   if len(output) > 0:
       # process exists -> kill it
       ql_pid = output.split()[0]
       os.kill(ql_pid,signal.SIGKILL)
       killedpid, stat = os.waitpid(ql_pid, os.WNOHANG)
       if killedpid == 0:
          # failed in killing old IDL process
          # print error message here
          logging.warn("Unable to kill existing apogeeql_IDL process %s" % ql_pid)

   sys.exit()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def main():
   import os, signal

   apogeeql = Apogeeql('apogeeql', 'apogeeql')
   apogeeql.connectQuickLook()
   apogeeql.startQuickLook()
   apogeeql.connectQuickReduce()
   apogeeql.startQuickReduce()
   reactor.callLater(3, apogeeql.periodicStatus)
   reactor.callLater(30, apogeeql.periodicDisksStatus)
   signal.signal(signal.SIGTERM, kill_handler)
   apogeeql.run()

#-------------------------------------------------------------
if __name__ == '__main__':
   try:
       main()
   except Exception,e:
       traceback.print_exc()


