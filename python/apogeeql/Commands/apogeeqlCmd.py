#!/usr/bin/env python

import pdb
import logging
import pprint
import re, sys, time, os
from time import sleep

import opscore.protocols.keys as keys
import opscore.protocols.types as types

from opscore.utility.qstr import qstr
import actorcore.help as help

class apogeeqlCmd(object):

   def __init__(self, actor):
      self.actor = actor

      # Define some typed command arguemetnts
      self.keys = keys.KeysDictionary("apogeeql_apogeeql", (1, 1),
                              keys.Key("actor", types.String(), help="Another actor to command"),
                              keys.Key("cmd", types.String(), help="A command string"),
                              keys.Key("count", types.Int(), help="A count of things to do"))
      #
      # Declare commands
      #
      self.vocab = [
         ('ping', '', self.ping),
         ('status', '', self.status),
         ('update', '', self.update),
         ('checkdisks', '', self.checkDisks),
         ('stopidl', '', self.stopIDL),
         ('startidl', '', self.startIDL),
         ('ql', '<cmd>', self.quicklook),
         ('doSomething', '<count>', self.doSomething),
         ('passAlong', 'actor <cmd>', self.passAlong),
      ]

   def ping(self, cmd):
      '''Query the actor for liveness/happiness.'''
      cmd.finish("text='Present and (probably) well'")

   def status(self, cmd):
      '''Report status and version; obtain and send current data'''
      self.actor.sendVersionKey(cmd)
      self.doStatus(cmd, flushCache=True)

   def update(self, cmd):
      '''Report status'''
      self.doStatus(cmd=cmd)

   def doStatus(self, cmd=None, flushCache=False, doFinish=True):
      '''Report full status'''
      if not cmd:
         cmd = self.actor.bcast

      cmd.inform('rootURL=%s' % (self.actor.rootURL))
      cmd.inform('snrAxisRange=%s,%s' % (self.actor.snrAxisRange[0],self.actor.snrAxisRange[1]))
      # keyStrings = ['text="nothing to say, really"']
      # keyMsg = '; '.join(keyStrings)

      # cmd.inform(keyMsg)
      # cmd.diag('text="still nothing to say"')
      cmd.finish()

   def _checkDisk(self,cmd,space,diskName):
      '''Send appropriate keyword message about current space remaining on diskName.'''
      if space > int(self.actor.warningDiskSpace):
          cmd.inform('%sDiskAlarm=Ok,%d' % (diskName,space))
      elif space > int(self.actor.seriousDiskSpace):
          cmd.warn('%sDiskAlarm=Warning,%d' % (diskName,space))
      elif space > int(self.actor.criticalDiskSpace):
          cmd.warn('%sDiskAlarm=Serious,%d' % (diskName,space))
      else:
          cmd.warn('%sDiskAlarm=Critical,%d' % (diskName,space))
   #...

   def checkDisks(self, cmd=None):
      '''Report disk space remaining on ICS, QL and Arch.'''
      if not cmd:
         cmd = self.actor.bcast

      s=os.statvfs(self.actor.ics_datadir)
      icsSpace = s.f_bsize * s.f_bavail / 1024 / 1024 / 1024  # free space in GB

      s=os.statvfs(self.actor.datadir)
      qlSpace = s.f_bsize * s.f_bavail / 1024 / 1024 / 1024  # free space in GB

      s=os.statvfs(self.actor.archive_dir)
      archSpace = s.f_bsize * s.f_bavail / 1024 / 1024 / 1024  # free space in GB

      cmd.inform('freeDiskSpace=%d,%d,%d' % (icsSpace, qlSpace, archSpace))

      self._checkDisk(cmd,icsSpace,'ics')
      self._checkDisk(cmd,qlSpace,'ql')
      self._checkDisk(cmd,archSpace,'arch')

      cmd.finish()

   def startIDL(self, cmd):
      '''Start a new IDL quicklook process.'''
      if not cmd:
         cmd = self.actor.bcast

      # the startQuickLook method will kill an existing IDL process before starting a new one
      keyStrings = ['text="starting new IDL process"']
      keyMsg = '; '.join(keyStrings)
      cmd.inform(keyMsg)
      self.actor.startQuickLook()
      cmd.finish()

   def stopIDL(self, cmd):
      '''Stop the currently running IDL quicklook process if it exists.'''
      if not cmd:
         cmd = self.actor.bcast

      if self.actor.ql_pid  != 0:
         keyStrings = ['text="stopping process %s"' % self.actor.ql_pid]
         keyMsg = '; '.join(keyStrings)
         cmd.inform(keyMsg)
         self.actor.stopQuickLook()
      else:
         keyStrings = ['text="no process to kill"']
         keyMsg = '; '.join(keyStrings)
         cmd.inform(keyMsg)
      cmd.finish()

   def doSomething(self, cmd):
      """ Do something pointless. (for testing) """
      cnt = cmd.cmd.keywords["count"].values[0]
      for i in range(cnt):
         cmd.inform('cnt=%d' % (i))
      cmd.finish()

   def quicklook(self, cmd):
      """ Send a command to the IDL quicklook process."""
      cmdString = cmd.cmd.keywords["cmd"].values[0]
      # pass the command string to the socket to quicklook_main.pro
      # print cmdString
      for s in self.actor.qlSources:
         s.sendLine(cmdString)
      cmd.finish()

   def passAlong(self, cmd):
      """ Pass a command along to another actor. """
      actor = cmd.cmd.keywords["actor"].values[0]
      cmdString = cmd.cmd.keywords["cmd"].values[0]

      cmdVar = self.actor.cmdr.call(actor=actor, cmdStr=cmdString, timeLim=30.0)
      if cmdVar.didFail:
         cmd.fail('text=%s' % (qstr('Failed to pass %s along to %s' % (cmdStr, actor))))
      else:
         cmd.finish()

