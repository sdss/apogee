from socket import *
import numpy as np
import pyfits as pyf
import sys
import os
import datetime
import CallEditFitsFunc
import StatFunc1m
#######################
# apogee-ql client 1  #
# to ICS              #
###############################################################################
# This will be the client that translates the commands from command1m for the 
# ICS, or allow the user to input commands in the command line and send them
# directly to the ICS

ICShostname = 'apogee-ics'         #Host name of ICS 
ICSport     = 33221                #Port through which to connect to ICS
ICSADRS = (ICShostname,ICSport)    #Combine IP and port into Address for server
ICSsoc = socket(AF_INET, SOCK_STREAM) #Create socket to server, named ICSsoc
ICSsoc.connect(ICSADRS)               #Connect to ICS socket at Address above
###############################################################################

#######################
# apogee-ql client 2  #
# to command1m        #
###############################################################################
# This client will send the status command to the 1-meter and receive the reply
IP1m   = '10.75.0.19'                #IP Address to connect to command1m 
port1m = 1053                        #Port to connect to command1m
ADRS1m = (IP1m,port1m)           #Combine IP and port into Address for server
soc1m = socket(AF_INET, SOCK_STREAM) #Create socket to ICS, named 1msoc
#soc1m.connect(ADRS1m)                #Connect to ICS socket at Address abov
###############################################################################

####################
# apgoee-ql server #
###############################################################################
serverIP   = ''
serverport = 1050
serverADRS = (serverIP,serverport)
svr = socket(AF_INET, SOCK_STREAM)
svr.bind(serverADRS)
svr.listen(5)
###############################################################################


user = 'apo.1m'	
def commandICS(command,cmdID,ICSsoc):
    if command.find('status')!=-1:
        cmd = user+' '+str(cmdID + 1)+' '+'status'
        ICSsoc.send(cmd+'\r\n')
        replystatICS = []
        while 1:
            ansICS = ICSsoc.recv(1024)
            print 'ansICS: '+ansICS
            replystatICS.append(ansICS)
            if ansICS.find(':')!= -1 : break
        print replystatICS
        return 0

    if command.find('shutter')!=-1:
        exp = command.split()
        cmd = user+' '+str(cmdID + 1)+' '+'shutter'+' '+exp[1]
        print 'cmd: '+cmd
        ICSsoc.send(cmd+'\r\n')
        replystatICS = []
        while 1:
            ansICS = ICSsoc.recv(1024)
            replystatICS.append(ansICS)
            if ansICS.find(':')!= -1 : break
        print replystatICS
        return 0

    if command.find('dither')!=-1:
        exp = command.split()
        cmd = user+' '+str(cmdID + 1)+' '+'dither namedpos='+exp[1]
        print 'cmd: '+cmd
        ICSsoc.send(cmd+'\r\n')
        replystatICS = []
        while 1:
            ansICS = ICSsoc.recv(1024)
            replystatICS.append(ansICS)
            if ansICS.find(':')!= -1 : break
        print replystatICS
        return 0
    
    if command.find('expose')!=-1:
        exp = command.split()
        cmd = user+' '+str(cmdID + 1)+' '+'expose'+' '+'object=Object'+' '+'nreads='+exp[1]
        print cmd
        ICSsoc.send(cmd+'\r\n') # Send the full expose command to ICS

    if (command.find('dark')== 0):
        exp = command.split()
        if (int(exp[1])>1) and (int(exp[1])<100):
            cmd = user+' '+str(cmdID+1)+' '+'expose'+' '+'object=Dark'+' '+'nreads='+exp[1]
            print cmd
        #exposecmd = raw_input("Enter full expose command(i.e. apo.1m 0001 expose object=Dark nreads=3)\r\n") 
            ICSsoc.send(cmd+'\r\n') # Send the full expose command to ICS
    #ICSsoc.send(cmd+'\r\n') # Send the full expose command to ICS

    while 1:
        ICSans = ICSsoc.recv(1024) # The replies will be received in socket ICS.soc
        print 'ICSans: ' + ICSans
        if ICSans.find('dayNumber=')!= -1:
            DayStart  = ICSans.find('dayNumber=')+len('dayNumber=')
            DayEnd    = len(ICSans)
            num=int(ICSans[DayStart:DayEnd].strip('\n'))
            DayNumber = str(num).zfill(4)
            #DayNumber = '0'+ICSans[DayStart:DayEnd].strip('\n')
        if ICSans.find('exposureWroteSummary')!=-1:
            statstring = StatFunc1m.stat1m(soc1m,ADRS1m)
            #STATUS = 'STATUS:  6.1720166  1.0863463 0.000000  6.1720166  1.0863463 1950.000000 56133.735195 6.845499 331.954057 20.567121 -94.000000 1 0 0 0 0 1 1 1 1 0.000000 0.000000 0.000000 0 0.000000 0 0 0 0 0 0 0.000000 0.000000 0.000000 -7351 -5940 -7414 18.333333 20.555556 11.111111 23 0 1 512 -0.450213 0.450213 0.000000 0.000000 3.128354 0.422420 0.422420 -1500.000000 -100.000000 -1.532399 0.000000 -846777 -5235000 0 11 1 c:\tocc\071027.mo'
            #statstring = STATUS.split()
            current_obs_ra  = float('%8.7f'%(float((statstring[1]))*(180/np.pi))) # Current 1m RA
            current_obs_dec = float('%8.7f'%(float((statstring[2]))*(180/np.pi)))  # Current 1m DEC
            current_obs_az  = float(statstring[9])                         # Current 1m az
            current_obs_alt = float(statstring[10])                        # Current 1m alt
            current_epoch = float(statstring[6])
            StartFileName=ICSans.find('exposureWroteFile=')+len('exposureWroteFile=')
            EndFileName=ICSans.find('";')
            CurrentFileName = ICSans[StartFileName:EndFileName].strip('"')
            # CallEditFitsFunc is a python module I created that moves the FITS files 
            # and edits their Headers as they become available.
            FITSTAT=CallEditFitsFunc.movefits(CurrentFileName,DayNumber,current_obs_ra,current_obs_dec,current_obs_az,current_obs_alt,current_epoch)
        if ICSans.find('exposureState=done')!= -1 :
            #ICSCOMMADSTAT = 'DONE!'
            break
        #else:
            #ICSCOMMANDSTAT = 'WORKING...'
    #return ICSCOMMANDSTAT
    name=CurrentFileName.split('-')
    return int(name[1])
# Main program block here to accept commands!
# Define a list of acceptable commands for the ICS, for now just status and expose
CMDLIST=['status','expose','dark','shutter','dither','done']
cmdID=0
nread=0
while 1:
  client,addr = svr.accept()
  command = client.recv(1024)
  #command = raw_input("Enter command for ICS(i.e. expose Dark 3):\r\n")
  testcmd = command.split()
  print 'received command' + command , nread
  print len(testcmd)
  if len(testcmd)>0:
    if testcmd[0] == 'done':
        ICSsoc.close()
        client.close()
        break
    if testcmd[0] in CMDLIST:
        print 'commandICS: '+command
        ret=commandICS(command,cmdID,ICSsoc)
        print 'returned: ',ret
        print 'sending done'
        client.send('done '+str(ret)+'\n')
        cmdID = cmdID + 1
    else:
        print 'Not Acceptable command'
  nread = nread + 1
  sys.stdout.flush()
svr.close()
print 'No longer taking commands'
