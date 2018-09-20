#!/usr/bin/env python
#
# server to communicate with the apql_wrapper on port 10033
#
# This is used to test apql_wrapper.
#    Start this program first in a window (first making sure the actual apogeeql actor is NOT running)
#    Start apql_wrapper in a separate window with "idl -e apql_wrapper -args localhost 1033"
#    You can now type in messages in the format that apql_wrapper expects.
#    You can change the default value for the environment variables pointing to other directories
#       (APQLDATA_DIR, /data-ql/data/)
#       (APQLSPECTRO_DIR, /data-ql/spectro/current/)
#       (APQLARCHIVE_DIR, /data/apogee/archive/)
#       (APQLQUICKRED_DIR, /data/apogee/quickred/)
#
#   To prevent additions to the database, comment out the exec_sql statement in apql_dbinsert.
#   To exit from this program, type in QUIT.
#   Make sure Xvfb is setup and that you "export DISPLAY=:1"
#
#   Expected messages:
#      PING
#      STARTING
#      plugMapInfo=4810,55685,2,/data-ql/plugmaps/plPlugMapA-4810-55685-02.par
#      UTR=/data-ql/data/55810/apRaw-02480036-001.fits,1,6
#      UTR=/data-ql/data/55810/apRaw-02480036-001.fits,2,6
#      UTR=/data-ql/data/55810/apRaw-02480036-001.fits,3,6
#      UTR=/data-ql/data/55810/apRaw-02480036-001.fits,4,6
#      UTR=/data-ql/data/55810/apRaw-02480036-001.fits,5,6
#      UTR=/data-ql/data/55810/apRaw-02480036-001.fits,6,6
#      UTR=DONE
#

import socket
import sys

# create a socket
s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s.setsockopt(socket.SOL_SOCKET,socket.SO_REUSEADDR,1)

# associate the socket with a port
host = 'localhost'
# port = int(sys.argv[1])
port = 10033
s.bind((host, port))

# accept "call" from client
s.listen(1)
conn, addr = s.accept()
print 'client is at', addr

while 1:
   s = raw_input('apql >>> ')
  
   if s.upper() == 'QUIT' or s.upper() == 'EXIT':
      break

   conn.sendall(s+'\n')

   # wait for response from apql_wrapper
   data = conn.recv(1024)
   print 'Received: %s' % data



# close the connection
conn.close()
