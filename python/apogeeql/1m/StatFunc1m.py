from socket import *
# Module to send/receive the status command from the 1-meter 
def stat1m(soc1m,ADRS1m):
    cmd1m = 'STATUS'+'\r\n' # Status command for 1meter
    soc1m=socket(AF_INET,SOCK_STREAM)
    soc1m.connect(ADRS1m)
    soc1m.send(cmd1m)       # Send status command to 1meter through socket soc1m
    replystat1m = []        # array to store the Status returned by the 1meter
    while 1:
        ans1m = soc1m.recv(1024)  # recveive replies from the 1meter through soc1m 
        if ans1m[0:7]=='STATUS:':
            replystat1m.append(ans1m) # store the reply in the array declared above 
            if ans1m.find('DONE:')!=-1: break # stop receiving replies from the 1-meter when DONE: is received
        else:
            start1m = ans1m.find('STATUS:')
            end1m   = len(ans1m)
            replystat1m.append(ans1m[start1m:end1m])
            if ans1m.find('DONE:')!=-1: break

    replystring = replystat1m[0].split()
    
    soc1m.close()
    return replystring
