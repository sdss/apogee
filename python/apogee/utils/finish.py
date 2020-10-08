import time
import datetime

def done(done,host,start) :

    now = datetime.datetime.now()
    print("End: ",now.strftime("%Y-%m-%d %H:%M:%S"))
    print("elapsed: ",time.time()-start)

    if done is not None :
        subprocess.call(['setdone',done])
        print('host', host)
        if host is not None :
            try: os.remove(done+'.'+host)
            except: pass

 
