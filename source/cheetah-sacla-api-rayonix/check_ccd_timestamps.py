import os
import glob
import datetime
import marccd

rayonix_root = "/xustrg0/rayonix/2018A/SFX"

"""
for(f in Sys.glob("*/*.dat")) {
 d=read.table(f,h=T)
 d$time.acq = as.POSIXct(d$time.acq)
 d$time.header = as.POSIXct(d$time.header)
 d$time.save = as.POSIXct(d$time.save)
 d$time.ctime = as.POSIXct(d$time.ctime)
 #print(summary(as.numeric(diff(d$time.acq))))
 df=diff(d$time.acq)
 df=df[diff(d$num)==1]
 print(sprintf("%s %5d %6.2f %.2f %6.2f %6.2f", f, length(df), 1000*mean(df), 1000*sd(df), 1000*max(df), 1000*min(df)))
}

"""

def run(runid):
    img_files = sorted(glob.glob(os.path.join(rayonix_root, str(runid), "data_*.img")))
    img_num = map(lambda f: int(f[f.rindex("_")+1:f.rindex(".")]), img_files)

    acq_time, header_time, save_time, ctime = [], [], [], []
    for f in img_files:
        img = marccd.MarCCD(f)
        acq_time.append(img.acquire_time)
        header_time.append(img.header_time)
        save_time.append(img.save_time)
        ctime.append(datetime.datetime.fromtimestamp(os.path.getctime(f)))

    #for i in xrange(1, len(img_files)):
    #    print img_num[i]-img_num[i-1], acq_time[i]-acq_time[i-1]
        
        
    print "run num file time.acq time.header time.save time.ctime"
    for i in xrange(len(img_files)):
        print "%d %5d %s" %(runid, img_num[i], img_files[i]),
        print acq_time[i].strftime('"%Y-%m-%d %H:%M:%S.%f"'),
        print header_time[i].strftime('"%Y-%m-%d %H:%M:%S.%f"'),
        print save_time[i].strftime('"%Y-%m-%d %H:%M:%S.%f"'),
        print ctime[i].strftime('"%Y-%m-%d %H:%M:%S.%f"')


if __name__ == "__main__":
    import sys
    run(int(sys.argv[1]))
