import os
import shutil
import subprocess
import h5py

def run(tag_lst, h5_in, h5_out):
    if os.path.exists(h5_out):
        print "Already exists:", h5_out
        return

    shutil.copy(h5_in, h5_out+".tmp")

    h = h5py.File(h5_out+".tmp", "a")
    
    for l in open(tag_lst):
        l = l.strip()
        if not l: continue
        path = "/tag-%d"%int(l)
        if path in h:
            print "removing", path
            del h[path]
        else:
            print "not found:", path

    h.close()
    p = subprocess.Popen(["h5repack", h5_out+".tmp",h5_out], shell=False)
    p.wait()
    os.remove(h5_out+".tmp")
        
if __name__ == "__main__":
    import sys
    tag_lst, h5_in, h5_out = sys.argv[1:]
    run(tag_lst, h5_in, h5_out)
