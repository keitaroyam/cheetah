import cPickle as pickle
import time
import numpy
from multiprocessing import Process
import zmq
import json
from yamtbx.dataproc import cbf
import pycbf
import tempfile
import os

nproc = 2
host = "127.0.0.1"

def get_pilatus_data(cbf_data):
    fd, cbftmp = tempfile.mkstemp(suffix="cheetah.cbf", dir="/dev/shm")
    os.write(fd, cbf_data)
    os.close(fd)
    h = pycbf.cbf_handle_struct()
    h.read_file(cbftmp, pycbf.MSG_DIGEST)
    h.require_category("array_data")
    h.find_column("header_contents")
    header = h.get_value()
    #h.require_category("array_data")
    h.find_column("data")
    compression, binary_id, elsize, elsigned, elunsigned, elements, minelement, maxelement, bo, ndimfast, ndimmid, ndimslow, padding = h.get_integerarrayparameters_wdims()
    assert elsize == 4 or elsize == 8
    assert elsigned == 1
    assert ndimslow <= 1
    arr = numpy.fromstring(h.get_integerarray_as_string(), dtype=numpy.int32 if elsize==4 else numpy.int64)
    arr = arr.reshape(ndimmid, ndimfast)
    os.remove(cbftmp)
    hdict = {}
    for l in header.splitlines():
        if l.startswith("# Wavelength"):
            hdict["wavelength"] = float(l.split()[2])
        elif l.startswith("# Beam_xy"):
            tmp = l[l.index("(")+1:l.rindex(")")].split(",")
            hdict["beam_center_x"] = float(tmp[0])
            hdict["beam_center_y"] = float(tmp[1])
        elif l.startswith("# Detector_distance"):
            hdict["distance"] = float(l.split()[2])*1000.
        elif l.startswith("# Pixel_size"):
            hdict["pixel_size_x"] = float(l.split()[2])*1000.
    
    return arr, hdict

def worker(widx):
    context = zmq.Context()
    receiver = context.socket(zmq.PULL)
    receiver.connect("tcp://{0}:9999".format(host))
    print "starting", widx
    times = []
    for i in xrange(10000):
        #frames = receiver.recv_multipart(copy = False)
        all_data = receiver.recv_pyobj()
        data, header = get_pilatus_data(all_data["cbf_data"])
        #all_data = map(lambda x: x.bytes, frames)
        print header
        times.append(time.time())
        print "recv:", i, numpy.mean(numpy.diff(times))

if __name__ == "__main__":
    for i in xrange(nproc):
        Process(target=worker, args=(i,)).start()

