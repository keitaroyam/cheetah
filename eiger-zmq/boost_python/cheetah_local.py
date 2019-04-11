import os
import h5py
import json
import struct
import numpy
import pyublas
import math
import time
import Queue
from multiprocessing import Process

from cheetah_ext import CheetahSpot, CheetahSinglePanel
from yamtbx.dataproc import eiger
import imageconv_ext 
import config_manager

import libtbx.phil

def make_result(retin, f_trans=None):
    ret = dict(spots=[], lookup={})

    for r in retin:
        x, y = r.peak_com_x, r.peak_com_y
        if f_trans: x, y = f_trans(x,y)
        ret["spots"].append((y, x, r.peak_snr, r.peak_com_res))

    ret["lookup"]["cheetah"] = tuple(range(len(retin)))

    return ret
# make_distl_stats_dict()

def cheetah_worker(header, data, algorithm, cut_roi, cheetah, binning=1):
    beamx, beamy = header["beam_center_x"], header["beam_center_y"]
    pixel_size = header["pixel_size_x"]
    # Cut ROI
    if cut_roi:
        r_max = header["distance"] * math.tan(2.*math.asin(header["wavelength"]/2./params.distl.res.outer)) / header["pixel_size_x"]
        roidata = data[max(0,beamy-r_max):min(beamy+r_max,data.shape[0]-1), max(0,beamx-r_max):min(beamx+r_max,data.shape[1]-1)]
        if roidata.shape[0] < r_max and beamy-r_max < 0: beamy = roidata.shape[0] - r_max
        else: beamy = r_max
        if roidata.shape[1] < r_max and beamx-r_max < 0: beamx = roidata.shape[1] - r_max
        else: beamx = r_max
        print "ROI cut:", roidata.shape
    else:
        roidata = data

    # Binning
    if binning > 1:
        roidata, f_xyconv, f_xyconvinv = software_binning(data, binning)
        (beamx, beamy) = f_xyconv(beamx, beamy)
        pixel_size *= binning
    else:
        f_xyconvinv = None

    if roidata.dtype == numpy.uint32: cheetah_run = cheetah.run_uint32 
    elif roidata.dtype == numpy.uint16: cheetah_run = cheetah.run_uint16

    ret = cheetah_run(roidata.shape[1], roidata.shape[0], roidata,
                      float(beamx), float(beamy),
                      float(header["wavelength"]), float(header["distance"]), float(pixel_size),
                      algorithm=algorithm)


    result = make_result(ret, f_xyconvinv)
    return result
# cheetah_worker()

def worker(wrk_num, h5in, header, framerange, queue, cut_roi, algorithm, params_str):
    """
    The code taken from yamtbx/dataproc/myspotfinder/command_line/spot_finder_backend.py (New BSD License)
    """

    master_params = libtbx.phil.parse(config_manager.master_params_str)
    working_params = master_params.fetch(sources=[libtbx.phil.parse(params_str)])
    params = working_params.extract()

    cheetah = CheetahSinglePanel()
    cheetah.set_params(ADCthresh=params.cheetah.ADCthresh,
                       MinSNR=params.cheetah.MinSNR,
                       MinPixCount=params.cheetah.MinPixCount,
                       MaxPixCount=params.cheetah.MaxPixCount,
                       LocalBGRadius=params.cheetah.LocalBGRadius,
                       MinPeakSeparation=params.cheetah.MinPeakSeparation,
                       Dmin=params.distl.res.outer,
                       Dmax=params.distl.res.inner)

    for frame in framerange:
        startt = time.time()
        data = eiger.extract_data(h5in, frame, apply_pixel_mask=False)
        data[data<0] = 0
        data = data.astype(numpy.uint32)
        result = cheetah_worker(header, data, algorithm, cut_roi, cheetah,
                                params.cheetah.binning)
        result["frame"] = frame
        eltime = time.time()-startt
        #print "Wrkr%3d %6d done in %.2f msec " % (wrk_num, frame, eltime*1.e3)
        print "%6d %6d %.2f" %(frame, len(result["spots"]), eltime*1.e3)
        queue.put(result)
    #print "Worker %d finished." % wrk_num
# worker()

def run(opts, args):
    assert opts.algorithm in (6, 8)

    h5in = args[0]
    h5 = h5py.File(h5in, "r")
    nframes = 0
    for k in sorted(h5["/entry/data"].keys()):
        try:
            nframes += h5["/entry/data"][k].shape[0]
        except KeyError:
            break

    header = dict(beam_center_x=h5["/entry/instrument/detector/beam_center_x"].value,
                  beam_center_y=h5["/entry/instrument/detector/beam_center_y"].value,
                  pixel_size_x=h5["/entry/instrument/detector/x_pixel_size"].value*1.e3,
                  distance=h5["entry/instrument/detector/detector_distance"].value*1.e3,
                  wavelength=h5["entry/instrument/beam/incident_wavelength"].value)
    print "#%s"%header

    ranges = map(tuple, numpy.array_split(numpy.arange(nframes)+1, opts.nproc))
    queue = Queue.Queue()
    pp = []
    print "frame nspots eltime"
    for i, rr in enumerate(ranges):
        if not rr: continue
        p = Process(target=worker, args=(i, h5in, header, rr, queue, opts.cut_roi, opts.algorithm, opts.params))
        p.start()
        pp.append(p) 

    for p in pp: p.join()
    
    n_spots = []
    while not queue.empty():
        ret = queue.get()
        n_spots.append((ret["frame"], len(ret["spots"])))
    
    n_spots.sort(key=lambda x:x[0])
    #print "frame nspots"
    #for frame, nsp in n_spots: print "%6d %6d"%(frame, nsp)

    
if __name__ == "__main__":
    import sys
    import optparse

    parser = optparse.OptionParser(usage="usage: %prog [options] coordinates...")

    parser.add_option("--nproc", action="store", dest="nproc", type=int, default=16)
    parser.add_option("--cut-roi", action="store_true", dest="cut_roi")
    parser.add_option("--algorithm", action="store", dest="algorithm", type=int, default=8)
    parser.add_option("--params", action="store", dest="params", type=str, default="")


    opts, args = parser.parse_args(sys.argv[1:])

    run(opts, args)
