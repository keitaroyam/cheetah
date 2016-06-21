import os
import h5py
import json
import struct
import numpy
import pyublas
import math
import time
from multiprocessing import Process, Queue

from cheetah_ext import CheetahSpot, CheetahSinglePanel
from yamtbx.dataproc import cbf
from yamtbx.dataproc.XIO import XIO

import libtbx.phil
from libtbx import easy_mp

def make_result(retin, f_trans=None):
    ret = dict(spots=[], lookup={})

    for r in retin:
        x, y = r.peak_com_x, r.peak_com_y
        if f_trans: x, y = f_trans(x,y)
        ret["spots"].append((x, y, r.peak_snr, r.peak_com_res))

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

def worker(wrk_num, header, imgfiles, queue, cut_roi, algorithm, binning,
           ADCthresh, MinSNR, MinPixCount, MaxPixCount, LocalBGRadius, MinPeakSeparation, d_min, d_max):
    """
    The code taken from yamtbx/dataproc/myspotfinder/command_line/spot_finder_backend.py (New BSD License)
    """

    cheetah = CheetahSinglePanel()
    cheetah.set_params(ADCthresh=ADCthresh,
                       MinSNR=MinSNR,
                       MinPixCount=MinPixCount,
                       MaxPixCount=MaxPixCount,
                       LocalBGRadius=LocalBGRadius,
                       MinPeakSeparation=MinPeakSeparation,
                       Dmin=d_min,
                       Dmax=d_max)

    for f in imgfiles:
        startt = time.time()
        if f.endswith(".cbf"):
            data, ndimfast, ndimmid = cbf.load_minicbf_as_numpy(f)
            data[data<0] = 0
            data = data.reshape(ndimmid, ndimfast).astype(numpy.uint32)
        else:
            im = XIO.Image(f)
            data = numpy.array(im.getData(), dtype=numpy.uint16).reshape(im.header["Height"], im.header["Width"])

        result = cheetah_worker(header, data, algorithm, cut_roi, cheetah, binning)
        result["frame"] = f
        eltime = time.time()-startt
        #print "Wrkr%3d %6d done in %.2f msec " % (wrk_num, frame, eltime*1.e3)
        print "%s %6d %.2f" %(f, len(result["spots"]), eltime*1.e3)
        queue.put(result)
    #print "Worker %d finished." % wrk_num
# worker()

def run(opts, args):
    assert opts.algorithm in (6, 8)

    nframes = len(args)
    first_img = args[0] # Assumes all images are with the same conditions
    im = XIO.Image(first_img)

    header = dict(beam_center_x=im.header["BeamX"]/im.header["PixelX"],
                  beam_center_y=im.header["BeamY"]/im.header["PixelX"],
                  pixel_size_x=im.header["PixelX"],
                  distance=im.header["Distance"],
                  wavelength=im.header["Wavelength"])
    print "#%s"%header

    ofs = open(opts.output, "w")

    ranges = map(tuple, numpy.array_split(numpy.arange(nframes), opts.nproc))
    queue = Queue()
    pp = []
    #print "frame nspots eltime"
    for i, rr in enumerate(ranges):
        if not rr: continue
        p = Process(target=worker, args=(i, header, map(lambda x: args[x], rr), queue, opts.cut_roi, opts.algorithm, opts.binning,
                                         opts.ADCthresh, opts.MinSNR, opts.MinPixCount, opts.MaxPixCount,
                                         opts.LocalBGRadius, opts.MinPeakSeparation, opts.d_min, opts.d_max))
        p.start()
        pp.append(p) 

    results = {}

    while any(map(lambda p: p.is_alive(), pp)):
        while not queue.empty():
            ret = queue.get()
            results[ret["frame"]] = ret

    for p in pp: p.join()

    if opts.gen_adx:
        for frame in sorted(results):
            ret = results[frame]
            adx_out = open(os.path.basename(frame)+".adx", "w")
            for x,y,snr,d in ret["spots"]: adx_out.write("%6d %6d %.2e\n" % (x,y,snr))
            adx_out.close()


    ofs.write("file nspots total_snr\n")
    for frame in sorted(results):
        ret = results[frame]
        n_spots = len(ret["spots"])
        total_snr = sum(map(lambda x: x[2], ret["spots"]))
        ofs.write("%s %6d %.3e\n"%(frame, n_spots, total_snr))
    
if __name__ == "__main__":
    import sys
    import optparse

    parser = optparse.OptionParser(usage="usage: %prog [options] .img or .cbf files...")

    parser.add_option("--nproc", action="store", dest="nproc", type=int, default=1)
    parser.add_option("--cut-roi", action="store_true", dest="cut_roi")
    parser.add_option("--output", action="store", dest="output", default="cheetah.dat")
    parser.add_option("--gen-adx", action="store_true", dest="gen_adx")

    parser.add_option("--dmin", action="store", dest="d_min", type=float, default=5)
    parser.add_option("--dmax", action="store", dest="d_max", type=float, default=30)
    parser.add_option("--adc-threshold", action="store", dest="ADCthresh", type=float, default=5)
    parser.add_option("--min-snr", action="store", dest="MinSNR", type=float, default=8)
    parser.add_option("--min-pixcount", action="store", dest="MinPixCount", type=int, default=3)
    parser.add_option("--max-pixcount", action="store", dest="MaxPixCount", type=int, default=40)
    parser.add_option("--local-bgradius", action="store", dest="LocalBGRadius", type=float, default=2)
    parser.add_option("--min-peaksep", action="store", dest="MinPeakSeparation", type=float, default=0)
    parser.add_option("--algorithm", action="store", dest="algorithm", type=int, default=8)
    parser.add_option("--binning", action="store", dest="binning", type=int, default=1)


    opts, args = parser.parse_args(sys.argv[1:])
    run(opts, args)
