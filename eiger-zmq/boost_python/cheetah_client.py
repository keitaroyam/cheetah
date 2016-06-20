import os
import zmq
import lz4
import bitshuffle
import json
import struct
import numpy
import pyublas
import math
import time
import traceback
#from multiprocessing import Process
import subprocess

from cheetah_ext import CheetahSpot, CheetahSinglePanel
from yamtbx.dataproc.myspotfinder import shikalog
import imageconv_ext 
import config_manager

import libtbx.phil

def read_eiger_stream_data(frames, bss_job_mode=4):
    """
    The code taken from yamtbx/dataproc/eiger.py (New BSD License)
    """
    if len(frames) != 5:
        return None, None

    header = json.loads(frames[0].bytes)
    for i in (1,3,4): header.update(json.loads(frames[i].bytes))

    if header.get("bss_job_mode", 4) != bss_job_mode:
        return None, None

    dtype = header["type"]
    shape = header["shape"][::-1]

    if dtype in ("int32","uint32"): byte = 4
    elif dtype in ("int16","uint16"): byte = 2
    else: raise RuntimeError("Unknown dtype (%s)"%dtype)

    size = byte*shape[0]*shape[1]

    if header["encoding"] == "lz4<":
        data = lz4.loads(struct.pack('<I', size) + frames[2].bytes)
        data = numpy.fromstring(data, dtype=dtype).reshape(shape)
        assert data.size * data.dtype.itemsize == size
    elif header["encoding"] == "bs32-lz4<":
        data = frames[2].bytes
        blob = numpy.fromstring(data[12:],dtype=numpy.uint8)
        # blocksize is big endian uint32 starting at byte 8, divided by element size
        blocksize = numpy.ndarray(shape=(),dtype=">u4", buffer=data[8:12])/4
        data = bitshuffle.decompress_lz4(blob, shape, numpy.dtype(dtype), blocksize)
        data = data.reshape(shape)
    elif header["encoding"] == "bs16-lz4<":
        data = frames[2].bytes
        blob = numpy.fromstring(data[12:],dtype=numpy.uint8)
        data = bitshuffle.decompress_lz4(blob, shape, numpy.dtype(dtype))
        data = data.reshape(shape)
    else:
        RuntimeError("Unknown encoding (%s)"%header["encoding"])

    data[data==2**(byte*8)-1] = 0
    return header, data
# read_eiger_stream_data()

def software_binning(data, binning):
    u = l = 0
    b = r = None
    newshape = [data.shape[0]//binning, data.shape[1]//binning]
    if data.shape[0]%binning > 0:
        res = data.shape[0]%binning
        u = res//2 # upper
        b = res - u # bottom
    if data.shape[1]%binning > 0:
        res = data.shape[1]%binning
        l = res//2
        r = res - l

    newdata = data[u:-b, l:-r]

    # Reference: http://stackoverflow.com/a/8090605
    newdata = newdata.reshape(newshape[0], binning, newshape[1], -1).sum(3).sum(1).astype(numpy.uint32)

    # xy transformation
    f_xyconv = lambda x,y: ((x-u)/binning, (y-l)/binning)
    f_xyconvinv = lambda x,y: (x*binning+u, y*binning+l)

    return newdata, f_xyconv, f_xyconvinv
# software_binning()

def make_result(retin, f_trans=None):
    ret = dict(spots=[], lookup={})

    for r in retin:
        x, y = r.peak_com_x, r.peak_com_y
        if f_trans: x, y = f_trans(x,y)
        ret["spots"].append((y, x, r.peak_snr, r.peak_com_res))

    ret["lookup"]["cheetah"] = tuple(range(len(retin)))

    return ret
# make_distl_stats_dict()

def write_zoomed_image(data, vendortype, width, height, beamx, beamy,
                       distance, wavelength, pixel_size, 
                       thumb_width, d_min,
                       jpgout, brightness=100, color_scheme=0):
    """
    The code taken from yamtbx/dataproc/myspotfinder/spot_finder_for_grid_scan.py (New BSD License)
    """

    def get_zoom_box(img_w, img_h, x, y, boxsize=400, mag=16): # from viewer.image.get_zoom_box()
        n_pixels = int(math.ceil(boxsize / mag)+.5)
        x0 = min(img_w - n_pixels, int(math.floor(x - (n_pixels / 2))+.5))
        y0 = min(img_h - n_pixels, int(math.floor(y - (n_pixels / 2))+.5))
        return (x0, y0, n_pixels, n_pixels)
    # get_zoom_box()

    tw = thumb_width
    mag = tw / width # Default
    if d_min is not None:
        try:
            clipsize = distance * math.tan(2.*math.asin(wavelength/2./d_min)) / pixel_size*2
            if clipsize > 0:
                mag = tw / min(clipsize, width)
        except ValueError:
            pass # Use default

    x0, y0, w, h = get_zoom_box(width, height, beamx, beamy, boxsize=tw, mag=mag)
    assert (w == h)
    x0, y0 = max(x0, 0), max(y0, 0)

    if data.dtype == numpy.uint32: MyImage = imageconv_ext.MyImage_uint32
    elif data.dtype == numpy.uint16: MyImage = imageconv_ext.MyImage_uint16

    try:
        img = MyImage(rawdata=data, width=width, height=height,
                      vendortype=vendortype,
                      brightness=brightness/100.,
                      saturation=65535)

        img.prep_string_cropped_and_scaled(x0,y0,w,h,tw,tw)
        img.save_as_jpeg(str(jpgout), tw, tw, quality=90)

    except ImportError:
        pass

    return x0, y0, mag
# write_zoomed_image()

def cheetah_worker(header, data, work_dir, imgfile, algorithm, cut_roi, cheetah, binning=1, jpgdir=None):
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
    else: shikalog.error("Unsupported dtype: %s"%roidata.dtype)
        
    ret = cheetah_run(roidata.shape[1], roidata.shape[0], roidata,
                      beamx, beamy,
                      header["wavelength"], header["distance"], pixel_size,
                      algorithm=algorithm)

    if jpgdir is None: jpgdir = work_dir
    x0, y0, mag = write_zoomed_image(data, vendortype="EIGER", width=data.shape[1], height=data.shape[0],
                                     beamx=header["beam_center_x"], beamy=header["beam_center_y"],
                                     distance=header["distance"], wavelength=header["wavelength"],
                                     pixel_size=header["pixel_size_x"],
                                     thumb_width=600, d_min=5,
                                     jpgout=os.path.join(jpgdir, os.path.basename(imgfile)+".jpg"),
                                     brightness=150, color_scheme=0)

    result = make_result(ret, f_xyconvinv)
    result["thumb_posmag"] = (x0, y0, mag)
    result["work_dir"] = work_dir
    result["imgfile"] = imgfile
    result["template"] = "%s_%s.img"%(str(header["file_prefix"]), "?"*6)
    result["idx"] = header["frame"]+1
    return result
# cheetah_worker()

def worker(wrk_num, ventilator_hosts, eiger_host, result_host, pub_host, mode, cut_roi, algorithm, bl):
    """
    The code taken from yamtbx/dataproc/myspotfinder/command_line/spot_finder_backend.py (New BSD License)
    """
    context = zmq.Context()
    shikalog.config(bl, "cheetah")

    # Set up a channel to receive work from the ventilator
    work_receivers = []
    for vh in ventilator_hosts.split(","):
        work_receivers.append(context.socket(zmq.PULL))
        work_receivers[-1].connect("tcp://%s"%vh)

    eiger_receiver = context.socket(zmq.PULL)
    eiger_receiver.connect("tcp://%s"%eiger_host)
 
    # Set up a channel to send result of work to the results reporter
    results_sender = context.socket(zmq.PUSH)
    results_sender.connect("tcp://%s"%result_host)
 
    # Set up a channel to receive control messages over
    control_receiver = context.socket(zmq.SUB)
    control_receiver.connect("tcp://%s"%pub_host)
    control_receiver.setsockopt(zmq.SUBSCRIBE, "")
 
    # Set up a poller to multiplex the work receiver and control receiver channels
    poller = zmq.Poller()
    for wr in work_receivers: poller.register(wr, zmq.POLLIN)
    poller.register(control_receiver, zmq.POLLIN)
    if mode == "eiger_streaming":
        poller.register(eiger_receiver, zmq.POLLIN)

    master_params = libtbx.phil.parse(config_manager.master_params_str)
    working_params = master_params.fetch(sources=[libtbx.phil.parse("")])
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

    shikalog.info("worker %d ready" % wrk_num)

    while True:
        socks = dict(poller.poll())

        # the message from ventilator
        if any(map(lambda x: socks.get(x)==zmq.POLLIN, work_receivers)):
            work_receiver = filter(lambda x: socks.get(x)==zmq.POLLIN, work_receivers)[0]
            msg = work_receiver.recv_json()
            imgfile = str(msg["imgfile"])
            header = msg["header"]
            startt = time.time()

            if "h5master" in msg: #imgfile.endswith(".h5"):
                from yamtbx.dataproc import eiger
                data = eiger.extract_data(str(msg["h5master"]), msg["idx"])
                data[data<0] = 0
                data = data.astype(numpy.uint32)
            else:
                raise "Unsupported type"

            work_dir = os.path.join(os.path.dirname(imgfile), "_spotfinder")
            params.work_dir = work_dir
            
            if os.path.exists(work_dir): assert os.path.isdir(work_dir)
            else:
                try: os.mkdir(work_dir)
                except: pass

            try:
                result = cheetah_worker(header, data, work_dir, imgfile, algorithm, cut_roi, cheetah,
                                        params.cheetah.binning)
                result["starttime"] = startt
                result["endtime"] = time.time()
                result["params"] = params
                eltime = time.time()-startt
                shikalog.info("Wrkr%3d %s done in %.2f msec " % (wrk_num, imgfile, eltime*1.e3))
                results_sender.send_pyobj(result)
            except:
                shikalog.error(traceback.format_exc())


        # the message from EIGER
        if socks.get(eiger_receiver) == zmq.POLLIN:
            frames = eiger_receiver.recv_multipart(copy = False)
            startt = time.time()
            header, data = read_eiger_stream_data(frames)
            if None in (header, data): continue

            work_dir = os.path.join(str(header["data_directory"]), "_spotfinder")
            params.work_dir = work_dir
            if os.path.exists(work_dir): assert os.path.isdir(work_dir)
            else:
                try: os.mkdir(work_dir)
                except: pass

            #print "Got data:", header

            imgfile = os.path.join(header["data_directory"],
                                   "%s_%.6d.img"%(str(header["file_prefix"]), header["frame"]+1))

            jpgdir = os.path.join(work_dir, "thumb_%s_%.3d" % (str(header["file_prefix"]), (header["frame"]+1)//1000))
            try: os.mkdir(jpgdir)
            except: pass
            result = cheetah_worker(header, data, work_dir, imgfile, algorithm, cut_roi, cheetah,
                                    params.cheetah.binning, jpgdir=jpgdir)
            result["starttime"] = startt
            result["endtime"] = time.time()
            result["params"] = params
            eltime = time.time()-startt
            shikalog.info("Wrkr%3d %s Done in %.2f msec " % (wrk_num, imgfile, eltime*1.e3))
            results_sender.send_pyobj(result)

        if socks.get(control_receiver) == zmq.POLLIN:
            msg = control_receiver.recv_pyobj()
            if "params" in msg:
                params = msg["params"][("BL32XU", "EIGER9M", None, None)]
                cheetah.set_params(ADCthresh=params.cheetah.ADCthresh,
                                   MinSNR=params.cheetah.MinSNR,
                                   MinPixCount=params.cheetah.MinPixCount,
                                   MaxPixCount=params.cheetah.MaxPixCount,
                                   LocalBGRadius=params.cheetah.LocalBGRadius,
                                   MinPeakSeparation=params.cheetah.MinPeakSeparation,
                                   Dmin=params.distl.res.outer,
                                   Dmax=params.distl.res.inner)
                algorithm = params.cheetah.algorithm
                shikalog.info("worker %d: Parameters updated" % wrk_num)


    shikalog.info("Worker %d finished." % wrk_num)
# worker()

def run(opts):
    assert opts.algorithm in (6, 8)

    pp = []
    for i in xrange(opts.nproc):
        p = subprocess.Popen(["%s"%sys.executable, "-"], shell=True, stdin=subprocess.PIPE)
        p.stdin.write("from cheetah_client import worker\nworker(*%s)\n"%((i, opts.ventilator_hosts, opts.eiger_host, opts.result_host, opts.pub_host, opts.mode, opts.cut_roi, opts.algorithm, opts.bl),))
        print "from cheetah_client import worker\nworker(*%s)\n"%((i, opts.ventilator_hosts, opts.eiger_host, opts.result_host, opts.pub_host, opts.mode, opts.cut_roi, opts.algorithm, opts.bl),)
        p.stdin.close()
        pp.append(p)
    
    for p in pp: p.wait()

if __name__ == "__main__":
    import sys
    import optparse

    parser = optparse.OptionParser(usage="usage: %prog [options] coordinates...")

    parser.add_option("--eiger-host", action="store", dest="eiger_host", default="192.168.163.204:9999")
    parser.add_option("--vent-hosts", action="store", dest="ventilator_hosts", default="192.168.163.10:5556,127.0.0.1:5557")
    parser.add_option("--result-host", action="store", dest="result_host", default="127.0.0.1:5558")
    parser.add_option("--pub-host", action="store", dest="pub_host", default="192.168.163.10:5559")
    parser.add_option("--mode", action="store", dest="mode", default="eiger_streaming")
    parser.add_option("--nproc", action="store", dest="nproc", type=int, default=16)
    parser.add_option("--cut-roi", action="store_true", dest="cut_roi")
    parser.add_option("--algorithm", action="store", dest="algorithm", type=int, default=8)
    parser.add_option("--bl", action="store", dest="bl", default="32xu")


    opts, args = parser.parse_args(sys.argv)

    run(opts)
