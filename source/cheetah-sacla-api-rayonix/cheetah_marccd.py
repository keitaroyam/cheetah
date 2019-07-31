"""
Reference: biochem-fan's code from https://github.com/biochem-fan/cheetah/blob/master/source/cheetah-sacla-api2/main-sacla-api2.cpp
"""
import os
import re
import h5py
import json
import struct
import numpy
import math
import time
import datetime
import glob
import traceback
from multiprocessing import Process, Queue

import dbpy
import marccd
from cheetah_ext import CheetahSinglePanel, CheetahSpot

import libtbx.phil
from libtbx import easy_mp

import bitshuffle.h5
import dxtbx.format

# Constants
PD_DARK_ANY = -2
PD_ANY = -1
PD_LIGHT = 0
parallel_size = 3

def make_geom(f, geom_out, beam_x=None, beam_y=None, clen=None):
    im = marccd.MarCCD(f)

    if beam_x and beam_x==beam_x: im.beam_x = beam_x # in px
    if beam_y and beam_y==beam_y: im.beam_y = beam_y # in px
    if clen and clen==clen: im.distance = clen # in mm

    # according to xds, gain of /xustrg0/rayonix/2018A/LRE1/lre-1343-1/001/sample_000????.img is ~0.2.
    # usually ccd's gain is underetimated (that is, photon is actually less) by a factor of 4 (a=4 in xds)
    # so 0.2*4=0.8 is appropriate value for adu_per_photon?
    
    s = """\
clen = %(clen)f
res = %(res).4f
data = /%%/data
photon_energy = /%%/photon_energy_ev

rigid_group_q0 = q0
rigid_group_collection_connected = q0
rigid_group_collection_independent = q0

q0/adu_per_eV = 0.8e-4 ; dummy for old hdfsee
q0/adu_per_photon = 0.8
q0/max_adu = %(max_adu)d
q0/min_fs = 0
q0/max_fs = %(fsmax)d
q0/min_ss = 0
q0/max_ss = %(ssmax)d
q0/corner_x = %(cornerx).2f
q0/corner_y = %(cornery).2f
q0/fs = -x
q0/ss = -y
""" % dict(fsmax=im.nfast-1, ssmax=im.nslow-1,
           cornerx=im.beam_x, cornery=im.beam_y,
           clen=im.distance/1000.,
           res=1./(im.pixel_x*1.e-3),
           max_adu=im.saturated_value)

    open(geom_out, "w").write(s)
# make_geom()

def make_result(retin, f_trans=None):
    ret = dict(spots=[], lookup={})

    for r in retin:
        x, y = r.peak_com_x, r.peak_com_y
        if f_trans: x, y = f_trans(x,y)
        ret["spots"].append((x, y, r.peak_snr, r.peak_com_res))

    return ret
# make_distl_stats_dict()

def cheetah_worker(header, roidata, algorithm, cheetah):
    beamx, beamy = header["beam_center_x"], header["beam_center_y"]
    pixel_size = header["pixel_size_x"]

    if roidata.dtype == numpy.uint32: cheetah_run = cheetah.run_uint32 
    elif roidata.dtype == numpy.uint16: cheetah_run = cheetah.run_uint16

    ret = cheetah_run(roidata.shape[1], roidata.shape[0], roidata,
                      float(beamx), float(beamy),
                      float(header["wavelength"]), float(header["distance"]), float(pixel_size),
                      algorithm=algorithm)

    result = make_result(ret, None)
    return result
# cheetah_worker()

def worker(wrk_num, header, imgfiles, queue, algorithm, 
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
        data = marccd.MarCCD(f).read_data()
        result = cheetah_worker(header, data, algorithm, cheetah)
        result["frame"] = f
        eltime = time.time()-startt
        #print "Wrkr%3d %6d done in %.2f msec " % (wrk_num, frame, eltime*1.e3)
        print "%s %6d %.2f" %(f, len(result["spots"]), eltime*1.e3)
        queue.put(result)
    #print "Worker %d finished." % wrk_num
# worker()

def make_h5(out, file_tag_ene, comment, default_energy=None, compression="shuf+gz"):
    startt = time.time()

    from multiprocessing.dummy import Pool as ThreadPool
    import threading
    import zlib

    lock = threading.Lock()
    of = h5py.File(out, "w")

    of["/metadata/detector"] = "Rayonix MX300HS"
    #of["/metadata/distance_in_mm"] = opts.clen_mm
    if comment: of["/metadata/run_comment"] = comment
    if file_tag_ene:
        tmp = marccd.MarCCD(file_tag_ene[0][0])
        of["/metadata/pixelsize_in_um"] = tmp.pixel_x*1000.
        if not default_energy: default_energy = 12.3984 / tmp.wavelength

    def write_par(i):
        f, tag, ene = file_tag_ene[i]
        tmp = time.time()
        data = marccd.MarCCD(f).read_data()
        grp = of.create_group("tag-%d"%tag)
        if ene!=ene: ene = default_energy

        if compression=="shuf+gz":
            as_uint8 = data.view(dtype=numpy.uint8)
            shuffled = as_uint8.reshape((-1, data.dtype.itemsize)).transpose().reshape(-1)
            shuffled_compressed = zlib.compress(shuffled.tobytes(), 4)

        lock.acquire()
        grp["photon_energy_ev"] = ene*1000.
        grp["photon_wavelength_A"] = 12.3984/ene
        grp["original_file"] = f

        if not compression:
            grp.create_dataset("data", data.shape,
                               dtype=data.dtype, data=data)        
        elif compression=="bslz4":
            grp.create_dataset("data", data.shape,
                               compression=bitshuffle.h5.H5FILTER,
                               compression_opts=(0, bitshuffle.h5.H5_COMPRESS_LZ4),
                               dtype=data.dtype, data=data)
        elif compression=="shuf+gz":
            dataset = grp.create_dataset("data", data.shape,
                                         compression="gzip", shuffle=True,
                                         dtype=data.dtype, chunks=data.shape)
            dataset.id.write_direct_chunk(offsets=(0, 0), data=shuffled_compressed, filter_mask=0)
        else:
            raise "Unknwon compression name (%s)" % compression

        print "# converted: %s %d %.4f %.2f" %(f, tag, ene, (time.time()-tmp)*1.e3)
        lock.release()

    pool = ThreadPool(opts.nproc)
    pool.map(write_par, xrange(len(file_tag_ene)))

    of.close()
    
    eltime = time.time()-startt
    print "# Processed: %s (%.2f sec for %d images)" % (out, eltime, len(file_tag_ene))
# make_h5()

def str2float(s):
    m = re.match("-?\d+(.\d+)?(e[+-]?\d+)?", s)
    if m is not None:
        return float(m.group(0))
    else:
        return float("nan")

def process_images(img_files, mean_photon_energy, opts):
    startt = time.time()

    nframes = len(img_files)
    first_img = img_files[0] # Assumes all images are with the same conditions
    im = marccd.MarCCD(first_img)

    header = dict(beam_center_x=im.beam_x,
                  beam_center_y=im.beam_y,
                  pixel_size_x=im.pixel_x,
                  distance=im.distance,
                  wavelength=im.wavelength)

    if opts.clen:
        print "# overriding camera distance = %f (header value: %f)" % (opts.clen, header["distance"])
        header["distance"] = opts.clen

    if mean_photon_energy and mean_photon_energy==mean_photon_energy:
        print "# overriding wavelength = %f (header value: %f)" % (12.3984/mean_photon_energy, header["wavelength"])
        header["wavelength"] = 12.3984/mean_photon_energy

    if opts.beam_x and opts.beam_x==opts.beam_x:
        print "# overriding beam_x = %f (header value: %f)" % (opts.beam_x, header["beam_center_x"])
        header["beam_center_x"] = opts.beam_x

    if opts.beam_y and opts.beam_y==opts.beam_y:
        print "# overriding beam_y = %f (header value: %f)" % (opts.beam_y, header["beam_center_y"])
        header["beam_center_y"] = opts.beam_y
        
    
    print "#%s"%header
    print "#frames= %d" % nframes

    ranges = map(tuple, numpy.array_split(numpy.arange(nframes), opts.nproc))
    queue = Queue()
    pp = []
    print "frame nspots eltime"
    for i, rr in enumerate(ranges):
        if not rr: continue
        p = Process(target=worker, args=(i, header, map(lambda x: img_files[x], rr), queue, opts.algorithm, 
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

    print "# hit finding for %d images finished in %.2f sec." % (len(img_files), time.time()-startt)

    return results
# process_images()

def error_status(err_str):
    open("status.txt", "w").write("""\
# Cheetah status
Update time: %(ctime)s
Status: Status=Error-%(err_str)s
""" % dict(ctime=time.ctime(), err_str=err_str))
# error_status()


def run(opts):
    eltime_from = time.time()
    print("#\n#Configurations:")
    print("# runNumber (-r/--run):         %d" % opts.runid)
    print("# output H5 file (-o/--output): %s (default = run######.h5)" % opts.outputH5)
    print("# beamline (--bl):              %d (default = 3)" % opts.bl)
    print("# img root (--rayonix-root):    %s" % opts.rayonix_root)
    print("# distance (--clen):            %s" % opts.clen)
    print("# beam center (--beam-x/y):     %s,%s" % (opts.beam_x,opts.beam_y))
    print("# Cheetah settings")
    print("#  --dmin, --dmax:              %s,%s" % (opts.d_min, opts.d_max))
    print("#  --adc-threshold:             %s" % opts.ADCthresh)
    print("#  --min-snr:                   %s" % opts.MinSNR)
    print("#  --min/max-pixcount:          %s,%s" % (opts.MinPixCount, opts.MaxPixCount))
    print("#  --local-bgradius:            %s" % opts.LocalBGRadius)
    print("#  --min-peaksep:               %s" % opts.MinPeakSeparation)
    print("#  --min-spots:                 %s" % opts.min_spots)
    print("#  --algorithm:                 %s" % opts.algorithm)
    print("# PD1 threshold (--pd1_thresh): %.3f (default = 0; ignore.)" % opts.pd1_threshold)
    print("# PD2 threshold (--pd2_thresh): %.3f (default = 0; ignore.)" % opts.pd2_threshold)
    print("# PD3 threshold (--pd3_thresh): %.3f (default = 0; ignore.)" % opts.pd3_threshold)
    print("# PD1 sensor name (--pd1_name): %s)" % opts.pd1_sensor_name)
    print("# PD2 sensor name (--pd2_name): %s)" % opts.pd2_sensor_name)
    print("# PD3 sensor name (--pd3_name): %s)" % opts.pd3_sensor_name)
    print("# nFrame after light:           %d (default = -1; accept all image. -2; accept all dark images)" % opts.light_dark)
    print("# parallel_block:               %d (default = -1; no parallelization)" % opts.parallel_block)
    print("# nproc:                        %d (default = 1)" % opts.nproc)
    print("")
    
    assert opts.algorithm in (6, 8)
    assert opts.runid is not None
    assert opts.bl is not None

    # Beamline specific constants
    if opts.bl == 2:
        sensor_spec = "xfel_bl_2_tc_spec_1/energy"
        sensor_shutter = "xfel_bl_2_shutter_1_open_valid/status"
    elif opts.bl == 3:
        sensor_spec = "xfel_bl_3_tc_spec_1/energy"
        sensor_shutter = "xfel_bl_3_shutter_1_open_valid/status"
    else:
        error_status("BadBeamline")
        return -1

    # Get run info
    try:
        run_info = dbpy.read_runinfo(opts.bl, opts.runid)
    except:
        error_status("BadRunID")
        return -1

    high_tag = dbpy.read_hightagnumber(opts.bl, opts.runid)
    start_tag = run_info['start_tagnumber']
    end_tag = run_info['end_tagnumber']

    tag_list = numpy.array(dbpy.read_taglist_byrun(opts.bl, opts.runid))
    print "# Run %d: HighTag %d, Tags %d (inclusive) to %d (exclusive), thus %d tags" % (opts.runid, high_tag, start_tag, end_tag, len(tag_list))
    comment = dbpy.read_comment(opts.bl, opts.runid)
    print "# Comment: %s" % comment
    print

    # Get shutter status and find images
    try:
        shutter = numpy.array(map(str2float, dbpy.read_syncdatalist(sensor_shutter, high_tag, tuple(tag_list))))
    except:
        print traceback.format_exc()
        error_status("NoShutterStatus")
        return -1

    # XXX
    valid_tags = tag_list[shutter==1] # [tag for tag, is_open in zip(tag_list, shutter) if is_open == 1]
    print "# DEBUG:: shutter=", shutter
    print "# DEBUG:: valid_tags=", valid_tags
    if 0:
        tag_offset = 3
        tag_list = tag_list[tag_offset:]
        valid_tags = tag_list[numpy.arange(1, len(tag_list)+1)%6==0]
        
    if valid_tags.size == 0:
        error_status("NoValidTags")
        return -1

    # Find images
    img_files = sorted(glob.glob(os.path.join(opts.rayonix_root, str(opts.runid), "data_*.img")))
    print "# DEBUG:: img_files=%d valid_tags=%d" % (len(img_files), len(valid_tags))
    if len(img_files)+1 != len(valid_tags): # last valid tag is not saved.
        print "# WARNING!! img_files and valid_tag number mismatch"

        img_numbers = map(lambda x: int(x[x.rindex("_")+1:-4]), img_files)
        dropped_frames = sorted(set(range(1, len(valid_tags))).difference(img_numbers))
        print "# Unsaved frame numbers =", tuple(dropped_frames)
        print "# DEBUG::", len(img_files)-len(dropped_frames)+1, len(valid_tags)
        if len(img_files)+len(dropped_frames)+1 == len(valid_tags):
            print "#  %d unsaved img files found, which explains number mismatch" % len(dropped_frames)
            valid_tags = numpy.delete(valid_tags, numpy.array(dropped_frames)-1)
            assert len(img_files)+1 == len(valid_tags)
        else:
            print "# Assuming last %d img files are generated after stopping run.." % (len(img_files)-len(valid_tags)+1)
            img_files = img_files[:len(valid_tags)-1]
            assert len(img_files)+1 == len(valid_tags)
    
    # Get photon energies
    photon_energies_in_keV  = numpy.array([str2float(s) for s in dbpy.read_syncdatalist(sensor_spec, high_tag, tuple(valid_tags))])
    mean_photon_energy = numpy.mean(photon_energies_in_keV[photon_energies_in_keV==photon_energies_in_keV]) # XXX if no valid data?
    print "# Photon energies obtained: %d valid numbers, %d invalid, average=%f sd=%f" % (len(photon_energies_in_keV),
                                                                                          sum(photon_energies_in_keV!=photon_energies_in_keV),
                                                                                          mean_photon_energy, numpy.std(photon_energies_in_keV[photon_energies_in_keV==photon_energies_in_keV]))
    photon_energies_in_keV[photon_energies_in_keV!=photon_energies_in_keV] = mean_photon_energy
    
    # Get PD values
    pd1_values, pd2_values, pd3_values = [], [], []
    if opts.pd1_threshold != 0:
        pd1_values = numpy.array(map(str2float, dbpy.read_syncdatalist(opts.pd1_sensor_name, high_tag, tuple(valid_tags))))
    if opts.pd2_threshold != 0:
        pd2_values = numpy.array(map(str2float, dbpy.read_syncdatalist(opts.pd2_sensor_name, high_tag, tuple(valid_tags))))
    if opts.pd3_threshold != 0:
        pd3_values = numpy.array(map(str2float, dbpy.read_syncdatalist(opts.pd3_sensor_name, high_tag, tuple(valid_tags))))

    # Identify bad tags
    # XXX not tested!! this feature must not be used. tags with bad PDs must be detected after experiment.
    frame_after_light = 9999
    bad_tag_idxes = []
    for i in xrange(len(valid_tags)):
        light = True
        if (opts.pd1_threshold != 0 and
            not (opts.pd1_threshold > 0 and opts.pd1_threshold <= pd1_values[i]) and
            not (opts.pd1_threshold < 0 and -opts.pd1_threshold > pd1_values[i])): light = False
        if (opts.pd2_threshold != 0 and
            not (opts.pd2_threshold > 0 and opts.pd2_threshold <= pd2_values[i]) and
            not (opts.pd2_threshold < 0 and -opts.pd2_threshold > pd2_values[i])): light = False
        if (opts.pd3_threshold != 0 and
            not (opts.pd3_threshold > 0 and opts.pd3_threshold <= pd3_values[i]) and
            not (opts.pd3_threshold < 0 and -opts.pd3_threshold > pd3_values[i])): light = False

        if light:
            frame_after_light = 0
        else:
            frame_after_light += 1

        if ((opts.light_dark >= 0 and frame_after_light != opts.light_dark) or
            (opts.light_dark == PD_DARK_ANY and frame_after_light == 0)):
            print "# PD check: %d is bad tag!" % valid_tags[i]
            bad_tag_idxes.append(i)

    if bad_tag_idxes:
        valid_tags = numpy.delete(valid_tags, numpy.array(bad_tag_idxes))
        for i in reversed(bad_tag_idxes): del img_files[i]


    # Debug code; this takes too much time!
    try:
        if 0 and opts.parallel_block == 0:
            tag_timestamp = map(lambda x: datetime.datetime.fromtimestamp(dbpy.read_timestamp_fromtag(high_tag, x, sensor_shutter)).strftime('%Y-%m-%d %H:%M:%S.%f'), valid_tags)
            img_timestamp = map(lambda x: marccd.MarCCD(x).acquire_time.strftime('%Y-%m-%d %H:%M:%S.%f'), img_files)
            ofs = open("tag_file_time.dat", "w")
            ofs.write("run tag file tag.time file.time\n")
            for i in xrange(len(img_files)):
                ofs.write('%d %d %s "%s" "%s"\n'%(opts.runid, valid_tags[i], img_files[i], tag_timestamp[i], img_timestamp[i]))
            ofs.close()
    except:
        pass

    # block spliting
    # TODO db query may be slow, which may need to be done only once?
    if opts.parallel_block >= 0:
        width = len(valid_tags)//parallel_size
        i_start = opts.parallel_block*width
        i_end = (opts.parallel_block+1)*width if opts.parallel_block < parallel_size-1 else None
        valid_tags = valid_tags[i_start:i_end]
        photon_energies_in_keV = photon_energies_in_keV[i_start:i_end]
        img_files = img_files[i_start:i_end]
        print "# parallel_block=%d: %d tags will be processed (%d..%d)" % (opts.parallel_block, len(valid_tags), valid_tags[0], valid_tags[-1])


    make_geom(img_files[0], opts.output_geom, beam_x=opts.beam_x, beam_y=opts.beam_y, clen=opts.clen)

    # Hit-finding
    results = process_images(img_files, mean_photon_energy, opts)
    file_tag_ene = []
    for frame, tag, ene in zip(sorted(results), valid_tags, photon_energies_in_keV):
        if len(results[frame]["spots"]) < opts.min_spots:
            continue
        file_tag_ene.append((frame, tag, ene))

    # TODO on-the-fly status updating
    open("status.txt", "w").write("""\
# Cheetah status
Update time: %(ctime)s
Elapsed time: %(eltime)f sec
Status: Total=%(ntotal)d,Processed=%(ntotal)d,LLFpassed=%(ntotal)d,Hits=%(nhits)d,Status=WritingH5
Frames processed: %(ntotal)d
Number of hits: %(nhits)d
""" % dict(ctime=time.ctime(), eltime=time.time()-eltime_from, ntotal=len(img_files), nhits=len(file_tag_ene)))

    # Save h5
    # TODO implement on-the-fly h5 file writing in hit-finding to avoid reading img file twice.
    make_h5(out=opts.outputH5,
            file_tag_ene=file_tag_ene,
            comment=comment)
    
    open("status.txt", "w").write("""\
# Cheetah status
Update time: %(ctime)s
Elapsed time: %(eltime)f sec
Status: Total=%(ntotal)d,Processed=%(ntotal)d,LLFpassed=%(ntotal)d,Hits=%(nhits)d,Status=Finished
Frames processed: %(ntotal)d
Number of hits: %(nhits)d
""" % dict(ctime=time.ctime(), eltime=time.time()-eltime_from, ntotal=len(img_files), nhits=len(file_tag_ene)))


    ofs = open("cheetah.dat", "w")
    ofs.write("file tag nspots total_snr\n")
    for frame, tag in zip(sorted(results), valid_tags):
        ret = results[frame]
        n_spots = len(ret["spots"])
        total_snr = sum(map(lambda x: x[2], ret["spots"]))
        ofs.write("%s %d %6d %.3e\n"%(frame, tag, n_spots, total_snr))
    ofs.close()

    if opts.gen_adx:
        for frame in sorted(results):
            ret = results[frame]
            adx_out = open(os.path.basename(frame)+".adx", "w")
            for x,y,snr,d in ret["spots"]: adx_out.write("%6d %6d %.2e\n" % (x,y,snr))
            adx_out.close()
# run()
    
if __name__ == "__main__":
    import sys
    import argparse
    
    parser = argparse.ArgumentParser(usage="usage: %s [options] @config.txt" % __file__, fromfile_prefix_chars="@")

    parser.add_argument("--nproc", action="store", dest="nproc", type=int, default=1)
    parser.add_argument("--gen-adx", action="store_true", dest="gen_adx")
    parser.add_argument("--clen", dest="clen", type=float, help="camera length in mm")
    parser.add_argument("--beam-x", dest="beam_x", type=float, help="beam center x in px")
    parser.add_argument("--beam-y", dest="beam_y", type=float, help="beam center y in px")

    parser.add_argument("--dmin", action="store", dest="d_min", type=float, default=0)
    parser.add_argument("--dmax", action="store", dest="d_max", type=float, default=30)
    parser.add_argument("--adc-threshold", action="store", dest="ADCthresh", type=float, default=5)
    parser.add_argument("--min-snr", action="store", dest="MinSNR", type=float, default=8)
    parser.add_argument("--min-pixcount", action="store", dest="MinPixCount", type=int, default=3)
    parser.add_argument("--max-pixcount", action="store", dest="MaxPixCount", type=int, default=40)
    parser.add_argument("--local-bgradius", action="store", dest="LocalBGRadius", type=float, default=2)
    parser.add_argument("--min-peaksep", action="store", dest="MinPeakSeparation", type=float, default=0)
    parser.add_argument("--algorithm", action="store", dest="algorithm", type=int, default=8)
    parser.add_argument("--min-spots", action="store", dest="min_spots", type=int, default=20)

    parser.add_argument("-r","--run", action="store", dest="runid", type=int)
    parser.add_argument("--pd1-thresh", action="store", dest="pd1_threshold", type=float, default=0)
    parser.add_argument("--pd2-thresh", action="store", dest="pd2_threshold", type=float, default=0)
    parser.add_argument("--pd3-thresh", action="store", dest="pd3_threshold", type=float, default=0)
    parser.add_argument("--pd1-name", action="store", dest="pd1_sensor_name", type=str)
    parser.add_argument("--pd2-name", action="store", dest="pd2_sensor_name", type=str)
    parser.add_argument("--pd3-name", action="store", dest="pd3_sensor_name", type=str)
    parser.add_argument("--type", action="store", dest="type", type=str)
    parser.add_argument("--bl", action="store", dest="bl", type=int)
    parser.add_argument("--rayonix-root", action="store", dest="rayonix_root", type=str, default="/xustrg0/SFX")
    parser.add_argument("-o","--output", action="store", dest="outputH5", type=str)
    parser.add_argument("--geom-out", action="store", dest="output_geom", type=str)

    #opts, args = parser.parse_args(sys.argv[1:])
    opts = parser.parse_args()

    opts.light_dark = PD_ANY
    opts.parallel_block = -1


    if opts.type:
        if opts.type == "light":
            opts.light_dark = PD_LIGHT
        elif opts.type == "dark":
            opts.light_dark = PD_DARK_ANY
        elif len(opts.type)==5 and opts.startswith("dark"):
            opts.light_dark = int(opts.type[4])
            if opts.light_dark < 1 or opts.light_dark > 9:
                print "ERROR: wrong type."
                sys.exit(1)
        else:
            opts.parallel_block = int(opts.type)
            if opts.parallel_block < -1 or opts.parallel_block >= parallel_size:
                print "ERROR: wrong type or parallel_block."
                sys.exit(1)

    if not opts.outputH5: opts.outputH5 = "run%d.h5" % opts.runid
    if not opts.output_geom: opts.output_geom = "%d.geom" % opts.runid

    
    run(opts)
    quit()

    """
    # Debug
    results = process_images(args, 10., opts)
    file_tag_ene = []
    for i, frame in enumerate(sorted(results)):
        if len(results[frame]["spots"]) < opts.min_spots:
            continue
        file_tag_ene.append((frame, i, float("nan")))
    
    # XXX
    make_h5(out=opts.outputH5,
            file_tag_ene=file_tag_ene,
            comment="hogehoge", compression="bslz4")

    if opts.gen_adx:
        for frame in sorted(results):
            ret = results[frame]
            adx_out = open(os.path.basename(frame)+".adx", "w")
            for x,y,snr,d in ret["spots"]: adx_out.write("%6d %6d %.2e\n" % (x,y,snr))
            adx_out.close()

    """
