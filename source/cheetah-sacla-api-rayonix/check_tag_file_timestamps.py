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
import dbpy
import marccd

def str2float(s):
    m = re.match("-?\d+(.\d+)?(e[+-]?\d+)?", s)
    if m is not None:
        return float(m.group(0))
    else:
        return float("nan")

def error_status(s):
    print s
    
def run(opts):
    eltime_from = time.time()

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
    print "# Run %d: HighTag %d, Tags %d (inclusive) to %d (exclusive), thus %d images" % (opts.runid, high_tag, start_tag, end_tag, len(tag_list))
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
    print "DEBUG:: shutter=", shutter
    print "DEBUG:: valid_tags=", valid_tags
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
    
    tag_timestamp = map(lambda x: datetime.datetime.fromtimestamp(dbpy.read_timestamp_fromtag(high_tag, x, sensor_shutter)).strftime('%Y-%m-%d %H:%M:%S.%f'), valid_tags)
    img_timestamp = map(lambda x: marccd.MarCCD(x).acquire_time.strftime('%Y-%m-%d %H:%M:%S.%f'), img_files)
    ofs = open("tag_file_time.dat", "w")
    ofs.write("run tag file tag.time file.time\n")
    for i in xrange(len(img_files)):
        ofs.write('%d %d %s "%s" "%s"\n'%(opts.runid, valid_tags[i], img_files[i], tag_timestamp[i], img_timestamp[i]))
    ofs.close()
    
# run()
    
if __name__ == "__main__":
    import sys
    import optparse

    parser = optparse.OptionParser(usage="usage: %prog [options]")

    parser.add_option("-r","--run", action="store", dest="runid", type=int)
    parser.add_option("--bl", action="store", dest="bl", type=int)
    parser.add_option("--rayonix-root", action="store", dest="rayonix_root", type=str, default="/xustrg0/rayonix/2018A/SFX")

    opts, args = parser.parse_args(sys.argv[1:])
    run(opts)
