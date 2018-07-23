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

import libtbx.phil
from libtbx import easy_mp

def str2float(s):
    m = re.match("-?\d+(.\d+)?(e[+-]?\d+)?", s)
    if m is not None:
        return float(m.group(0))
    else:
        return float("nan")

def run(opts):
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

    # Get shutter status and find images
    try:
        shutter = numpy.array(map(str2float, dbpy.read_syncdatalist(sensor_shutter, high_tag, tuple(tag_list))))
    except:
        print traceback.format_exc()
        error_status("NoShutterStatus")
        return -1

    # XXX
    valid_tags = tag_list[shutter==1] # [tag for tag, is_open in zip(tag_list, shutter) if is_open == 1]
    if 0:
        tag_offset = 3
        tag_list = tag_list[tag_offset:]
        valid_tags = tag_list[numpy.arange(1, len(tag_list)+1)%6==0]
        
    if valid_tags.size == 0:
        error_status("NoValidTags")
        return -1

    # Get PD values
    pd1_values, pd2_values, pd3_values = map(lambda y: map(lambda x: float("nan"), xrange(len(valid_tags))), xrange(3))
    if opts.pd1_sensor_name:
        pd1_values = numpy.array(map(str2float, dbpy.read_syncdatalist(opts.pd1_sensor_name, high_tag, tuple(valid_tags))))
    if opts.pd2_sensor_name:
        pd2_values = numpy.array(map(str2float, dbpy.read_syncdatalist(opts.pd2_sensor_name, high_tag, tuple(valid_tags))))
    if opts.pd3_sensor_name:
        pd3_values = numpy.array(map(str2float, dbpy.read_syncdatalist(opts.pd3_sensor_name, high_tag, tuple(valid_tags))))

    print "tag pd1 pd2 pd3"
    for tag, pd1, pd2, pd3 in zip(valid_tags, pd1_values, pd2_values, pd3_values):
        print tag, pd1, pd2, pd3

    for i in xrange(len(valid_tags)):
        bad = []
        if (opts.pd1_threshold != 0 and
            not (opts.pd1_threshold > 0 and opts.pd1_threshold <= pd1_values[i]) and
            not (opts.pd1_threshold < 0 and -opts.pd1_threshold > pd1_values[i])): bad.append(1)
        if (opts.pd2_threshold != 0 and
            not (opts.pd2_threshold > 0 and opts.pd2_threshold <= pd2_values[i]) and
            not (opts.pd2_threshold < 0 and -opts.pd2_threshold > pd2_values[i])): bad.append(2)
        if (opts.pd3_threshold != 0 and
            not (opts.pd3_threshold > 0 and opts.pd3_threshold <= pd3_values[i]) and
            not (opts.pd3_threshold < 0 and -opts.pd3_threshold > pd3_values[i])): bad.append(3)

        if bad:
            print "# Bad tag=%d BadPD=%s" %(valid_tags[i], bad)

# run()
    
if __name__ == "__main__":
    import sys
    import argparse
    
    parser = argparse.ArgumentParser(usage="usage: %s [options] @config.txt" % __file__, fromfile_prefix_chars="@")

    parser.add_argument("-r","--run", action="store", dest="runid", type=int)
    parser.add_argument("--pd1-name", action="store", dest="pd1_sensor_name", type=str)
    parser.add_argument("--pd2-name", action="store", dest="pd2_sensor_name", type=str)
    parser.add_argument("--pd3-name", action="store", dest="pd3_sensor_name", type=str)
    parser.add_argument("--pd1-thresh", action="store", dest="pd1_threshold", type=float, default=0)
    parser.add_argument("--pd2-thresh", action="store", dest="pd2_threshold", type=float, default=0)
    parser.add_argument("--pd3-thresh", action="store", dest="pd3_threshold", type=float, default=0)
    parser.add_argument("--bl", action="store", dest="bl", type=int)

    opts = parser.parse_args()
    run(opts)
