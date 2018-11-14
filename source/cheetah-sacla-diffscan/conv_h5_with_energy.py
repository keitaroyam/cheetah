import h5py
import numpy
import math
import os
import time
import dbpy
from yamtbx.dataproc.XIO import XIO
from yamtbx.dataproc import cbf
from yamtbx.dataproc import bl_logfiles
from libtbx import easy_mp

"""
GAIN (0.294*4?) * h["wavelength"]/12398.4
"""

def make_geom(f, geom_out, corner_x=None, corner_y=None, rotated=False):
    h, junk = read_image(f, read_data=False)

    if corner_x is None or corner_x!=corner_x: corner_x = h["orgx"]
    if corner_y is None or corner_y!=corner_y: corner_y = h["orgy"]

    if rotated:
        fs, ss = "y", "-x"
        corner_x, corner_y = corner_y, -corner_x
    else:
        fs, ss = "-x", "-y"

    print h
    s = """\
data = /%%/data
photon_energy = /%%/photon_energy_ev

rigid_group_q0 = q0
rigid_group_collection_connected = q0
rigid_group_collection_independent = q0

q0/res = %(res).4f
q0/min_fs = 0
q0/max_fs = %(fsmax)d
q0/min_ss = 0
q0/max_ss = %(ssmax)d
q0/corner_x = %(cornerx).2f
q0/corner_y = %(cornery).2f
q0/fs = %(fs)s
q0/ss = %(ss)s
q0/clen = %(clen)f
q0/adu_per_eV = /LCLS/adu_per_eV
q0/max_adu = %(max_adu)d
""" % dict(fsmax=h["size1"]-1, ssmax=h["size2"]-1,
           cornerx=corner_x, cornery=corner_y, fs=fs, ss=ss,
           clen=h["distance"]/1000., res=1./(h["pixel_size"]*1.e-3), max_adu=65535)

    open(geom_out, "w").write(s)
# make_geom()

def read_cbf(f, read_data=True):
    h = {}
    h["size1"] = 2880
    h["size2"] = 2880
    h["pixel_size"] = 0.0782
    h["distance"] = 94.
    h["wavelength"] = 1.23984
    h["beamx"] = 112.608
    h["beamy"] = 112.608
    h["orgx"], h["orgy"] = h["beamx"]/h["pixel_size"], h["beamy"]/h["pixel_size"]

    if read_data:
        data, nfast, nslow = cbf.load_cbf_as_numpy(f)
        data = numpy.array(data, dtype=numpy.uint16).reshape(h["size2"],h["size1"])

    return h, data

def read_image(f, read_data=True):
    if f.endswith(".cbf"): return read_cbf(f, read_data)

    h = {}
    data = None

    im = XIO.Image(f)
    h["size1"] = im.header["Width"]
    h["size2"] = im.header["Height"]
    h["pixel_size"] = im.header["PixelX"]
    h["distance"] = im.header["Distance"]
    h["wavelength"] = im.header["Wavelength"]
    h["beamx"] = im.header["BeamX"]
    h["beamy"] = im.header["BeamY"]
    h["orgx"], h["orgy"] = h["beamx"]/h["pixel_size"], h["beamy"]/h["pixel_size"]

    if read_data:
        data = numpy.array(im.getData(), dtype=numpy.uint16).reshape(h["size2"],h["size1"])

    return h, data
# read_cmos_image()

def get_energies(tags, taghi, bl):
    if not tags: return []

    valid_tags = filter(lambda x: type(x) is int, tags)

    if valid_tags:
        ene = dbpy.read_syncdatalist("xfel_bl_%d_tc_spec_1/energy"%bl, taghi, tuple(valid_tags))
        assert len(ene) == len(valid_tags)
    else:
        ene = ()

    tag_ene = dict(zip(valid_tags, ene))

    ret = []
    for tag in tags:
        e = tag_ene.get(tag)
        try: e = float(e)
        except: e = float("nan")
        ret.append(e)
    return ret            
# get_energies()

def read_cheetah_dat(datin, scanlog, datout):
    if not os.path.isfile(datin): return []

    slog = bl_logfiles.BssDiffscanLog(scanlog)
    slog.remove_overwritten_scans()
    filename_tag = {}
    for scan in slog.scans:
        for f, tag in scan.filename_tags: filename_tag[os.path.basename(f)] = tag

    out = open(datout, "w")
    out.write("file tag x y nspots\n")

    ret = []
    fin = open(datin)
    fin.readline()
    for l in fin:
        sp = l.split()
        if len(sp) == 3:
            f = sp[0]
            nspots = int(sp[1])

            xy = slog.calc_grid_coord(filename=f)
            if xy[0]!=xy[0] or xy[1]!=xy[1]: continue # nan if dummy

            tag = filename_tag.get(os.path.basename(f))
            ret.append((f, nspots, tag))
            out.write("%s %s %f %f %d\n" % (os.path.basename(f), tag, xy[0], xy[1], nspots))

    return ret
# read_cheetah_dat()

def make_h5(out, file_spots_tag, energies, comment, default_energy):
    of = h5py.File(out, "w")

    of["/metadata/detector"] = "Rayonix MX300HS"
    if comment: of["/metadata/run_comment"] = comment
    if file_spots_tag: of["/metadata/pixelsize_in_um"] = read_image(file_spots_tag[0][0], False)[0]["pixel_size"]*1000.
    #of["/metadata/distance_in_mm"] = opts.clen_mm

    for (f, s, tag), ene in zip(file_spots_tag, energies):
        print "converting", f, s, tag, ene
        h, data = read_image(f)
        group_name = "tag-%d"%tag if tag is not None else "file-%s" % os.path.splitext(os.path.basename(f))[0]
        grp = of.create_group(group_name)
        if ene!=ene: ene = default_energy
        dset = grp.create_dataset("photon_energy_ev", (1,), dtype=numpy.float)
        dset[...] = ene * 1000.
        dset = grp.create_dataset("photon_wavelength_A", (1,), dtype=numpy.float)
        dset[...] = 12.3984 / ene
        dset = grp.create_dataset("original_file", (1,), "S%d"%len(f))
        dset[...] = f

        dset = grp.create_dataset("data", data.shape, dtype=data.dtype, compression="gzip", shuffle=True)
        dset[...] = data

    of.close()
    print "Processed:", out
# make_h5()

def run(opts, cheetah_dat, scanlog):#, ene_dat):

    #energies = read_ene_dat(ene_dat, opts.default_energy)
    file_spots_tag = read_cheetah_dat(cheetah_dat, scanlog, os.path.splitext(opts.out)[0]+"_forplot.dat")
    n_proccessed_images = len(file_spots_tag)

    if opts.geom_out and file_spots_tag:
        make_geom(file_spots_tag[0][0], opts.geom_out, corner_x=opts.corner_x, corner_y=opts.corner_y, rotated=opts.rotated)

    file_spots_tag = filter(lambda x: x[1] >= opts.min_spots, file_spots_tag)
    energies = get_energies(map(lambda x: x[2], file_spots_tag), opts.taghi, opts.bl)

    lst_out = open(os.path.splitext(opts.out)[0]+".dat", "w")
    print >>lst_out, "# cheetah_dat=", cheetah_dat
    print >>lst_out, "# diffscan_log=", scanlog
    print >>lst_out, "# --min-spots=", opts.min_spots
    print >>lst_out, "# --default-energy=", opts.default_energy
    print >>lst_out, "# --comment=", opts.comment
    print >>lst_out, "file spots tag energy"
    for (f, s, t), e in zip(file_spots_tag, energies):
        print >>lst_out, f, s, t, e
    lst_out.close()

    n_h5files = int(math.ceil(float(len(file_spots_tag))/float(opts.max_in_h5)))
    idxes = range(len(file_spots_tag))
    splitter = lambda l: map(lambda x: l[x:x+opts.max_in_h5], xrange(0, len(idxes), opts.max_in_h5))
    file_spots_tag_sp, energies_sp = splitter(file_spots_tag), splitter(energies)
    
    if n_h5files > 0:
        easy_mp.pool_map(fixed_func=lambda x: make_h5("%s_%.4d.h5" % (opts.out[:-3], x+1), 
                                                      file_spots_tag_sp[x],
                                                      energies_sp[x],
                                                      opts.comment, opts.default_energy),
                         args=range(n_h5files),
                         processes=opts.nproc)

    """
    of = h5py.File(opts.out, "w")

    of["/metadata/detector"] = "Rayonix MX225HS"
    if opts.comment: of["/metadata/run_comment"] = opts.comment
    if file_spots_tag: of["/metadata/pixelsize_in_um"] = read_image(file_spots_tag[0][0], False)[0]["pixel_size"]*1000.
    #of["/metadata/distance_in_mm"] = opts.clen_mm

    for (f, s, tag), ene in zip(file_spots_tag, energies):
        print "converting", f, s, tag, ene
        h, data = read_image(f)
        grp = of.create_group("tag-%d"%tag)
        if ene!=ene: ene = opts.default_energy
        dset = grp.create_dataset("photon_energy_ev", (1,), dtype=numpy.float)
        dset[...] = ene * 1000.
        dset = grp.create_dataset("photon_wavelength_A", (1,), dtype=numpy.float)
        dset[...] = 12.3984 / ene
        dset = grp.create_dataset("original_file", (1,), "S%d"%len(f))
        dset[...] = f

        dset = grp.create_dataset("data", data.shape, dtype=data.dtype, compression="gzip", shuffle=True)
        dset[...] = data

    of.close()
    print "Processed:", opts.out
    """

    if opts.status_out:
        open(opts.status_out, "w").write("""\
# Cheetah status
Update time: %(ctime)s
Elapsed time: %(eltime)f sec
Status: Total=%(ntotal)d,Processed=%(ntotal)d,LLFpassed=%(ntotal)d,Hits=%(nhits)d,Status=Finished
Frames processed: %(ntotal)d
Number of hits: %(nhits)d
""" % dict(ctime=time.ctime(), eltime=time.time()-opts.eltime_from if opts.eltime_from is not None else -1, ntotal=n_proccessed_images, nhits=len(file_spots_tag)))

# run()

if __name__ == "__main__":
    #print "To make calculation faster, exiting! Bye."
    #quit()
    import optparse

    parser = optparse.OptionParser()
    parser.add_option("--min-spots", dest="min_spots", type=int,
                      help="Minimum number of spots")
    parser.add_option("--out", dest="out", type=str,
                      help="Output h5 file")
    parser.add_option("--default-energy", dest="default_energy", type=float)
    parser.add_option("--status-out", dest="status_out", type=str, default="status.txt")
    parser.add_option("--max-in-h5", dest="max_in_h5", type=int, default=200)
    parser.add_option("--nproc", dest="nproc", type=int, default=1)
    parser.add_option("--eltime-from", dest="eltime_from", type=float)
    parser.add_option("--corner-x", dest="corner_x", type=float)
    parser.add_option("--corner-y", dest="corner_y", type=float)
    parser.add_option("--rotated", dest="rotated", action="store_true", help="rotated MX300HS")
    parser.add_option("--geom-out", dest="geom_out", type=str)
    parser.add_option("--comment", dest="comment", type=str)
    parser.add_option("--taghi", dest="taghi", type=int)
    parser.add_option("--bl", dest="bl", type=int)

    opts, args = parser.parse_args()

    assert None not in (opts.default_energy, opts.min_spots, opts.out, opts.taghi, opts.bl)
    assert opts.out.endswith(".h5")

    cheetah_dat, diffscan_log = args
    
    run(opts, cheetah_dat, diffscan_log)
