import os
import glob
import zmq
import json
import pickle
import time
from yamtbx.dataproc import cbf

if __name__ == "__main__":
    context = zmq.Context()
    sender = context.socket(zmq.PUSH)
    sender.bind("tcp://127.0.0.1:9999")

    #x = pickle.load(open("test_bslz4.pkl"))
    files = glob.glob("/dev/shm/shika-pilatus-test/*.cbf")
    t0 = time.time()
    #header = dict(beam_center_x=1230.00, beam_center_y=1296.00, pixel_size_x=0.172, wavelength=1, distance=400,
    #              file_prefix="scan02", frame=0,
    #              raster_horizontal_number=100, raster_vertical_number=100,
    #              raster_horizontal_step=0.01, raster_vertical_step=0.015,
    #              raster_scan_direction="horizontal", raster_scan_path="zig-zag")
    header = dict(file_prefix="scan02", frame=0,
                  raster_horizontal_number=100, raster_vertical_number=100,
                  raster_horizontal_step=0.01, raster_vertical_step=0.015,
                  raster_scan_direction="horizontal", raster_scan_path="zig-zag")
    for i, f in enumerate(files):
        header["frame"] = i
        #data = dict(imgfile=f[:f.rindex(".cbf")]+".img", header=header, cbf_data=open(f, "rb").read()) # GUI only accepts .img
        data = dict(imgfile=f, header=header, cbf_data=open(f, "rb").read()) # GUI only accepts .img
        #data, ndimfast, ndimmid = cbf.load_minicbf_as_numpy(f)
        #data = data.reshape(ndimmid, ndimfast)
        #data = dict(f=f, data=data)
        print time.time(), "sending", os.path.basename(f)
        sender.send_pyobj(data)

    print (time.time()-t0)/len(files), "sec/frame"


    quit()

    for f in glob.glob("lyz/*.pkl"):
        data = pickle.load(open(f))
        print "sending", f
        sender.send_multipart(data, copy=False)
        #time.sleep(2)
	#break

