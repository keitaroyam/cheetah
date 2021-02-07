# Cheetah dispatcher by Takanori Nakane
# modified for diffraction scan by Keitaro Yamashita

# This script must be edited for your environment,
# especially (1) values quoted by @@ and (2) job submission commands.
# At SACLA, the configured script is installed.

# If you have font problems with cctbx.python,
# try "export PHENIX_GUI_ENVIRONMENT=1"

"""
Modifications are being made for Rayonix detectors
Currently time-resolved stuff are not considered.

PHENIX_GUI_ENVIRONMENT=1 phenix.python /home/hirata/program/cheetah-biochem-fan/source/cheetah-rayonix/cheetah_dispatcher.py --data-root=/work/hirata/Mar2017 --clen=140


Raw images are uploaded into /work/hirata/** with serially-incremented directory names.
This script watches the latest serial number.

** run.info should look like:

Status   :Stopped (Ready to Read)
BeamLine :3
RunNumber:405432
Start    :2016/03/17 18:30:03.815567
Stop     :2016/03/17 18:32:55.459369
TriggerHi:201601
Trigger  :10300trigger 437433684-437443984
Comment  :2F_A_021

, where Comment is only used in this script.

** status.txt should look like:

# Cheetah status
Update time: Thu Mar 17 03:09:47 2016
Elapsed time: 0hr 10min 54sec
Status: Total=5000,Processed=5000,LLFpassed=1666,Hits=443,Status=Finished
Frames processed: 1666
Number of hits: 443

, where "Status:" line is used in this script.

TO be resubmitted:
337 342 354 357
371 373 376 378 379 391 394-
"""

import glob
import os
import os.path
import re
import sys
import time
import datetime
import threading
import traceback
#from subprocess import *
import subprocess
import commands
import optparse

import wx
import wx.grid
import wx.lib.newevent

from yamtbx.dataproc import bl_logfiles

PARALLEL_SIZE = 3 # MUST match hard-coded value in Cheetah

re_filename = re.compile("^[0-9]+(-dark[0-9]?|-light|-\d)?$")
re_status = re.compile("^Status:")
(ThreadEvent, EVT_THREAD) = wx.lib.newevent.NewEvent()
(SlogFoundEvent, EVT_SLOG_FOUND) = wx.lib.newevent.NewEvent()

job_script = '''\
#!/bin/bash
#PBS -l nodes=1:ppn=14
#PBS -e {runname}/{prefix}_cheetah.stderr
#PBS -o {runname}/{prefix}_cheetah.stdout
#PBS -N {runname}/{prefix}
#PBS -q {queuename}

cd $PBS_O_WORKDIR/{runname}

# This is the master job for runid, runid-light, runid-0.
# Subjobs must be submitted separatedly.

#if [ -e job.id ]; then
#   exit
#fi

echo $PBS_JOBID > {prefix}_job.id
# TODO run.info stuff
#cp -pv {datadir}/info.log .
#title=`grep Title: info.log | cut -b 8-`
#echo "Comment  :$title" > run.info
#echo "Comment  :" > run.info

eltime_from=`date +"%s"`
# XXX only prefix_[0-9]+.img should be processed
/usr/bin/time -o {prefix}_eltime.log --append -f "cheetah %E" /home/sacla_sfx_app/packages/rayonix/cheetah-diffscan/eiger-zmq/bin/cheetah.local_singles {datadir}/{prefix}_*.img --nproc=$NCPUS --output={prefix}_cheetah_spots.dat --min-snr=8 --dmin=5
/usr/bin/time -o {prefix}_eltime.log --append -f "convh5 %E" /home/sacla_sfx_app/packages/rayonix/dials-v1-10-2/build/bin/dials.python /home/sacla_sfx_app/packages/rayonix/cheetah-diffscan/conv_h5_with_energy.py \
        {prefix}_cheetah_spots.dat {diffscan_log} --min-spots={min_spots} --default-energy={default_energy} --out {prefix}_hits.h5 --status-out={prefix}_status.txt --eltime-from=$eltime_from --geom-out={prefix}.geom --bl={beamline} --nproc=$NCPUS --max-in-h5=100 --corner-x={corner_x} --corner-y={corner_y} {rotated_option}

R --vanilla <<+
d=read.table("{prefix}_hits_forplot.dat",h=T)
library(ggplot2)
p=ggplot(d, aes(x=x*1000,y=y*1000,fill=nspots)) +geom_tile() +scale_fill_gradient(low="white",high="red") +scale_x_continuous(expand=c(0,0)) +scale_y_continuous(expand=c(0,0)) +coord_fixed() +theme_bw() +labs(x="",y="",title="{runname}/{prefix}")
asp=(max(d\$x)-min(d\$x)) / (max(d\$y)-min(d\$y))
ggsave("{prefix}_nspots.pdf", p, width=3*asp, height=3, limitsize=F)
ggsave("{prefix}_nspots.png", p, width=3*asp, height=3, limitsize=F)
+

#exit

# th 100 gr 5000000 for > 10 keV
#/usr/bin/time -o {prefix}_eltime.log --append -f "crystfel %E" /home/sacla_sfx_app/local/bin/indexamajig -g {prefix}.geom --indexing=dirax-raw --peaks=zaef --threshold=400 --min-gradient=10000 --min-snr=5 --int-radius=3,4,7 -o {prefix}_out.stream -j $NCPUS -i - {crystfel_args} <<EOF
\ls {prefix}_hits*.h5 > {prefix}_hits_h5.lst
/usr/bin/time -o {prefix}_eltime.log --append -f "crystfel %E" /home/sacla_sfx_app/local/bin/indexamajig -g {prefix}.geom --indexing=dirax-raw --peaks=peakfinder8 --threshold=100 --int-radius=3,4,7 -o {prefix}_out.stream -j $NCPUS -i {prefix}_hits_h5.lst {crystfel_args}
rm -fr indexamajig.*
grep -c Cell {prefix}_out.stream > {prefix}_indexed.cnt

# compress files in /work
#find {datadir} -name {prefix}_??????.img | xargs -P$NCPUS -t -n1 bzip2 

##ruby @@SCRIPT_PATH@@/parse_stream.rb < {runname}.stream > {runname}.csv
'''

class LogWatcher(threading.Thread):
    #def __init__(self, filename, window, runid):
    def __init__(self, datadir, window, runid, cv):
        super(LogWatcher, self).__init__()
        #self.filename = filename
        self.datadir = datadir
        self.prefix = runid[1]
        #self.index_cnt = re.sub("status.txt", "indexed.cnt", self.filename)
        #self.run_info = re.sub("status.txt", "run.info", self.filename)
        self.status_txt = os.path.join(runid[0], "%s_status.txt"%self.prefix)
        self.jobid_file = os.path.join(runid[0], "%s_job.id"%self.prefix)
        self.index_cnt_file = os.path.join(runid[0], "%s_indexed.cnt"%self.prefix)
        self.window = window
        self.running = True
        self.cv = cv#threading.Condition()
        self.runid = runid
        self.comment = None
        self.update_scan()

        #print "HOGEEEEE",self.status_txt, os.path.isfile(self.status_txt)
        if os.path.isfile(self.jobid_file) and os.path.getsize(self.jobid_file):
            self.status = "running"
        else:
            self.status = "waiting"

        try:
            self.start()
        except: # To avoid "too many theads" error
            time.sleep(1) # wait "Finished" threads to be stopped
            self.start()

    def stop(self):
        self.running = False
        self.cv.acquire()
        self.cv.notify()
        self.cv.release()

    def update_scan(self):
        slog = bl_logfiles.BssDiffscanLog(os.path.join(self.datadir, "diffscan.log"))
        slog.remove_overwritten_scans()
        fltr = filter(lambda x: x.get_prefix()[:-1]==self.runid[1], slog.scans)
        self.scan = fltr[-1]

    def can_hit_find_start(self):
        # Count actually-exists
        exist = filter(lambda x: os.path.isfile(x),
                       map(lambda x: os.path.join(self.datadir, os.path.basename(x[0])),
                           self.scan.filename_coords))

        # Completed
        if len(exist) == self.scan.vpoints * self.scan.hpoints:
            print "All data exist!"
            return True

        # Check time-out
        timediff = (datetime.datetime.now() - self.scan.date).total_seconds()
        s = self.scan
        _hscan = s.scan_direction in ("horizontal",None) # True if horizontal-scan
        _readout = 0 # detector readout time
        #_turn = 10*(s.hpoints*s.hstep if _hscan else s.vpoints*s.vstep) # time to go to next line; assuming 10 mm/s
        _turn = 5 # time to go to next line (should be same when zig-zag..)
        time_est = ((_readout+1./s.frame_rate)*(s.hpoints if _hscan else s.vpoints)+_turn) * (s.vpoints if _hscan else s.hpoints)
        print "timediff=%.1f, time_est=%.1f" %(timediff, time_est)
        if timediff > time_est + 10*60 : # 10 minutes buffer
            return True

        return False

    def prepare_job(self):
        run_dir = self.runid[0]
        f = open("%s/%s_run.sh" % (run_dir, self.prefix), "w")
        opts = self.window.opts
        f.write(job_script.format(runname=run_dir, clen=opts.clen, queuename=opts.queue, prefix=self.prefix, diffscan_log=os.path.join(self.datadir, "diffscan.log"),
                                  datadir=self.datadir, crystfel_args=opts.crystfel_args, beamline=opts.bl, corner_x=opts.corner_x, corner_y=opts.corner_y,
                                  rotated_option="--rotated" if opts.rotated else "", default_energy=opts.default_energy, min_spots=opts.min_spots))
        f.close()

    def start_job(self):
        run_dir = self.runid[0]
        print "Starting job", self.runid
        os.system("qsub {rundir}/{prefix}_run.sh > {rundir}/{prefix}_job.id".format(rundir=run_dir, prefix=self.prefix))
        self.status = "running"

    def get_log(self):
        ret = None
        try:
            with open(self.status_txt) as f:
                for line in f:
                    if re_status.match(line):
                        ret = line[8:].rstrip() #  len("Status: ")
                        break
        except:
#            print "File not ready: %s" % self.filename
            return None
                
        return ret

    def parse_log(self, str):
        dict = {}
        try:
            key_vals = str.split(",")
            for kv in key_vals:
                key, val = kv.split("=")
                dict[key] = val
        except:
            return None
        return dict

    def parse_info(self):
        try:
            with open(self.run_info) as f:
                for line in f:
                    if line.startswith("Comment"):
                        self.comment = line[10:-1]
        except:
            return None

    def qstat(self):
        if not os.path.isfile(self.jobid_file): return None
        job_id = open(self.jobid_file).read().strip()
        p = subprocess.Popen(["qstat", job_id], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        outs, errs = p.communicate()
        #p.wait()
        if p.returncode != 0:
            return "finished"

    def count_actually_collected(self):
        slog = bl_logfiles.BssDiffscanLog(os.path.join(self.datadir, "diffscan.log"))
        slog.remove_overwritten_scans()
        fltr = filter(lambda x: x.get_prefix()[:-1]==self.runid[1], slog.scans)
        n = 0
        for f in glob.glob(os.path.join(self.datadir, "%s[0-9]*.img"%fltr[-1].get_prefix())):
            xy = slog.calc_grid_coord(filename=f)
            if xy[0]==xy[0] and xy[1]==xy[1]:
                n += 1

        return n

    def run(self):
        while self.running:
            self.cv.acquire()
            self.cv.wait(1.5)
            #self.cv.release()
            #print "doing", self.runid

            self.update_scan()

            tmp = {}
            tmp['runid'] = self.runid
            tmp["Total"] = self.scan.vpoints * self.scan.hpoints
            tmp["Collected"] = self.count_actually_collected() # len(self.scan.filename_coords)

            if self.status == "waiting":
                tmp["Status"] = "waiting"
                self.prepare_job()
                if self.can_hit_find_start():
                    self.start_job()

            elif self.status == "running":
                tmp2 = self.parse_log(self.get_log())
                if tmp2:
                    tmp["Status"] = tmp2["Status"]
                    tmp["Processed"] = tmp2["Processed"]
                    tmp["Hits"] = tmp2["Hits"]

                    if tmp['Status'].startswith("Finished"):
                        if os.path.isfile(self.index_cnt_file):
                            tmp['Status'] = "Finished"
                            try:
                                n = int(open(self.index_cnt_file).read().strip())
                                tmp['indexed'] = n
                            except:
                                tmp['indexed'] = 0
                        else:
                            tmp['Status'] = "Indexing"

                elif self.qstat()=="finished":
                    tmp["Status"] = "Stopped"
                else:
                    tmp["Status"] = "Running"

            event = ThreadEvent(msg=tmp)
            wx.PostEvent(self.window, event)
            if tmp['Status'].startswith(("Stopped", "Error", "Finished")):
                self.running = False # FIXME


            self.cv.release()
            continue

            #if self.comment is None:
            #    self.parse_info()

            tmp = self.parse_log(self.get_log())
            if tmp != None:
                tmp['runid'] = self.runid
                tmp['Comment'] = self.comment
                try:
                    f = open(self.index_cnt)
                    tmp['indexed'] = f.read()
                    f.close()
                except:
                    tmp['indexed'] = "NA"
                    if tmp['Status'].startswith("Finished"):
                        tmp['Status'] = "Indexing"
                event = ThreadEvent(msg=tmp)
                wx.PostEvent(self.window, event)
                if (tmp['indexed'] != "NA" and tmp['indexed'] != "") or tmp['Status'].startswith("Error"):
                    self.running = False # FIXME

        self.stop()

class ScanlogFinder:
    def __init__(self, parent, data_root):
        self.thread = None
        self.parent = parent
        self.data_root = data_root
        self.scanlog_cache = {}
    # __init__()

    def start(self, interval=5):
        self.stop()

        self.keep_going = True
        self.running = True
        if interval is not None:
            self.interval = interval

        self.thread = threading.Thread(None, self.run)
        self.thread.daemon = True
        self.thread.start()
        print "Started."

    def stop(self):
        if self.is_running():
            self.keep_going = False
            self.thread.join()
            print "Stopped!"

    def is_running(self):
        return self.thread is not None and self.thread.is_alive()

    def run(self):
        while self.keep_going:

            for root, dirnames, filenames in os.walk(self.data_root, followlinks=True):
                dirnames.sort()
                #print "looking at::", root
                if "diffscan.log" in filenames:
                    slogfile = os.path.join(root, "diffscan.log")
                    stat = os.stat(slogfile)
                    stat = stat.st_mtime, stat.st_size
                    if self.scanlog_cache.get(slogfile) == stat:
                        #print "Not updated:", slogfile
                        continue
                    self.scanlog_cache[slogfile] = stat

                    rootrel = os.path.relpath(root, self.data_root)
                    try:
                        slog = bl_logfiles.BssDiffscanLog(slogfile)
                        slog.remove_overwritten_scans()
                        if not os.path.exists(rootrel): os.makedirs(rootrel)
                        for scan in slog.scans:
                            prefix = scan.get_prefix()[:-1]

                            #self.addRun((rootrel, prefix))
                            print "notifying", (rootrel, prefix)
                            ev = SlogFoundEvent(run_id=(rootrel, prefix))
                            wx.PostEvent(self.parent, ev)
                    except:
                        print "Ignoring following error occurred with", slogfile
                        print traceback.format_exc()

            time.sleep(self.interval)


class AutoQsub(threading.Thread):
    LONG_WAIT = 20
    SHORT_WAIT = 1

    def __init__(self, queue="serial", maxjobs=14):
        super(AutoQsub, self).__init__()

        self.maxjobs = maxjobs
        self.cv = threading.Condition()
        self.running = True
        self.start()
        self.queue = queue

    def stop(self):
        self.running = False
        self.cv.acquire()
        self.cv.notify()
        self.cv.release()
        self.join()
        print "AutoQsub stopped."

    def checkAndSubmit(self, forceSubmit=False):
#        print "started directory scan ---"
        for filename in sorted(glob.glob("*"), reverse=True):
            if os.path.isdir(filename) and re_filename.match(filename) != None:
                dir = filename
                if os.path.exists(dir + "/run.sh") and not os.path.exists(dir + "/job.id"):
                    self.cv.acquire()
                    self.cv.wait(self.SHORT_WAIT) # check again
                    self.cv.release()

                    if os.path.exists(dir + "/job.id"):
                        print "WARNING: Another instance of auto_qsub is runnning?"
                        break
                    print "Unsubmitted job found in " + dir
                    njobs = int(commands.getoutput("qstat -w -u $USER | grep -c %s" % self.queue)) # FIXME: deprecated method 
                    if forceSubmit or njobs <= self.maxjobs:
                        print " submitted this job. %d jobs in the queue" % (njobs + 1)
                        os.system("qsub {dir}/run.sh > {dir}/job.id".format(dir=dir))
                    else:
                        print " the queue is full"
                        break

    def run(self):
        while self.running:
            self.checkAndSubmit(False)
            self.cv.acquire()
            self.cv.wait(self.LONG_WAIT)
            self.cv.release()

class MainWindow(wx.Frame):
    #COL_RUNID = 0
    COL_DIR = 0
    COL_PREFIX= 1
    COL_STATUS = 2
    COL_TOTAL = 3
    COL_COLLECTED = 4
    COL_PROCESSED = 5
    #COL_LLF_PASSED = 5
    COL_HITS = 6
    COL_INDEXED = 7
    COL_COMMENT = 8

    MENU_KILLJOB = 0
    MENU_HDFSEE = 1
    MENU_CELLEXPLORER = 2
    MENU_PLOT = 3
    MENU_COUNTSUMS = 4

    def __init__(self, parent, opts):
        self.opts = opts
        self.subthreads = []
        self.rows = {}
        self.scanlog_cache = {}
        self.waitFor = None
        self.autoqsub = None
        self.cv = threading.Condition()
        self.slogfinder = ScanlogFinder(self, opts.data_root)
        title = "Cheetah Dispatcher for *diffscan* on " + os.getcwd()
        wx.Frame.__init__(self, parent, title=title, size=(900,800))
        self.table = wx.grid.Grid(self, size=(900, -1))
        self.table.CreateGrid(0, 9, wx.grid.Grid.SelectRows) # row, column
        self.table.EnableEditing(False)
        self.table.SetRowLabelSize(0)
        # FIXME: better col width
        self.table.SetColLabelValue(MainWindow.COL_DIR,     "Directory     ")
        self.table.SetColLabelValue(MainWindow.COL_PREFIX,  "Prefix     ")
        self.table.SetColLabelValue(MainWindow.COL_STATUS,  "Status       ")
        self.table.SetColLabelValue(MainWindow.COL_TOTAL, "Total    ")
        self.table.SetColLabelValue(MainWindow.COL_COLLECTED, "Collected      ")
        self.table.SetColLabelValue(MainWindow.COL_PROCESSED, "Processed     ")
        self.table.SetColLabelValue(MainWindow.COL_HITS, "Hits          ")
        self.table.SetColLabelValue(MainWindow.COL_INDEXED, "Indexed       ")
        self.table.SetColLabelValue(MainWindow.COL_COMMENT, "Comment  ")
        attr = wx.grid.GridCellAttr()
        attr.SetRenderer(ProgressCellRenderer())
        self.table.SetColAttr(MainWindow.COL_PROCESSED, attr)
        self.table.SetColAttr(MainWindow.COL_COLLECTED, attr)
        self.table.SetColAttr(MainWindow.COL_HITS, attr)
        self.table.SetColAttr(MainWindow.COL_INDEXED, attr)
        for i in xrange(8):
            self.table.AutoSizeColLabelSize(i)
        self.table.Bind(wx.grid.EVT_GRID_CELL_RIGHT_CLICK, self.OnGridRightClick)
        self.table.Bind(wx.EVT_SIZE, self.OnGridResize)

        self.start_button = wx.Button(self, wx.ID_ANY, "  Start job(s)  ")
        self.start_button.Bind(wx.EVT_BUTTON, self.onPush)
        self.text_runid = wx.TextCtrl(self)
        self.text_runid.Disable()
        self.text_runid.SetValue(self.opts.data_root)
        self.label_runid = wx.StaticText(self, wx.ID_ANY, "Data root:")
        self.text_maxI = wx.TextCtrl(self, wx.ID_ANY, "0")
        self.text_maxI.Disable()
        self.label_maxI = wx.StaticText(self, wx.ID_ANY, "MaxI threshold:")
        self.combo_station = wx.ComboBox(self, wx.ID_ANY, value="ST4", choices=["ST2", "ST3", "ST4"],
                                         size=(80, -1), style=wx.CB_READONLY)
        self.combo_station.SetSelection(2) # somehow value= is ignored...
        self.combo_station.Disable()

        self.hsizer = wx.BoxSizer(wx.HORIZONTAL)
        self.hsizer.Add(self.label_runid, 0, wx.ALIGN_CENTER_VERTICAL)
        self.hsizer.Add(self.text_runid, 1, wx.EXPAND | wx.ALL, 5)
        self.hsizer.Add(self.label_maxI, 0, wx.ALIGN_CENTER_VERTICAL)
        self.hsizer.Add(self.text_maxI, 0, wx.ALL, 5)
        self.hsizer.Add(self.combo_station, 0, wx.ALIGN_CENTER_VERTICAL)
        self.hsizer.Add(self.start_button, 0, wx.ALL, 3)

        self.label_pd = wx.StaticText(self, wx.ID_ANY, "Dark/Light threshold (0 to disable)")
        self.text_pd1 = wx.TextCtrl(self, wx.ID_ANY, "0")
        self.text_pd1.Disable()
        self.label_pd1 = wx.StaticText(self, wx.ID_ANY, "")#"%s :" % self.opts.pd1_name)
        self.text_pd2 = wx.TextCtrl(self, wx.ID_ANY, "0")
        self.text_pd2.Disable()
        self.label_pd2 = wx.StaticText(self, wx.ID_ANY, "")#"%s:" % self.opts.pd2_name)
        self.text_pd3 = wx.TextCtrl(self, wx.ID_ANY, "0")
        self.text_pd3.Disable()
        self.label_pd3 = wx.StaticText(self, wx.ID_ANY, "")#"%S:" % self.opts.pd3_name)

        self.hsizer_pd = wx.BoxSizer(wx.HORIZONTAL)
        self.hsizer_pd.Add(self.label_pd, 0, wx.ALIGN_CENTER_VERTICAL)
        self.hsizer_pd1 = wx.BoxSizer(wx.HORIZONTAL)
        self.hsizer_pd1.Add(self.label_pd1, 0, wx.ALIGN_CENTER_VERTICAL)
        self.hsizer_pd1.Add(self.text_pd1, 0, wx.EXPAND | wx.ALL, 5)
        self.hsizer_pd2 = wx.BoxSizer(wx.HORIZONTAL)
        self.hsizer_pd2.Add(self.label_pd2, 0, wx.ALIGN_CENTER_VERTICAL)
        self.hsizer_pd2.Add(self.text_pd2, 0, wx.EXPAND | wx.ALL, 5)
        self.hsizer_pd3 = wx.BoxSizer(wx.HORIZONTAL)
        self.hsizer_pd3.Add(self.label_pd3, 0, wx.ALIGN_CENTER_VERTICAL)
        self.hsizer_pd3.Add(self.text_pd3, 0, wx.EXPAND | wx.ALL, 5)
        
        self.vsizer = wx.BoxSizer(wx.VERTICAL)
        self.vsizer.Add(self.hsizer, 0, wx.EXPAND | wx.RIGHT)
        self.vsizer.Add(self.hsizer_pd, 0, wx.EXPAND | wx.RIGHT)
        self.vsizer.Add(self.hsizer_pd1, 0, wx.EXPAND | wx.RIGHT)
        self.vsizer.Add(self.hsizer_pd2, 0, wx.EXPAND | wx.RIGHT)
        self.vsizer.Add(self.hsizer_pd3, 0, wx.EXPAND | wx.RIGHT)
        
        self.vsizer.Add(self.table, 1, wx.EXPAND | wx.ALL)

        self.Bind(EVT_THREAD, self.onUpdate)
        self.Bind(EVT_SLOG_FOUND, lambda x: self.addRun(x.run_id))
        self.Bind(wx.EVT_CLOSE, self.onClose)
        
        self.SetSizer(self.vsizer)
        self.SetAutoLayout(1)
#        self.vsizer.Fit(self)
        self.Show()

        # Monitor mode shows only the table
        if self.opts.monitor is not False:
            self.vsizer.Hide(self.hsizer)
            self.vsizer.Hide(self.hsizer_pd)
            self.vsizer.Hide(self.hsizer_pd1)
            self.vsizer.Hide(self.hsizer_pd1)
            self.vsizer.Hide(self.hsizer_pd1)
        # Show PD threshold only when specified
        if 1:#self.opts.pd1_name is None and self.opts.pd2_name is None and self.opts.pd3_name is None:
            self.vsizer.Hide(self.hsizer_pd)
        if 1:#self.opts.pd1_name is None:
            self.vsizer.Hide(self.hsizer_pd1)
        if 1:#self.opts.pd2_name is None:
            self.vsizer.Hide(self.hsizer_pd2)
        if 1:#self.opts.pd3_name is None:
            self.vsizer.Hide(self.hsizer_pd3)
        self.Layout()
        
        self.timer = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, self.OnTimer)
        #self.timer.Start(5000)

        self.scanDirectory()

        if self.opts.quick == 1: 
            self.startAutoQsub()

    def get_runid(self, row):
        wdir = self.table.GetCellValue(row, self.COL_DIR)
        prefix = self.table.GetCellValue(row, self.COL_PREFIX)
        return (wdir, prefix)

    def startAutoQsub(self):
        self.autoqsub = AutoQsub(self.opts.queue, self.opts.max_jobs)
        
    def OnGridResize(self, event):
        # Reference: https://forums.wxwidgets.org/viewtopic.php?t=90
        sum_width = self.table.GetRowLabelSize()
        ncols = self.table.GetNumberCols()
        for i in xrange(ncols - 1):
            sum_width += self.table.GetColSize(i)

        target = self.table.GetClientSize().GetWidth() - sum_width
        self.table.SetColSize(ncols - 1, target)
        self.table.ForceRefresh()
        event.Skip()

    def OnGridRightClick(self, event):
        row = event.GetRow()
        runid = self.get_runid(row) # self.table.GetCellValue(row, self.COL_RUNID)
        
        point = event.GetPosition()
        popupmenu = wx.Menu()
        killjob = popupmenu.Append(MainWindow.MENU_KILLJOB, "Kill this job")
        hdfsee = popupmenu.Append(MainWindow.MENU_HDFSEE, "View hits")
        hdfsee = popupmenu.Append(MainWindow.MENU_CELLEXPLORER, "Check cell")
        plot = popupmenu.Append(MainWindow.MENU_PLOT, "View grid")
        hdfsee = popupmenu.Append(MainWindow.MENU_COUNTSUMS, "Count sums")
        self.table.Bind(wx.EVT_MENU, lambda event: self.OnPopupmenuSelected(event, runid))
#        saturation_hist = popupmenu.Append(runid, "Saturation histogram")
#        self.table.Bind(wx.EVT_MENU, self.OnSaturationHistSelected, saturation_hist)
        self.table.PopupMenu(popupmenu, point)

    def OnPopupmenuSelected(self, event, runid):
        id = event.GetId()
    
        if id == MainWindow.MENU_KILLJOB:
            self.KillJob(runid)
        elif id == MainWindow.MENU_HDFSEE:
            self.HDFsee(runid)
        elif id == MainWindow.MENU_CELLEXPLORER:
            self.CellExplorer(runid)
        elif id == MainWindow.MENU_PLOT:
            plot_file = os.path.join(self.get_data_dir(runid, False), runid[1]+"_nspots.png")
            if os.path.isfile(plot_file):
                subprocess.Popen(["eog %s"%plot_file], shell=True)
        elif id == MainWindow.MENU_COUNTSUMS:
#            GetSelectionBlock no longer works on newer wxPython as before...
#             https://forums.wxwidgets.org/viewtopic.php?t=41607
#            row1 = self.table.GetSelectionBlockTopLeft()[0][0]
#            row2 = self.table.GetSelectionBlockBottomRight()[0][0]
            self.CountSums(self.table.GetSelectedRows())

    def CountSums(self, rows): # XXX
        total = {}
        processed = {}
        accepted = {}
        hits = {}
        indexed = {}
        types = ("normal", "light", "dark") + tuple("dark%d" % i for i in xrange(1, 10))
        for t in types:
            for x in (total, processed, accepted, hits, indexed):
                x[t] = 0
        for row in rows:
            text = self.table.GetCellValue(row, 0)
            typ = "normal"
            for t in types:
                if text.endswith(t):
                    typ = t
            try:
                total[typ] += int(self.table.GetCellValue(row, self.COL_TOTAL))
                processed[typ] += int(self.table.GetCellValue(row, self.COL_PROCESSED).rsplit(" ")[0])
                accepted[typ] += 0# XXX #  int(self.table.GetCellValue(row, self.COL_LLF_PASSED).rsplit(" ")[0])
                hits[typ] += int(self.table.GetCellValue(row, self.COL_HITS).rsplit(" ")[0])
                indexed[typ] += int(self.table.GetCellValue(row, self.COL_INDEXED).rsplit(" ")[0])
            except:
                pass
        message = ""
        for t in types:
            if total[t] != 0:
                if message != "":
                    message += "\n"
                message += "Type: %s\nTotal: %d\nProcessed: %d\nAccepted: %d\nHits: %d" % (t, total[t], processed[t], accepted[t], hits[t])
                if accepted[t] != 0:
                    message += " (%.1f%% of accepted)" % (100.0 * hits[t] / accepted[t])
                message += "\nIndexed: %d " % indexed[t]
                if hits[t] != 0:
                    message += " (%.1f%% of hits)" % (100.0 * indexed[t] / hits[t])
                message += "\n"                
        dlg = wx.MessageDialog(None, message, "Summary", wx.OK)
        dlg.ShowModal()
        dlg.Destroy()

    def HDFsee(self, runid):
        # FIXME: warn when multiple rows are selected
        def launchHDFsee():
            ddir = self.get_data_dir(runid, False)
            if os.path.exists(os.path.join(ddir, "%s_out.stream"%runid[1])):
                command = "cd {ddir}; /home/sacla_sfx_app/local/bin/hdfsee -g {prefix}.geom {prefix}_out.stream -c invmono -i 10 &".format(ddir=ddir,prefix=runid[1])
            else:
                command = "/home/sacla_sfx_app/local/bin/hdfsee -g {prefix}.geom {ddir}/{prefix}_hits.h5 -c invmono -i 10 &".format(ddir=ddir,prefix=runid[1])
            os.system(command)

        threading.Thread(target=launchHDFsee).start()

    def CellExplorer(self, runid):
        ddir = self.get_data_dir(runid, False)
        # FIXME: warn when multiple rows are selected
        def launchCellExplorer():
            command = "cd {ddir}; /home/sacla_sfx_app/local/bin/cell_explorer {prefix}_out.stream &".format(ddir=ddir,prefix=runid[1])
            os.system(command)
            
        if os.path.exists("%s/%s_out.stream" % (ddir, runid[1])):
            t = threading.Thread(target=launchCellExplorer).start()
        else:
            self.showError("Indexing result for %s/%s is not available (yet)." % runid)

    def KillJob(self, runid):
        message = "Are you sure to kill job %s/%s?" % runid
        dlg = wx.MessageDialog(None, message, "Cheetah dispatcher", wx.YES_NO | wx.NO_DEFAULT)
        ret = dlg.ShowModal()
        dlg.Destroy()
        if ret == wx.ID_NO:
            return
        try:
            jobid = open("%s/%s_job.id" % (self.get_data_dir(runid, False), runid[1])).read()
            m = re.match("([0-9]+)", jobid)
            jobid = m.group(1)
            os.system("qdel %s" % jobid)
#            f = open("d/status.txt", "w")
#            f.write("Status: Status=Killed") # TODO: should better change status only
            import shutil
#            shutil.rmtree("%s/" % runid)
        except:
            self.showError("Failed to delete run %s." % runid)
            print traceback.format_exc()

    def scanDirectory(self):
        target = "."
        for root, dirnames, filenames in os.walk(target, followlinks=True):
            dirnames.sort()
            sel = filter(lambda x: x.endswith("_job.id"), filenames)
            for f in sel:
                prefix = f[:-len("_job.id")]
                rootrel = os.path.relpath(root, target)
                self.addRun((rootrel, prefix))
        #for dirname in sorted(glob.glob("*/*/")):
        #    runid = dirname.split("/")[1]
        #    self.addRun(runid)
        ## TODO: Support deleting existing rows. This is more complex since the row idx changes.

    def scanDataDirectory(self):
        for root, dirnames, filenames in os.walk(self.opts.data_root, followlinks=True):
            dirnames.sort()
            #print "looking at::", root
            if "diffscan.log" in filenames:
                slogfile = os.path.join(root, "diffscan.log")
                stat = os.stat(slogfile)
                stat = stat.st_mtime, stat.st_size
                if self.scanlog_cache.get(slogfile) == stat:
                    #print "Not updated:", slogfile
                    continue
                self.scanlog_cache[slogfile] = stat

                rootrel = os.path.relpath(root, self.opts.data_root)
                try:
                    slog = bl_logfiles.BssDiffscanLog(slogfile)
                    slog.remove_overwritten_scans()
                    if not os.path.exists(rootrel): os.makedirs(rootrel)
                    for scan in slog.scans:
                        prefix = scan.get_prefix()[:-1]
                        self.addRun((rootrel, prefix))
                except:
                    print "Ignoring following error occurred with", slogfile
                    print traceback.format_exc()


    def onPush(self, event):
        if self.timer.IsRunning():
            self.stopWatch()
        else:
            self.startWatch()

        #if self.waitFor != None:
        #    self.stopWatch()
        #else:
        #    self.prepareSubmitJob(self.text_runid.GetValue())

    def stopWatch(self):
        self.waitFor = None
        self.start_button.SetLabel("  Start job(s)  ")
        #self.text_runid.Enable()
        ##self.text_pd1.Enable()
        ##self.text_pd2.Enable()
        ##self.text_pd3.Enable()
        ##self.text_maxI.Enable()
        ##self.combo_station.Enable()
        #self.timer.Stop()
        self.slogfinder.stop()

    def startWatch(self):#From(self, runid):
        #self.waitFor = runid
        self.start_button.SetLabel("Stop automode")
        self.text_runid.Disable()
        self.text_pd1.Disable()
        self.text_pd2.Disable()
        self.text_pd3.Disable()
        self.text_maxI.Disable()
        self.combo_station.Disable()
        #self.OnTimer(None)
        #self.timer.Start(5000)
        self.slogfinder.start()

    def get_data_dir(self, n, start_with_root=True):
        """
        if start_with_root==True, return the original data location
        otherwise, rern the working directory location
        """
        wdir, prefix = n
        if start_with_root:
            return os.path.join(self.opts.data_root, wdir)
        else:
            return wdir
        #reld = os.path.join("%.3d"%(n//500+1), "%.6d"%n)
        #return os.path.join(self.opts.data_root, reld) if start_with_root else reld

    def writeCSV(self):
        import csv
        with open('cheetah.csv', 'w') as f:
            writer = csv.writer(f, lineterminator='\n')
            writer.writerow(map(lambda x: self.table.GetColLabelValue(x), xrange(self.table.GetNumberCols())))
            for i in xrange(self.table.GetNumberRows()):
                writer.writerow(map(lambda x: x.split()[0] if x else "",
                                    map(lambda x: self.table.GetCellValue(i, x), xrange(self.table.GetNumberCols()))))

    def OnTimer(self, event):
        #self.scanDirectory()
        #print "onTimer!!"
        self.scanDataDirectory()
        #if (self.waitFor == None):
        #    return

        ## TODO check the latest directory and find ready ones
        ##out = Popen(["ShowRunInfo", "-b", "%d" % self.opts.bl, "-r", "%d" % self.waitFor], stdout=PIPE).communicate()[0]
        #data_dir = self.get_data_dir(self.waitFor)
        #if os.path.isdir(self.get_data_dir(self.waitFor+1)) or len(glob.glob(os.path.join(data_dir, "*.img"))) >= 1000:
        #    print "\rRun %d became ready." % self.waitFor
        #    self.prepareSubmitJob("%d" % self.waitFor)
        #    self.waitFor = self.waitFor + 1
        #    self.text_runid.SetValue("%d-" % self.waitFor)
        self.writeCSV()

    def showError(self, message):
        dlg = wx.MessageDialog(None, message, "Cheetah dispatcher", wx.OK)
        dlg.ShowModal()
        dlg.Destroy()

    def addRun(self, runid):
        if self.rows.get(runid) is not None:
            return
        if not os.path.isfile(os.path.join(self.get_data_dir(runid, True), "diffscan.log")):
            return
        #try: int(runid)
        #except ValueError: return

        row = self.table.GetNumberRows()
        self.subthreads.append(LogWatcher(self.get_data_dir(runid, True), self, runid, self.cv))
        self.rows[runid] = row
        self.table.AppendRows(1)
        #self.table.SetCellValue(row, MainWindow.COL_RUNID, runid)
        self.table.SetCellValue(row, MainWindow.COL_DIR, runid[0])
        self.table.SetCellValue(row, MainWindow.COL_PREFIX, runid[1])
        self.table.SetCellValue(row, MainWindow.COL_STATUS, "waiting")

        #self.writeCSV()

    def onClose(self, event):
        self.writeCSV()
        dlg = wx.MessageDialog(self, "Do you really want to exit? Submitted jobs are kept running.",
                               "Exit", wx.YES_NO|wx.ICON_QUESTION)
        result = dlg.ShowModal()
        dlg.Destroy()
        if result != wx.ID_YES:
            return
        
        # WaitFor is implemented using wxTimer, so we don't have to worry about it.
        for th in self.subthreads:
            th.stop()
            th.join()
        if self.autoqsub != None:
            self.autoqsub.stop()
            dlg = wx.MessageDialog(None, "Do you want to submit all remaining jobs (if any) now?\nWarning: This may delay execution of newer jobs.", "Cheetah dispatcher", wx.YES_NO | wx.NO_DEFAULT)
            ret = dlg.ShowModal()
            dlg.Destroy()
            if ret == wx.ID_YES:
                self.autoqsub.checkAndSubmit(True)

        self.Destroy()
        return True

    def onUpdate(self, event):
        # TODO fix this for rayonix
        row = self.rows.get(event.msg['runid'])
        if row == None:
            return

        comment = event.msg.get('Comment')
        if comment:
            self.table.SetCellValue(row, MainWindow.COL_COMMENT, comment)

        status = event.msg['Status']
        self.table.SetCellValue(row, MainWindow.COL_STATUS, status)
        if (status.startswith("Error")):
            return 

        total = int(event.msg['Total'])
        collected = int(event.msg.get('Collected', 0))
        processed = int(event.msg.get('Processed', 0))

        self.table.SetCellValue(row, MainWindow.COL_TOTAL, 
                                "%d" % total)
        self.table.SetCellValue(row, MainWindow.COL_COLLECTED, 
                                "%d (%.1f%%)" % (collected, 100.*collected/total))
        if collected > 0:
            self.table.SetCellValue(row, MainWindow.COL_PROCESSED, 
                                    "%d (%.1f%%)" % (processed, 100.0 * processed / collected))

        if (status == "DarkAveraging" or status == "waiting"):
            return
        #LLFpassed = int(event.msg['LLFpassed'])
        hits = int(event.msg.get('Hits', 0))
        #self.table.SetCellValue(row, MainWindow.COL_LLF_PASSED,
        #                        "%d (%.1f%%)" % (LLFpassed, 100.0 * LLFpassed / processed))
        if processed > 0:
            self.table.SetCellValue(row, MainWindow.COL_HITS,
                                    "%d (%.1f%%)" % (hits, 100.0 * hits / processed))

        try:
            indexed = int(event.msg['indexed'])
            if hits == 0:
                self.table.SetCellValue(row, MainWindow.COL_INDEXED, "0 (0.0%)")
            else:
                self.table.SetCellValue(row, MainWindow.COL_INDEXED,
                                            "%d (%.1f%%)" % (indexed, 100.0 * indexed / hits))
        except:
            self.table.SetCellValue(row, MainWindow.COL_INDEXED, "NA")

# Based on https://groups.google.com/forum/#!topic/wxpython-dev/MRVAPULwG70
#          http://screeniqsys.com/blog/2009/11/01/progresscellrenderer-gauges-for-grid-cells/
#       by David Watson
class ProgressCellRenderer(wx.grid.PyGridCellRenderer):
    re_percent = re.compile('([0-9.]+)%')
    def __init__(self):
        wx.grid.PyGridCellRenderer.__init__(self)
        self.progColor = wx.Colour(124, 252, 0)

    def Draw(self, grid, attr, dc, rect, row, col, isSelected):
        text = grid.GetCellValue(row, col)
        prog = 0
        try:
            prog = float(ProgressCellRenderer.re_percent.search(text).group(1))
        except:
            pass

        if prog > 100:
            prog = 100

        total_width = rect.width
        units = int(float(prog / 100.0) * rect.width)

        prog_rect = wx.Rect(rect.x, rect.y, units, rect.height)
        other_rect = wx.Rect(rect.x + units, rect.y, rect.width - prog_rect.width, rect.height)
        hAlign, vAlign = attr.GetAlignment()

        dc.SetFont(attr.GetFont())

        if isSelected:
            bg = grid.GetSelectionBackground()
        else:
            bg = wx.Colour(255, 255, 255)

        dc.SetBrush(wx.Brush(self.progColor, wx.SOLID))
        dc.SetPen(wx.BLACK_PEN)
        dc.DrawRectangleRect(prog_rect)

        # This fills in the rest of the cell background so it doesn't shear
        dc.SetBrush(wx.Brush(bg, wx.SOLID))
        dc.SetPen(wx.TRANSPARENT_PEN)
        dc.DrawRectangleRect(other_rect)

        dc.SetTextForeground(wx.Colour(0, 0, 0))
        grid.DrawTextRectangle(dc, text, rect, hAlign, vAlign)

    def GetBestSize(self, grid, attr, dc, row, col):
        text = grid.GetCellValue(row, col)
        dc.SetFont(attr.GetFont())
        w, h = dc.GetTextExtent(text)
        return wx.Size(w, h)

    def Clone(self):
        return ProgressCellRenderer() 

print
print "Cheetah dispatcher GUI version 2017/02/18"
print "   by Takanori Nakane (takanori.nakane@bs.s.u-tokyo.ac.jp)"
print "with dirty and quick modifications for Rayonix MX225HS"
print "   by Keitaro Yamashita"
print
print "Please cite the following paper when you use this software."
print " \"Data processing pipeline for serial femtosecond crystallography at SACLA\""
print " Nakane et al., J. Appl. Cryst. (2016). 49"
print

##if not os.path.exists("sacla-photon.ini"):
##    sys.stderr.write("ERROR: Configuration file was not found!\n\n")
##    sys.stderr.write("You should copy @@TEMPLATE_FILE@@ into this directory\n")
##    sys.stderr.write("and confirm the settings.\n")
##    sys.exit(-1)

parser = optparse.OptionParser()
parser.add_option("--monitor", dest="monitor", type=int, default=False, help="Monitor only")
parser.add_option("--bl", dest="bl", type=int, default=3, help="Beamline")
parser.add_option("--clen", dest="clen", type=float, default=None, help="camera length in mm")
parser.add_option("--quick", dest="quick", type=int, default=False, help="enable quick mode")
parser.add_option("--queue", dest="queue", type=str, default="serial", help="queue name")
parser.add_option("--max_jobs", dest="max_jobs", type=int, default=14, help="maximum number of jobs to submit when --quick is enabled")
#parser.add_option("--pd1_name", dest="pd1_name", type=str, default=None, help="PD1 sensor name e.g. xfel_bl_3_st_4_pd_laser_fitting_peak/voltage")
# e.g. xfel_bl_3_st_4_pd_user_10_fitting_peak/voltage
#parser.add_option("--pd2_name", dest="pd2_name", type=str, default=None, help="PD2 sensor name")
#parser.add_option("--pd3_name", dest="pd3_name", type=str, default=None, help="PD3 sensor name")
#parser.add_option("--submit_dark2", dest="submit_dark2", type=int, default=False, help="(DEPRECATED) accepts second darks (Ln-D2) and divide into light, dark1 and dark2")
#parser.add_option("--submit_dark_to", dest="submit_dark_to", type=int, default=False, help="accepts up to M (<=9) darks (Ln-Dm) and divide into light, dark1, ..., darkM")
#parser.add_option("--submit_dark_any", dest="submit_dark_any", type=int, default=False, help="accepts any lights and darks (Ln-Dm) and divide into light and dark")
parser.add_option("--crystfel_args", dest="crystfel_args", type=str, default="", help="optional arguments to CrystFEL")
parser.add_option("--data-root", dest="data_root", type=str, help="Root directory to find img files")
parser.add_option("--corner-x", dest="corner_x", type=float, help="override beam center", default=float("nan"))
parser.add_option("--corner-y", dest="corner_y", type=float, help="override beam center", default=float("nan"))
parser.add_option("--rotated", dest="rotated", action="store_true", help="rotated MX300HS")
parser.add_option("--default-energy", dest="default_energy", type=float, help="default energy for h5 (keV)", default=10.)
parser.add_option("--min-spots", dest="min_spots", type=int, help="Minimum spots for hit", default=10)

opts, args = parser.parse_args()

#if opts.submit_dark2 is not False:
#    sys.stderr.write("WARNING: submit_dark2 has been deprecated. Use --submit-dark-to=2 instead.\n")
#    sys.stderr.write("         --submit-dark-to was set to 2 for you.\n")
#    opts.submit_dark_to = 2

#if opts.submit_dark_to is not False and opts.submit_dark_any == 1:
#    sys.stderr.write("ERROR: You cannot enable submit_dark_any and set submit_dark_to simultaneously!\n")
#    sys.exit(-1)

#if opts.submit_dark_to is False:
#    opts.submit_dark_to = 1

#if opts.submit_dark_to > 9 or opts.submit_dark_to < 1:
#    sys.stderr.write("ERROR: submit_dark_to must be within 1 to 9.\n")
#    sys.exit(-1)

if opts.monitor is not False and opts.quick is not False:
    sys.stderr.write("ERROR: You cannot enable 'quick' when 'monitor' mode is enabled.\n")
    sys.exit(-1)

if opts.max_jobs < 1:
    sys.stderr.write("ERROR: max_jobs must be a positive integer.\n")
    sys.exit(-1)

if opts.bl != 2 and opts.bl != 3:
    sys.stderr.write("ERROR: beamline must be 2 or 3.\n")
    sys.exit(-1)

if not opts.data_root or not os.path.isdir(opts.data_root):
    sys.stderr.write("ERROR: invalid or unspecified --data-root directory name\n")
    sys.exit(-1)

opts.data_root = os.path.abspath(opts.data_root)

#print "Option: monitor          = %s" % opts.monitor
print "Option: bl               = %d" % opts.bl
if opts.clen: print "Option: clen             = %f mm" % opts.clen
print "Option: quick            = %s" % opts.quick
print "Option: max_jobs         = %s" % opts.max_jobs
print "Option: queue            = %s" % opts.queue
#print "Option: pd1_name         = %s" % opts.pd1_name
#print "Option: pd2_name         = %s" % opts.pd2_name
#print "Option: pd3_name         = %s" % opts.pd3_name
#print "Option: submit_dark_to   = %s" % opts.submit_dark_to
#print "Option: submit_dark_any  = %s" % opts.submit_dark_any
print "Option: crystfel_args    = %s" % opts.crystfel_args

app = wx.App(False)
frame = MainWindow(None, opts)
app.MainLoop()

