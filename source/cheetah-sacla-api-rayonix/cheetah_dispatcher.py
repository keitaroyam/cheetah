# Cheetah dispatcher by Takanori Nakane
# modified for Rayonix detector with SACLA API by Keitaro Yamashita

# This script must be edited for your environment,
# especially (1) values quoted by @@ and (2) job submission commands.
# At SACLA, the configured script is installed.

# If you have font problems with cctbx.python,
# try "export PHENIX_GUI_ENVIRONMENT=1"

import glob
import os
import os.path
import re
import sys
import time
import threading
import traceback
from subprocess import *
import commands
import optparse

import wx
import wx.grid
import wx.lib.newevent

PARALLEL_SIZE = 3 # MUST match hard-coded value in Cheetah

re_filename = re.compile("^[0-9]+(-dark[0-9]?|-light|-\d)?$")
re_status = re.compile("^Status:")
(ThreadEvent, EVT_THREAD) = wx.lib.newevent.NewEvent()

job_script = '''#!/bin/bash
#PBS -l nodes=1:ppn=14
#PBS -e {runname}/cheetah.stderr
#PBS -o {runname}/cheetah.stdout
#PBS -N {runname}
#PBS -q {queuename}

cd $PBS_O_WORKDIR/{runname}

# This is the master job for runid, runid-light, runid-0.
# Subjobs must be submitted separatedly.

#if [ -e job.id ]; then
#   exit
#fi

echo $PBS_JOBID > job.id
hostname > job.host
source @@SETUP_SCRIPT@@
ShowRunInfo -b {beamline} -r {runid} > run.info
@@CHEETAH_PATH@@/prepare-cheetah-sacla-api2.py {runid} --bl={beamline} --clen={clen}
grep Error status.txt
if [ $? -eq 0 ]; then # Found
   for i in {subjobs}; do
      ln -s ../{runname}/status.txt ../{runid}-$i/
   done
   exit
fi

ln -s {runid}-geom.h5 sacla-geom.h5
ln -s {runid}-dark.h5 sacla-dark.h5

for i in {subjobs}; do
   if [ ! -e ../{runid}-$i/metadata.h5 ]; then
      cp {runid}.h5 ../{runid}-$i/
   fi
   ln -s ../{runname}/{runid}-geom.h5 ../{runid}-$i/sacla-geom.h5
   ln -s ../{runname}/{runid}.geom ../{runid}-$i/{runid}-$i.geom
   ln -s ../{runname}/run.info ../{runid}-$i/run.info
   ln -s ../{runname}/{runid}-dark.h5 ../{runid}-$i/sacla-dark.h5
done

if [ ! -e run{runname}.h5 ]; then
   cp {runid}.h5 run{runname}.h5
fi

@@CHEETAH_PATH@@/cheetah-sacla-api2 --ini ../sacla-photon.ini --run {runid} -o run{runname}.h5 --bl {beamline} {arguments}
rm {runid}.h5

# th 100 gr 5000000 for > 10 keV
@@INDEXAMAJIG_PATH@@/indexamajig -g {runid}.geom --indexing=dirax --peaks=zaef --threshold=400 --min-gradient=10000 --min-snr=5 --int-radius=3,4,7 -o {runname}.stream -j 14 -i - {crystfel_args} <<EOF
run{runname}.h5
EOF
rm -fr indexamajig.*
grep Cell {runname}.stream | wc -l > indexed.cnt
ruby @@SCRIPT_PATH@@/parse_stream.rb < {runname}.stream > {runname}.csv

rm job.id job.host
'''

job_script_dark = '''#!/bin/bash
#PBS -l nodes=1:ppn=14
#PBS -e {runname}/cheetah.stderr
#PBS -o {runname}/cheetah.stdout
#PBS -N {runname}
#PBS -q {queuename}

cd $PBS_O_WORKDIR/{runname}/

#if [ -e job.id ]; then
#   exit
#fi

echo $PBS_JOBID > job.id
hostname > job.host
source @@SETUP_SCRIPT@@

i=0
while :; do
   let i=i+1
   if [ -e sacla-dark.h5 ]; then
      break
   fi
   grep Error status.txt
   if [ $? -eq 0 ]; then
      exit
   fi

   if [ $i -gt 500 ]; then
      echo "Status: Status=Error-TimeoutWaitingDarkAverage" > status.txt
      exit -1
   fi

   sleep 2
done

cp {runid}.h5 run{runname}.h5
@@CHEETAH_PATH@@/cheetah-sacla-api2 --ini ../sacla-photon.ini --run {runid} -o run{runname}.h5 --bl {beamline} {arguments}
rm {runid}.h5

# th 100 gr 5000000 for > 10 keV
@@INDEXAMAJIG_PATH@@/indexamajig -g {runname}.geom --indexing=dirax --peaks=zaef --threshold=400 --min-gradient=10000 --min-snr=5 --int-radius=3,4,7 -o {runname}.stream -j 14 -i - {crystfel_args} <<EOF
run{runname}.h5
EOF
rm -fr indexamajig.*
grep Cell {runname}.stream | wc -l > indexed.cnt
ruby @@SCRIPT_PATH@@/parse_stream.rb < {runname}.stream > {runname}.csv

rm job.id job.host
'''

# JOB SCRIPTS FOR RAYONIX
job_script_hs = '''#!/bin/bash
#PBS -l nodes=1:ppn=14
#PBS -e {runname}/cheetah.stderr
#PBS -o {runname}/cheetah.stdout
#PBS -N {runname}
#PBS -q {queuename}

cd $PBS_O_WORKDIR/{runname}

# This is the master job for runid, runid-light, runid-0.
# Subjobs must be submitted separatedly.

#if [ -e job.id ]; then
#   exit
#fi

echo "PBS_JOBID=$PBS_JOBID\nhostname=`hostname`"
echo $PBS_JOBID > job.id
hostname > job.host
source @@SETUP_SCRIPT@@
ShowRunInfo -b {beamline} -r {runid} > run.info

dials.python @@CHEETAH_PATH@@/cheetah_marccd.py --run {runid} -o run{runname}.h5 --bl {beamline} --nproc $NCPUS {arguments}


# th 100 gr 5000000 for > 10 keV
@@INDEXAMAJIG_PATH@@/indexamajig -g {runid}.geom --indexing=dirax --peaks=zaef --threshold=400 --min-gradient=10000 --min-snr=5 --int-radius=3,4,7 -o {runname}.stream -j 14 -i - {crystfel_args} <<EOF
run{runname}.h5
EOF
rm -fr indexamajig.*
grep Cell {runname}.stream | wc -l > indexed.cnt
ruby @@SCRIPT_PATH@@/parse_stream.rb < {runname}.stream > {runname}.csv

rm job.id job.host
'''

job_script_dark_hs = '''#!/bin/bash
#PBS -l nodes=1:ppn=14
#PBS -e {runname}/cheetah.stderr
#PBS -o {runname}/cheetah.stdout
#PBS -N {runname}
#PBS -q {queuename}

cd $PBS_O_WORKDIR/{runname}/

#if [ -e job.id ]; then
#   exit
#fi

echo "PBS_JOBID=$PBS_JOBID\nhostname=`hostname`"
echo $PBS_JOBID > job.id
hostname > job.host
source @@SETUP_SCRIPT@@

dials.python @@CHEETAH_PATH@@/cheetah_marccd.py --run {runid} -o run{runname}.h5 --bl {beamline} --nproc $NCPUS {arguments}

# th 100 gr 5000000 for > 10 keV
@@INDEXAMAJIG_PATH@@/indexamajig -g {runid}.geom --indexing=dirax --peaks=zaef --threshold=400 --min-gradient=10000 --min-snr=5 --int-radius=3,4,7 -o {runname}.stream -j 14 -i - {crystfel_args} <<EOF
run{runname}.h5
EOF
rm -fr indexamajig.*
grep Cell {runname}.stream | wc -l > indexed.cnt
ruby @@SCRIPT_PATH@@/parse_stream.rb < {runname}.stream > {runname}.csv

rm job.id job.host
'''


class LogWatcher(threading.Thread):
    def __init__(self, filename, window, runid):
        super(LogWatcher, self).__init__()
        self.filename = filename
        self.index_cnt = re.sub("status.txt", "indexed.cnt", self.filename)
        self.run_info = re.sub("status.txt", "run.info", self.filename)
        self.window = window
        self.running = True
        self.cv = threading.Condition()
        self.runid = runid
        self.comment = None
        try:
            self.start()
        except: # To avoid "too many theads" error
            time.sleep(1) # wait "Finished" threads to be stopped
            self.start

    def stop(self):
        self.running = False
        self.cv.acquire()
        self.cv.notify()
        self.cv.release()

    def get_log(self):
        ret = None
        try:
            with open(self.filename) as f:
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

    def run(self):
        while self.running:
            self.cv.acquire()
            self.cv.wait(1.5)
            self.cv.release()

            if self.comment is None:
                self.parse_info()
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
    COL_RUNID = 0
    COL_STATUS = 1
    COL_TOTAL = 2
    COL_PROCESSED = 3
    COL_LLF_PASSED = 4
    COL_HITS = 5
    COL_INDEXED = 6
    COL_COMMENT = 7

    MENU_KILLJOB = 0
    MENU_HDFSEE = 1
    MENU_CELLEXPLORER = 2
    MENU_COUNTSUMS = 3

    def __init__(self, parent, opts, nogui=False):
        self.opts = opts
        self.subthreads = []
        self.rows = {}
        self.waitFor = None
        self.autoqsub = None
        self.noGUI = nogui
        if nogui:
            return

        title = "Cheetah Dispatcher on " + os.getcwd()
        wx.Frame.__init__(self, parent, title=title, size=(900,800))
        self.table = wx.grid.Grid(self, size=(900, -1))
        self.table.CreateGrid(0, 8, wx.grid.Grid.SelectRows) # row, column
        self.table.EnableEditing(False)
        self.table.SetRowLabelSize(0)
        # FIXME: better col width
        self.table.SetColLabelValue(MainWindow.COL_RUNID, "Run ID         ")
        self.table.SetColLabelValue(MainWindow.COL_STATUS, "Status       ")
        self.table.SetColLabelValue(MainWindow.COL_TOTAL, "Total    ")
        self.table.SetColLabelValue(MainWindow.COL_PROCESSED, "Processed     ")
        self.table.SetColLabelValue(MainWindow.COL_LLF_PASSED, "Accepted      ")
        self.table.SetColLabelValue(MainWindow.COL_HITS, "Hits          ")
        self.table.SetColLabelValue(MainWindow.COL_INDEXED, "Indexed       ")
        self.table.SetColLabelValue(MainWindow.COL_COMMENT, "Comment  ")
        attr = wx.grid.GridCellAttr()
        attr.SetRenderer(ProgressCellRenderer())
        self.table.SetColAttr(MainWindow.COL_PROCESSED, attr)
        self.table.SetColAttr(MainWindow.COL_LLF_PASSED, attr)
        self.table.SetColAttr(MainWindow.COL_HITS, attr)
        self.table.SetColAttr(MainWindow.COL_INDEXED, attr)
        for i in xrange(7):
            self.table.AutoSizeColLabelSize(i)
        self.table.Bind(wx.grid.EVT_GRID_CELL_RIGHT_CLICK, self.OnGridRightClick)
        self.table.Bind(wx.EVT_SIZE, self.OnGridResize)

        self.start_button = wx.Button(self, wx.ID_ANY, "  Submit job(s)  ")
        self.start_button.Bind(wx.EVT_BUTTON, self.onPush)
        self.text_runid = wx.TextCtrl(self)
        self.label_runid = wx.StaticText(self, wx.ID_ANY, "Run ID:")
	self.combo_bl = wx.ComboBox(self, wx.ID_ANY, value="BL2", choices=["BL2", "BL3"],
                                            size=(80, -1), style=wx.CB_READONLY)
        self.combo_bl.SetSelection(self.opts.bl - 2) # TODO: refactor 

        self.hsizer = wx.BoxSizer(wx.HORIZONTAL)
        self.hsizer.Add(self.label_runid, 0, wx.ALIGN_CENTER_VERTICAL)
        self.hsizer.Add(self.text_runid, 1, wx.EXPAND | wx.ALL, 5)
        self.hsizer.Add(self.combo_bl, 0, wx.ALIGN_CENTER_VERTICAL)
        self.hsizer.Add(self.start_button, 0, wx.ALL, 3)

        self.label_pd = wx.StaticText(self, wx.ID_ANY, "Dark/Light threshold (0 to disable)")
        self.text_pd1 = wx.TextCtrl(self, wx.ID_ANY, "%f" % opts.pd1_thresh)
        self.label_pd1 = wx.StaticText(self, wx.ID_ANY, "%s :" % self.opts.pd1_name)
        self.text_pd2 = wx.TextCtrl(self, wx.ID_ANY, "%f" % opts.pd2_thresh)
        self.label_pd2 = wx.StaticText(self, wx.ID_ANY, "%s:" % self.opts.pd2_name)
        self.text_pd3 = wx.TextCtrl(self, wx.ID_ANY, "%f" % opts.pd3_thresh)
        self.label_pd3 = wx.StaticText(self, wx.ID_ANY, "%s:" % self.opts.pd3_name)

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
        if self.opts.pd1_name is None and self.opts.pd2_name is None and self.opts.pd3_name is None:
            self.vsizer.Hide(self.hsizer_pd)
        if self.opts.pd1_name is None:
            self.vsizer.Hide(self.hsizer_pd1)
        if self.opts.pd2_name is None:
            self.vsizer.Hide(self.hsizer_pd2)
        if self.opts.pd3_name is None:
            self.vsizer.Hide(self.hsizer_pd3)
        self.Layout()
        
        self.timer = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, self.OnTimer)
        self.timer.Start(5000)

        self.scanDirectory()

        if self.opts.quick == 1: 
            self.startAutoQsub()

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
        runid = self.table.GetCellValue(row, self.COL_RUNID)
        
        point = event.GetPosition()
        popupmenu = wx.Menu()
        killjob = popupmenu.Append(MainWindow.MENU_KILLJOB, "Kill this job")
        hdfsee = popupmenu.Append(MainWindow.MENU_HDFSEE, "View hits")
        hdfsee = popupmenu.Append(MainWindow.MENU_CELLEXPLORER, "Check cell")
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
        elif id == MainWindow.MENU_COUNTSUMS:
#            GetSelectionBlock no longer works on newer wxPython as before...
#             https://forums.wxwidgets.org/viewtopic.php?t=41607
#            row1 = self.table.GetSelectionBlockTopLeft()[0][0]
#            row2 = self.table.GetSelectionBlockBottomRight()[0][0]
            self.CountSums(self.table.GetSelectedRows())

    def CountSums(self, rows):
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
                accepted[typ] += int(self.table.GetCellValue(row, self.COL_LLF_PASSED).rsplit(" ")[0])
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
        runnum = re.sub("-light|-0", "", runid)
        def launchHDFsee():
            if os.path.exists("%s/%s.stream" % (runid, runid)):
                command = "cd {runid}; hdfsee -g {runnum}.geom {runid}.stream -c invmono -i 10 &".format(runid=runid, runnum=runnum)
            else:
                command = "hdfsee -g {runid}/{runnum}.geom {runid}/run{runid}.h5 -c invmono -i 10 &".format(runid=runid, runnum=runnum)
            os.system(command)

        threading.Thread(target=launchHDFsee).start()

    def CellExplorer(self, runid):
        # FIXME: warn when multiple rows are selected
        def launchCellExplorer():
            command = "cd {runid}; cell_explorer {runid}.stream &".format(runid=runid)
            os.system(command)
            
        runnum = re.sub("-light|-0", "", runid)
        if os.path.exists("%s/%s.stream" % (runid, runid)):
            t = threading.Thread(target=launchCellExplorer).start()
        else:
            self.showError("Indexing result for {runid} is not available (yet).".format(runid=runid))

    def KillJob(self, runid):
        message = "Are you sure to kill job %s?" % runid
        dlg = wx.MessageDialog(None, message, "Cheetah dispatcher", wx.YES_NO | wx.NO_DEFAULT)
        ret = dlg.ShowModal()
        dlg.Destroy()
        if ret == wx.ID_NO:
            return
        try:
            jobid = open("%s/job.id" % runid).read()
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
        for filename in sorted(glob.glob("*")):
            if os.path.isdir(filename) and re_filename.match(filename) != None:
                self.addRun(filename)
        # TODO: Support deleting existing rows. This is more complex since the row idx changes.

    def onPush(self, event):
        if self.waitFor != None:
            self.stopWatch()
        else:
            self.prepareSubmitJob(self.text_runid.GetValue())

    def prepareSubmitJob(self, run_str):
        runids = None
        try:
            bl = self.combo_bl.GetSelection() + 2 # TODO: refactor
            pd1_thresh = float(self.text_pd1.GetValue())
            pd2_thresh = float(self.text_pd2.GetValue())
            pd3_thresh = float(self.text_pd3.GetValue())
            re_single = re.match('^([0-9]+)$', run_str)
            re_range = re.match('^([0-9]+)-([0-9]+)$', run_str)
            re_autofollow = re.match('^([0-9]+)-$', run_str)
            if pd1_thresh < 0 or pd2_thresh < 0 or pd3_thresh < 0:
                raise

            if re_single != None: # I wish I could use = within if !!
                runids = [int(re_single.group(1))]
            elif re_autofollow != None:
                self.startWatchFrom(int(re_autofollow.group(1)))
                return
            elif re_range != None:
                run_start = int(re_range.group(1))
                run_end = int(re_range.group(2))
                if (run_end - run_start > 100):
                    self.showError("Too many runs to submit")
                    return
                runids = range(run_start, run_end + 1)
            else:
                raise
        except:
            print traceback.format_exc()
            self.showError("Invalid run ID or paramters were specified.")
            return

        for runid in runids:
            self.startRun("%d" % runid, bl, pd1_thresh, pd2_thresh, pd3_thresh)

    def stopWatch(self):
        self.waitFor = None
        self.start_button.SetLabel("  Submit job(s)  ")
        self.text_runid.Enable()
        self.text_pd1.Enable()
        self.text_pd2.Enable()
        self.text_pd3.Enable()
        self.combo_bl.Enable()

    def startWatchFrom(self, runid):
        self.waitFor = runid
        self.start_button.SetLabel("Stop automode")
        self.text_runid.Disable()
        self.text_pd1.Disable()
        self.text_pd2.Disable()
        self.text_pd3.Disable()
        self.combo_bl.Disable()

    def OnTimer(self, event):
        self.scanDirectory()
        if (self.waitFor == None):
            return

        out = Popen(["ShowRunInfo", "-b", "%d" % self.opts.bl, "-r", "%d" % self.waitFor], stdout=PIPE).communicate()[0]
        lines = out.split("\n")
        if lines[0].find("Ready to Read") != -1:
                print "\rRun %d became ready." % self.waitFor
                self.prepareSubmitJob("%d" % self.waitFor)
                self.waitFor = self.waitFor + 1
                self.text_runid.SetValue("%d-" % self.waitFor)

    def showError(self, message):
        dlg = wx.MessageDialog(None, message, "Cheetah dispatcher", wx.OK)
        dlg.ShowModal()
        dlg.Destroy()

    def loadCrystFELArgs(self):
        ret = " "
        try:
            f = open("crystfel.args", "r")
            for line in f:
                line = line.strip()
                pos = line.find("#")
                if pos >= 0:
                    line = line[:pos]
                if len(line) > 0:
                    ret = ret + line + " "
        except:
            return ""
        return ret

    def startRun(self, runid, bl, pd1_thresh=0, pd2_thresh=0, pd3_thresh=0):
        run_dir = runid
        arguments = ""
        master_arguments = ""
        subjobs = []
        
        if opts.detector == "rayonix":
            if opts.clen and opts.clen==opts.clen: arguments += " --clen=%f " % opts.clen
            if opts.beam_x and opts.beam_x==opts.beam_x: arguments += " --beam-x=%f " % opts.beam_x
            if opts.beam_y and opts.beam_y==opts.beam_y: arguments += " --beam-y=%f " % opts.beam_y

        if (pd1_thresh != 0 or pd2_thresh != 0 or pd3_thresh != 0):
            run_dir += "-light"
            if (pd1_thresh != 0 and self.opts.pd1_name is not None):
                arguments += " --pd1_thresh=%.3f --pd1_name=%s " % (pd1_thresh, self.opts.pd1_name)
            if (pd2_thresh != 0 and self.opts.pd2_name is not None):
                arguments += " --pd2_thresh=%.3f --pd2_name=%s " % (pd2_thresh, self.opts.pd2_name)
            if (pd3_thresh != 0 and self.opts.pd3_name is not None):
                arguments += " --pd3_thresh=%.3f --pd3_name=%s " % (pd3_thresh, self.opts.pd3_name)
            master_arguments = arguments + " --type=light "
            if self.opts.submit_dark_any == 1:
                subjobs.append("dark")
            else:
                for i in xrange(1, self.opts.submit_dark_to + 1):
                    subjobs.append("dark%d" % i)

        else:
            master_arguments = arguments + " --type=0 "
            run_dir += "-0"
            for i in xrange(1, PARALLEL_SIZE):
                subjobs.append("%d" % i)

        crystfel_args = self.opts.crystfel_args + " " + self.loadCrystFELArgs()

        # Master
        if os.path.exists(run_dir):
            self.showError("You told me to process run %s, but a directory for the run already exists. Please remove it before re-processing." % run_dir)
            return
        os.mkdir(run_dir)
        f = open("%s/run.sh" % run_dir, "w")
        job_scr = job_script_hs if opts.detector == "rayonix" else job_script
        f.write(job_scr.format(runid=runid, runname=run_dir, clen=self.opts.clen, queuename=self.opts.queue,
                               subjobs=" ".join(subjobs), arguments=master_arguments,
                               crystfel_args=crystfel_args, beamline=bl))
        f.close()
        os.system("qsub {rundir}/run.sh > {rundir}/job.id".format(rundir=run_dir))
        self.addRun(run_dir)

        # Children
        for subjob in subjobs:
            run_dir = runid + "-" + subjob
            child_arguments = arguments + " --type=" + subjob
            
            if os.path.exists(run_dir):
                self.showError("You told me to process run %s, but a directory for the run already exists. Please remove it before re-processing." % run_dir)
                return
            os.mkdir(run_dir)
            f = open("%s/run.sh" % run_dir, "w")
            job_scr = job_script_dark_hs if opts.detector == "rayonix" else job_script_dark
            f.write(job_scr.format(runid=runid, runname=run_dir, 
                                   queuename=self.opts.queue, arguments=child_arguments,
                                   crystfel_args=crystfel_args, beamline=bl))
            f.close()
            if self.opts.quick != 1:
                os.system("qsub {runid}/run.sh > {runid}/job.id".format(runid=run_dir))
            self.addRun(run_dir)

    def addRun(self, runid):
        if self.noGUI or self.rows.get(runid) is not None:
            return

        row = self.table.GetNumberRows()
        self.subthreads.append(LogWatcher("%s/status.txt" % runid, self, runid))
        self.rows[runid] = row
        self.table.AppendRows(1)
        self.table.SetCellValue(row, MainWindow.COL_RUNID, runid)
        self.table.SetCellValue(row, MainWindow.COL_STATUS, "waiting")

    def onClose(self, event):
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
        row = self.rows.get(event.msg['runid'])
        if row == None:
            return

        comment = event.msg['Comment']
        if comment is not None:
            self.table.SetCellValue(row, MainWindow.COL_COMMENT, comment)

        status = event.msg['Status']
        self.table.SetCellValue(row, MainWindow.COL_STATUS, status)
        if (status == "Finished"):
            self.table.SetCellBackgroundColour(row, MainWindow.COL_STATUS, (124, 252, 0))
        elif (status.startswith("Error")):
            self.table.SetCellBackgroundColour(row, MainWindow.COL_STATUS, (252, 124, 0))
            return 

        total = int(event.msg['Total'])
        processed = int(event.msg['Processed'])
        self.table.SetCellValue(row, MainWindow.COL_TOTAL, 
                                "%d" % total)
        self.table.SetCellValue(row, MainWindow.COL_PROCESSED, 
                                "%d (%.1f%%)" % (processed, 100.0 * processed / total))

        if (status == "DarkAveraging" or status == "waiting"):
            return
        LLFpassed = int(event.msg['LLFpassed'])
        hits = int(event.msg['Hits'])
        acceptance_rate = 0
        if processed != 0: acceptance_rate = 100.0 * LLFpassed / processed
        hit_rate = 0
        if LLFpassed != 0: hit_rate = 100.0 * hits / LLFpassed
        self.table.SetCellValue(row, MainWindow.COL_LLF_PASSED,
                                "%d (%.1f%%)" % (LLFpassed, acceptance_rate))
        self.table.SetCellValue(row, MainWindow.COL_HITS,
                                "%d (%.1f%%)" % (hits, hit_rate))

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
        #dc.SetPen(wx.BLACK_PEN)
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
print "Cheetah dispatcher GUI version 20171123"
print "   by Takanori Nakane (takanori.nakane@bs.s.u-tokyo.ac.jp)"
print
print "Please cite the following paper when you use this software."
print " \"Data processing pipeline for serial femtosecond crystallography at SACLA\""
print " Nakane et al., J. Appl. Cryst. (2016). 49"
print

parser = optparse.OptionParser()
parser.add_option("--monitor", dest="monitor", type=int, default=False, help="Monitor only")
parser.add_option("--bl", dest="bl", type=int, default=2, help="Beamline")
parser.add_option("--clen", dest="clen", type=float, default=50.0, help="camera length in mm")
parser.add_option("--beam-x", dest="beam_x", type=float, help="beam center x in px")
parser.add_option("--beam-y", dest="beam_y", type=float, help="beam center y in px")
parser.add_option("--quick", dest="quick", type=int, default=False, help="enable quick mode")
parser.add_option("--queue", dest="queue", type=str, default="serial", help="queue name")
parser.add_option("--max_jobs", dest="max_jobs", type=int, default=14, help="maximum number of jobs to submit when --quick is enabled")
parser.add_option("--pd1_name", dest="pd1_name", type=str, default=None, help="PD1 sensor name e.g. xfel_bl_3_st_4_pd_laser_fitting_peak/voltage")
# e.g. xfel_bl_3_st_4_pd_user_10_fitting_peak/voltage
parser.add_option("--pd2_name", dest="pd2_name", type=str, default=None, help="PD2 sensor name")
parser.add_option("--pd3_name", dest="pd3_name", type=str, default=None, help="PD3 sensor name")
parser.add_option("--submit_dark2", dest="submit_dark2", type=int, default=False, help="(DEPRECATED) accepts second darks (Ln-D2) and divide into light, dark1 and dark2")
parser.add_option("--submit_dark_to", dest="submit_dark_to", type=int, default=False, help="accepts up to M (<=9) darks (Ln-Dm) and divide into light, dark1, ..., darkM")
parser.add_option("--submit_dark_any", dest="submit_dark_any", type=int, default=False, help="accepts any lights and darks (Ln-Dm) and divide into light and dark")
parser.add_option("--crystfel_args", dest="crystfel_args", type=str, default="", help="optional arguments to CrystFEL")
parser.add_option("--submit", dest="submit", type=int, default=-1, help="submit run (noGUI)")
parser.add_option("--pd1_thresh", dest="pd1_thresh", type=float, default=0, help="PD1 threshold")
parser.add_option("--pd2_thresh", dest="pd2_thresh", type=float, default=0, help="PD2 threshold")
parser.add_option("--pd3_thresh", dest="pd3_thresh", type=float, default=0, help="PD3 threshold")
parser.add_option("--detector", dest="detector", type=str, default="mpccd", help="mpccd or rayonix")

opts, args = parser.parse_args()

if opts.detector not in ("mpccd", "rayonix"):
    sys.stderr.write("ERROR: Invalid --detector parameter (%s was given)\n" % opts.detector)
    sys.exit(-1)

if opts.detector=="mpccd" and not os.path.exists("sacla-photon.ini"):
    sys.stderr.write("ERROR: Configuration file was not found!\n\n")
    sys.stderr.write("You should copy @@TEMPLATE_FILE@@ into this directory\n")
    sys.stderr.write("and confirm the settings.\n")
    sys.exit(-1)

if opts.submit_dark2 is not False:
    sys.stderr.write("WARNING: submit_dark2 has been deprecated. Use --submit-dark-to=2 instead.\n")
    sys.stderr.write("         --submit-dark-to was set to 2 for you.\n")
    opts.submit_dark_to = 2

if opts.submit_dark_to is not False and opts.submit_dark_any == 1:
    sys.stderr.write("ERROR: You cannot enable submit_dark_any and set submit_dark_to simultaneously!\n")
    sys.exit(-1)

if opts.submit_dark_to is False:
    opts.submit_dark_to = 1

if opts.submit_dark_to > 9 or opts.submit_dark_to < 1:
    sys.stderr.write("ERROR: submit_dark_to must be within 1 to 9.\n")
    sys.exit(-1)

if opts.monitor is not False and opts.quick is not False:
    sys.stderr.write("ERROR: You cannot enable 'quick' when 'monitor' mode is enabled.\n")
    sys.exit(-1)

if opts.max_jobs < 1:
    sys.stderr.write("ERROR: max_jobs must be a positive integer.\n")
    sys.exit(-1)

if opts.bl != 2 and opts.bl != 3:
    sys.stderr.write("ERROR: beamline must be 2 or 3.\n")
    sys.exit(-1)

print "Option: detector         = %s" % opts.detector
print "Option: monitor          = %s" % opts.monitor
print "Option: bl               = %d" % opts.bl
print "Option: clen             = %f mm" % opts.clen
print "Option: quick            = %s" % opts.quick
print "Option: max_jobs         = %s" % opts.max_jobs
print "Option: queue            = %s" % opts.queue
print "Option: pd1_name         = %s" % opts.pd1_name
print "Option: pd2_name         = %s" % opts.pd2_name
print "Option: pd3_name         = %s" % opts.pd3_name
print "Option: pd1_thresh       = %f" % opts.pd1_thresh
print "Option: pd2_thresh       = %f" % opts.pd2_thresh
print "Option: pd3_thresh       = %f" % opts.pd3_thresh
print "Option: submit_dark_to   = %s" % opts.submit_dark_to
print "Option: submit_dark_any  = %s" % opts.submit_dark_any
print "Option: crystfel_args    = %s" % opts.crystfel_args

if opts.submit > 0: # headless
    mw = MainWindow(None, opts, True)
    print mw.loadCrystFELArgs()
    mw.startRun("%d" % opts.submit, opts.bl, pd1_thresh=opts.pd1_thresh, pd2_thresh=opts.pd2_thresh, pd3_thresh=opts.pd3_thresh)
else:
    app = wx.App(False)
    frame = MainWindow(None, opts)
    app.MainLoop()

