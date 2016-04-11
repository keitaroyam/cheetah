from collections import OrderedDict

master_params_str = """\
work_dir = None
  .type = path
  .help = working directory

distl {
  res {
    outer = 5.
      .type=float
      .help="High resolution limit in angstroms"
    inner = 30.
      .type=float
      .help="Low resolution limit in angstroms"
  }
}
cheetah {
 ADCthresh = 5
  .type = float
 MinSNR = 8
  .type = float
 MinPixCount = 2
  .type = int
 MaxPixCount = 40
  .type = int
 LocalBGRadius = 2
  .type = float
 MinPeakSeparation = 0
  .type = float
}
"""

sp_params_strs =  OrderedDict(((("BL32XU", "EIGER9M", None, None), """\
"""),
                               ))
