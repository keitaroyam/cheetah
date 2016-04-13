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
 MinPixCount = 3
  .type = int
 MaxPixCount = 40
  .type = int
 LocalBGRadius = 2
  .type = int
 MinPeakSeparation = 0
  .type = float
 algorithm = 8
  .type = int
  .help = 6 or 8
}
"""

sp_params_strs =  OrderedDict(((("BL32XU", "EIGER9M", None, None), """\
"""),
                               ))
