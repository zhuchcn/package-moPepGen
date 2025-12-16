""" moPepGen util module """
# Import all util modules to make them accessible as util.module_name
# This is required for __main__.py which accesses them as:
#   util.brute_force.parse_args(subparsers)
#   util.fuzz_test.parse_args(subparsers)
#   etc.
from . import brute_force
from . import brute_force_novel_orf
from . import downsample_reference
from . import validate_variant_calling
from . import fuzz_test
from . import extract_gvf
from . import validate_novel_orf_calling
from . import add_fuzz_test_log
from . import common
from . import ResourcesMonitor

PROG_NAME = 'moPepGen-util'
