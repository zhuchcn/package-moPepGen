""" moPepGen util module """

PROG_NAME = 'moPepGen-util'

# Import brute_force module to maintain backward compatibility
# This allows: from moPepGen import util; util.brute_force.main()
from . import brute_force
