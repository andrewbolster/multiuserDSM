"""
The main module for the multiuserdsm framework
"""

from optparse import OptionParser
from bundle import Bundle
from utility import *
from osb import OSB
import cProfile

log.info("Starting Up...")
profiling=True
parser = OptionParser()

parser.add_option("-n","--network", dest="network", help="read network configuration from FILE",metavar="NETWORK",default="test.net")
parser.add_option("-K","--tones", dest="K", help="specify number of DMTCHANNELS",metavar="K",default=224)
parser.add_option("-S","--scenario", dest="scenarioname", help="specify scenario name (for file IO)",metavar="SCENARIONAME",default="OSB-3k_5k-near_far")

(options,args) = parser.parse_args()
bundle = Bundle(options.network,options.K,options.scenarioname)

if __name__ == "__main__":
    """
    Perform algorithm selection and option passing
    """
    algo = OSB(bundle)
    if profiling:
        cProfile.run('algo.run()','profiles/'+options.scenarioname)
    else:
        algo.run()
    algo.tofile(options.scenarioname)