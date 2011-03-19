"""
The main module for the multiuserdsm framework
"""

from optparse import OptionParser
from bundle import Bundle
from utility import *
from osb import OSB

log.info("Starting Up...")

parser = OptionParser()

parser.add_option("-n","--network", dest="network", help="read network configuration from FILE",metavar="NETWORK",default="test.net")
parser.add_option("-K","--tones", dest="K", help="specify number of DMTCHANNELS",metavar="K",default=224)

(options,args) = parser.parse_args()
bundle = Bundle(options.network,options.K)

if __name__ == "__main__":
    """
    Perform algorithm selection and option passing
    """
    algo = OSB(bundle)
    algo.run()