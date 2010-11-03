"""
The main module for the multiuserdsm framework
"""

from optparse import OptionParser
from bundle import Bundle


parser = OptionParser()

parser.add_option("-n","--network", dest="network", help="read network configuration from FILE",metavar="NETWORK")
parser.add_option("-K","--tones", dest="K", help="specify number of DMTCHANNELS",metavar="K")

(options,args) = parser.parse_args()

if __name__ == "__main__":
	bundle = Bundle(options.network,options.K)
