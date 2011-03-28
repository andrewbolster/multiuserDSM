'''
The main module for the multiuserdsm framework
'''

from optparse import OptionParser
from bundle import Bundle
from utility import *
from osb import OSB
from mipb import MIPB
import cProfile,os

log.info("Starting Up...")
profiling=True
parser = OptionParser()

parser.add_option("-n","--network", dest="network", help="read network configuration from FILE",metavar="NETWORK",default="test.net")
parser.add_option("-K","--tones", dest="K", help="specify number of DMTCHANNELS",metavar="K",default=224)
parser.add_option("-S","--scenario", dest="scenarioname", help="specify scenario name (for file IO)",metavar="SCENARIONAME",default="MIPB-3k_5k-near_far")

(options,args) = parser.parse_args()
bundle = Bundle(options.network,options.K,options.scenarioname)

if __name__ == "__main__":
    '''
    Perform algorithm selection and option passing
    '''
    algo = MIPB(bundle)
    if profiling:
        profiledir="profiles/"
        if not os.path.isdir(profiledir):
            os.makedirs(profiledir)
        cProfile.run('algo.run()',profiledir+options.scenarioname)
    else:
        algo.run()
    algo.tofile(options.scenarioname)