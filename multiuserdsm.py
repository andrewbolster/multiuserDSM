'''
The main module for the multiuserdsm framework
'''
import matplotlib
matplotlib.use("Agg")

from optparse import OptionParser
from bundle import Bundle
from utility import *
from osb import OSB
from mipb import MIPB
from graphgenerator import graphgen
import cProfile,os

log.info("Starting Up...")
profiling=True
parser = OptionParser()

parser.add_option("-n","--network", dest="network", help="read network configuration from FILE",metavar="NETWORK",default="test.net")
parser.add_option("-K","--tones", dest="K", help="specify number of DMTCHANNELS",metavar="K",default=224)
parser.add_option("-S","--scenario", dest="scenarioname", help="specify scenario name (for file IO)",metavar="SCENARIONAME",default="defaultscenario")
parser.add_option("-A","--altscenario", dest="altscenario", help="specify an existing scenario of check against", metavar="ALTSCENARIO")
parser.add_option("-G","--nographing",dest="graphing", action="store_false",help="disable graphing", default=True)

(options,args) = parser.parse_args()
bundle = Bundle(options.network,options.K,options.scenarioname)

if __name__ == "__main__":
    '''
    Perform algorithm selection and option passing
    '''
    algo = MIPB(bundle)
    if profiling:
        if not os.path.isdir(profdir):
            os.makedirs(profdir)
        cProfile.run('algo.run()',profdir+options.scenarioname)
    else:
        algo.run()
    algo.tofile(options.scenarioname)
    
    if options.altscenario:
        algo.test_compliance(options.altscenario)
    
    if options.graphing:
        g=graphgen(options.scenarioname)
        g.graph_all()
    