#!/usr/bin/env python
'''
The main module for the multiuserdsm framework
'''
import matplotlib
matplotlib.use("Agg")

from optparse import OptionParser

from graphgenerator import graphgen
import cProfile,os

parser = OptionParser()

parser.add_option("-n","--network", dest="network", help="read network configuration from FILE",metavar="NETWORK",default="")
parser.add_option("-a","--algo", dest="algo", help="specify algo to be used",metavar="ALGO",default="")
parser.add_option("-K","--tones", dest="K", help="specify number of DMTCHANNELS",metavar="K",default=224)
parser.add_option("-S","--scenario", dest="scenarioname", help="specify scenario name (for file IO)",metavar="SCENARIONAME",default="defaultscenario")
parser.add_option("-A","--altscenario", dest="altscenario", help="specify an existing scenario of check against", metavar="ALTSCENARIO")
parser.add_option("-G","--nographing",dest="graphing", action="store_false",help="disable graphing", default=True)
parser.add_option("-p","--profiling",dest="profiling", action="store_true",help="enable profiling", default=False)
parser.add_option("-c","--cache",dest="cache", action="store_true",help="attempt to read psd_cache from scenariofile", default=False)
parser.add_option("-g","--gpu",dest="gpu",action="store_true",help="attempt to use gpu", default=False)
parser.add_option("-r","--rate",dest="rate_search", action="store_true",help="enable rate searching", default=False)
parser.add_option("-l","--logfile", dest="logfile", help="specify logfile",metavar="LOGFILE",default="multiuserdsm.log")
parser.add_option("-N","--ngpu",dest="gpu",action="store", type="int", help="attempt to use N gpus", default=False)

from bundle import Bundle
from osb import OSB
from mipb import MIPB
from isb import ISB
import utility as util

(options,args) = parser.parse_args()

h = util.logging.FileHandler(options.logfile)
h.setFormatter(util.f)
util.log.addHandler(h)

if os.path.isfile(util.rawdir+options.scenarioname+'-lines.npy') and options.network=="":
    options.network=util.rawdir+options.scenarioname+'-lines.npy'
if options.cache and os.path.isfile(util.rawdir+options.scenarioname+'-cache.npy'):
    log.info("Using Cached PSD values. This is very dangerous")
    options.cache=util.rawdir+options.scenarioname+'-cache.npy'

util.log.info("Starting Up... logging to %s"%options.logfile)

bundle = Bundle(network_file=options.network,K=options.K,scenarioname=options.scenarioname, cachefile=options.cache, useGPU=options.gpu)
algos={"OSB":OSB,"MIPB":MIPB,"ISB":ISB}



if __name__ == "__main__":
    
    if options.algo=="":
        util.log.info("Generating scenario files %s and exiting"%options.scenarioname)
        sys.exit(0)
    
    '''
    Perform algorithm selection and option passing
    '''
    algo = algos[options.algo](bundle,options.gpu, options.rate_search)
    if options.profiling:
        log.info("Profiling Run")
        if not os.path.isdir(profdir):
            os.makedirs(profdir)
        cProfile.run('algo.run()',profdir+options.scenarioname)
    else:
        algo.run()
    algo.tofile(options.scenarioname)
    bundle.tofile(options.scenarioname)
    
    if options.altscenario:
        algo.test_compliance(options.altscenario)
    
    if options.graphing:
        g=graphgen(options.scenarioname)
        g.graph_all()
    
