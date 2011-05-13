'''
Algorithm Parent module 
'''

import sys
from utility import *
import math
import numpy as np
import pylab as pl
import time
import os

from bundle import Bundle
from line import Line



class Algorithm(object):
    '''
    Algorithm Constants; usually set within child classes
    '''
    MAXRATE = None #Need to set this for load_fm and ISB
    MAXPOWER = None
    MAXBITSPERTONE = 16 #max 15 + 0
    name="Default Algorithm Name; You should never see this!"
    stats={}


    def __init__(self,bundle,useGPU=False, rate_search=True):
        assert isinstance(bundle,Bundle), "Did you give me a bundle to use?"
        self.bundle=bundle
        self.useGPU=useGPU
        self.rate_search=rate_search
        self.stats['ratesearch']=rate_search
        self.stats['gpu']=useGPU
        #Default values
        self.defaults={'maxval':sys.maxint*1.0,
                'l':1.0,             #From osb_bb.c::bisect_l() (modified from 0)
                'l_min':0.0,
                'l_max':1.0,
                'w':1.0, #same behaviour for any real value
                'w_max':100,
                'w_min':0.001,
                'p_budget':0.110,    #watts #from scenarios.c (0.110
                'rate_target':False,
                'min_step':500,      #was min_sl
                'p_tol':0.015,      #should be 0.015
                'rate_tol':10,
                'GAMMA':self.bundle.GAMMA,
                'spectral_mask':dbmhz_to_watts(50)
                }



    def preamble(self):
        '''
        Stuff that needs to be done for any algorithm, including 
        *initialisation of power and bitloading matrices
        *timing information
        '''
        #Power and bitloading are dim(KxN)
        self.p=np.zeros((self.bundle.K,self.bundle.N)) #per tone per user power in watts
        self.b=np.asmatrix(np.zeros((self.bundle.K,self.bundle.N)))
        self.rate_targets = np.asarray([line.rate_target for line in self.bundle.lines]) 
        log.info("Starting Spectrum Balancing using %s"%self.name)
        self.stats['algo']=self.name
        self.stats['start']=time.time()

    def postscript(self):
        '''
        Stuff that should be done for any algorithm, including
        *downtyping of b -> ndarray
        *reassignment of relevent p's and b's to line values
        *SNR et al calculations
        *Timing statistics
        '''
        #Recalculate P from Bitloads to be careful
        self.b=mat2arr(self.b)

        if self.useGPU:
            self.p=self.bundle.recalcpsd(self.b)
        if False:
            for k in range(self.bundle.K):
                self.p[k]=self.bundle.calc_psd(self.b[k],k,gpu=False)#disable manual psd for the time being


        self.bundle.calculate_snr()
        self.stats['end']=time.time()
        self.stats['duration']=self.stats['end']-self.stats['start']
        try:
            self.stats['hitratio']=self.bundle._psd_cache['hits']/(self.bundle._psd_cache['hits']+self.bundle._psd_cache['misses'])
        except ZeroDivisionError:
            self.stats['hitratio']=0.0

        self.stats['linerates']=np.asarray([line.b_total for line in self.bundle.lines])
        self.stats['lines']=np.asarray([(line.nt,line.lt) for line in self.bundle.lines])
        self.stats['linetargets']=np.asarray([line.rate_target for line in self.bundle.lines])
        self.stats['linepowers']=np.asarray([line.p_total for line in self.bundle.lines])
        self.stats['bundlerate']=mat2arr(sum([line.rate() for line in self.bundle.lines]))

        log.info("All Done Here, took %f, hit ratio of %.4f"%(self.stats['duration'],self.stats['hitratio']))


    def rate_bisect(self,bisect):
        '''
    Common Rate bisection algorithm. Operates using lambda bisection function passed to it
    '''
        log.info("Rate Region Search: Warning: This doesn't work!")
        while not self.bundle.rates_converged(self.defaults['rate_tol']):
            for lineid,line in enumerate(self.bundle.lines):
                if self.rate_targets[lineid] == self.defaults['rate_target']:
                    continue
                bisect()
                log.info("Line:%d,Current:%d,Target:%d,Weight:%.3f "%(lineid,line.rate(),line.rate_target,self.w[lineid]))
                if line.rate_converged(self.defaults['rate_tol']):
                    log.info("Converged very early on line:%s"%lineid)
                    continue
                if line.rate() > line.rate_target:
                    log.info("Line:%d,Current:%d,Target:%d,Weight:%.3f Overshot"%(lineid,line.rate(),line.rate_target,self.w[lineid]))
                    self.w_max[lineid]=self.w[lineid]
                    self.w[lineid]/=2
                    while True:
                        bisect()
                        if line.rate()>line.rate_target:
                            log.info("Line:%d,Current:%d,Target:%d,Weight:%.3f Overshot"%(lineid,line.rate(),line.rate_target,self.w[lineid]))
                            self.w_max[lineid]=self.w[lineid]
                            self.w[lineid]/=2
                        else:
                            log.info("Line:%d,Current:%d,Target:%d,Weight:%.3f Settled"%(lineid,line.rate(),line.rate_target,self.w[lineid]))
                            self.w_min[lineid]=self.w[lineid]
                            break
                else: #Rate less than target
                    log.info("Line:%d,Current:%d,Target:%d,Weight:%.3f Undershot"%(lineid,line.rate(),line.rate_target,self.w[lineid]))
                    self.w_min[lineid]=self.w[lineid]
                    self.w[lineid]*=2
                    while True:
                        bisect()
                        if line.rate()<line.rate_target:
                            log.info("Line:%d,Current:%d,Target:%d,Weight:%.3f Undershot"%(lineid,line.rate(),line.rate_target,self.w[lineid]))
                            self.w_min[lineid]=self.w[lineid]
                            self.w[lineid]*=2
                        else:
                            log.info("Line:%d,Current:%d,Target:%d,Weight:%.3f Settled"%(lineid,line.rate(),line.rate_target,self.w[lineid]))
                            self.w_max[lineid]=self.w[lineid]
                            break
                if  line.rate_converged(self.defaults['rate_tol']):
                    log.info("Line:%d,Current:%d,Target:%d,Weight:%.3f Max/Min Converged[%.3f,%.3f]"%(lineid,line.rate(),line.rate_target,self.w[lineid], self.w_max[lineid],self.w_min[lineid]))
                    continue
                else:
                    self.w[lineid]=(self.w_max[lineid]+self.w_min[lineid])/2
                    log.info("Line:%d,Current:%d,Target:%d,Weight:%.3f Full Bisection"%(lineid,line.rate(),line.rate_target,self.w[lineid]))
                    while True:
                        bisect()
                        if line.rate>line.rate_target:
                            self.w_max[lineid]=self.w[lineid]
                        else:
                            self.w_min[lineid]=self.w[lineid]

                        if line.rate_converged(self.defaults['rate_tol']):
                            log.info("Line:%d,Current:%d,Target:%d,Weight:%.3f Full Converged"%(lineid,line.rate(),line.rate_target,self.w[lineid]))
                            break
                        else:
                            self.w[lineid]=(self.w_max[lineid]+self.w_min[lineid])/2
                            log.info("Line:%d,Current:%d,Target:%d,Weight:%.3f New Weight"%(lineid,line.rate(),line.rate_target,self.w[lineid]))
        log.info("Rates converged")

    def _l_converged(self,line=False,last=False):
        '''
    l converged
    Decides whether the line/bundle is done
    '''
        if isinstance(line,Line): #if called with a line
            if last == False: 
                return False #Force the first loop through
            else:
                thispower = self.total_power(line)
                howfar = abs(self.power_budget[line.id]-thispower)
                return howfar < self.defaults['p_tol'] or thispower==last
        else: #if called without a line, assume operation on the bundle
            if (self.l_last == self.l).all():
                #Optimisation done since all values the same as last time
                assert (self.l > 0).all()
                return True
            else:
                return False

    def _l_k(self,bitload,lambdas,weights,k):
        '''
    L_k; Return the Lagrangian given a bit-loading combination and tone
    Effectively the cost function
    '''
        P=self.bundle.calc_psd(bitload,k)
        #If anything's broken, this run is screwed anyway so feed optimise_p a bogus value
        #THIS CHECK HAS BEEN GPU CHECKED 19/4
        if (P < 0).any() or (P>self.defaults['spectral_mask']).any():
            return -self.defaults['maxval']

        # p is a matrix of power values required on each line[k] to support bitload[line] bits
        # l is the current lambda array
        # bitload is the current bitload array
        # w is the omega-weight array

        #After profiling, np.reduce is faster than looping
        bw=np.add.reduce(bitload*weights)
        lp=np.add.reduce(lambdas*P)

        lk=bw-lp

        return lk

    def total_power(self,line=False):
        '''
    Total Power: return a line's planned total power
    '''
        if isinstance(line,Line):
            #swap the tone and line dimensions for cleanliness
            power = np.add.reduce(self.p)[line.id]
        else: 
            power = np.add.reduce(np.add.reduce(self.p))[0,0]
        assert(isinstance(power,np.float64))
        return power

    def update_b_p(self):
        '''
    Update bundle b/p values with respect to algorithm progression
    '''   
        for line in self.bundle.lines:
            line.b=self.b[:,line.id]
            line.p=self.p[:,line.id]

            line.p_total = sum(map(dbmhz_to_watts,line.p))
            line.b_total = sum(line.b)
            #log.info("Line:%d,Power:%fW,Rate:%dbpf"%(line.id,line.p_total,line.b_total))

    def eta(self,how_done,start=False):
        '''
        Takes a 'work done' value from the child algo in the range 0-1 and log.info's the estimated time left
        Can also work out eta's for algorithm sub-processes when giving a non-global start time 
        '''
        if start == False:
            #Not given any start time, work of algo time
            start = self.stats['start']
        time_done=time.time()-start
        try:
            eta=time_done/how_done
            log.info("ETA:%s,%f"%(str(eta),how_done))
        except ZeroDivisionError:
            pass

    def _calc_delta_p(self,line,tone_full):
        '''
    Calculate cost of additional bit on all the tones in a line
    :from am_load.c (implicit)
    Optimised from original version (Thank you wolframalpha)
    '''
        delta_p=np.zeros(self.bundle.K)
        for tone in self.bundle.K:
            if not tone_full[tone]:
                delta_p[tone] = (pow(2,(line.b[tone]-1)) * 3 - 2 )* self.bundle.gamma_hat/line.cnr[tone]
            else:
                delta_p[tone] = 100
        return delta_p

    def update_psds(self,line):
        '''
    Update PSD's for this line
    :from am_load.c (implicit)
    Potentially integratable into calculatesnr?
    '''
        for tone in self.bundle.K:
            line.p[tone] = watts_to_dbmhz((pow(2,line.b[tone])-1)*(self.gamma_hat/line.cnr[tone]))
            if (line.p[tone] < line.MINPSD) and (line.b[tone] != 0):
                log.debug("Changing p from %e to MINPSD"%line.p[tone])
                line.p[tone] = line.MINPSD

    def tofile(self,filename):
        '''
    Print Statistics, Power and Bitrates to file after generation
    '''
        assert os.path.isdir(rawdir)==True, "No "
        np.save(rawdir+filename+'-power', self.p)
        np.save(rawdir+filename+'-stats', self.stats)
        np.save(rawdir+filename+'-bitrate',self.b)

    def test_compliance(self,alt_scenario):
        '''
    Test Compliance checks the calculated bit and power values against a previous scenario.
    Will not work for non-determinitic algorithms (duh)
    '''
        log.info("Testing against %s"%alt_scenario)
        assert os.path.isdir(rawdir), "No Path, Dunno wut happened there...%s"%graphdir
        alt_p=np.load(rawdir+alt_scenario+"-power.npy")
        alt_b=np.load(rawdir+alt_scenario+"-bitrate.npy")
        alt_cm=np.load(rawdir+alt_scenario+"-channelmatrix.npy")

        if not np.equal(alt_cm,self.bundle.xtalk_gain).all():
            log.critical("Channel Matrix is broken")
        if not np.equal(alt_p,self.p).all():
            log.critical("P is broken")
        if not np.equal(alt_b,self.b).all():
            log.critical("B is broken")
