'''
Created on 22 Mar 2011

@author: bolster
'''
#Global Imports
import numpy as np
import pylab as pl
import sys
from itertools import product

#Local Imports
from bundle import Bundle
from algorithm import Algorithm
from line import Line
import utility as util

try:
    from gpu import GPU
    using_gpu=True
except:
    util.log.info("Not using GPU")
    using_gpu=False
    
if util.mp:
    #Thread Stuffs
    from multiprocessing import Process, Queue, cpu_count
    from Queue import Empty
 

class MIPB(Algorithm):
    '''
    Multiuser Incremental Power Balancing
    Implementation of greedy loading algorithm with adaptive power
    penalty function.
    constructor overloaded by Algorithm
    '''

    def run(self):
        '''
        Execute MIPB on the pre-initialised bundle
        #FIXME Not Done

        '''
        self.name = "MIPB"
        self.defaults = {'maxval':sys.maxint*1.0,
                         'bit_inc':1,
                         'rate_target':False,
                         'p_budget':0.110,    #watts #from mipb.c _p_budget[user]=
                         'w':1.0,
                         'w_min':0.0,
                         'min_sw':1,
                         'rate_tol':10,
                         'osc_tol':5
                         }
        
        self.spectral_mask=util.dbmhz_to_watts(-30) #from mipb.c _spectral_mask
        self.p_ave=0.0
        #Finished[tone,user]=False
        self.finished=np.tile(False,(self.bundle.K,self.bundle.N))
        self.delta_p=np.zeros((self.bundle.K,self.bundle.N,self.bundle.N))
        self.w=np.tile(self.defaults['w'],self.bundle.N)
        self.w_last=np.tile(self.defaults['w_min'],self.bundle.N)
        
        self.cost=np.zeros((self.bundle.K,self.bundle.N))

        self.power_budget=np.tile(self.defaults['p_budget'], self.bundle.N)
        self.line_p=np.zeros(self.bundle.N) #aka _p_used
        
        self.line_b=np.zeros(self.bundle.N)
        self.line_b_last=np.copy(self.line_b)
        
        self.preamble()
        hold_counter=0
        self.stepsize=self.defaults['min_sw']
        if any(self.rate_targets):
            util.log.info("Running with targets:%s"%str(self.rate_targets))

            while not self.converged():
                self.load_bundle(self.w)
                self.update_totals()
                util.log.info("After LoadBundle, b:%s:%d"%(self.line_b, sum(self.line_b)))
                if self.oscillating():
                    self.stepsize/=2
                    util.log.error("Oscillation Detected, Decreased Stepsize to %lf"%self.stepsize)
                    util.log.info("Regressing W and B")
                    self.w=np.copy(self.w_last)
                    self.line_b=np.copy(self.line_b_last)
                else:
                    self.line_b_last=np.copy(self.line_b)
                    self.w_last = np.copy(self.w)
                    if (hold_counter>0):
                        util.log.info("Holding...")
                        hold_counter-=1
                    else:
                        self.stepsize+=1
                        util.log.info("Increased stepsize to :%lf"%self.stepsize)
                self.w=self.update_w(self.w)
        else: #No Rate Targets Given
            util.log.info("No Target Given")
            self.load_bundle(self.w) #w defaults to 1 anyway, so no effect
        
        self.postscript()
    
    def converged(self):
        '''
        Checks if we are within the rate-margin
        :I can see this being a problem if attempting parallelisation; non singular entry exit.
        '''
        rate_tol = self.defaults['rate_tol']
        for line in range(self.bundle.N):
            if not self.rate_targets[line] == False : #If rate has been set
                if (self.rate_targets[line]-self.line_b[line])>0:
                    return False
        return True
    
    def oscillating(self):
        '''
        Simple enough, return true if the results are oscillating
        #FIXME Parallelize
        '''
        osc_tol=self.defaults['osc_tol']
        delta_b=(self.line_b-self.rate_targets)
        delta_b_last=(self.line_b-self.rate_targets)
        osc=0
        osc += sum(True for n in range(self.bundle.N) if (delta_b[n]>osc_tol and -delta_b_last[n]>osc_tol) ) #Loads are spiking down
        osc += sum(True for n in range(self.bundle.N) if (delta_b_last[n]>osc_tol and -delta_b[n]>osc_tol) ) #Loads are spiking up        
        util.log.info("Found %d oscillating channels"%osc)
        return (osc>=2)
    
    def load_bundle(self,weights):
        '''
        Apply the current weights and 'pour' power into the bundle
        '''
        #mipb.c::load
        
        #Bring some variables into local scope to speed things up
        MAXBITSPERTONE = self.MAXBITSPERTONE
        bit_inc=self.defaults['bit_inc']
        
        #reset data for new weight run
        self.delta_p=np.zeros((self.bundle.K,self.bundle.N,self.bundle.N))
        self.finished=np.tile(False,(self.bundle.K,self.bundle.N))                
        self.cost=np.zeros((self.bundle.K,self.bundle.N),util.CostValue)
        self.p=np.zeros((self.bundle.K,self.bundle.N))
        self.b=np.zeros((self.bundle.K,self.bundle.N))
        self.line_p=np.zeros(self.bundle.N) #aka _p_used
        self.line_b=np.zeros(self.bundle.N)
        
        #Initialise DeltaP (and cost matrix?)
        for k in range(self.bundle.K):
            self.update_delta_p(k)


        self.update_cost_matrix(weights)

        while not (self.finished == True).all():
            util.log.debug("LineP:\n%s"%self.line_p)
            #and return mins = (k_min,user_min)
            mins=self.min_tone_cost()
            if mins: #If a minimum is found
                (k_min,n_min)=mins
                self.b[mins]+=bit_inc
                assert (self.b[mins] <= MAXBITSPERTONE), "You cannot break the laws of physics! %d:%d,%d"%(self.b[mins],mins[0],mins[1])
                #according to mipb.c, greedy invokes update_wp(), not implemented here as per AMK
                
                #Recalculate delta-p's for all lines wrt this bit/power change
                #NB mins[0] -> min tone
                #FIXME PARALLELISE
                self.update_delta_p(k_min)
                #Update the powers 
                self.update_power(mins)
                self.update_cost_matrix(weights, k_min)  
            else:
                util.log.info("%s"%self.cost)
                #assert (self.finished == True).all(), "No min, but not finished"
                pass #Completed this run, all tones full
        #end while tones not full loop               
                    
    def _constraints_broken(self,tone,line):
        '''
        Check power, mask, constraints and return false if broken
        '''
        #Check Power
        for xline in range(self.bundle.N):
            if (self.line_p[xline] + self.delta_p[tone,line,xline] > self.power_budget[xline]):
                #util.log.info("Power Broken: %g + %g > %g,(%d,%d),b=%d"%(self.line_p[xline],self.delta_p[tone,line,xline],self.power_budget[xline],tone,line,self.b[tone,line]))
                return True
            
        #Check Spectral Mask
        if (self.p[tone,line] + self.delta_p[tone,line,line] > self.spectral_mask):
            #util.log.info("Spectrum Broken: %g + %g > %g,(%d,%d),b=%d"%(self.p[tone,line],self.delta_p[tone,line,line],self.spectral_mask,tone,line,self.b[tone,line]))
            return True
        else: #If we got this far, all is well.
            return False
        
    def update_delta_p(self,tone):
        if util.mp:
            self.update_delta_p_MP(tone)
        else:
            for line in range(self.bundle.N):
                self.delta_p[tone,line]=self._calc_delta_p(line,k=tone)
    
    def update_delta_p_MP(self,tone):
        queue=Queue()
        for line in range(self.bundle.N):
            queue.put(line)
        threads = [Process(target=self._calc_delta_p, args=(queue,), kwargs={'k':tone}) for i in range(cpu_count())]
        for p in threads:
            p.start()
        for p in threads:
            p.join()
    
    def update_cost_matrix(self,weights,tone=False):
        '''
        Update cost matrix and return a tuple of the minimum bit addition (tone,user)
        if not given a tone, assume operation of (min_cost/init_cost_matrix)
        if given a tone, assume operation of (recalc_costs)
        Raises NameError on no min found/all tones full

        '''
        
        if not isinstance(tone,bool):
            #In this case, update was called with a tone to update, so
            self.update_tone_cost(weights, tone)
        else:
            #recalculate costs matrix
            for k in range(self.bundle.K):
                self.update_tone_cost(weights,k)
            #util.log.info("%s"%str(self.cost))
    
    def update_tone_cost(self,weights,tone):
        '''
        Experimental Single-shot cost generation
        '''
        N=self.bundle.N
        tonecost=np.zeros(N)
        try:
            for line in range(N):
                #Check some blocking conditions first
                if self._constraints_broken(tone, line) or self.b[tone,line]>=self.MAXBITSPERTONE:
                    self.finished[tone,line]=True
                
                if self.finished[tone,line]:
                    tonecost[line]=self.defaults['maxval']
                else:
                    for xline in range(N):
                        if line==xline:
                            tonecost[line]+=weights[line]*self.delta_p[tone,line,line]
                        elif weights[xline]<1:
                            tonecost[line]+=(1-weights[xline])*self.delta_p[tone,line,xline]
            self.cost[tone]=tonecost
            #util.log.info("Cost on tone %d\n%s"%(tone,str(self.cost[tone])))
            
        except:
            util.log.info("W:%s,DP:%s"%(str(weights),str(self.cost[tone])))
            raise
        #Fancy Functional stuff won't work with ratetarget==false
        
    def min_tone_cost(self):
        #TODO This could be replaced with a few ndarray operations
        min_cost=float(sys.maxint)

        mindex=self.cost.argmin()
        min=np.unravel_index(mindex,self.cost.shape)
        if self.cost[min]<self.defaults['maxval']:
            return min
        else:
            util.log.info("No Cost Minimum Found")
            return False 
        
    def _calc_delta_p(self,line,k=False):
        '''
        (re)Calculate the delta_p matrix for powers in the bundle
        Since this implementation os Bundle.calc_psd does not update
        a global power setting, there is no need for a local old_p
        #FIXME Not Tested        
        '''
        _b=util.mat2arr(self.b[k])
        _b[line]+=1
        new_p=self.bundle.calc_psd(_b,k)
        delta_p=np.zeros(self.bundle.N)
        for xline in range(self.bundle.N):
            if new_p[xline]<0:
                #Eliminate anyone else from talking to this line.
                ##I Don't Think This Makes Any Sense...
                delta_p[xline]=self.defaults['maxval']
            else:
                delta_p[xline]=new_p[xline]-self.p[k,xline]
        return delta_p
        
    def update_power(self,(tone,line)):
        '''
        Update power totals for all lines for this line change
        #FIXME Not Tested
        '''
        this_delta=0.0
        for xline in range(self.bundle.N):
            this_delta=self.delta_p[tone,line,xline]
            self.p[tone,xline]+=this_delta
            self.line_p[xline]+=this_delta
            self.p_ave+=this_delta/self.bundle.N
        
    def update_totals(self):
        '''
        Swapout current total values to last, and update current totals
        '''
        self.line_b_last=np.copy(self.line_b)
        self.w_last = np.copy(self.w)
        self.line_b=np.sum(self.b,axis=0) #FIXME test this
    
    def update_w(self,weights):
        '''
        Update the weights assigned to each line
        '''
        weights=map(self._update_single_w,range(self.bundle.N))
        return weights
    
    def _update_single_w(self,line):
        current=self.w[line]
        ratio=0.0
        if not self.rate_targets[line]==False:
            rate_tol = self.defaults['rate_tol']
            stepsize=self.stepsize
            diff = self.rate_targets[line]-self.line_b[line]
            ratio= diff/self.rate_targets[line]
            #If we're within tolerance, do nothing, otherwise...
            if abs(diff)>rate_tol:
                #TODO Need To Talk To AMK about this; I think its better. And its certainly faster.
                #Need to deal with if b -> 0, just reset the weight and hope it recovers
                if ratio == 1:
                    new=1
                elif ratio>0:
                    new= current*(ratio)
                else:
                    new=current*(1+-ratio)
            else:
                new= current
        else:
            new= current#Keep the line unweighted
        util.log.info("Weight on Line:%d,%f,ratio:%.3f"%(line,new,ratio))

        return new
                