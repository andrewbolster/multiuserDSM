'''
Created on 22 Mar 2011

@author: bolster
'''
#Global Imports
import numpy as np
import sys
from itertools import product

#Local Imports
from bundle import Bundle
from algorithm import Algorithm
from line import Line
import utility as util

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
                         'bit_inc':1.0,
                         'rate_target':False,
                         'p_budget':0.110,    #watts #from scenarios.c (0.110
                         'w':1.0,
                         'w_min':0.0,
                         'min_sw':1e-4
                         }
        dim2=(self.bundle.K,self.bundle.N)
        
        self.spectral_mask=0.0 #FIXME Need a value for this.
        self.p_ave=0.0
        #Finished[tone,user]=False
        self.finished=np.tile(False,dim2)
        self.delta_p=np.zeros((self.bundle.K,self.bundle.N,self.bundle.N))
        self.w=np.tile(self.defaults['w'],self.bundle.N)
        self.w_last=np.tile(self.defaults['w_min'],self.bundle.N)
        
        self.cost=np.zeros((self.bundle.K,self.bundle.N))

        self.power_budget=np.tile(self.defaults['p_budget'], self.bundle.N)
        self.line_p=np.zeros(self.bundle.N) #aka _p_used
        
        self.line_b=np.zeros(self.bundle.N)
        self.line_b_last=np.copy(self.line_b)
        self.rate_targets=np.tile(self.defaults['rate_target'], self.bundle.N)

        
        self.preamble()
        hold_counter=0
        self.stepsize=self.defaults['min_sw']
        while not self.converged():
            self.load_bundle(self.w)
            self.update_totals()
            util.log.debug("After LoadBundle, lineb:%s"%self.line_b)
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
                    self.stepsize+=self.stepsize #*=2
                    util.log.debug("Increased stepsize to :%lf"%self.stepsize)
            self.w=self.update_w(self.w)
    
    def converged(self):
        '''
        Checks if we are within the rate-margin
        since rate targets not implemented yet, returns false if the stepsize
        the default, and true if its been changed
        '''
        return (self.stepsize > 1)
    
    
    def oscillating(self):
        '''
        Simple enough, return true if the results are oscillating
        #FIXME Not implemented since rate_targets arn't implemented, re AMK
        '''
        return False
    
    def load_bundle(self,weights):
        #Algorithm Specific Initialisation
        #mipb.c::load
        
        #Bring some variables into local scope to speed things up
        MAXBITSPERTONE = self.MAXBITSPERTONE
        util.log.info("load_bundle:%s"%str(weights))
        
        #Initialise DeltaP
        for kn in product(range(self.bundle.K),range(self.bundle.N)):
            self._calc_delta_p(*kn)

        while not (self.finished == True).all():
            util.log.info("Still not finished")

            #Generate Cost Calculations and return mins = (k_min,user_min)
            mins=self.update_cost_matrix(weights)
            (k_min,n_min)=mins
            util.log.debug("Made it through update with min cost of %lf"%(self.cost[mins]))
            self.b[mins]+=self.defaults['bit_inc']
            #Update the powers 
            self.update_power(mins)
            if self.b[mins] >= MAXBITSPERTONE:
                self.finished[mins]=True
            #according to mipb.c, greedy invokes update_wp(), not implemented here as per AMK
            
            #Recalculate delta-p's for all lines wrt this bit/power change
            #NB mins[0] -> min tone
            #FIXME PARALLELISE
            for line in range(self.bundle.N):
                self._calc_delta_p(k_min,line)
            self.update_cost_matrix(weights,k_min) #I think this is duplication wrt update_cost_matrix
        #end while tones not full loop
    
    def update_cost_matrix(self,weights,tone=False):
        '''
        Update cost matrix and return a tuple of the minimum bit addition (tone,user)
        if not given a tone, assume operation of (min_cost/init_cost_matrix)
        if given a tone, assume operation of (recalc_costs)
        #FIXME Not Done

        '''
        min_cost=float(sys.maxint)
        
        #recalculate costs matrix
        for kn in product(range(self.bundle.K),range(self.bundle.N)):
            self.cost[kn]=self._cost_function(kn,weights)
            
        '''
        Since the 'cost function' is simply the sum of delta_p's, does it not 
        make more sense to do this as all together?
        '''    
        
        if not isinstance(tone,bool):
            #In this case, update was called with a tone to update, so
            return
        else:
            #TODO This could be replaced with a few ndarray operations
            for kn in product(range(self.bundle.K),range(self.bundle.N)):
                if not self.finished[kn]:
                    assert self.cost[kn] > 0, "Non-positive cost value for (k,n):(%d,%d)"%kn
                    if self.cost[kn]<min_cost:
                        if not self._constraints_ok(*kn):
                            self.finished[kn]=True
                        else:
                            min_cost=self.cost[kn]
                            min=kn
            try:
                return kn
            except NameError:
                util.log.critical("No Cost Minimum Found")
                raise StandardError                    
                    
    def _constraints_ok(self,tone,line):
        '''
        Check power, mask, constraints and return false if broken
        '''
        if (self.line_p + self.delta_p[tone,line] > self.power_budget).any():
            return False
        elif (self.p[tone,line] + self.delta_p[tone,line,line] > self.spectral_mask):
            return False
        else:
            return True
        
        
    def _cost_function(self,(tone,line),weights):
        '''
        Cost function; returns the the current cost of adding a bit to this tone
        line combination.
        Assumes that constraint test has already been applied to the bundle
        '''
        #Unless I'm being retarded, all that cost_function does is sum the deltap's along the last dimension...
        #util.log.debug("Hello, I'm the Cost Function, trying to sum this:%s"%self.delta_p[tone,line,:])
        return np.sum(self.delta_p[tone,line,:])*weights[line]
                
        
    def _calc_delta_p(self,tone,line):
        '''
        (re)Calculate the delta_p matrix for powers in the bundle
        Since this implementation os Bundle.calc_psd does not update
        a global power setting, there is no need for a local old_p
        #FIXME Not Tested        
        '''
        _b=util.mat2arr(self.b[tone])
        _b[line]+=1
        new_p=self.bundle.calc_psd(_b,tone)
        
        #This could be very wrong
        self.delta_p[tone,line]=new_p-self.p[tone,:] #update delta_p in a slice rather than looping
        
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
        #FIXME not implemented yet
        '''
        return weights
        