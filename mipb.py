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
                         'bit_inc':1,
                         'rate_target':False,
                         'p_budget':0.110,    #watts #from mipb.c _p_budget[user]=
                         'w':1.0,
                         'w_min':0.0,
                         'min_sw':1e-4
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
        self.rate_targets=np.tile(self.defaults['rate_target'], self.bundle.N)

        
        self.preamble()
        hold_counter=0
        self.stepsize=self.defaults['min_sw']
        if any(self.rate_targets):
            util.log.info("Running with targets:%s"%str(self.rate_targets))

            while not self.converged():
                self.load_bundle(self.w)
                self.update_totals()
                util.log.info("After LoadBundle, lineb:%s"%self.line_b)
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
                if abs(self.rate_targets[line]-self.line_b[line])>rate_tol:
                    return False
        return True
    
    def oscillating(self):
        '''
        Simple enough, return true if the results are oscillating
        '''
        osc_tol=self.defaults['osc_tol']
        delta_b=(self.line_b-self.rate_targets)
        delta_b_last=(self.line_b-self.rate_targets)
        osc=0
        osc += sum(delta_b>osc_tol and -delta_b_last>osc_tol) #Loads are spiking down
        osc += sum(delta_b_last>osc_tol and -delta_b>osc_tol) #Loads are spiking up
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
        util.log.info("load_bundle:%s"%str(weights))
        
        #reset data for new weight run
        self.delta_p=np.zeros((self.bundle.K,self.bundle.N,self.bundle.N))
        self.finished=np.tile(False,(self.bundle.K,self.bundle.N))                
        self.cost=np.zeros((self.bundle.K,self.bundle.N))
        self.p=np.zeros((self.bundle.K,self.bundle.N))
        self.b=np.zeros((self.bundle.K,self.bundle.N))
        self.line_p=np.zeros(self.bundle.N) #aka _p_used
        self.line_b=np.zeros(self.bundle.N)
        
        #Initialise DeltaP (and cost matrix?)
        for kn in product(range(self.bundle.K),range(self.bundle.N)):
            self._calc_delta_p(*kn)

        while not (self.finished == True).all():
            util.log.debug("LineP:\n%s"%self.line_p)

            try:
                #Generate Cost Calculations and return mins = (k_min,user_min)
                mins=self.update_cost_matrix_experimental(weights)
                (k_min,n_min)=mins
                util.log.debug("Made it through update with min cost of %lf"%(self.cost[mins]))
                util.log.info("%f%% There"%((util.watts_to_dbmhz(self.defaults['p_budget'])*100)/util.watts_to_dbmhz(self.line_p[mins[1]])))
                self.b[mins]+=bit_inc
                #Update the powers 
                self.update_power(mins)                
                if self.b[mins] == MAXBITSPERTONE:
                    self.finished[mins]=True
                #assert (self.b[mins] <= MAXBITSPERTONE), "You cannot break the laws of physics! %d:%d,%d"%(self.b[mins],mins[0:2])
                #according to mipb.c, greedy invokes update_wp(), not implemented here as per AMK
                
                #Recalculate delta-p's for all lines wrt this bit/power change
                #NB mins[0] -> min tone
                #FIXME PARALLELISE
                for line in range(self.bundle.N):
                    self._calc_delta_p(k_min,line)
            except NameError:
                util.log.info("All Done (Hopefully)")
        #end while tones not full loop
    
    def update_cost_matrix(self,weights,tone=False):
        '''
        Update cost matrix and return a tuple of the minimum bit addition (tone,user)
        if not given a tone, assume operation of (min_cost/init_cost_matrix)
        if given a tone, assume operation of (recalc_costs)
        Raises NameError on no min found/all tones full

        '''
        min_cost=float(sys.maxint)
        
        '''
        Since the 'cost function' is simply the sum of delta_p's, does it not 
        make more sense to do this as all together?
        '''    
        
        if not isinstance(tone,bool):
            #In this case, update was called with a tone to update, so
            self.update_cost_function(weights,tone)
            for n in range(self.bundle.N):
                self.cost[tone,n]=self._cost_function((tone,n),weights)
            return
        else:
            #recalculate costs matrix
            for kn in product(range(self.bundle.K),range(self.bundle.N)):
                self.cost[kn]=self._cost_function(kn,weights)
            #TODO This could be replaced with a few ndarray operations
            for kn in product(range(self.bundle.K),range(self.bundle.N)):
                if not self.finished[kn]:
                    #assert self.cost[kn] > 0, "Non-positive cost value for (k,n):(%d,%d)"%kn
                    if self.cost[kn]<min_cost:
                        util.log.debug("Testing: %d,%d"%kn)
                        if self._constraints_broken(*kn):
                            util.log.debug("And this one is done: %f<%f"%(self.cost[kn],min_cost))
                            self.finished[kn]=True
                        else:
                            util.log.debug("New Min: %f<%f"%(self.cost[kn],min_cost))
                            min_cost=self.cost[kn]
                            min=kn
            try:
                return min
            except NameError:
                util.log.info("No Cost Minimum Found, raising NameError to notify parent")
                raise NameError                    
                    
    def _constraints_broken(self,tone,line):
        '''
        Check power, mask, constraints and return false if broken
        '''
        #Check Power
        for xline in range(self.bundle.N):
            if (self.line_p[xline] + self.delta_p[tone,line,xline] > self.power_budget[xline]):
                util.log.info("Power Broken: %g + %g > %g,(%d,%d)"%(self.line_p[xline],self.delta_p[tone,line,xline],self.power_budget[xline],tone,line))
                return True
            
        #Check Spectral Mask
        if (self.p[tone,line] + self.delta_p[tone,line,line] > self.spectral_mask):
            util.log.info("Spectrum Broken")
            return True
        else: #If we got this far, all is well.
            return False
          
    def _cost_function(self,(tone,line),weights):
        '''
        Cost function; returns the the current cost of adding a bit to this tone
        line combination.
        Assumes that constraint test has already been applied to the bundle
        #FIXME Once DP>0 issue fixed, remove all this down to single return statement
        '''
        #Unless I'm being retarded, all that cost_function does is sum the deltap's along the last dimension...
        #util.log.info("Hello, I'm the Cost Function, trying to sum this:%s"%self.delta_p[tone,line,:])
        dp_sum=np.sum(self.delta_p[tone,line,:])
        if dp_sum <=0:
            #Implies no solution on this tone
            util.log.critical("No solution on tone,line:(%d,%d)"%(tone,line))
            self.finished[tone,line]=True
            return self.defaults['maxval']
        return dp_sum*weights[line]
    
    def update_cost_matrix_experimental(self,weights,tone=False):
        '''
        Update cost matrix and return a tuple of the minimum bit addition (tone,user)
        if not given a tone, assume operation of (min_cost/init_cost_matrix)
        if given a tone, assume operation of (recalc_costs)
        Raises NameError on no min found/all tones full

        '''
        min_cost=float(sys.maxint)
        
        '''
        Since the 'cost function' is simply the sum of delta_p's, does it not 
        make more sense to do this as all together?
        '''    
        if not isinstance(tone,bool):
            #In this case, update was called with a tone to update, so
            self.update_cost_function(weights, tone)
            return
        else:
            #recalculate costs matrix
            for k in range(self.bundle.K):
                self.update_cost_function(weights,k)
            #TODO This could be replaced with a few ndarray operations
            for kn in product(range(self.bundle.K),range(self.bundle.N)):
                if not self.finished[kn]:
                    #assert self.cost[kn] > 0, "Non-positive cost value for (k,n):(%d,%d)"%kn
                    if self.cost[kn]<min_cost:
                        util.log.debug("Testing: %d,%d"%kn)
                        if self._constraints_broken(*kn):
                            util.log.debug("And this one is done: %f<%f"%(self.cost[kn],min_cost))
                            self.finished[kn]=True
                        else:
                            util.log.debug("New Min: %f<%f"%(self.cost[kn],min_cost))
                            min_cost=self.cost[kn]
                            min=kn
            try:
                return min
            except NameError:
                util.log.info("No Cost Minimum Found, raising NameError to notify parent")
                raise NameError 
    
    def update_cost_function(self,weights,tone=False):
        '''
        Experimental Single-shot cost generation
        #NOT TESTED
        '''
        dp_sum=np.sum(self.delta_p[tone],axis=1)
        self.cost[tone]=weights*dp_sum
        for n in range(self.bundle.N):
            if dp_sum[n] <=0:
                #Implies no solution on this tone
                util.log.critical("No solution on tone,line:(%d,%d)"%(tone,n))
                self.finished[tone,n]=True
                self.cost[tone,n]=self.defaults['maxval']

        return
        
        
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
        '''
        weights=map(self._update_single_w,range(self.bundle.N))
        
        return weights
    
    def _update_single_w(self,line):
        if not self.rate_targets[line]==False:
            rate_tol = self.defaults['rate_tol']
            stepsize=self.stepsize
            diff = self.line_b[line]-self.rate_targets[line]
            #If we're within tolerance, do nothing, otherwise...
            if abs(diff)>rate_tol:
                if diff < 0: #if b<r
                    #I don't think this makes sense but implementing it for graphical test
                    if (stepsize*(-diff))>self.w[line]: #if updated weight is going to be non positive
                        return self.w[line]*0.9
                #if b>r or update !>w
                return self.w[line]+(stepsize*(-diff))
                #I reckon the entire thing could be replaced by
                return max(self.w[line]+(stepsize*(-diff)),self.w[line]*0.9)
            else:
                return self.w[line]
                