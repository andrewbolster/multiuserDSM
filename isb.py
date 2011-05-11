'''
Algorithm Modules
'''

#Global Imports
import numpy as np
import math
import sys
import logging
import hashlib
import multiprocessing


#Local Imports
from bundle import Bundle
from algorithm import Algorithm
from line import Line
from utility import *

class ISB(Algorithm):
    '''
    Iterative Spectrum Balancing
    '''
    def run(self):
        self.name = "ISB"
        #Aimed-for rates per-line
        self.rate_targets = np.tile(100,self.bundle.N)
        #Tolerance of rate-targets (how close is close enough)
        self.target_tolerance=1
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
                       'p_tol':0.0015,      #should be 0.015
                       'rate_tol':10,
                       'GAMMA':self.bundle.GAMMA
                       }
        
        #rate targeting #TODO
        self.rate_targets=np.tile(self.defaults['rate_target'], self.bundle.N)
        self.power_budget=np.tile(self.defaults['p_budget'], self.bundle.N)
        
        if not self.useGPU:
            #Create multithreading pool
            pass
            #self.threadpool=multiprocessing.Pool(multiprocessing.cpu_count())
            
        self.preamble()
        #lambda values and weights are dim(N)
        self.l=np.tile(self.defaults['l'],(self.bundle.N))
        self.w=np.tile(self.defaults['w'],(self.bundle.N))
        #status tracking values
        self.l_best=np.zeros((self.bundle.N))
        self.l_last=np.tile(-self.defaults['maxval'],self.bundle.N)
        self.l_last_last=np.tile(+self.defaults['maxval'],self.bundle.N)

        self.w_min=np.zeros((self.bundle.N))
        self.w_max=np.zeros((self.bundle.N))
                              
        #TODO Implement rate region searching
        if self.rate_search and (all(line.rate_target != self.defaults['rate_target'] for line in self.bundle.lines)):
            self.rate_bisect(self._bisect_l)
        else:
            log.info("Bisection")
            self._bisect_l();
        #init_lines() What the hell is this?
        self.postscript()
        return
    '''
    Lambda Bisection
    :from isb3g.c
    '''
    def _bisect_l(self):
        logit=log.info
        logit("Beginning Bisection")
        self.l=np.tile(self.defaults['l'],(self.bundle.N))
        self.p_total = np.zeros(self.bundle.N)
        self.p_total_last = np.tile(1.0,self.bundle.N)
        
        while not self._l_converged():
            #Send across last-run's values as last
            self.l_last = self.l
            self.p_total_last=self.p_total
            for lineid,line in enumerate(self.bundle.lines):
                l_min=self.defaults['l_min']
                l_max=self.defaults['l_max']
                self.l[lineid]=self.defaults['l']
                lastpower=self.defaults['maxval']                
                #L-range hunting
                logit("Beginning l-range hunt;line:%d"%lineid)
                self.optimise_p(self.l, self.w)
                linepower=self.total_power(line)
                if (linepower > self.power_budget[lineid]): #current is lmin
                    while True: #Until power budget satisfied
                        self.optimise_p(self.l, self.w)
                        linepower=self.total_power(line)
                        if ( linepower > self.power_budget[lineid]): 
                            logit("Overshot power budget:linepower:%.3f,lastdrop:%.0f%%,budget:%s"%((linepower),(100*(lastpower-linepower)/lastpower),str(self.power_budget[lineid])))
                            lastpower=linepower                      
                            if (self.l[lineid] < 1):
                                self.l[lineid]=1 #0*2=0
                            else:
                                self.l[lineid]*=2
                        else:
                            logit("Satisfied power budget:linepower:%.5f,budget:%s"%((linepower),str(self.power_budget[lineid])))
                            break
                    l_max=self.l[lineid]
                    l_min=l_max/2
                else: #current is lmax
                    while True: #Until power budget satisfied
                        self.optimise_p(self.l, self.w)
                        linepower=self.total_power(line)
                        if ( linepower < self.power_budget[lineid]): 
                            lastpower=linepower                      
                            if (self.l[lineid] < 1):
                                #Too low, let above take care of it
                                logit("Zeroed lambda on %d:linepower:%.5f,budget:%s"%(lineid,(linepower),str(self.power_budget[lineid])))
                                break
                            else:
                                self.l[lineid]/=2
                        else:
                            logit("Undershot power budget:linepower:%.5f,budget:%s"%((linepower),str(self.power_budget[lineid])))
                            break
                    l_min=self.l[lineid]
                    l_max=2*l_min
                logit("Completed l-range hunt; max:%f,min:%f"%(l_max,l_min))

                #Actual optimisation
                last=False #force _l_converged to do first loop
                logit("Beginning optimisation run;line:%d"%lineid)           
                while not self._l_converged(line,last):
                    last=self.total_power(line)
                    self.l[lineid]=(l_max+l_min)/2
                    self.optimise_p(self.l, self.w)
                    #assert 1==0, "WIN"
                    linepower=self.total_power(line)
                    if linepower > self.power_budget[lineid]:
                        l_min=self.l[lineid]
                    else:
                        l_max=self.l[lineid]
                logit("Completed optimisation run;line:%d, l:%f"%(lineid,self.l[lineid]))         

            #End line loop
            
        #End while loop
        self.update_b_p()   
    
                
    '''
    Optimise Power (aka optimise_s)
    :from OSB_original.pdf paper
    '''   
    def optimise_p(self,lambdas,weights):
        #Maybe try this? http://sites.google.com/site/sachinkagarwal/home/code-snippets/parallelizing-multiprocessing-commands-using-python
        if (self.useGPU):
            #Multideviceable GPU execution
            (self.p,self.b)=self.bundle.gpu.isb_optimise_p(lambdas,weights,self.bundle.xtalk_gain)
            #log.info("%s"%self.b)
        else:
            self.optimise_p_k(lambdas,weights)            
            #Now we have b hopefully optimised
            
    def optimise_p_k(self,lambdas,weights):
        for k in range(self.bundle.K):
            #Convergence value check
            b_this=np.tile(0,self.bundle.N)
            b_last=np.tile(-1,self.bundle.N)
            log.debug("Launched channel %d search"%k)
            
            #Until convergence of this channels bitload
            while not (b_last==b_this).all():
                b_last=b_this.copy()       
                for line in xrange(self.bundle.N):
                    lk_max=-self.defaults['maxval']
                    b_max=[]
                    #for each bit modification
                    b_this[line]=0
                    while b_this[line] <= self.MAXBITSPERTONE:
                        #The lagrangian value for this bit combination
                        lk=self._l_k(b_this,lambdas,weights,k)
                        if lk >= lk_max:
                            lk_max=lk
                            b_max=b_this[line]
                        b_this[line]+=1
                            
                    #By now we have b_max for this user on this channel
                    b_this[line]=b_max
                #at this point we are hopefully maximised
                #print "CPU LKmax %d:%s:%s:%s"%(k,str(lk_max),str(b_max),str(self.p[k]))
            self.b[k]=b_this
            self.p[k]=self.bundle.calc_psd(b_this,k)
        #end while
        
    def optimise_p_k_alt(self,lambdas, weights):
        #Convergence value check
        b_this=np.tile(0,(self.bundle.K,self.bundle.N))
        b_last=np.tile(-1,(self.bundle.K,self.bundle.N))
        #Until convergence of these channels bitload
        while (b_last!=b_this).any():
            b_last=b_this.copy()       
            for k in range(self.bundle.K):
                log.debug("Launched channel %d search"%k)
                for line in xrange(self.bundle.N):
                    lk_max=-self.defaults['maxval']
                    b_max=[]
                    #for each bit modification
                    b_this[k,line]=0
                    while b_this[k,line] <= self.MAXBITSPERTONE:
                        #The lagrangian value for this bit combination
                        lk=self._l_k(b_this[k],lambdas,weights,k)
                        if lk >= lk_max:
                            lk_max=lk
                            b_max=b_this[k,line]
                        b_this[k,line]+=1
                            
                    #By now we have b_max for this user on this channel
                    b_this[k,line]=b_max
                #at this point we are hopefully maximised
                #print "CPU LKmax %d:%s:%s:%s"%(k,str(lk_max),str(b_max),str(self.p[k]))
        self.b=b_this
        for k in range(self.bundle.K):
            self.p[k]=self.bundle.calc_psd(b_this[k],k)
        #end while
        
