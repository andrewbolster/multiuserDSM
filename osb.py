"""
Algorithm Modules
"""

#Global Imports
import numpy as np
import math
import sys

#Local Imports
from bundle import Bundle
from algorithm import Algorithm
from line import Line
import utility

class OSB(Algorithm):
    """
    Optimum Spectrum Balancing
    """
    def run(self):
        self.name = "OSB"
        #Aimed-for rates per-line
        self.rate_targets = np.tile(100,self.bundle.N)
        #Tolerance of rate-targets (how close is close enough)
        self.target_tolerance=1
        #Default values
        self.defaults={'maxval':sys.maxint*1.0,
                       'l':0.0,
                       'l_min':1.0,
                       'l_max':0.0,
                       'w':1.0/self.bundle.N,
                       'w_max':100,
                       'w_min':0,
                       'p_budget':0.02818, #watts #from scenarios.c
                       'rate_target':False,
                       'min_step':500, #was min_sl
                       'p_tol':0.015
                       }
        
        #rate targeting #TODO
        self.rate_targets=np.tile(self.defaults['rate_target'], self.bundle.N)
        self.power_budget=np.tile(self.defaults['p_budget'], self.bundle.N)

        
        self.preamble
        #Power and bitloading are dim(KxN)
        self.p=np.asmatrix(np.zeros((self.bundle.K,self.bundle.N)))
        self.b=np.asmatrix(np.zeros((self.bundle.K,self.bundle.N)))
        #lambda values and weights are dim(N)
        self.l=np.zeros((self.bundle.N))
        self.w=np.tile(self.defaults['w'],(self.bundle.N))
        #status tracking values
        self.l_best=np.zeros((self.bundle.N))
        self.l_last=np.tile(-self.defaults['maxval'],self.bundle.N)

        self.w_min=np.zeros((self.bundle.N))
        self.w_max=np.zeros((self.bundle.N))
                              
        #TODO Implement rate region searching
        if (all(rate != self.defaults['rate_target'] for rate in self.rate_targets)):
            pass
        else:
            self._bisect_l();
        #init_lines() What the hell is this?
        self.bundle.calculate_snr() #move into postscript?
        self.postscript
        return
    
    """
    Lambda Bisection
    :from osb_bb.c
    """
    def _bisect_l(self):
        utility.log.info("Beginning Bisection")
        p_now = np.zeros(self.bundle.N)
        while not self._l_converged():
            for id,line in enumerate(self.bundle.lines):
                l_min=self.defaults['l_min']
                l_max=self.defaults['l_max']
                self.l[id]=self.defaults['l']
                
                #L-range hunting
                evertested=False
                utility.log.info("Beginning l-range hunt;line:%d"%id)
                while True: #FIXME there must be a better way of doing this
                    self.optimise_p()
                    utility.log.debug("After optimise_p(), p:%s"%str(self.total_power(line)))
                    if self.total_power(line) > self.power_budget[id]:
                        evertested=True
                        if (self.l[id] == self.defaults['l']):
                            self.l[id]+=1 #0*2=0
                        else:
                            self.l[id]<<1 #*=2
                    else:
                        assert(evertested==True)
                        break
                l_max=self.l[id]
                l_min=self.l[id]/2.0
                utility.log.info("Completed l-range hunt; max:%f,min:%f"%(l_max,l_min))

                
                #Actual optimisation
                last=False #force _l_converged to do first loop
                utility.log.info("Beginning optimisation run")            
                while not self._l_converged(line,last):
                    self.l[id]=(l_max+l_min)/2
                    self.optimise_p()
                    power=self.total_power(line)
                    if power > self.power_budget[id]:
                        l_min=self.l[id]
                    else:
                        l_max=self.l[id]
                    last=power
            #End line loop
            p_now=np.asmatrix(map(self.total_power,self.bundle.lines))
        #End while loop
            
                    
                    
    """
    l converged
    Decides whether the line/bundle is done
    """
    def _l_converged(self,line=False,last=False):
        if isinstance(line,Line): #if called with a line
            if last == False: 
                return False #Force the first loop through
            else:
                howfar = abs(self.power_budget[line.id]-self.total_power(line))
                return (self.total_power(line) == last or howfar < self.defaults('p_tolerance'))
        else:
            #if called without a line, assume operation on the bundle
            if (self.l_last == self.l).all(): 
                #Optimisation done since all values the same as last time
                return True
            else:
                #TODO Need to add rate checking in here for rate mode
                return False
            
    """
    Total Power: return a line's planned total power
    """
    def total_power(self,line=False):
        if isinstance(line,Line):
            #swap the tone and line dimensions for cleanliness
            #FIXME Make sure this works when N>MAXBITSPERTONE, not sure of the shape of this.
            power = np.add.reduce(self.p)[0,line.id]
        else: 
            power = np.add.reduce(np.add.reduce(self.p),axis=1)[0,0]
        assert(isinstance(power,np.float64))
        return power
    
    
                
    """
    Optimise Power (aka optimise_s)
    :from OSB_original.pdf paper
    """   
    def optimise_p(self):
        #for each subchannel
        for k in range(self.bundle.K): #Loop in osb_bb.c:optimise_p
            #utility.log.info("optimise_p(k:%d)"%k)
            lk_max=-self.defaults['maxval']
            b_max=[]
            #for each bit combination
            b_combinator=utility.combinations(range(self.MAXBITSPERTONE), self.bundle.N)
            for b_combo in b_combinator:
                b_combo=np.asarray(b_combo)
                lk=self._l_k(b_combo,k)
                #utility.log.debug("LK/LK_max/combo/b_max:%s %s %s %s"%(lk,lk_max,b_combo,b_max))
                if lk >= lk_max:
                    lk_max=lk
                    b_max=b_combo
            #By now we have b_max[k]
            self.p[k]=self.bundle.calc_psd(b_max,self.w,k,self.p)
            utility.log.info("Max[k=%d][lk=%s][bcombo=%s]%s"%(k,str(lk),str(b_max),str(self.p[k])))
            self.b[k]=b_max
        #Now we have b hopefully optimised
            
    """
    L_k; Return the Lagrangian given a bit-loading combination and tone
    """
    def _l_k(self,bitload,k):
        #use a local p for later parallelism
        #utility.log.debug("bitload,w,k,p:%s,%s,%s,%s"%(str(bitload),str(self.w),str(k),str(self.total_power())))
        p=self.bundle.calc_psd(bitload,self.w,k,self.p)
        #If anything's broken, this run is screwed anyway so feed optimise_p a bogus value
        if (p < 0).any(): #TODO Spectral Mask
            return -self.defaults['maxval']
        
        # p is a matrix of power values required on each line[k] to support bitload[line] bits
        # l is the current lambda array #TODO parallelise?
        # bitload is the current bitload array
        # w is the gamma-weight array
        lp=np.add.reduce(self.l*utility.mat2arr(p))
        #lp=\sum_i^N{l_i*p_i}
        bw=np.add.reduce(bitload*self.w)
        #bw=\sum_i^N{bitload_i*w_i}

        utility.log.debug("LP:%s,BW:%s,p:%s,l:%s"%(str(lp),str(bw),str(p),str(self.l)))
        lk=bw-lp
        
        assert(isinstance(lk,np.float64))
        return lk
       
