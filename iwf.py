'''
Algorithm Modules
'''

#Global Imports
import numpy
import math

#Local Imports
from bundle import Bundle
from algorithm import Algorithm
import utility

class IWF(Algorithm):
    '''
    Iterative Water Filling
    '''
    
    def __init__(self):
        self.name = "IWF"
        
        assert self.bundle is Bundle # "You dun goofed the Bundle!"
        
        '''
        Set initial values for power
        '''
        if self.bundle.type == "ADSL_DOWNSTREAM":
            p_initial = 0.111
        else:
            p_initial = 0.02818
        
        self.p = numpy.tile(p_initial, self.bundle.K) 		#size bundle.K
        self.bundle.calculate_snr()
        
        self.rate_targets = numpy.tile()
        
        #TODO start timer system here
        self.preamble()

        iterations = rate_iterations = 0
        '''
        Iteratively Execute am_load_X and check rate targets
        '''
        while iterations < 30 and rate_iterations < 10:
            while self.bundle.check_all_margins(2):
                for i,line in enumerate(self.bundle.lines):
                    if line.rate_target == None:    #TODO Where does line.rate_target[] get set?
                        self.am_load_ra(line)       #am_load_X functions in parent class
                    else:
                        self.am_load_fm(line)
                    self.bundle.calculate_snr()
                iterations+=1
            
            #assume the best
            targets_met=True
            
            '''
            Target Checking
            '''
            for i,line in enumerate(self.bundle.lines): #TODO List Comprehension/Reduce?
                if line.rate_target == None:
                    continue
                if line.rate < line.rate_target:
                    targets_met=False
            
            if targets_met: #we're done here
                utility.log.info("Rate Targets satisfied after %d rate iterations"%rate_iterations)
                break
            
            rate_iterations+=1 #start a new attempt
            
        else: #will only run after 30 failed iterations or 10 failed 
            utility.log.warning("Convergence Failure after "+iterations+" iterations, "+rate_iterations+" rate iterations")
        
        self.postscript()
