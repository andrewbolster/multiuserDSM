"""
Algorithm Modules
"""

import bundle
import math

class OSB(Algorithm):
    """
    Iterative Water Filling
    """
    
    def __init__(self):
        self.name = "OSB"
        self.preamble
        for w in xrange(0,1,0.0001):
            self.optimise_L1(w)
        self.postscript
        return   
    
    def _calc_w_distance(self,w):
        #sum of square of differences in rate targets
        sum=reduce(lambda x,y: x+y, map(lambda t,c: math.pow(t-c,2),target_rates,current_rates))
        return math.sqrt(sum)