'''
Algorithm Modules
'''

#Global Imports
import numpy
import math

#Local Imports
from bundle import Bundle
from algorithm import Algorithm

class ISB(Algorithm):
    '''
    Iterative Spectrum Balancing
    '''
    def __init__(self):
        self.name = "ISB"
        self.preamble
        
        self.current_rate = numpy.zeros(self.bundle.K)
        