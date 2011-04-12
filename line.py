'''
Class that describes a DSL Line object
'''
#Global Imports
import numpy

#Local Imports
import utility

class Line(object):
    def __init__(self,nt,lt,rate,id,bundle):
        self.MINPSD=-60     #from am_load.c
        
        # nt and lt are network termination (CO or RT) and line termination (CPE)
        self.bundle = bundle
        self.nt = int(nt)
        self.lt = int(lt)
        self.length = abs(self.nt - self.lt)
        self.rate_target=rate
        self.gain = numpy.zeros(bundle.K) 
        self.p = numpy.tile(-36.5,bundle.K) #PSD of this line, initial value
        self.id = id #could this be removed as an array of lines?
        self.noise = utility.dbmhz_to_watts(-140) #Assuming standard background noise
        self.type = 3   #I'm assuming this declares the material of the line
                        #so in theory it can be removed from transfer_fn

        self.p_total = 0    #configured in bundle.calculate_snr
        self.b_total = 0    #ditto
        
        self.p_max = 1 #P=Total AFE power / delta_Freq #TODO
        
        self.snr = numpy.zeros(bundle.K) 
        self.cnr = numpy.zeros(bundle.K) 
        self.gamma_m = numpy.zeros(bundle.K) 
        self.b = numpy.tile(2,bundle.K) #bit loading on each subchannel
        
        
        
    def __str__(self):
        '''
        Super cool graphics
        '''
        s = ""
        s += "-" * (self.nt / 250)
        s += "|"
        s += "-" * ((self.lt-self.nt) / 250)
        s += "|"
        return s
    
    '''
    Return the RLCG Parameterised Transfer Function
    '''   
    def transfer_fn(self, freq):
        #C Version uses bundle as a linked list of lines so most of this function is redundant.
        return utility.do_transfer_function(abs(self.length), freq ,measure="m", type=self.type)
    
    '''
    Line Sanity Check
    :from snr.c check_line_sanity
    Called from bundle
    '''
    def sanity(self):
        assert (self.b >= 0).all(), utility.log.error("Tone less than zero")
        assert (self.gain >= 0).all(), utility.log.error("Gain less than zero")
        utility.log.debug("Line %d is sane"%self.id)
        return
    
    '''
    Calculate Far End XT noise on a channel in watts
    :from snr.c
    '''
    def calc_fext_noise(self,k):
        noise=0
        for xtalker in self.bundle.lines:
            if xtalker == self : 
                continue   #stop hitting yourself!
            else:
                noise += utility.dbmhz_to_watts(self.p[k])*self.bundle._h2(self,xtalker,k)
        return noise
    
    '''
    Return AlienXT noise
    :from snr.c
    '''
    def alien_xtalk(self,channel): #TODO
        return 0
    
    '''
    Return current rate on line
    '''
    def rate(self):
        return sum(self.b)