"""
Class that describes a DSL Line object
"""
#Global Imports
import numpy

#Local Imports
import utility

class Line(object):
    def __init__(self,nt,lt,id,bundle):
        # nt and lt are network termination (CO or RT) and line termination (CPE)
        self.bundle = bundle
        self.nt = int(nt)
        self.lt = int(lt)
        self.length = self.nt - self.lt
        self.gain = numpy.zeros(bundle.K) 
        self.p = numpy.zeros(bundle.K) #PSD of this line
        self.id = id #could this be removed as an array of lines?
        self.type = 2   #I'm assuming this declares the material of the line
                        #so in theory it can be removed from transfer_fn

        self.p_total = 0    #configured in bundle.calculate_snr
        self.b_total = 0    #ditto
        
        self.p_max = 1 #P=Total AFE power / delta_Freq #TODO
        
        self.snr = []
        self.cnr = []
        self.gamma_m = []
        self.b = numpy.zeros(bundle.K) #bit loading on each subchannel
        
    def __str__(self):
        """
        Super cool graphics
        """
        s = ""
        s += "-" * (self.nt / 250)
        s += "|"
        s += "-" * ((self.lt-self.nt) / 250)
        s += "|"
        return s
    
    """
    Return the RLCG Parameterised Transfer Function
    """   
    def transfer_fn(self, freq):
        #C Version uses bundle as a linked list of lines so most of this function is redundant.
        return utility.do_transfer_function(abs(self.length), freq ,type=self.type)
    
    """
    Line Sanity Check
    :from snr.c check_line_sanity
    Called from bundle
    """
    def sanity(self):
        for t in range(self.bundle.K):
            assert self.b[t] >= 0, "Tone "+t+" less than zero:"+self.b[t]
            assert self.gain[t] >= 0, "Gain "+t+" less than zero"+self.gain[t]
        assert self.background_noise >=0, "Background noise is not set on line "+self.id
        utility.log.debug("Line %d is sane"%self.id)
        return
    
    """
    Calculate Far End XT noise on a channel
    :from snr.c
    """
    def calc_fext_noise(self,channel):
        noise=0
        for xtalker in self.bundle.lines:
            if xtalker == self : continue   #stop hitting yourself!
            else:
                noise += utility.dbmhz_to_watts(self.p[channel])*self.bundle.get_xtalk_gain(xtalker,self,channel)
        return noise
    
    """
    Return AlienXT noise
    :from snr.c
    """
    def alien_xtalk(self,channel): #TODO
        return 0
    
    """
    Return current rate on line
    """
    def rate(self):
        return sum(self.b)