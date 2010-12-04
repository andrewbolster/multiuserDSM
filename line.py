
import utility

class Line(Bundle):
    def __init__(self,nt,lt):
        # nt and lt are network termination (CO or RT) and line termination (CPE)
        self.nt = int(nt)
        self.lt = int(lt)
        self.length = self.nt - self.lt
        self.gain = None
        self.id = none #could this be removed as an array of lines?
        self.type = 2   #I'm assuming this declares the material of the line
                        #so in theory it can be removed from transfer_fn

        self.p_total = 0    #configured in bundle.calculate_snr
        self.b_total = 0    #ditto
        
        self.rate = []
        self.snr = []
        self.cnr = []
        self.gamma_m = []
        self.b = []
        

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

    def transfer_fn(self, freq): #TODO
        """
        Return the RLCG Parameterised Transfer Function
        """
        #C Version uses bundle as a linked list of lines so most of this function is redundant.
        return utility.do_transfer_fn(type = self.type, length=self.length/1000, freq )
    
    def sanity(self):
        for t in range(super(Bundle,K)):
            assert self.b[t] >= 0, "Tone "+t+" less than zero:"+self.b[t]
            assert self.gain[t] >= 0, "Gain "+t+" less than zero"+self.gain[t]
        
        assert self.background_noise >=0, "Background noise is not set on line "+self.id
        return
    
        
        