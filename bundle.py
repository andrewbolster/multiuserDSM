"""
Class that describes the DSL bundle object
"""

import sys
from math import all
from utility import all
import cmath
from line import Line

class Bundle(object):
    def __init__(self,network_file,K):
        self.K = K                # the number of DMT channels
        self.lines = []            # the DSL line objects
        self.xtalk_gain = None    # XT Matrix (NB Must be FULLY instantiated before operating)

        #Assuming that each bundle is going to be one medium
        """
        Material Parameters
        """
        self.material={
            "r_0c"    :0,
            "a_c"    :0,
            "r_0s"    :0,
            "a_s"    :0,
        
            "l_0"    :0,
            "l_inf"    :0,
        
            "b"        :0,
            "f_m"    :0,
        
            "c_inf"    :0,
            "c_0"    :0,
            "c_e"    :0,
        
            "g_0"    :0,
            "g_e"    :0
            }

        """
        Try to open and parse the network configuration file
        """
        try:
            with open(network_file,"r") as nf:
                for line in nf:
                    # nt and lt are network termination (CO or RT) and line termination (CPE)
                    nt,lt = line.split(",")
                    self.lines.append(Line(nt,lt))

            self.N = len(self.lines)

        except IOError:
            print "Cannot open the network file",network_file
            sys.exit(1)


        """
        Calculate the channel matrix
        """
        self.xtalk_gain = zeroes(self.N,self.N,self.K)
        self.channel_matrix = self.calc_channel_matrix()


    def calc_channel_matrix(self): #TODO
        for K in range(self.K):
            for i,li in enumerate(self.lines):
                for j,lj in enumerate(self.lines):
                    if i == j:
                        """
                        TODO:
                        Could the li.gain assignment be moved inside the object?
                        """
                        li.gain = self.xtalk_gain[i,j,k] = li.transfer_fn(freq_on_tone(K))
                    else:
                        self.xtalk_gain[i,j,k] = self.calc_fext_xtalk_gain(li,lj,K,DOWN) #This makes more sense in passing line objects instead of id's



    def freq_on_tone(self,K): #TODO
        """
        Assume ADSL downstream for now
        """
        return K * 4312.5 + 140156.25;

    def calc_fext_xtalk_gain(self,victim,xtalker,freq,dir): #TODO Currently does UPSTREAM only, swap lt/nt for victim and xtalker for downstream
        """
        Calculate Far End XT Gain
        
        TODO: Look into "Impact of Crosstalk Channel Estimation on the DSM
        Performance for DSL Networks" for an alternative
        
        Uses:
            transfer_func
            insertion_loss
            fext
        (Upstream/Downstream) Cases:
            victim and xtalker == lt
                victim length < xtalker length (4/8)
                victim length > xtalker length (6/2)
                victim and xtalker == length (5)
            victim lt < xtalker lt
                victim nt < xtalker nt (1/9)
                victim nt < xtalker nt (3)
            victim lt > xtalker lt
                victim nt > xtalker nt (9/1)
                victim nt < xtalker nt (7)
            victim and xtalker == nt
                victim length < xtalker length (8/4)
                victim length > xtalker length (2/6)
            victim nt <= xtalker lt (0)
            victim lt >= xtalker nt (0)
        
        What do the cases mean?
        http://www.nordstrom.nu/tomas/publications/BalEtAl_2003_ETSITM6_033w07.pdf
        "Proposed Method on Crosstalk Calculations in a Distributed Environment"
        In Section 2: Top line is 'victim', A->B is nt->lt
        
        
        4 bundle segment types for 2 lines (xtalker and victim);
        H1:xtalker before victim (1,4,7)
        H2:xtalker and victim (all)
        H3:victim after xtalker(1,2,3)
        H4:victim before xtalker(3,6)
        
        """
        
        if dir == "DOWNSTREAM":                 #If DOWNSTREAM, swap nt/lt's
            (dnt,dlt)=(xtalker.lt,xtalker.nt)
            (vnt,vlt)=(victim.lt,victim.nt)
        else:
            (dnt,dlt)=(xtalker.nt,xtalker.lt)
            (vnt,vlt)=(victim.nt,victim.lt)
            
        
        h1 = self.insertion_loss(abs(dnt-vnt), freq)
        h2 = self.insertion_loss(max(vnt,dnt)-dlt, freq)
        h3 = self.insertion_loss(dlt-vlt, freq)
        """
        This _H takes away the biggest pain of the fext gain cases; 
        H1 > 0 in cases 1,3,4,6,7 (also acting as H4 using abs())
        H2 active in all cases, but using max(v.nt,x.lt)
        H3 > 0 in cases 1,2,3
        """
        H = h1*h2*h3
        
        gain = H * self.fext(freq, max(vnt,dnt)-max(vlt,dlt))
        return gain
    
    
    
    
    def fext(self,freq,length):
        """
        Calculate FEXT for given Freq and Length
        Uses:
            UndB
        """
        
        """
        Constants taken from
            “Transmission and Multiplexing (TM); Access transmission systems on metallic
            access cables; Very high speed Digital Subscriber Line (VDSL); Part 1: Functional
            requirements,”
        """
        K_xn=0.0032 #NOTUSED
        K_xf=0.0056 #NOTUSED
        
        """
        Model from BT's simulation parameters document cp38-2
        """
        x=-55
        n=(lines-1)/(max_bundle_size-1) #n disturbers REFACTOR lines max_bundle_size
        x+= 20*log10(freq/90e3)         #f in Hz
        x+= 6*log10(n)                  #shared length in km
        x+= 10*log10(length)
        
        try:
            return UndB(x)
        except ZeroDivisionError:
            return 0
        
        
    """
    insertion loss makes more sense in the bundle as most of the time,
    it is looking at a common length between two lines
    """
    def insertion_loss(self,length,freq): #TODO: Reintroduce factor removed from C version?
        if length > 0:
            return utility.do_transfer_fn(length,freq )
        else: 
            return 0
    
    """
    Calculate snr statistics for all lines in bundle
    """
    def calculate_snr(self):
        for line in self.lines:
            line.sanity()
            
            noise = line.calc_fext_noise(k) + line.alien_xtalk(k) + dbmhz_to_watts(line.noise) #TODO alien
            
            line.cnr = [ line.gain[x]/noise for x in range(self.K)]
            line.snr = [ dbmhz_to_watts(line.psd[x])*line.cnr[x] for x in range(self.K)]
            line.gamma_m = [ 10*math.log10(line.snr[x]/math.pow(2,line.b[x]-1)) for x in range(self.K)]
            line.symerr = [ self._calc_sym_err(line,xtalker) for xtalker in range(self.lines) ] #TODO symerr.c
            
            line.p_total = sum(map(dbmhz_to_watts(line.psd)))
            line.b_total = sum(line.b)
            line.rate[line.service] = "sum of line.b[k]" #TODO How to vectorise this? Where Does Service Come From?
            if line.is_frac():
                line.rate[line.service] =  "sum of line._b[k]"
            """
            With This whole Fractional thing; 
            Are b/_b ever different? If its a fractional line is b[k] ever needed?
            """
    """
    Calculate Symbol Error Rate between two lines
    THERE MUST BE A CLEARER WAY!
    """
    def _calc_sym_err(self,xtalker): #TODO
        return 0
            