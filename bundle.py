"""
Class that describes the DSL bundle object
"""

#Global imports
import sys
import cmath
import numpy
import pydot
import pyparsing
import scipy.special as ss
import pylab
from math import pow,log10,sqrt

#Local Imports
from utility import log,freq_on_tone,dbmhz_to_watts,do_transfer_function, UndB
from line import Line

class Bundle(object):
    def __init__(self,network_file,K=512,rates_file=None, weights_file=None):
        self.K = int(K)                # the number of DMT channels
        self.N = 0                  #Filled in on loading file, here to remind me
        self.lines = []            # the DSL line objects
        self.xtalk_gain = []    # XT Matrix (NB Must be FULLY instantiated before operating)

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
                for n,line in enumerate(nf):
                    # nt and lt are network termination (CO or RT) and line termination (CPE)
                    nt,lt = line.split(",")
                    self.lines.append(Line(nt,lt,n,self))

            self.N = len(self.lines)
            log.info("Read %d Lines",self.N)
            

        except IOError:
            log.error("Cannot open the network file: %s",network_file)
            sys.exit(1)


        log.info("Calculating the channel matrix")
        #Initialise the gains
        self.xtalk_gain = numpy.zeros((self.N,self.N,self.K))
        #The real work begins
        self.calc_channel_matrix()
        
        log.info("Printing xtalk_gain")
        print self.xtalk_gain
        #self.graph_channel_matrix()
        log.info("Running self check:")
        self.check_xtalk_gains() #This is only ever used once; could be sent into calc_channel_matrix?
        


    """
    Calculates the bundle channel gain matrix, generating xtalk_gain[][][]
    :from channel_matrix.c
    """
    def calc_channel_matrix(self): #TODO
        for k in range(self.K):                         #For Each Tone
            for i,li in enumerate(self.lines):          #Between Every Victim
                for j,lj in enumerate(self.lines):      # and every xtalker
                    if i == j:                          #If you're talking to yourself, do lazy transfer_fn
                        li.gain[k] = self.xtalk_gain[i,j,k] = li.transfer_fn(freq_on_tone(k))
                    else:                               #Otherwise look at XT
                        self.xtalk_gain[i,j,k] = self.calc_fext_xtalk_gain(li,lj,freq_on_tone(k),"DOWNSTREAM") #This makes more sense in passing line objects instead of id's
                    log.debug("Set channel %d,%d,%d to %g",i,j,k,self.xtalk_gain[i,j,k])
    
    """
    Check Normalised XT Gains
    :from channel_matrix.c
    """
    def check_xtalk_gains(self):
        yeses=0
        """ Original Way to do it
        for k in self.K:            #tone
            for x,lineX in enumerate(self.lines):    #xtalker
                for v,lineV in enumerate(self.lines):#victim
                    if x==v: continue
                    if self.xtalk_gain[x,v,k]/lineV.gain[k] > 0.5:
                        yeses+=1
        """
        
        """Lets turn this on its head"""
        for v,victim in enumerate(self.lines):
            #listcomprehension for x,v,k on xtalk_gains and gain[k]
            gainratio=[self.xtalk_gain[x,v,k]/victim.gain[k] for x in range(self.N) for k in range(self.K)]
            yeses+=len([1 for i in gainratio if i>0.5])
        
        log.info("Total:%d,%%Yes:%f%%"%(len(gainratio),yeses/(1.0*len(gainratio))))
               
    
    """
    Get XT Gain for line objects
    """
    def get_xtalk_gain(self,line1,line2,channel):
        return self.xtalk_gain(line1.id, line2.id,channel)
    
    """
    Calculate Far End XT Gain
    :from channel_matrix.c 
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
    def calc_fext_xtalk_gain(self,victim,xtalker,freq,dir): #TODO Currently does UPSTREAM only, swap lt/nt for victim and xtalker for downstream
               
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
        log.debug("Returned xtalk_gain: %g",gain)

        return gain  
    
    """
    Calculate FEXT for given Freq and Length
    Uses:
        UndB
    Model from BT's simulation parameters document cp38-2 http://www.niccstandards.org.uk/files/current/NICC%20ND%201513%20(2010-01).pdf
    :from channel_matrix.c
    """
    def fext(self,freq,length):
       
        x=-55   #This is different from the NICC Model  #FUDGE
        x+= 20*log10(freq/90e3)         #f in Hz
        x+= 6*log10(len(self.lines)-1)/(49) #n disturbers 
        x+= 10*log10(length)            #shared length in km
        
        try:
            return UndB(x)
        except ZeroDivisionError:
            return 0
        
    """
    Calculate Insertion Loss
    :from transfer_fn.c
    """    
    def insertion_loss(self,length,freq): #TODO: Reintroduce factor removed from C version?
        """
        insertion loss makes more sense in the bundle as most of the time,
        it is looking at a common length between two lines
        """
        if length > 0:
            return do_transfer_function(length,freq )
        else: 
            return 1 #Leads straight into multiplication; so ends up x*0=0 otherwise
    
    """
    Calculate snr statistics for all lines in bundle
    :from snr.c
    """
    def calculate_snr(self):
        for line in self.lines:
            line.sanity()
            
            noise = [ line.calc_fext_noise(k) + line.alien_xtalk(k) + dbmhz_to_watts(line.noise) for k in range(self.K)] #TODO alien
            
            line.cnr = [ line.gain[k]/noise[k] for k in range(self.K)]
            line.snr = [ dbmhz_to_watts(line.psd[k])*line.cnr[k] for k in range(self.K)]
            line.gamma_m = [ 10*log10(line.snr[k]/pow(2,line.b[k]-1)) for x in range(self.K)]
            line.symerr = [ self._calc_sym_err(line,xtalker) for xtalker in range(self.lines) ] #TODO _calc_sym_err
            log.debug("Symbol Err for n:%d,%lf"%(line.id,line.symerr))
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
    Calculate Symbol Error Rate on line
    :from symerr.c
    """
    def _calc_sym_err(self,line,k): #TODO
        
        M=pow(2,line.b[k])
               
        return 1 - (1 - (1 - 1/sqrt(M))*ss.erf(sqrt((3*line.snr[k])/(2*(M-1)))))
    
    """
    Checks inter-service margins
    :from snr.c
    #I suspect this can be ignored
    """
    def check_all_margins(self,gap):
        pass
            
    
    """
    Pretty Print channel Matrix
    #TODO I've got no idea how to display this....
    """
    def graph_channel_matrix(self):
        
        pylab.contourf(self.xtalk_gain[...,0])
        pylab.figure()
        pylab.show

            