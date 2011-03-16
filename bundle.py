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
import pprint
import itertools
from math import pow,log10,sqrt

#Local Imports
from utility import *
from line import Line

class Bundle(object):
    N = 0                  #Filled in on loading file, here to remind me
    lines = []            # the DSL line objects
    xtalk_gain = []    # XT Matrix (NB Must be FULLY instantiated before operating)
    #Assuming that each bundle is going to be one medium

    MARGIN = 3.0   #Performance Margin db
    C_G = 0        #Coding Gain db
    def __init__(self,network_file,K=512,rates_file=None, weights_file=None):
        self.K = int(K)                # the number of DMT channels
        self.GAMMA= get_GAMMA(1e-7, 4)       #Uncoded SNR Gap (setup elsewere)
        self.gamma_hat = pow(10,(self.GAMMA+self.MARGIN-self.C_G)/10)


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
        self.print_xtalk()
        self.print_xtalk_onetone(self.K/2)

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
                    #log.debug("Set channel %d,%d,%d to %g",i,j,k,self.xtalk_gain[i,j,k])
        assert(not numpy.allclose(self.xtalk_gain[0].T,self.xtalk_gain[0]))
    
    """
    Check Normalised XT Gains and xtalk symmetry 
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
        
        #Check symmetry if all lines are the same
        ntlt=set([(line.nt,line.lt) for line in self.lines])
        log.info("%s"%str(ntlt))

        if len(ntlt)!=1:
            #Not all the lines are the same so the xtalk matrix cannot be symmetric
            log.info("Lines are different, xtalk should not be symmetric")
            assert (self.xtalk_gain.transpose(1, 0, 2) != self.xtalk_gain).all(), "Xtalk Symmtric"
        else:
            log.info("Lines are identical, xtalk should be symmetric")
            assert (self.xtalk_gain.transpose(1, 0, 2) == self.xtalk_gain).all(), "Xtalk not Symmtric"
              
            
        log.info("Total:%d,%%Yes:%f%%"%(len(gainratio),yeses/(1.0*len(gainratio))))
    
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
            victim nt > xtalker nt (3)
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
            (D1,D2)=(xtalker.lt,xtalker.nt)
            (A,B)=(victim.lt,victim.nt)
        else: #Upstream
            (D1,D2)=(xtalker.nt,xtalker.lt)
            (A,B)=(victim.nt,victim.lt)
            
        #Shared length is the H2 Span
        shared_length=abs(min(B,D2)-max(A,D1))
        #Head Length is the H1/H4 Span
        head_length=abs(A-D1)
        #Tail Length is the H3 Span
        tail_length=B-D2
        
        h1 = self.insertion_loss(head_length, freq)
        h2 = self.insertion_loss(shared_length, freq)
        h3 = self.insertion_loss(tail_length, freq)
        """
        This _H takes away the biggest pain of the fext gain cases; 
        H1 > 0 in cases 1,3,4,6,7 (also acting as H4 using abs())
        H2 active in all cases, but using max(v.nt,x.lt)
        H3 > 0 in cases 1,2,3
        
        NB, insertion_loss(<0,f)=1
        
        I think AMK's Case1/9 fext(length) is wrong. Should be the
        common length between the victim and xtalker, not the combined
        length
        """
        H = h1*h2*h3
        
        gain = H * self.fext(freq, shared_length)
        #log.debug("Returned xtalk_gain: %g",gain)

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
        x+= 6*log10((self.N-1.0)/49)    #n disturbers 
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
            noise = numpy.zeros(self.K)
            
            for tone in xrange(self.K):
                
                noise = line.calc_fext_noise(tone) + line.alien_xtalk(tone) + dbmhz_to_watts(line.noise)
                line.cnr[tone] = line.gain[tone]/noise
                #log.debug("b%d,b%e,g%e"%(line.b[tone],(pow(2,line.b[tone])-1),(self.gamma_hat/line.cnr[tone])))
                line.snr[tone] = dbmhz_to_watts(line.p[tone])*line.cnr[tone]
                print line.calc_fext_noise(tone)
                print line.gain
                print line.cnr
                print line.snr
                print line.p
                print line.b
                
                line.gamma_m[tone] = TodB(line.snr[tone]/pow(2,line.b[tone]-1))
                
            #line.symerr = [ self._calc_sym_err(line,xtalker) for xtalker in xrange(self.K)] #TODO _calc_sym_err
            line.p_total = sum(map(dbmhz_to_watts,line.p))
            line.b_total = sum(line.b)
            
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
    Generate PSD vector matrix between lines and return matrix
    :from psd_vector.c
    """
    def calc_psd(self,bitload,gamma,k,power):
        #Generate A Matrix (See Intro.pdf 2.23)
        A=numpy.asmatrix(numpy.zeros((self.N,self.N)))
        for i in range(self.N):
            for j in range(self.N):
                if i==j: 
                    A[i,j]=1
                else:
                    A[i,j]=self._psd_A_elem(i,j,bitload,gamma,k)

        #Generate B Matrix (same reference)
        B=numpy.asmatrix(numpy.zeros(self.N))
        for i,line in enumerate(self.lines):
            B[0,i]=self._psd_B_elem(line,bitload,gamma,k)
            
        #Transpose B to be single-column matrix (1xN.NxN)
        B=B.T               
        assert(B.shape==(self.N,1))
        #log.debug("A,B:%s,%s"%(str(A),str(B)))
        #Everyone loves linear algebra...dont they?
        if (False): #QR or regular way?
            q,r = numpy.linalg.qr(A)
            assert numpy.allclose(A,numpy.dot(q,r))
            p= numpy.dot(q.T,B)
            P=numpy.dot(numpy.linalg.inv(r),p)
        else:
            P=numpy.linalg.solve(A,B)
        
        #Because I'm paranoid
        assert numpy.allclose(B,(numpy.dot(A,P))), log.error("Turns out linear algebra is hard")
        log.info("A:\n%s"%str(A))
        log.info("B:\n%s"%str(B))
        log.info("P:\n%s"%str(P))
        P=P.T

        assert(P.shape==(1,self.N)),"Non-single-row P:%s"%str(P.shape)
        #assert (P>0).any(),"Zero or Negative P\n%s"%str(P) #TODO Spectral Mask?
        
        #FIXME This should just return [p0 p1 p2...] not [[px...]]
        return P[0]

        
    """
    PSD A-Element
    """
    def _psd_A_elem(self,i,j,bitload,gamma,k):
        return (-self._f(i,bitload,gamma)*self.xtalk_gain[j][i][k])/self.xtalk_gain[i][i][k]
    
    """
    PSD B-Element
    """
    def _psd_B_elem(self,line,bitload,gamma,k):
        return self._f(line.id,bitload,gamma)*(line.noise+line.alien_xtalk(k))/self.xtalk_gain[line.id][line.id][k]
                    
    """
    Transfer function lookup - |h_{i,j}|^2
    This more or less does nothing but return xtalk_gain but is closer to the math
    """
    def _h2(self,line1,line2,channel):
        return self.xtalk_gain[line1.id][line2.id][channel]
    
    """
    F- Gamma function - 
    f(line,k)=\Gamma(2^{b_n(k)} -1)
    :from psd_vector.c
    """
    def _f(self,lineid,bitload,gamma):
        assert isinstance(bitload,numpy.ndarray) and len(bitload)>0,"WTF is this bitload? %s,%s"%(str(bitload),type(bitload))
        assert(isinstance(gamma,float))
        g=gamma #TODO gamma values from symerr.c ? initially 9.95 for osb  
        b=bitload[lineid]
        return pow(10,(g+3)/10)*(pow(2,b)-1) #TODO initially, b=0, so this doesnt work
    """    
    Pretty Print channel Matrix
    #TODO I've got no idea how to display this....
    """
    def graph_channel_matrix(self):
        
        pylab.contourf(self.xtalk_gain[...,0])
        pylab.figure()
        pylab.show
    
    """
    Print Channel Matrix to file after generation
    """
    def print_cm(self,filename):
        f=open(filename,'w')
        
    """
    Print Channel Matrix to screen
    """
    def print_xtalk(self):
        for tone in range(self.K):
            print("\n%d:"%(tone)),
            for i in range(self.N):
                for j in range(self.N):
                    print("%e"%self.xtalk_gain[i][j][tone]), 
    
    """
    Print Channel Matrix to screen
    """
    def print_xtalk_onetone(self,tone):
        print("\n%d:"%(tone))
        for i in range(self.N):
            print("")
            for j in range(self.N):
                print("%e"%self.xtalk_gain[i][j][tone]), 
        