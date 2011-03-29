'''
Class that describes the DSL bundle object
'''

#Global imports
import sys
import re
import cmath
import numpy as np
import pyparsing
import scipy.special as ss
import pylab as pl
import pprint
import itertools
import hashlib
import os
from tempfile import mkdtemp
from math import pow,log10,sqrt

#Local Imports
from utility import *
from line import Line


class Bundle(object):
    MARGIN = 3.0   #Performance Margin db
    C_G = 0        #Coding Gain db
    GAMMA= 9.95
    UGAMMA= get_GAMMA(1e-7, 4)       #Uncoded SNR Gap (setup elsewere)
    gamma_hat = pow(10,(UGAMMA+MARGIN-C_G)/10)
    
    def __init__(self,network_file,K=512,scenarioname="NoScenario",rates_file=None, weights_file=None):
        self.N = 0                  #Filled in on loading file, here to remind me
        self.lines = []            # the DSL line objects
        self.xtalk_gain = []    # XT Matrix (NB Must be FULLY instantiated before operating)
        #Assuming that each bundle is going to be one medium
        self.K = int(K)                # the number of DMT channels
        self._psd_cache = {'hits':0,'misses':0}
        self.freq=np.asarray([140156.25 + 4312.5 * i for i in range(self.K)])
        


        '''
        Try to open and parse the network configuration file
        '''
        try:
            with open(network_file,"r") as nf:
                for n,line in enumerate(nf):
                    if '#' in line: #not used yet but could be handy later
                        scenarioline=re.split('#|:|,',line)
                        for key in range(1,len(scenarioline),2): 
                            self.scenariooptions[scenarioline[key]]=scenarioline[key+1]
                    else:
                        # nt and lt are network termination (CO or RT) and line termination (CPE)
                        nt,lt = line.split(",")
                    self.lines.append(Line(nt,lt,n,self))

            self.N = len(self.lines)
            log.info("Read %d Lines",self.N)
            

        except IOError:
            log.error("Cannot open the network file: %s",network_file)
            sys.exit(1)


        log.info("Calculating the channel matrix")
        #The real work begins
        self.calc_channel_matrix()

        #self.graph_channel_matrix()
        log.info("Running self check:")
        self.check_xtalk_gains() #This is only ever used once; could be sent into calc_channel_matrix?
        self.tofile(scenarioname)
        


    '''
    Calculates the bundle channel gain matrix, generating xtalk_gain[][][]
    :from channel_matrix.c
    '''
    def calc_channel_matrix(self):
        #Initialise the gains
        self.xtalk_gain = np.zeros((self.N,self.N,self.K))
        for k in range(self.K):                         #For Each Tone
            for x,lx in enumerate(self.lines):          #Between Every xtalker
                for v,lv in enumerate(self.lines):      # and every victim
                    if x == v:                          #If you're talking to yourself, do lazy transfer_fn
                        lx.gain[k] = self.xtalk_gain[x][v][k] = lx.transfer_fn(self.freq[k])
                    else:                               #Otherwise look at XT
                        self.xtalk_gain[x][v][k] = self.calc_fext_xtalk_gain(lx,lv,self.freq[k],"DOWNSTREAM") #This makes more sense in passing line objects instead of id's
                   
    '''
    Check Normalised XT Gains and xtalk symmetry 
    :from channel_matrix.c
    Normalised Gains arn't really needed but its just as easy to keep it in
    XTalk symmetry needs to be checked because if all the lines are NOT the same
    the xtalk_gain matrix can NOT be symmetric in the [i,j] axis
    '''
    def check_xtalk_gains(self):
        yeses=0
        ''' Original Way to do it
        for k in self.K:            #tone
            for x,lineX in enumerate(self.lines):    #xtalker
                for v,lineV in enumerate(self.lines):#victim
                    if x==v: continue
                    if self.xtalk_gain[x,v,k]/lineV.gain[k] > 0.5:
                        yeses+=1
        '''
        
        '''Lets turn this on its head'''
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
            for k in range(self.K):
                assert (self.xtalk_gain[:,:,k].T!=self.xtalk_gain[:,:,k]).all()==False, "Xtalk Symmtric on tone %d"%k

        else:
            log.info("Lines are identical, xtalk should be symmetric")
            for k in range(self.K):
                assert (self.xtalk_gain[:,:,k].T!=self.xtalk_gain[:,:,k]).any()==True, "Xtalk Not Symmtric on tone %d"%k
              
            
        log.info("Total:%d,%%Yes:%f%%"%(len(gainratio),yeses/(1.0*len(gainratio))))
    
    '''
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
    
    '''
    def calc_fext_xtalk_gain(self,victim,xtalker,freq,dir): #TODO Currently does UPSTREAM only, swap lt/nt for victim and xtalker for downstream
               
        if dir == "DOWNSTREAM":                 #If DOWNSTREAM, swap  and negate nt/lt's
            (D1,D2)=(-xtalker.lt,-xtalker.nt)
            (A,B)=(-victim.lt,-victim.nt)
        else: #Upstream
            (D1,D2)=(xtalker.nt,xtalker.lt)
            (A,B)=(victim.nt,victim.lt)
            
        #Shared length is the H2 Span
        shared_length=abs(max(A,D1)-min(B,D2))/1000
        #Head Length is the H1/H4 Span
        head_length=(A-D1)/1000
        #Tail Length is the H3 Span
        tail_length=(B-D2)/1000
        
        h1 = self.insertion_loss(head_length, freq)
        h2 = self.insertion_loss(shared_length, freq)
        h3 = self.insertion_loss(tail_length, freq)
        '''
        This _H takes away the biggest pain of the fext gain cases; 
        H1 > 0 in cases 1,3,4,6,7
        H2 active in all cases, but using max(v.nt,x.lt)
        H3 > 0 in cases 1,2,3
        H4 Not relevent for FEXT
        
        NB, insertion_loss(<0,f)=1
        
        I think AMK's Case1/9 fext(length) is wrong. Should be the
        common length between the victim and xtalker, not the combined
        length
        '''
        H = h1*h2*h3
        
        gain = H * self.fext(freq, shared_length)
        return gain  
    
    '''
    Calculate FEXT for given Freq and Length in m
    Uses:
        UndB
    Model from BT's simulation parameters document cp38-2 http://www.niccstandards.org.uk/files/current/NICC%20ND%201513%20(2010-01).pdf
    :from channel_matrix.c
    '''
    def fext(self,freq,length):
       
        x=-55   #This is different from the NICC Model  #FUDGE
        x+= 20*log10(freq/90e3)         #f in Hz
        x+= 6*log10(float(self.N-1)/49)    #n disturbers 
        x+= 10*log10(length)       #shared length in km
        
        try:
            return UndB(x)
        except ZeroDivisionError:
            return 0
        
    '''
    Calculate Insertion Loss
    :from transfer_fn.c
    '''    
    def insertion_loss(self,length,freq): #TODO: Reintroduce factor removed from C version?
        '''
        insertion loss makes more sense in the bundle as most of the time,
        it is looking at a common length between two lines
        '''
        if length > 0:
            return do_transfer_function(length,freq,measure="km")
        else: 
            return 1 #Leads straight into multiplication; so ends up x*0=0 otherwise
    
    '''
    Calculate snr statistics for all lines in bundle
    :from snr.c
    '''
    def calculate_snr(self):
        for line in self.lines:
            line.sanity()
            noise = np.zeros(self.K)
            
            log.debug("fext\tgain\tcnr\tsnr")
            for tone in xrange(self.K):
                
                noise = line.calc_fext_noise(tone) + line.alien_xtalk(tone) + dbmhz_to_watts(line.noise)
                line.cnr[tone] = line.gain[tone]/noise #gain calculated from xtalk_gain generation
                line.snr[tone] = dbmhz_to_watts(line.p[tone])*noise
                
                line.gamma_m[tone] = TodB(line.snr[tone]/pow(2,line.b[tone]-1))
                
            #line.symerr = [ self._calc_sym_err(line,xtalker) for xtalker in xrange(self.K)] #TODO _calc_sym_err
            line.p_total = sum(map(dbmhz_to_watts,line.p))
            line.b_total = sum(line.b)
            log.info("Line:%d,Power:%fW,Rate:%dbps"%(line.id,line.p_total,line.b_total))
            
            '''
            With This whole Fractional thing; 
            Are b/_b ever different? If its a fractional line is b[k] ever needed?
            '''
    '''
    Calculate Symbol Error Rate on line
    :from symerr.c
    '''
    def _calc_sym_err(self,line,k): #TODO
        
        M=pow(2,line.b[k])
               
        return 1 - (1 - (1 - 1/sqrt(M))*ss.erf(sqrt((3*line.snr[k])/(2*(M-1)))))
    
    '''
    Checks inter-service margins
    :from snr.c
    #I suspect this can be ignored
    '''
    def check_all_margins(self,gap):
        pass
            
    
    '''
    Generate PSD vector matrix between lines and return matrix
    :from psd_vector.c
    '''
    def calc_psd(self,bitload,k):
        #Caching implementation (AMK)
        key = hashlib.sha1(bitload.view()).hexdigest()+str(k)
        try:
            ret = self._psd_cache[key]
            self._psd_cache['hits']+=1
            return ret
        except:
            pass
        
        #Generate Matrices (See Intro.pdf 2.23)
        A=np.asmatrix(np.zeros((self.N,self.N)))
        B=np.asmatrix(np.zeros(self.N))
        
        log.debug("Channel:%d,Bitload:%s"%(k,str(bitload)))
        
        for v in range(self.N): #victims
            B[0,v]=pow(10,(self.GAMMA+3)/10)*(pow(2,bitload[v])-1)*(dbmhz_to_watts(-140)/self.xtalk_gain[v][v][k])                
            for x in range(self.N): #xtalkers
                if v==x: 
                    A[v,x]=1
                else:
                    A[v,x]=(-1.0*pow(10,(self.GAMMA+3)/10)*(pow(2,bitload[v])-1)*self.xtalk_gain[x][v][k])/self.xtalk_gain[v][v][k]
        
        B=B.T

        #Everyone loves linear algebra...dont they?
        P=np.linalg.solve(A,B)
        
        #Useful debugging
        '''
        log.debug("A:\n%s"%str(A))
        log.debug("B:\n%s"%str(B))
        log.debug("P:\n%s"%str(P))
        '''
        
        
        P=P.T

        P=mat2arr(P[0])

        self._psd_cache[key]=P   
        self._psd_cache['misses']+=1
     
        return P 


    '''
    Transfer function lookup - |h_{i,j}|^2
    This more or less does nothing but return xtalk_gain but is closer to the math
    '''
    def _h2(self,line1,line2,channel):
        return self.xtalk_gain[line1.id][line2.id][channel]
    
    '''
    F- Gamma function - 
    f(line,k)=\Gamma(2^{b_n(k)} -1)
    :from psd_vector.c
    #TODO Memoize
    '''
    def _f(self,bitload,gamma=GAMMA):
        key = "%d-%f"%(bitload,gamma)
        if key in self._f_cache:
            return self._f_cache[key]
        else:
            result=pow(10,(gamma+3)/10)*(pow(2,bitload)-1)
            self._f_cache[key]=result
        return result #TODO initially, b=0, so this doesnt work
    '''    
    Pretty Print channel Matrix
    #TODO I've got no idea how to display this....
    '''
    def graph_channel_matrix(self):
        
        pl.contourf(self.xtalk_gain[...,0])
        pl.figure()
        pl.show
    
    '''
    Print CM to file
    '''
    def tofile(self,filename):
        resultsdir=rawdir
        if not os.path.isdir(resultsdir):
            os.makedirs(resultsdir)
        np.save(resultsdir+filename+'-channelmatrix', self.xtalk_gain)
        
    '''
    Print Channel Matrix to screen
    '''
    def print_xtalk(self):
        for tone in range(self.K):
            print("\n%d:"%(tone)),
            for x in range(self.N):
                for v in range(self.N):
                    print("%e"%self.xtalk_gain[x][v][tone]), 
    
    '''
    Print Channel Matrix to screen
    '''
    def print_xtalk_onetone(self,tone):
        print("\n%d:"%(tone))
        for x in range(self.N):
            print("")
            for v in range(self.N):
                print("%e"%self.xtalk_gain[x][v][tone]),   