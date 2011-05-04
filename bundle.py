'''
Class that describes the DSL bundle object
'''

#Global imports
import sys
import re
import cmath
import numpy as np
import pyparsing
#import scipy.special as ss
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
from gpu import GPU

class Bundle(object):
    MARGIN = 3.0   #Performance Margin db
    C_G = 0        #Coding Gain db
    GAMMA= 9.95
    UGAMMA= get_GAMMA(1e-7, 4)       #Uncoded SNR Gap (setup elsewere)
    gamma_hat = pow(10,(UGAMMA+MARGIN-C_G)/10)
    NOISE=dbmhz_to_watts(-140)
    MAXBITSPERTONE=15
    
    def __init__(self,network_file="",K=224,scenarioname="NoScenario",cachefile=False,useGPU=False):
        self.N = 0                  #Filled in on loading file, here to remind me
        self.lines = []            # the DSL line objects
        self.xtalk_gain = []    # XT Matrix (NB Must be FULLY instantiated before operating)
        #Assuming that each bundle is going to be one medium
        self.K = int(K)                # the number of DMT channels
        self._psd_cache = {'hits':0,'misses':0}
                
        assert cachefile==False, "Someone doesn't know what they're doing!"
        
        self.freq=np.asarray([140156.25 + 4312.5 * i for i in range(self.K)])
        
        #Try if network file is a pre-generated npy pickle
        try:
            self.lines=np.load(network_file)
        #Will throw IOError if file is not pickle; so assume is regular network file
        except IOError:
            log.info("Reading fresh Network File %s"%network_file)
            try:
                with open(network_file,"r") as nf:
                    for n,line in enumerate(nf):
                        if '#' in line: #not used yet but could be handy later
                            scenarioline=re.split('#|:|,',line)
                            for key in range(1,len(scenarioline),2): 
                                self.scenariooptions[scenarioline[key]]=scenarioline[key+1]
                        else:
                            try:
                                # nt and lt are network termination (CO or RT) and line termination (CPE)
                                nt,lt,rate= line.split(",")
                                self.lines.append(Line(nt,lt,int(rate),n,self))
    
                            except ValueError:
                                nt,lt=line.split(",")
                                self.lines.append(Line(nt,lt,False,n,self))               
            except IOError:
                log.error("Cannot open the network file: %s",network_file)
                sys.exit(1)
        finally:
            self.N = len(self.lines)
            log.info("Successfully read %d lines from %s"%(self.N,network_file))

        if(useGPU):
            self.gpu=GPU(self)


        log.info("Calculating the channel matrix for %d channels"%self.K)
        #The real work begins
        self.calc_channel_matrix()

        #self.graph_channel_matrix()
        log.info("Running self check:")
        self.check_xtalk_gains() #This is only ever used once; could be sent into calc_channel_matrix?
        self.tofile(scenarioname)
        

    '''
    Return GAMMA and NOISE (for GPU stuff)
    '''
    def get_GAMMA(self,gamma=GAMMA):
        return (gamma)

    def get_NOISE(self,noise=NOISE):
        return (noise)
    
    def get_MBPT(self,mbpt=MAXBITSPERTONE):
        return mbpt
    
    '''
    Calculates the bundle channel gain matrix, generating xtalk_gain[][][]
    :from channel_matrix.c
    '''
    def calc_channel_matrix(self):
        #Initialise the gains
        self.xtalk_gain = np.zeros((self.K,self.N,self.N))
        for k in range(self.K):                         #For Each Tone
            for x,lx in enumerate(self.lines):          #Between Every xtalker
                for v,lv in enumerate(self.lines):      # and every victim
                    if x == v:                          #If you're talking to yourself, do lazy transfer_fn
                        lx.gain[k] = self.xtalk_gain[k][x][v] = lx.transfer_fn(self.freq[k])
                    else:                               #Otherwise look at XT
                        self.xtalk_gain[k][x][v] = self.calc_fext_xtalk_gain(lx,lv,self.freq[k],"DOWNSTREAM") #This makes more sense in passing line objects instead of id's
                   
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
            gainratio=[self.xtalk_gain[k,x,v]/victim.gain[k] for x in range(self.N) for k in range(self.K)]
            yeses+=len([1 for i in gainratio if i>0.5])
        
        #Check symmetry if all lines are the same
        ntlt=set([(line.nt,line.lt) for line in self.lines])
        log.info("%s"%str(ntlt))
        
        if len(ntlt)!=1:
            #Not all the lines are the same so the xtalk matrix cannot be symmetric
            log.info("Lines are different, xtalk should not be symmetric")
            for k in range(self.K):
                assert (self.xtalk_gain[k,:,:].T!=self.xtalk_gain[k,:,:]).all()==False, "Xtalk Symmtric on tone %d"%k

        else:
            log.info("Lines are identical, xtalk should be symmetric")
            for k in range(self.K):
                assert (self.xtalk_gain[k,:,:].T!=self.xtalk_gain[k,:,:]).any()==True, "Xtalk Not Symmtric on tone %d"%k
              
            
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
        shared_length=abs(max(A,D1)-min(B,D2))/1000.0
        #Head Length is the H1/H4 Span
        head_length=(A-D1)/1000.0
        #Tail Length is the H3 Span
        tail_length=(B-D2)/1000.0
        
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
        
        try:
            gain = H * self.fext(freq, shared_length)
        except ValueError:
            log.error("Failed on (A,B,D1,D2)=(%d,%d,%d,%d):Shared:%f"%(A,B,D1,D2,shared_length))
            raise
        
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
        #Cache Stats
        log.info("Cache Hits:%d"%self._psd_cache['hits'])
        log.info("Cache Mises:%d"%self._psd_cache['misses'])
        
        for line in self.lines:
            line.sanity()
            noise = np.zeros(self.K)
            
            log.debug("fext\tgain\tcnr\tsnr")
            for tone in xrange(self.K):
                try:
                    noise = line.calc_fext_noise(tone) + line.alien_xtalk(tone) + (line.noise)
                except OverflowError:
                    log.error("Overflowed on tone:%d:%d"%(line.id,tone))
                line.cnr[tone] = line.gain[tone]/noise #gain calculated from xtalk_gain generation
                line.snr[tone] = dbmhz_to_watts(line.p[tone])*noise
                
                line.gamma_m[tone] = TodB(line.snr[tone]/pow(2,line.b[tone]-1))
                
            #line.symerr = [ self._calc_sym_err(line,xtalker) for xtalker in xrange(self.K)] #TODO _calc_sym_err
            line.p_total = sum(map(dbmhz_to_watts,line.p))
            line.b_total = sum(line.b)
            log.info("Line:%d,Power:%fW,Rate:%dbpf"%(line.id,line.p_total,line.b_total))
            
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
              
        log.error("Tried to execute broken _calc_cym_error")
        raise ValueError
        return 1 - (1 - (1 - 1/sqrt(M)))#*ss.erf(sqrt((3*line.snr[k])/(2*(M-1)))))
    
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
    def calc_psd(self,bitload,k,gamma=GAMMA, precompute=True, noise=NOISE, gpu=False):
        #Caching implementation (AMK)
        key = hashlib.sha1(bitload.view()).hexdigest()+str(k)
        
        try:
            ret = self._psd_cache[key]
            self._psd_cache['hits']+=1
            return ret
        except:
            (not precompute) and log.info("Calculating PSD Outside of PreCompute")
            pass
        
        #Generate Matrices (See Intro.pdf 2.23)
        A=np.zeros((self.N,self.N))
        B=np.zeros((self.N,1))
        
        XTG=self.xtalk_gain[k,:,:]
        
        channelgap=pow(10,(gamma+3)/10)
        
        log.debug("Channel:%d,Bitload:%s"%(k,str(bitload)))
        for v in range(self.N): #victims
            B[v,0]=(noise*channelgap*(pow(2,bitload[v])-1)/XTG[v,v])
            A[v]=-(channelgap*(pow(2,bitload[v])-1)*XTG[:,v]/XTG[v,v])
            A[v,v]=1
        '''
            B[v,0]=((channelgap*pow(2,bitload[v]-1)*noise)/XTG[v,v]) #equiv f(b_v(k))*o / H_vv
            for x in range(self.N):
                A[v,x]=(-1.0*channelgap*(pow(2,bitload[v]-1))*XTG[x,v])/XTG[v,v]  #equiv -f(b_v(k))*H_xv/H_vv
            A[v,v]=1          #equiv if i==j=>Aij=1
        '''
    
        #Everyone loves linear algebra...dont they?
        if (gpu):
            P=self.gpu.solve(A,B,224) #FIXME currently useless and disabled
            log.debug("GPU Outside Scope:%d:%s:%s"%(k,str(P),str(bitload)))
        else:
            P=np.linalg.solve(A,B)
        
        #Useful debugging
        
        log.debug("A:\n%s"%str(A))
        log.debug("B:\n%s"%str(B))
        log.debug("P:\n%s"%str(P))
        
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
        return self.xtalk_gain[channel][line1.id][line2.id]
    
    '''    
    Pretty Print channel Matrix
    #TODO I've got no idea how to display this....
    '''
    def graph_channel_matrix(self):
        pl.contourf(self.xtalk_gain[...,0])
        pl.figure()
        pl.show
    
    '''
    Print CM and lines to file
    '''
    def tofile(self,filename):
        resultsdir=rawdir
        if not os.path.isdir(resultsdir):
            os.makedirs(resultsdir)
        np.save(resultsdir+filename+'-channelmatrix', self.xtalk_gain)
        #np.save(resultsdir+filename+'-lines', self.lines)
        
    '''
    Print Channel Matrix to screen
    '''
    def print_xtalk(self):
        for tone in range(self.K):
            print("\n%d:"%(tone)),
            for x in range(self.N):
                for v in range(self.N):
                    print("%e"%self.xtalk_gain[tone][x][v]), 
    
    '''
    Print Channel Matrix to screen
    '''
    def print_xtalk_onetone(self,tone):
        print("\n%d:"%(tone))
        for x in range(self.N):
            print("")
            for v in range(self.N):
                print("%e"%self.xtalk_gain[tone][x][v]),   
                
    
