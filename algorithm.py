'''
Algorithm Parent module 
'''

import sys
import utility as util
import math
import numpy as np
import pylab as pl
import time
import os

from bundle import Bundle



class Algorithm(object):
    '''
    Algorithm Constants; usually set within child classes
    '''
    MAXRATE = None #Need to set this for load_fm and ISB
    MAXPOWER = None
    MAXBITSPERTONE = 15
    name="Default Algorithm Name; You should never see this!"
    stats={}

        
    def __init__(self,bundle):
        assert isinstance(bundle,Bundle), "Did you give me a bundle to use?"
        self.bundle=bundle

    def preamble(self):
        '''
        Stuff that needs to be done for any algorithm, including 
        *initialisation of power and bitloading matrices
        *timing information
        '''
        #Power and bitloading are dim(KxN)
        self.p=np.zeros((self.bundle.K,self.bundle.N)) #per tone per user power in watts
        self.b=np.asmatrix(np.zeros((self.bundle.K,self.bundle.N)))     
        util.log.info("Lets get this party started!")
        self.stats['start']=time.time()
    
    def postscript(self):
        '''
        Stuff that should be done for any algorithm, including
        *downtyping of b -> ndarray
        *reassignment of relevent p's and b's to line values
        *SNR et al calculations
        *Timing statistics
        '''
        self.b=util.mat2arr(self.b)
        for line in self.bundle.lines:
            line.p=self.p[:,line.id]
            line.b=self.b[:,line.id]
        self.bundle.calculate_snr() #move into postscript?
        self.stats['end']=time.time()
        self.stats['duration']=self.stats['end']-self.stats['start']

        util.log.info("All Done Here, took %f"%self.stats['duration'])
        
    #FIXME I'm fairly sure nothing below here is layed out properly yet; am_load_ra is used in SOME algos...
    def am_load_ra(self,line):  #am_load.c used in IWF, Multiuser_greedy, RR_IWF, RR_OSB
        line.b = np.zeros(self.bundle.K)
        tone_full = np.tile(False,self.bundle.K)
        p = np.zeros(self.bundle.K)
        p_total=0
        b_total=0
        delta_p = []
        
        self._calc_delta_p(line)
        
        while (0 <= min(delta_p)): #TODO I know this doesn't work, Don't understand the breaking state see am_load.c find_min()
            tone = min(delta_p)
            line.b[tone]+=1
            p_total += delta_p[tone]
            p[tone] += delta_p[tone]
            
            if (line.b[tone] >= self.MAXBITSPERTONE):
                tone_full[tone] = True
            
            if (p_total > self.MAXPOWER):
                line.b[tone]-=1
                p_total -=delta_p[tone]
                break
            elif (p[tone] > util.dbmhz_to_watts(-20)): #if this tone is greater than the masking power
                line.b[tone]-=1
                p[tone]-=delta_p[tone]
                p_total -= delta_p[tone]
                tone_full[tone] = True
            
            self._calc_delta_p(line,tone_full)
        else:
            util.log.info("All Tones are full!") #Breaking statement where min(delta_p) < 0
        
        self.update_psds(line)
        
        
        line.service = np.zeros(self.bundle.K) #Still have no idea what this does, but reset it to zero anyway
        
        return line.rate()
    
    def am_load_fm(self,line,half=0):  #am_load.c used in IWF, Multiuser_greedy, RR_IWF, Single_IWF
        #half:(-1:bottom,0:both,1:top)
        '''
        Initialise Everything
        '''
        line.b = np.zeros(self.bundle.K)
        tone_full = np.tile(False,self.bundle.K)     #set all tones to full
        p_total=0
        b_total=0
               
        '''
        Calculate Initial bit/tone costs
        '''
        delta_p=self._calc_delta_p(line,tone_full)
        
        '''
        Deal with Halves
        ''' #Why is there a check on K==512?
        if ( half == 1 ): #Top Half
            for tone in range(self.bundle.K/2):
                tone_full[tone] = True
                util.log.debug("Top Half only")
        elif (half == -1): #Bottom Half
            for tone in range(self.bundle.K/2,self.bundle.K):
                tone_full[tone] = True        
                util.log.debug("Bottom Half only")
        
        while (0 < min(delta_p)):
            tone = np.argmin(delta_p) #return the tone with the lowest bit-adding cost
            line.b[tone]+=1
            b_total += 1            #Keep track of yer bits!
            p_total += delta_p[tone]
            
            '''Rate/Power/Bit Checking'''
            if (b_total == self.MAXRATE):
                util.log.debug("n:%d,k:%d - rate-satisfied"%(line.id,tone))
                break #I'm done          
            if (line.b[tone] >= self.MAXBITSPERTONE):
                tone_full[tone] = True
                util.log.info("n:%d,k:%d - tone full"%(line.id,tone))
            if (p_total > self.MAXPOWER):
                util.log.info("n:%d,k:%d - exceeded power budget, rolling back"%(line.id,tone))
                line.b[tone]-=1
                b_total -=1
                p_total -=delta_p[tone]
                break
            
            #Recalculate bit/tone costs
            delta_p=self._calc_delta_p(line,tone_full)
        else:
            util.log.info("All Tones are full!") #Breaking statement where min(delta_p) <= 0
        
        #Update powers
        self.update_psds(line)
                        
        line.service = np.zeros(self.bundle.K) #Still have no idea what this does, but reset it to zero anyway
        
        if b_total != self.MAXRATE: #This should be right as b_total is incrementally added
            util.log.error("Could not reach target data rate. Desired:",self.MAXRATE," Achieved:",b_total)                
        return b_total
    
    '''
    Calculate cost of additional bit on all the tones in a line
    :from am_load.c (implicit)
    Optimised from original version (Thank you wolframalpha)
    '''
    def _calc_delta_p(self,line,tone_full):
        delta_p=np.zeros(self.bundle.K)
        for tone in self.bundle.K:
            if not tone_full[tone]:
                delta_p[tone] = (pow(2,(line.b[tone]-1)) * 3 - 2 )* self.bundle.gamma_hat/line.cnr[tone]
            else:
                delta_p[tone] = 100
        return delta_p
                
    '''
    Update PSD's for this line
    :from am_load.c (implicit)
    Potentially integratable into calculatesnr?
    '''
    def update_psds(self,line):
        for tone in self.bundle.K:
            line.p[tone] = util.watts_to_dbmhz((pow(2,line.b[tone])-1)*(self.gamma_hat/line.cnr[tone]))
            if (line.p[tone] < line.MINPSD) and (line.b[tone] != 0):
                util.log.debug("Changing p from %e to MINPSD"%line.p[tone])
                line.p[tone] = line.MINPSD
                
    '''
    Print Power and Bitrate settings to file
    '''
    '''
    Print Channel Matrix, Power and Bitrates to file after generation
    '''
    def tofile(self,filename):
        assert os.path.isdir(util.rawdir)==True, "No "
        np.save(util.rawdir+filename+'-power', self.p)
        np.save(util.rawdir+filename+'-bitrate',self.b)
