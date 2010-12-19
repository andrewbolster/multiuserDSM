"""
Algorithm Parent module 
"""

import sys
import utility
import math
import numpy


class Algorithm(object):
    def __init__(self,bundle):
        self.name="Default Algorithm Name; You should never see this!"
        self.bundle=bundle
        
        """
        Algorithm Constants; usually set within child classes
        """
        self.MARGIN = 3.0
        self.C_G = 0
        self.MAXBITSPERTONE = 15
        self.MINPSD=60
        self.target_rate = None #Need to set this for load_fm and ISB
        pass
    
    def preamble(self):
        #Stuff that needs to be done for any algorithm
        pass
    
    def am_load_ra(self,line):  #am_load.c used in IWF, Multiuser_greedy, RR_IWF, RR_OSB
        line.b = numpy.zeros(self.bundle.K)
        tone_full = numpy.tile(False,self.bundle.K)
        gamma_hat = pow(10,(self.GAMMA+self.MARGIN-self.C_G)/10)
        p = numpy.zeros(self.bundle.K)
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
            
            if (p_total > self.p_budget):
                line.b[tone]-=1
                p_total -=delta_p[tone]
                break
            elif (p[tone] > utility.dbmhz_to_watts(-20)): #if this tone is greater than the masking power
                line.b[tone]-=1
                p[tone]-=delta_p[tone]
                p_total -= delta_p[tone]
                tone_full[tone] = True
            
            self._calc_delta_p(line)
        else:
            print("All Tones are full!") #Breaking statement where min(delta_p) < 0
        
        for tone in xrange(self.bundle.K):
            line.psd[tone] = utility.watts_to_dbmhz((math.pow(2,line.b[tone])-1)*(gamma_hat/line.cnr[tone])) #TODO gamma_hat/cnr[k] functionalise
            if (line.psd[tone] < self.MINPSD) and (line.b[tone] != 0):
                line.psd[tone] = self.MINPSD
        
        b_total = sum(line.b)
        
        line.service = numpy.zeros(self.bundle.K) #Still have no idea what this does, but reset it to zero anyway
        
        return b_total
    
    def am_load_fm(self,line,half=0):  #am_load.c used in IWF, Multiuser_greedy, RR_IWF, Single_IWF
        #half:(-1:bottom,0:both,1:top)
        
        """
        Need to discuss with AM merging am_load_* and passing in a 'satisfaction' function;
            am_load_fm == am_load(sumofline.b[]=b_target
            am_load_ra == am_load(sumofline.p[]=p_target
        """
        line.b = numpy.zeros(self.bundle.K)
        tone_full = numpy.tile(False,self.bundle.K)     #set all tones to full
        gamma_hat = pow(10,(self.GAMMA+self.MARGIN-self.C_G)/10)
        p_total=0
        b_total=0
        delta_p = []
        
        self._calc_delta_p(line)
        
        """Deal with Halves""" #Why is there a check on K==512?
        if ( half == 1 ): #Top Half
            for tone in range(self.bundle.K/2):
                tone_full[tone] = True
        elif (half == -1): #Bottom Half
            for tone in range(self.bundle.K/2,self.bundle.K):
                tone_full[tone] = True        
        
        
        while (0 <= min(delta_p)): #TODO I know this doesn't work, Don't understand the breaking state see am_load.c find_min()
            tone = min(delta_p)
            line.b[tone]+=1
            b_total += 1            #Keep track of yer bits!
            p_total += delta_p[tone]
            if (b_total == self.target_rate):
                break #I'm done          
            if (line.b[tone] >= self.MAXBITSPERTONE):
                tone_full[tone] = True
            
            if (p_total > self.p_budget):
                line.b[tone]-=1
                b_total -=1
                p_total -=delta_p[tone]
                break
            
            self._calc_delta_p(line)
        else:
            print("All Tones are full!") #Breaking statement where min(delta_p) < 0
        
        for tone in xrange(self.bundle.K):
            line.psd[tone] = utility.watts_to_dbmhz((math.pow(2,line.b[tone])-1)*(gamma_hat/line.cnr[tone])) #TODO gamma_hat/cnr[k] functionalise
                
        line.service = numpy.zeros(self.bundle.K) #Still have no idea what this does, but reset it to zero anyway
        
        if b_total != self.target_rate:
            print("Could not reach target data rate. Desired:",self.target_rate," Achieved:",b_total)                
        return b_total
        
    def _calc_delta_p(self,line): #TODO:Used in several am_load related functions; need to make variable safe
        for tone in self.bundle.K:
            if not self.tone_full[tone]:
                self.delta_p[tone] = (pow(2,(line.b[tone]-1)) * 3 - 2 )* self.gamma_hat/line.cnr[tone]
            else:
                self.delta_p[tone] = 100