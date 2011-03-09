"""
Algorithm Modules
"""

#Global Imports
import numpy
import math

#Local Imports
from bundle import Bundle
from algorithm import Algorithm
import utility

class OSB(Algorithm):
    """
    Optimum Spectrum Balancing
    """
    def run(self):
        self.name = "OSB"
        #Aimed-for rates per-line
        self.rate_targets = numpy.tile(100,self.bundle.N)
        #Tolerance of rate-targets (how close is close enough)
        self.target_tolerance=1
        
        self.preamble
        self.w=[]
        #How the fuck does this iterate over the bundle?
        for linea in self.bundle.lines:
            for lineb in self.bundle.lines:
                #Ignore comparing lines to themselves
                if linea.id != lineb.id:
                    self._optimise_w(linea,lineb)
        self.postscript
        return  
    
    """
    Optimise W -Currently not used
    """
    def optimise_weights(self,linea,lineb):
        
        #FIXME There is no way this makes sense for multi-user system >2
        if (self.rate_targets[lineb.id]==0):
            self._optimise_w(linea,lineb)
        elif (self.rate_targets[lineb.id]==0):
            #Swap and carry on regardless
            self._optimise_w(lineb,linea)
        else:
            utility.log.error("What the fuck am I supposed to do now?!")
        
    
    """
    Optimise W- working function
    
    """
    def _optimise_w(self,linea,lineb):
        w_max=1
        w_min=0
        while (abs(linea.rate()-self.rate_targets[linea.id]) > self.target_tolerance):
            w=(w_max+w_min)/2
            self.optimise_l1(w,linea,lineb)
            if linea.rate() > self.rate_targets[linea.id]:
                w_max=w
            else:
                w_min=w
    
    """
    Optimise Lambda 1
    :from OSB_original.pdf paper
    """
    def optimise_l1(self,w,linea,lineb):
        l1_max=1
        l1_min=0
        power=-1.0 #initialised so that converged will work
        #First need to find the max permissable value of lambda_n (l1_max)
        while sum(linea.p) > linea.p_max:
            l1_max = 2 * l1_max  #This could be replaced by a bitshift?
            self.optimise_l2(w,l1_max,linea,lineb)
        while not self._converged(linea,power): #TODO
            l1 = (l1_max + l1_min)/2
            self.optimise_l2(w,l1,linea,lineb)
            power=sum(linea.p)
            if power > linea.p_max:
                l1_min=l1
            else:
                l1_max=l1
            
                
    """
    Optimise Lambda 2
    :from OSB_original.pdf paper
    """
    def optimise_l2(self,w,l1,linea,lineb):
        l2_max=1
        l2_min=0
        power=-1.0 #initialised so that converged will work
        #First need to find the max permissable value of lambda_n (l1_max)
        while sum(lineb.p) > lineb.p_max:
            l2_max = 2 * l2_max  #This could be replaced by a bitshift?
            self.optimise_s(w,l1,l2_max,linea,lineb)
            
        while not self._converged(lineb,power): #TODO
            l2 = (l2_max + l2_min)/2
            self.optimise_p(w,l1,l2,linea,lineb)
            power=sum(lineb.p)
            if power > lineb.p_max:
                l2_min=l2
            else:
                l2_max=l2
                
    """
    Optimise Power (aka optimise_s)
    :from OSB_original.pdf paper
    """   
    def optimise_p(self,w,l1,l2,linea,lineb):
        L_max=-1.0 #max L_k seen
        b1_max = 0 #overall max for b1
        b2_max = 0 #overall max for b2
        
        #for each subchannel
        for k in range(self.bundle.K):
            #and each bit possibility
            for b1 in range(self.MAXBITSPERTONE):
                #for both lines
                for b2 in range(self.MAXBITSPERTONE):
                    L_k=self._L(w, l1, l2, linea, lineb, k, b1, b2)
                    if L_k > L_max:
                        L_max = L_k
                        #Store these new maxe's
                        linea.b[k] = b1_max = b1
                        lineb.b[k] = b2_max = b2
                        
                    #Powers need to be updated at some point. see PSD_vector.h
            #end of search double loop on b's
            #Update 'best' bit loads
            linea.b[k]=b1_max
            lineb.b[k]=b2_max
            #Update powers
            self.update_psds(linea)
            self.update_psds(lineb)
            #The reason why this seems to make sense HERE is that
            #in theory, the above double loop would be 'instant'.
            

    """
    Convergence check for both lambda optimisations
    """                 
    def _converged(self,line,lastpower):
        #How much convergence is satisfactory?
        e=0.01
        return (abs(sum(line.p)-lastpower)<e)
    
    """
    Calculate the Lagrangian
    :from osb_original.pdf and osb.c
    """
    def _L(self,w,l1,l2,linea,lineb,k,b1,b2):
        result=0
        #Weighted Sections
        result+= (w*b1+(1-w)*b2)
        #Lagrangian Sections
        result-= (l1*linea.p[k])*(self.bundle.xtalk_gain[b1,b2,k])#this isn't supposed to be xtalk, its p_matrix. No idea where the fuck it comes from tho.
        result-= (l2*lineb.p[k])*(self.bundle.xtalk_gain[b1,b2,k])
        return result
